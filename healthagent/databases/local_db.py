"""SQLite local database manager for genomics data.

Provides a thread-safe connection pool and helper methods for querying
the local copy of GWAS Catalog, ClinVar, PharmGKB, and wellness traits.
"""

import sqlite3
import threading
from pathlib import Path
from typing import Any, Optional

DB_PATH = Path(__file__).parent.parent.parent / "data" / "healthagent.db"
SCHEMA_PATH = Path(__file__).parent / "schema.sql"

_local = threading.local()


def _get_conn() -> sqlite3.Connection:
    """Return a per-thread SQLite connection, creating it if needed."""
    if not hasattr(_local, "conn") or _local.conn is None:
        DB_PATH.parent.mkdir(parents=True, exist_ok=True)
        conn = sqlite3.connect(str(DB_PATH), check_same_thread=False)
        conn.row_factory = sqlite3.Row
        conn.execute("PRAGMA journal_mode=WAL")
        conn.execute("PRAGMA foreign_keys=ON")
        _local.conn = conn
    return _local.conn


def init_db() -> None:
    """Create all tables from schema.sql if they don't exist."""
    schema = SCHEMA_PATH.read_text()
    conn = _get_conn()
    conn.executescript(schema)
    conn.commit()


def close() -> None:
    """Close the current thread's connection."""
    if hasattr(_local, "conn") and _local.conn:
        _local.conn.close()
        _local.conn = None


# ── Query helpers ─────────────────────────────────────────────────

def query(sql: str, params: tuple = ()) -> list[sqlite3.Row]:
    conn = _get_conn()
    cur = conn.execute(sql, params)
    return cur.fetchall()


def execute(sql: str, params: tuple = ()) -> int:
    """Execute a write statement; returns lastrowid."""
    conn = _get_conn()
    cur = conn.execute(sql, params)
    conn.commit()
    return cur.lastrowid


def executemany(sql: str, rows: list[tuple]) -> int:
    """Bulk insert; returns number of rows inserted."""
    conn = _get_conn()
    cur = conn.executemany(sql, rows)
    conn.commit()
    return cur.rowcount


# ── Domain queries ────────────────────────────────────────────────

_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def _geno_variants(geno: str) -> tuple[str, ...]:
    """Return all equivalent representations of a genotype.

    Handles both allele-order swaps (CT == TC) and strand flips
    (CC == GG on opposite strand). AncestryDNA and some other arrays
    report SNPs on the reverse complement strand for certain rsIDs.
    """
    comp = geno.translate(_COMPLEMENT)
    return (geno, geno[::-1], comp, comp[::-1])


def get_wellness_traits(rsids: list[str], genotypes: dict[str, str]) -> list[dict]:
    """
    Return wellness trait results for the given rsids, filtered to matching
    genotypes. genotypes maps rsid → called genotype (e.g. "AT").

    Handles strand flips: AncestryDNA reports some SNPs on the reverse strand,
    so GG in the file may match CC in the wellness table (complement).

    Queries wellness_trait first (small table, ~48 rows) then intersects with
    the profile genotype map — avoids hitting SQLite's 32k variable limit with
    large ancestry files.
    """
    if not rsids:
        return []

    # Build a set for fast membership check
    rsid_set = set(rsids)

    # Fetch all wellness trait rows (small table — never >500 rows)
    rows = query("SELECT * FROM wellness_trait")

    results = []
    for row in rows:
        if row["rsid"] not in rsid_set:
            continue
        geno = genotypes.get(row["rsid"], "")
        if not geno:
            continue
        if row["genotype"] in _geno_variants(geno):
            results.append(dict(row))
    return results


_SQLITE_VAR_LIMIT = 30_000   # SQLite default max_variable_number is 32766


def _load_rsids_temp(conn, rsids: list[str]) -> str:
    """Insert rsids into a temp table and return its name.

    Used when len(rsids) > _SQLITE_VAR_LIMIT to avoid the SQLite
    'too many SQL variables' error. The temp table is scoped to the
    current connection and is automatically dropped at session end.
    """
    conn.execute(
        "CREATE TEMP TABLE IF NOT EXISTS _profile_rsids (rsid TEXT PRIMARY KEY)"
    )
    conn.execute("DELETE FROM _profile_rsids")
    conn.executemany("INSERT OR IGNORE INTO _profile_rsids VALUES (?)", [(r,) for r in rsids])
    return "_profile_rsids"


def _rsid_clause(rsids: list[str]) -> tuple[str, tuple]:
    """Return (WHERE clause fragment, params) for an rsid list.

    Automatically switches to a JOIN on a temp table when the list is
    too large for SQLite's variable limit.
    Returns (clause, params) where clause is e.g. 'rsid IN (?,?,?)'.
    Caller must still do conn-level setup for the temp table variant.
    """
    if len(rsids) <= _SQLITE_VAR_LIMIT:
        ph = ",".join("?" * len(rsids))
        return f"rsid IN ({ph})", tuple(rsids)
    # Temp table path — caller must call _load_rsids_temp first
    return "rsid IN (SELECT rsid FROM _profile_rsids)", ()


def get_gwas_associations(rsids: list[str], limit: int = 5) -> list[dict]:
    """Return top GWAS associations for a list of rsids."""
    if not rsids:
        return []
    conn = _get_conn()
    if len(rsids) > _SQLITE_VAR_LIMIT:
        _load_rsids_temp(conn, rsids)
        clause, params = "rsid IN (SELECT rsid FROM _profile_rsids)", ()
    else:
        ph = ",".join("?" * len(rsids))
        clause, params = f"rsid IN ({ph})", tuple(rsids)
    cur = conn.execute(
        f"""SELECT rsid, trait, trait_category, p_value, odds_ratio, risk_allele,
                   study_title, pubmed_id
            FROM gwas_association
            WHERE {clause}
            ORDER BY p_value ASC LIMIT {limit * max(1, min(len(rsids), 500))}""",
        params,
    )
    return [dict(r) for r in cur.fetchall()]


def get_clinvar_variants(rsids: list[str]) -> list[dict]:
    """Return ClinVar clinical significance entries for a list of rsids."""
    if not rsids:
        return []
    conn = _get_conn()
    if len(rsids) > _SQLITE_VAR_LIMIT:
        _load_rsids_temp(conn, rsids)
        clause, params = "rsid IN (SELECT rsid FROM _profile_rsids)", ()
    else:
        ph = ",".join("?" * len(rsids))
        clause, params = f"rsid IN ({ph})", tuple(rsids)
    cur = conn.execute(
        f"""SELECT rsid, clinical_sig, condition, review_status, gene, molecular_consequence
            FROM clinvar_variant
            WHERE {clause}
            ORDER BY
                CASE clinical_sig
                    WHEN 'Pathogenic' THEN 1
                    WHEN 'Likely pathogenic' THEN 2
                    WHEN 'Uncertain significance' THEN 3
                    WHEN 'Likely benign' THEN 4
                    WHEN 'Benign' THEN 5
                    ELSE 6
                END""",
        params,
    )
    return [dict(r) for r in cur.fetchall()]


def get_drug_interactions(rsids: list[str] = None, genes: list[str] = None) -> list[dict]:
    """Return PharmGKB drug interactions for rsids or gene names."""
    conn = _get_conn()
    conditions, params = [], []

    if rsids:
        if len(rsids) > _SQLITE_VAR_LIMIT:
            _load_rsids_temp(conn, rsids)
            conditions.append("rsid IN (SELECT rsid FROM _profile_rsids)")
        else:
            placeholders = ",".join("?" * len(rsids))
            conditions.append(f"rsid IN ({placeholders})")
            params.extend(rsids)

    if genes:
        placeholders = ",".join("?" * len(genes))
        conditions.append(f"gene IN ({placeholders})")
        params.extend(genes)

    if not conditions:
        return []
    where = " OR ".join(conditions)
    cur = conn.execute(
        f"""SELECT gene, drug_name, phenotype, significance, plain_english,
                   category, rsid
            FROM drug_interaction
            WHERE {where}
            ORDER BY
                CASE significance
                    WHEN '1A' THEN 1 WHEN '1B' THEN 2
                    WHEN '2A' THEN 3 WHEN '2B' THEN 4
                    WHEN '3'  THEN 5 ELSE 6
                END""",
        tuple(params),
    )
    return [dict(r) for r in cur.fetchall()]


def get_disgenet_diseases(genes: list[str], min_score: float = 0.2) -> list[dict]:
    """Return top DisGeNET gene-disease associations for given genes."""
    if not genes:
        return []
    ph = ",".join("?" * len(genes))
    rows = query(
        f"""SELECT gene_symbol, disease_name, disease_class, score, ei
            FROM disgenet_association
            WHERE gene_symbol IN ({ph}) AND score >= ?
            ORDER BY score DESC LIMIT {len(genes) * 5}""",
        tuple(genes) + (min_score,),
    )
    return [dict(r) for r in rows]


def get_opentargets(genes: list[str]) -> list[dict]:
    """Return OpenTargets gene-disease associations for given genes."""
    if not genes:
        return []
    ph = ",".join("?" * len(genes))
    rows = query(
        f"""SELECT gene_symbol, disease_name, overall_score,
                   genetic_score, drug_score
            FROM opentargets_association
            WHERE gene_symbol IN ({ph})
            ORDER BY overall_score DESC LIMIT {len(genes) * 5}""",
        tuple(genes),
    )
    return [dict(r) for r in rows]


def get_ensembl_consequences(rsids: list[str]) -> list[dict]:
    """Return Ensembl VEP consequences for given rsids."""
    if not rsids:
        return []
    conn = _get_conn()
    if len(rsids) > _SQLITE_VAR_LIMIT:
        _load_rsids_temp(conn, rsids)
        clause, params = "rsid IN (SELECT rsid FROM _profile_rsids)", ()
    else:
        ph = ",".join("?" * len(rsids))
        clause, params = f"rsid IN ({ph})", tuple(rsids)
    cur = conn.execute(
        f"""SELECT rsid, gene_symbol, consequence, impact,
                   sift_prediction, polyphen_pred
            FROM ensembl_consequence
            WHERE {clause} AND impact IN ('HIGH','MODERATE')
            ORDER BY CASE impact WHEN 'HIGH' THEN 1 WHEN 'MODERATE' THEN 2 ELSE 3 END""",
        params,
    )
    return [dict(r) for r in cur.fetchall()]


def get_biobank_phenotypes(rsids: list[str]) -> list[dict]:
    """Return biobank phenotype associations for given rsids."""
    if not rsids:
        return []
    conn = _get_conn()
    if len(rsids) > _SQLITE_VAR_LIMIT:
        _load_rsids_temp(conn, rsids)
        clause, params = "rsid IN (SELECT rsid FROM _profile_rsids)", ()
    else:
        ph = ",".join("?" * len(rsids))
        clause, params = f"rsid IN ({ph})", tuple(rsids)
    cur = conn.execute(
        f"""SELECT rsid, phenotype, phenotype_category, study_name,
                   n_cases, p_value, beta
            FROM biobank_phenotype
            WHERE {clause}
            ORDER BY p_value ASC LIMIT 500""",
        params,
    )
    return [dict(r) for r in cur.fetchall()]


def get_db_stats() -> dict:
    """Return row counts per table for status reporting."""
    tables = ["snp", "gwas_association", "clinvar_variant", "drug_interaction",
              "wellness_trait", "disgenet_association", "opentargets_association",
              "ensembl_consequence", "biobank_phenotype", "download_log"]
    stats = {}
    for t in tables:
        try:
            row = query(f"SELECT COUNT(*) AS n FROM {t}")
            stats[t] = row[0]["n"] if row else 0
        except Exception:
            stats[t] = 0
    return stats


def save_profile(profile: "DNAProfile") -> int:
    """Persist a DNAProfile's SNPs to the snp table.

    Uses INSERT OR REPLACE so re-uploading the same file is safe.
    Returns number of SNPs saved.

    This makes the profile survive server restarts — on startup,
    load_profile() can restore the session from the database.
    """
    rows = [
        (s.rsid, s.chromosome, str(s.position), s.genotype)
        for s in profile.snps
    ]
    if not rows:
        return 0

    # Clear existing stored profile, insert new one
    conn = _get_conn()
    conn.execute("DELETE FROM snp")

    executemany(
        """INSERT OR REPLACE INTO snp (rsid, chromosome, position, alleles)
           VALUES (?, ?, ?, ?)""",
        rows,
    )

    # Store format and filename for reload
    conn.execute("DELETE FROM download_log WHERE source = 'profile_meta'")
    conn.commit()
    execute(
        "INSERT INTO download_log(source, status, records_added) VALUES (?,?,?)",
        ("profile_meta", profile.source_format.value, len(rows)),
    )
    return len(rows)


def load_profile() -> "Optional[DNAProfile]":
    """Reload the last uploaded DNA profile from the snp table.

    Returns None if no profile has been saved.
    """
    from healthagent.dna_importer import DNAProfile, SNP, DNAFormat

    meta = query(
        "SELECT status, records_added FROM download_log WHERE source='profile_meta' ORDER BY id DESC LIMIT 1"
    )
    if not meta:
        return None

    fmt_val  = meta[0]["status"]
    expected = meta[0]["records_added"]

    rows = query("SELECT rsid, chromosome, position, alleles FROM snp")
    if not rows:
        return None

    try:
        fmt = DNAFormat(fmt_val)
    except ValueError:
        fmt = DNAFormat.ANCESTRY_DNA

    snps = [
        SNP(
            rsid=r["rsid"],
            chromosome=r["chromosome"],
            position=int(r["position"] or 0),
            genotype=r["alleles"] or "",
        )
        for r in rows
    ]

    return DNAProfile(source_format=fmt, snps=snps, metadata={"restored": True, "snp_count": expected})


def log_download(source: str, status: str, records: int = 0, error: str = None) -> int:
    if status == "started":
        return execute(
            "INSERT INTO download_log(source, status) VALUES (?,?)",
            (source, "started"),
        )
    else:
        execute(
            """UPDATE download_log SET status=?, records_added=?, error_msg=?,
               finished_at=datetime('now')
               WHERE id=(SELECT MAX(id) FROM download_log WHERE source=?)""",
            (status, records, error, source),
        )
        return 0
