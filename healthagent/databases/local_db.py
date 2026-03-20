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

def get_wellness_traits(rsids: list[str], genotypes: dict[str, str]) -> list[dict]:
    """
    Return wellness trait results for the given rsids, filtered to matching
    genotypes. genotypes maps rsid → called genotype (e.g. "AT").
    """
    if not rsids:
        return []
    placeholders = ",".join("?" * len(rsids))
    rows = query(
        f"SELECT * FROM wellness_trait WHERE rsid IN ({placeholders})",
        tuple(rsids),
    )
    results = []
    for row in rows:
        geno = genotypes.get(row["rsid"], "")
        if row["genotype"] in (geno, geno[::-1]):
            results.append(dict(row))
    return results


def get_gwas_associations(rsids: list[str], limit: int = 5) -> list[dict]:
    """Return top GWAS associations for a list of rsids."""
    if not rsids:
        return []
    placeholders = ",".join("?" * len(rsids))
    rows = query(
        f"""
        SELECT rsid, trait, trait_category, p_value, odds_ratio, risk_allele,
               study_title, pubmed_id
        FROM gwas_association
        WHERE rsid IN ({placeholders})
        ORDER BY p_value ASC
        LIMIT {limit * len(rsids)}
        """,
        tuple(rsids),
    )
    return [dict(r) for r in rows]


def get_clinvar_variants(rsids: list[str]) -> list[dict]:
    """Return ClinVar clinical significance entries for a list of rsids."""
    if not rsids:
        return []
    placeholders = ",".join("?" * len(rsids))
    rows = query(
        f"""
        SELECT rsid, clinical_sig, condition, review_status, gene, molecular_consequence
        FROM clinvar_variant
        WHERE rsid IN ({placeholders})
        ORDER BY
            CASE clinical_sig
                WHEN 'Pathogenic' THEN 1
                WHEN 'Likely pathogenic' THEN 2
                WHEN 'Uncertain significance' THEN 3
                WHEN 'Likely benign' THEN 4
                WHEN 'Benign' THEN 5
                ELSE 6
            END
        """,
        tuple(rsids),
    )
    return [dict(r) for r in rows]


def get_drug_interactions(rsids: list[str] = None, genes: list[str] = None) -> list[dict]:
    """Return PharmGKB drug interactions for rsids or gene names."""
    conditions, params = [], []
    if rsids:
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
    rows = query(
        f"""
        SELECT gene, drug_name, phenotype, significance, plain_english,
               category, rsid
        FROM drug_interaction
        WHERE {where}
        ORDER BY
            CASE significance
                WHEN '1A' THEN 1 WHEN '1B' THEN 2
                WHEN '2A' THEN 3 WHEN '2B' THEN 4
                WHEN '3'  THEN 5 ELSE 6
            END
        """,
        tuple(params),
    )
    return [dict(r) for r in rows]


def get_db_stats() -> dict:
    """Return row counts per table for status reporting."""
    tables = ["snp", "gwas_association", "clinvar_variant",
              "drug_interaction", "wellness_trait", "download_log"]
    stats = {}
    for t in tables:
        try:
            row = query(f"SELECT COUNT(*) AS n FROM {t}")
            stats[t] = row[0]["n"] if row else 0
        except Exception:
            stats[t] = 0
    return stats


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
