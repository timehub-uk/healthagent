"""TCGA (The Cancer Genome Atlas) on-demand client.

Queries the NCI GDC (Genomic Data Commons) API for the user's specific
DNA markers only — no bulk downloads. Results are cached locally in SQLite.

What TCGA is:
    The Cancer Genome Atlas is a landmark cancer genomics program that
    molecularly characterized over 20,000 primary cancers across 33 cancer
    types. All somatic mutation, copy number, and expression data is freely
    available via the NCI GDC API.

Why it matters for users:
    If a user's SNP appears in TCGA somatic mutation data, it can tell them:
    - Which cancer types that variant has been observed in
    - How frequently it occurs (mutation frequency)
    - Whether it is a known driver mutation or a passenger
    This does NOT mean the user has cancer — germline SNPs from ancestry
    files are different from somatic mutations. We clearly explain this.

API documentation:
    https://docs.gdc.cancer.gov/API/Users_Guide/Getting_Started/
    Base URL: https://api.gdc.cancer.gov
    No API key required for open-access data.
"""

import json
import logging
import time
import urllib.request
from typing import Optional

from healthagent.databases import local_db

log = logging.getLogger(__name__)

GDC_BASE = "https://api.gdc.cancer.gov"

# ── Schema extension ──────────────────────────────────────────────

TCGA_SCHEMA = """
CREATE TABLE IF NOT EXISTS tcga_mutation (
    id                  INTEGER PRIMARY KEY AUTOINCREMENT,
    rsid                TEXT,
    gene_symbol         TEXT,
    mutation_id         TEXT,          -- GDC mutation UUID
    cancer_type         TEXT,          -- e.g. "Breast Invasive Carcinoma"
    cancer_abbr         TEXT,          -- e.g. "BRCA"
    mutation_type       TEXT,          -- Missense_Mutation, Nonsense_Mutation, etc.
    variant_class       TEXT,          -- SNP, DEL, INS
    amino_acid_change   TEXT,          -- e.g. "p.Val600Glu"
    case_count          INTEGER,       -- how many TCGA cases have this mutation
    frequency           REAL,          -- case_count / total_cases_in_project
    cosmic_id           TEXT,          -- COSMIC mutation ID if known
    is_known_driver     INTEGER DEFAULT 0,
    plain_english       TEXT,          -- consumer-friendly explanation
    queried_at          TEXT DEFAULT (datetime('now')),
    UNIQUE(rsid, cancer_abbr)
);
CREATE INDEX IF NOT EXISTS idx_tcga_rsid ON tcga_mutation(rsid);
CREATE INDEX IF NOT EXISTS idx_tcga_gene ON tcga_mutation(gene_symbol);
"""

# TCGA project codes → readable cancer names
CANCER_NAMES = {
    "TCGA-BRCA": "Breast Cancer",
    "TCGA-LUAD": "Lung Adenocarcinoma",
    "TCGA-LUSC": "Lung Squamous Cell Carcinoma",
    "TCGA-COAD": "Colon Adenocarcinoma",
    "TCGA-READ": "Rectal Adenocarcinoma",
    "TCGA-PRAD": "Prostate Cancer",
    "TCGA-THCA": "Thyroid Cancer",
    "TCGA-BLCA": "Bladder Urothelial Carcinoma",
    "TCGA-KIRC": "Kidney Renal Clear Cell Carcinoma",
    "TCGA-UCEC": "Uterine Corpus Endometrial Carcinoma",
    "TCGA-GBM":  "Glioblastoma",
    "TCGA-OV":   "Ovarian Serous Cystadenocarcinoma",
    "TCGA-HNSC": "Head and Neck Squamous Cell Carcinoma",
    "TCGA-LGG":  "Brain Lower Grade Glioma",
    "TCGA-STAD": "Stomach Adenocarcinoma",
    "TCGA-SKCM": "Skin Cutaneous Melanoma",
    "TCGA-CESC": "Cervical Squamous Cell Carcinoma",
    "TCGA-LIHC": "Liver Hepatocellular Carcinoma",
    "TCGA-KIRP": "Kidney Renal Papillary Cell Carcinoma",
    "TCGA-SARC": "Sarcoma",
    "TCGA-PAAD": "Pancreatic Adenocarcinoma",
    "TCGA-PCPG": "Pheochromocytoma and Paraganglioma",
    "TCGA-ESCA": "Esophageal Carcinoma",
    "TCGA-TGCT": "Testicular Germ Cell Tumors",
    "TCGA-THYM": "Thymoma",
    "TCGA-MESO": "Mesothelioma",
    "TCGA-UVM":  "Uveal Melanoma",
    "TCGA-UCS":  "Uterine Carcinosarcoma",
    "TCGA-DLBC": "Lymphoid Neoplasm Diffuse Large B-cell Lymphoma",
    "TCGA-CHOL": "Cholangiocarcinoma",
    "TCGA-ACC":  "Adrenocortical Carcinoma",
    "TCGA-KICH": "Kidney Chromophobe",
    "TCGA-LAML": "Acute Myeloid Leukaemia",
}

KNOWN_DRIVERS = {
    "TP53", "KRAS", "PIK3CA", "PTEN", "APC", "BRAF", "EGFR", "IDH1",
    "IDH2", "RB1", "CDKN2A", "MET", "ALK", "RET", "BRCA1", "BRCA2",
    "ATM", "ARID1A", "CDH1", "ERBB2", "FGFR1", "FGFR2", "FGFR3",
    "NF1", "NF2", "SMAD4", "STK11", "TSC1", "TSC2", "VHL",
}


def _gdc_post(endpoint: str, payload: dict, timeout: int = 20) -> dict:
    """POST to GDC API and return parsed JSON."""
    url  = f"{GDC_BASE}/{endpoint}"
    data = json.dumps(payload).encode()
    req  = urllib.request.Request(
        url, data=data,
        headers={"Content-Type": "application/json", "Accept": "application/json",
                 "User-Agent": "HealthAgent/1.0"},
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read())


def _gdc_get(endpoint: str, params: dict = None, timeout: int = 20) -> dict:
    """GET from GDC API and return parsed JSON."""
    from urllib.parse import urlencode
    url = f"{GDC_BASE}/{endpoint}"
    if params:
        url += "?" + urlencode(params)
    req = urllib.request.Request(
        url, headers={"Accept": "application/json", "User-Agent": "HealthAgent/1.0"}
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read())


def init_tcga_schema() -> None:
    """Create TCGA table if it doesn't exist."""
    local_db.init_db()
    local_db._get_conn().executescript(TCGA_SCHEMA)
    local_db._get_conn().commit()


def query_tcga_for_rsids(
    rsids: list[str],
    genes: list[str] = None,
    progress_cb=print,
) -> list[dict]:
    """Query TCGA/GDC for specific rsIDs and return cancer mutation context.

    Only fetches data for the provided rsIDs — no bulk downloads.
    Results are cached in SQLite; re-queried SNPs return cached data.

    Args:
        rsids:       List of rsIDs from the user's DNA profile.
        genes:       Optional list of gene names to supplement rsID queries.
        progress_cb: Progress reporting callback.

    Returns:
        List of TCGA mutation dicts with cancer context.
    """
    init_tcga_schema()

    # Check cache first (query all, filter in Python to avoid 32k variable limit)
    if rsids:
        rsid_set = set(rsids)
        all_cached = local_db.query("SELECT * FROM tcga_mutation")
        cached = [r for r in all_cached if r["rsid"] in rsid_set]
        if cached:
            progress_cb(f"[TCGA] Returning {len(cached)} cached results.")
            return [dict(r) for r in cached]

    progress_cb(f"[TCGA] Querying GDC for {len(rsids)} SNPs (no bulk download)...")

    rows_to_insert = []
    returned = []

    for rsid in rsids:
        try:
            result = _query_single_rsid(rsid, progress_cb)
            rows_to_insert.extend(result)
            returned.extend(result)
            time.sleep(0.3)   # respect GDC rate limits
        except Exception as e:
            log.debug(f"[TCGA] {rsid}: {e}")
            continue

    # Also query by gene if provided and rsID queries return nothing
    if not returned and genes:
        for gene in genes[:10]:   # limit gene queries
            try:
                result = _query_by_gene(gene, progress_cb)
                rows_to_insert.extend(result)
                returned.extend(result)
                time.sleep(0.3)
            except Exception as e:
                log.debug(f"[TCGA] gene {gene}: {e}")
                continue

    if rows_to_insert:
        local_db.executemany(
            """INSERT OR IGNORE INTO tcga_mutation
               (rsid, gene_symbol, mutation_id, cancer_type, cancer_abbr,
                mutation_type, variant_class, amino_acid_change,
                case_count, frequency, cosmic_id, is_known_driver, plain_english)
               VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)""",
            [_row_to_tuple(r) for r in rows_to_insert],
        )
        progress_cb(f"[TCGA] Cached {len(rows_to_insert)} results locally.")

    return returned


def _query_single_rsid(rsid: str, progress_cb) -> list[dict]:
    """Query GDC ssm (simple somatic mutations) for one rsID."""
    payload = {
        "filters": {
            "op": "=",
            "content": {"field": "ssm.consequence.transcript.annotation.dbsnp_rs", "value": rsid},
        },
        "fields": (
            "ssm_id,gene.symbol,consequence.transcript.consequence_type,"
            "consequence.transcript.annotation.amino_acids,"
            "consequence.transcript.annotation.dbsnp_rs,"
            "mutation_subtype,cosmic_id,"
            "occurrence.case.project.project_id,"
            "occurrence.case.project.name"
        ),
        "size": 50,
    }
    try:
        resp = _gdc_post("ssms", payload)
    except Exception:
        return []

    hits = resp.get("data", {}).get("hits", [])
    return _parse_ssm_hits(hits, rsid_hint=rsid)


def _query_by_gene(gene: str, progress_cb) -> list[dict]:
    """Query GDC for top mutations in a gene across TCGA projects."""
    payload = {
        "filters": {
            "op": "and",
            "content": [
                {"op": "=", "content": {"field": "ssm.consequence.transcript.gene.symbol", "value": gene}},
                {"op": "in",  "content": {"field": "ssm.consequence.transcript.consequence_type",
                                           "value": ["missense_variant", "stop_gained", "frameshift_variant"]}},
            ],
        },
        "fields": (
            "ssm_id,gene.symbol,consequence.transcript.consequence_type,"
            "consequence.transcript.annotation.amino_acids,"
            "consequence.transcript.annotation.dbsnp_rs,"
            "mutation_subtype,cosmic_id,"
            "occurrence.case.project.project_id,"
            "occurrence.case.project.name"
        ),
        "size": 20,
    }
    try:
        resp = _gdc_post("ssms", payload)
    except Exception:
        return []

    hits = resp.get("data", {}).get("hits", [])
    return _parse_ssm_hits(hits)


def _parse_ssm_hits(hits: list, rsid_hint: str = None) -> list[dict]:
    """Convert GDC SSM hits into our internal format."""
    results = []
    for hit in hits:
        gene_list = hit.get("gene", [])
        gene = gene_list[0].get("symbol", "") if gene_list else ""

        conseqs = hit.get("consequence", [])
        conseq_type = ""
        aa_change   = ""
        rsid        = rsid_hint or ""
        for c in conseqs:
            t = c.get("transcript", {})
            if t.get("consequence_type"):
                conseq_type = t["consequence_type"]
            ann = t.get("annotation", {})
            if ann.get("amino_acids"):
                aa_change = ann["amino_acids"]
            if ann.get("dbsnp_rs") and not rsid_hint:
                rsid = ann["dbsnp_rs"]

        # Count occurrences per cancer project
        occurrences = hit.get("occurrence", [])
        project_counts: dict[str, int] = {}
        for occ in occurrences:
            proj = occ.get("case", {}).get("project", {})
            proj_id = proj.get("project_id", "")
            project_counts[proj_id] = project_counts.get(proj_id, 0) + 1

        mutation_id  = hit.get("ssm_id", "")
        cosmic_id    = hit.get("cosmic_id", "")
        mut_subtype  = hit.get("mutation_subtype", "SNP")
        is_driver    = int(gene in KNOWN_DRIVERS)

        for proj_id, case_count in project_counts.items():
            cancer_name = CANCER_NAMES.get(proj_id, proj_id)
            plain = _make_plain_english(gene, conseq_type, aa_change, cancer_name, case_count, is_driver)
            results.append({
                "rsid":               rsid,
                "gene_symbol":        gene,
                "mutation_id":        mutation_id,
                "cancer_type":        cancer_name,
                "cancer_abbr":        proj_id.replace("TCGA-", ""),
                "mutation_type":      conseq_type,
                "variant_class":      mut_subtype,
                "amino_acid_change":  aa_change,
                "case_count":         case_count,
                "frequency":          0.0,
                "cosmic_id":          cosmic_id,
                "is_known_driver":    is_driver,
                "plain_english":      plain,
            })

    return results


def _make_plain_english(gene, conseq_type, aa_change, cancer_name, case_count, is_driver) -> str:
    """Generate a consumer-friendly TCGA finding explanation."""
    conseq_plain = {
        "missense_variant":   "an amino acid change",
        "stop_gained":        "a premature stop in the protein",
        "frameshift_variant": "a reading frame shift in the protein",
        "splice_region_variant": "a change near a splice site",
    }.get(conseq_type, "a genetic change")

    driver_note = " This is a well-known cancer driver gene." if is_driver else ""

    aa_note = f" (protein change: {aa_change})" if aa_change else ""

    return (
        f"This variant — {conseq_plain}{aa_note} in {gene} — has been observed "
        f"in {case_count} {cancer_name} tumour sample(s) in TCGA research data.{driver_note} "
        f"Important: this is somatic (tumour) data, not a personal cancer risk prediction. "
        f"Your germline variant may differ from these tumour mutations."
    )


def _row_to_tuple(r: dict) -> tuple:
    return (
        r["rsid"], r["gene_symbol"], r["mutation_id"],
        r["cancer_type"], r["cancer_abbr"], r["mutation_type"],
        r["variant_class"], r["amino_acid_change"], r["case_count"],
        r["frequency"], r["cosmic_id"], r["is_known_driver"], r["plain_english"],
    )


def get_cached_tcga(rsids: list[str]) -> list[dict]:
    """Return any cached TCGA results for the given rsids."""
    init_tcga_schema()
    if not rsids:
        return []
    # tcga_mutation is populated on-demand — query all rows then filter
    # to avoid SQLite 32k variable limit with large ancestry files
    rows = local_db.query(
        "SELECT * FROM tcga_mutation ORDER BY case_count DESC"
    )
    if not rows:
        return []
    rsid_set = set(rsids)
    return [dict(r) for r in rows if r["rsid"] in rsid_set]
