-- HealthAgent local genomics database schema
-- Stores data downloaded from open-source databases:
--   GWAS Catalog (EBI)    : trait/disease associations
--   ClinVar (NCBI)        : clinical significance
--   PharmGKB              : drug/medication interactions
--   dbSNP (NCBI)          : variant functional annotations

PRAGMA journal_mode=WAL;
PRAGMA foreign_keys=ON;

-- ── SNP core table ────────────────────────────────────────────────
CREATE TABLE IF NOT EXISTS snp (
    rsid            TEXT PRIMARY KEY,
    chromosome      TEXT,
    position        INTEGER,
    gene            TEXT,
    gene_id         TEXT,
    alleles         TEXT,       -- e.g. "A/G"
    functional_class TEXT,      -- missense, synonymous, intron, etc.
    maf             REAL,       -- minor allele frequency (global)
    updated_at      TEXT DEFAULT (datetime('now'))
);

-- ── GWAS Catalog associations ─────────────────────────────────────
CREATE TABLE IF NOT EXISTS gwas_association (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    rsid            TEXT NOT NULL,
    trait           TEXT NOT NULL,   -- plain English trait name
    trait_efo       TEXT,            -- EFO ontology ID (e.g. EFO_0000270)
    trait_category  TEXT,            -- e.g. "Immune system disease"
    p_value         REAL,
    odds_ratio      REAL,
    beta            REAL,
    risk_allele     TEXT,
    study_title     TEXT,
    pubmed_id       TEXT,
    source          TEXT DEFAULT 'gwas_catalog',
    downloaded_at   TEXT DEFAULT (datetime('now')),
    FOREIGN KEY (rsid) REFERENCES snp(rsid)
);

-- ── ClinVar clinical significance ─────────────────────────────────
CREATE TABLE IF NOT EXISTS clinvar_variant (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    rsid            TEXT,
    variation_id    TEXT,
    clinical_sig    TEXT,   -- Pathogenic, Likely pathogenic, VUS, Benign, etc.
    condition       TEXT,   -- disease/condition name
    condition_id    TEXT,   -- MedGen/OMIM ID
    review_status   TEXT,   -- criteria provided, practice guideline, etc.
    last_evaluated  TEXT,
    gene            TEXT,
    molecular_consequence TEXT,
    downloaded_at   TEXT DEFAULT (datetime('now'))
);

-- ── PharmGKB drug interactions ────────────────────────────────────
CREATE TABLE IF NOT EXISTS drug_interaction (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    rsid            TEXT,
    gene            TEXT,
    drug_name       TEXT NOT NULL,
    drug_rxnorm     TEXT,           -- RxNorm drug ID
    phenotype       TEXT,           -- e.g. "Decreased metabolism"
    significance    TEXT,           -- Level 1A, 1B, 2A, 2B, 3, 4
    plain_english   TEXT,           -- consumer-friendly explanation
    category        TEXT,           -- dosage, toxicity, efficacy, metabolism
    source          TEXT DEFAULT 'pharmgkb',
    downloaded_at   TEXT DEFAULT (datetime('now'))
);

-- ── Wellness trait annotations (curated) ─────────────────────────
CREATE TABLE IF NOT EXISTS wellness_trait (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    rsid            TEXT NOT NULL,
    gene            TEXT,
    category        TEXT NOT NULL,  -- nutrition, fitness, sleep, allergies, etc.
    trait_name      TEXT NOT NULL,
    icon            TEXT,           -- emoji
    genotype        TEXT NOT NULL,  -- e.g. "TT"
    risk_level      TEXT,           -- typical, variant, elevated, reduced
    your_result     TEXT NOT NULL,  -- plain English one-liner
    detail          TEXT NOT NULL,  -- 2-3 sentence explanation
    why_it_matters  TEXT,           -- mechanism in plain English
    what_to_do      TEXT,           -- actionable advice
    source_study    TEXT,           -- PubMed ID or GWAS study
    UNIQUE(rsid, genotype)
);

-- ── Download log ─────────────────────────────────────────────────
CREATE TABLE IF NOT EXISTS download_log (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    source          TEXT NOT NULL,   -- gwas_catalog, clinvar, pharmgkb, etc.
    status          TEXT NOT NULL,   -- started, completed, failed
    records_added   INTEGER DEFAULT 0,
    error_msg       TEXT,
    started_at      TEXT DEFAULT (datetime('now')),
    finished_at     TEXT
);

-- ── Indexes ───────────────────────────────────────────────────────
CREATE INDEX IF NOT EXISTS idx_gwas_rsid      ON gwas_association(rsid);
CREATE INDEX IF NOT EXISTS idx_gwas_trait     ON gwas_association(trait);
CREATE INDEX IF NOT EXISTS idx_clinvar_rsid   ON clinvar_variant(rsid);
CREATE INDEX IF NOT EXISTS idx_drug_rsid      ON drug_interaction(rsid);
CREATE INDEX IF NOT EXISTS idx_drug_gene      ON drug_interaction(gene);
CREATE INDEX IF NOT EXISTS idx_wellness_rsid  ON wellness_trait(rsid);
CREATE INDEX IF NOT EXISTS idx_wellness_cat   ON wellness_trait(category);
