"""Open-source genomics database downloader.

Downloads and ingests data from:
  - GWAS Catalog (EBI)      : SNP → trait/disease associations
  - ClinVar (NCBI)          : SNP → clinical significance
  - PharmGKB                : gene/SNP → drug interactions
  - DisGeNET                : gene → disease associations (scored)
  - OpenTargets Platform    : gene → disease evidence (multi-source)
  - Ensembl REST API        : variant functional consequences (per-SNP)
  - FinnGen / UK Biobank    : large-scale biobank phenome associations

All data is stored in the local SQLite database (data/healthagent.db).
Each source can be downloaded independently or all at once.

Usage:
    python -m healthagent.databases.downloader --all
    python -m healthagent.databases.downloader --gwas
    python -m healthagent.databases.downloader --clinvar
    python -m healthagent.databases.downloader --pharmgkb
    python -m healthagent.databases.downloader --disgenet
    python -m healthagent.databases.downloader --opentargets
    python -m healthagent.databases.downloader --ensembl
    python -m healthagent.databases.downloader --finngen
    python -m healthagent.databases.downloader --status
"""

import argparse
import csv
import gzip
import io
import json
import logging
import sys
import time
import urllib.request
from pathlib import Path
from typing import Callable, Generator, Iterator

from healthagent.databases import local_db

log = logging.getLogger(__name__)

# ── Download cache dir ────────────────────────────────────────────
CACHE_DIR = Path(__file__).parent.parent.parent / "data" / "cache"

# ── Source URLs ───────────────────────────────────────────────────
GWAS_ASSOCIATIONS_URL = (
    "https://www.ebi.ac.uk/gwas/api/search/downloads/alternative"
)
CLINVAR_VARIANT_SUMMARY_URL = (
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
)
PHARMGKB_RELATIONSHIPS_URL = (
    "https://api.pharmgkb.org/v1/download/file/data/relationships.zip"
)
PHARMGKB_VARIANTS_URL = (
    "https://api.pharmgkb.org/v1/download/file/data/clinicalAnnotations.zip"
)

# DisGeNET — gene-disease associations (all sources, CC BY-NC-SA 4.0)
DISGENET_ALL_URL = (
    "https://www.disgenet.org/static/disgenet_ap1/files/downloads/all_gene_disease_associations.tsv.gz"
)

# OpenTargets — overall association scores (Apache 2.0, Parquet/JSON)
# We use the smaller JSON summary file rather than full Parquet
OPENTARGETS_ASSOC_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/latest/output/etl/json/associationByOverallDirect/"
)

# FinnGen summary statistics endpoint (R10 public release)
FINNGEN_MANIFEST_URL = (
    "https://r10.finngen.fi/api/phenos"
)

# Ensembl REST API base (per-variant queries, no bulk download needed)
ENSEMBL_REST_BASE = "https://rest.ensembl.org"


def _fetch(url: str, label: str, cache_name: str) -> Path:
    """Download url to cache, returning path. Re-uses cached file if <24h old."""
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    dest = CACHE_DIR / cache_name
    if dest.exists():
        age_h = (time.time() - dest.stat().st_mtime) / 3600
        if age_h < 24:
            log.info(f"[{label}] Using cached file ({age_h:.1f}h old): {dest.name}")
            return dest
    log.info(f"[{label}] Downloading from {url}")
    req = urllib.request.Request(url, headers={"User-Agent": "HealthAgent/1.0"})
    with urllib.request.urlopen(req, timeout=120) as resp:
        data = resp.read()
    dest.write_bytes(data)
    log.info(f"[{label}] Saved {len(data):,} bytes → {dest}")
    return dest


def _open_maybe_gz(path: Path):
    """Return a text file handle, decompressing .gz if needed."""
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


# ═══════════════════════════════════════════════════════════════════
#  GWAS CATALOG
# ═══════════════════════════════════════════════════════════════════

def download_gwas(progress_cb: Callable[[str], None] = print) -> int:
    """Download GWAS Catalog associations and store in local DB.

    What it is:
        The GWAS Catalog is a curated collection of all published genome-wide
        association studies. It links specific SNPs (rsIDs) to traits and
        diseases, along with p-values, odds ratios, and source publications.

    Why it matters for users:
        It answers "has this SNP been scientifically associated with any
        health conditions?" — giving users real research context rather than
        marketing copy.

    Returns:
        Number of association records inserted.
    """
    log_id = local_db.log_download("gwas_catalog", "started")
    progress_cb("[GWAS] Starting download from EBI GWAS Catalog...")

    try:
        path = _fetch(GWAS_ASSOCIATIONS_URL, "GWAS", "gwas_associations.tsv")
    except Exception as e:
        local_db.log_download("gwas_catalog", "failed", error=str(e))
        progress_cb(f"[GWAS] Download failed: {e}")
        return 0

    progress_cb("[GWAS] Parsing associations...")

    # Columns in GWAS alternative download TSV:
    # DATE ADDED TO CATALOG, PUBMEDID, FIRST AUTHOR, DATE, JOURNAL, LINK,
    # STUDY, DISEASE/TRAIT, INITIAL SAMPLE SIZE, REPLICATION SAMPLE SIZE,
    # REGION, CHR_ID, CHR_POS, REPORTED GENE(S), MAPPED_GENE, UPSTREAM_GENE_ID,
    # DOWNSTREAM_GENE_ID, SNP_GENE_IDS, UPSTREAM_GENE_DISTANCE,
    # DOWNSTREAM_GENE_DISTANCE, STRONGEST SNP-RISK ALLELE, SNPS, MERGED,
    # SNP_ID_CURRENT, CONTEXT, INTERGENIC, RISK ALLELE FREQUENCY, P-VALUE,
    # PVALUE_MLOG, P-VALUE (TEXT), OR or BETA, 95% CI (TEXT), PLATFORM,
    # CNV, MAPPED_TRAIT, MAPPED_TRAIT_URI, STUDY ACCESSION, GENOTYPING TECHNOLOGY

    rows_to_insert = []
    snps_to_insert = []
    count = 0

    with _open_maybe_gz(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rsid_raw = row.get("SNPS", "").strip()
            if not rsid_raw or not rsid_raw.startswith("rs"):
                continue

            trait = row.get("DISEASE/TRAIT", "").strip()
            if not trait:
                continue

            try:
                p_val = float(row.get("P-VALUE", "") or 0)
            except ValueError:
                p_val = None

            try:
                or_beta = float(row.get("OR or BETA", "") or 0)
            except ValueError:
                or_beta = None

            risk_allele_full = row.get("STRONGEST SNP-RISK ALLELE", "")
            risk_allele = risk_allele_full.split("-")[-1] if "-" in risk_allele_full else risk_allele_full

            mapped_trait = row.get("MAPPED_TRAIT", "").strip()
            trait_uri    = row.get("MAPPED_TRAIT_URI", "").strip()
            efo_id       = trait_uri.split("/")[-1] if "/" in trait_uri else trait_uri
            pubmed_id    = row.get("PUBMEDID", "").strip()
            study_title  = row.get("STUDY", "").strip()
            gene         = row.get("MAPPED_GENE", "").strip()
            chrom        = row.get("CHR_ID", "").strip()

            try:
                pos = int(row.get("CHR_POS", "") or 0)
            except ValueError:
                pos = 0

            snps_to_insert.append((rsid_raw, chrom, pos, gene))
            rows_to_insert.append((
                rsid_raw, mapped_trait or trait, efo_id,
                row.get("CONTEXT", "").strip(),
                p_val, or_beta, None, risk_allele,
                study_title, pubmed_id,
            ))
            count += 1

            if count % 10000 == 0:
                progress_cb(f"[GWAS] Parsed {count:,} records...")

    progress_cb(f"[GWAS] Inserting {len(snps_to_insert):,} SNP records...")
    local_db.executemany(
        "INSERT OR IGNORE INTO snp(rsid, chromosome, position, gene) VALUES (?,?,?,?)",
        snps_to_insert,
    )

    progress_cb(f"[GWAS] Inserting {len(rows_to_insert):,} associations...")
    local_db.executemany(
        """INSERT OR IGNORE INTO gwas_association
           (rsid, trait, trait_efo, trait_category, p_value, odds_ratio, beta,
            risk_allele, study_title, pubmed_id)
           VALUES (?,?,?,?,?,?,?,?,?,?)""",
        rows_to_insert,
    )

    local_db.log_download("gwas_catalog", "completed", records=count)
    progress_cb(f"[GWAS] Done — {count:,} associations stored.")
    return count


# ═══════════════════════════════════════════════════════════════════
#  CLINVAR
# ═══════════════════════════════════════════════════════════════════

def download_clinvar(progress_cb: Callable[[str], None] = print) -> int:
    """Download ClinVar variant summary and store clinical significance.

    What it is:
        ClinVar is NCBI's database of genomic variants and their clinical
        significance. Expert panels and labs submit interpretations of whether
        a variant is Pathogenic, Benign, or of Uncertain Significance (VUS).

    Why it matters for users:
        It tells users whether their specific DNA variant has been reviewed
        by medical geneticists and what health conditions it's associated with.
        This is the most authoritative source for clinical variant interpretation.

    Returns:
        Number of variants inserted.
    """
    local_db.log_download("clinvar", "started")
    progress_cb("[ClinVar] Starting download from NCBI FTP...")

    try:
        path = _fetch(CLINVAR_VARIANT_SUMMARY_URL, "ClinVar", "variant_summary.txt.gz")
    except Exception as e:
        local_db.log_download("clinvar", "failed", error=str(e))
        progress_cb(f"[ClinVar] Download failed: {e}")
        return 0

    progress_cb("[ClinVar] Parsing variant summary (this may take a minute)...")

    # variant_summary.txt columns (tab-delimited):
    # AlleleID, Type, Name, GeneID, GeneSymbol, HGNC_ID, ClinicalSignificance,
    # ClinSigSimple, LastEvaluated, RS# (dbSNP), nsv/esv (dbVar), RCVaccession,
    # PhenotypeIDS, PhenotypeList, Origin, OriginSimple, Assembly, ChromosomeAccession,
    # Chromosome, Start, Stop, ReferenceAllele, AlternateAllele, Cytogenetic,
    # ReviewStatus, NumberSubmitters, Guidelines, TestedInGTR, OtherIDs, SubmitterCategories,
    # VariationID, PositionVCF, ReferenceAlleleVCF, AlternateAlleleVCF

    rows_to_insert = []
    count = 0

    KEEP_SIG = {
        "Pathogenic", "Likely pathogenic",
        "Pathogenic/Likely pathogenic",
        "Uncertain significance",
        "Benign", "Likely benign",
        "Benign/Likely benign",
    }

    with _open_maybe_gz(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rs_raw = row.get("RS# (dbSNP)", "").strip()
            if not rs_raw or rs_raw in (".", "-1", ""):
                continue
            rsid = f"rs{rs_raw}"

            clinical_sig = row.get("ClinicalSignificance", "").strip()
            # Only store clinically meaningful classifications
            if not any(k in clinical_sig for k in KEEP_SIG):
                continue

            condition = row.get("PhenotypeList", "").strip()
            condition_id = row.get("PhenotypeIDS", "").strip()
            gene = row.get("GeneSymbol", "").strip()
            review_status = row.get("ReviewStatus", "").strip()
            last_eval = row.get("LastEvaluated", "").strip()
            variation_id = row.get("VariationID", "").strip()
            mol_consequence = row.get("Type", "").strip()

            rows_to_insert.append((
                rsid, variation_id, clinical_sig, condition, condition_id,
                review_status, last_eval, gene, mol_consequence,
            ))
            count += 1

            if count % 10000 == 0:
                progress_cb(f"[ClinVar] Parsed {count:,} records...")

    progress_cb(f"[ClinVar] Inserting {len(rows_to_insert):,} variants...")
    local_db.executemany(
        """INSERT OR IGNORE INTO clinvar_variant
           (rsid, variation_id, clinical_sig, condition, condition_id,
            review_status, last_evaluated, gene, molecular_consequence)
           VALUES (?,?,?,?,?,?,?,?,?)""",
        rows_to_insert,
    )

    local_db.log_download("clinvar", "completed", records=count)
    progress_cb(f"[ClinVar] Done — {count:,} variants stored.")
    return count


# ═══════════════════════════════════════════════════════════════════
#  PHARMGKB
# ═══════════════════════════════════════════════════════════════════

def download_pharmgkb(progress_cb: Callable[[str], None] = print) -> int:
    """Download PharmGKB clinical annotations and store drug interactions.

    What it is:
        PharmGKB curates how genetic variants affect drug response —
        metabolism, efficacy, dosage requirements, and adverse reactions.
        Evidence levels range from 1A (FDA-approved labeling) to 4 (case reports).

    Why it matters for users:
        Some medications work very differently depending on your DNA. This data
        tells users which drugs they may metabolize faster/slower, which may
        require dose adjustment, and which carry higher risk of side effects.

    Returns:
        Number of drug interaction records inserted.
    """
    local_db.log_download("pharmgkb", "started")
    progress_cb("[PharmGKB] Downloading clinical annotations...")

    try:
        path = _fetch(PHARMGKB_VARIANTS_URL, "PharmGKB", "clinical_annotations.zip")
    except Exception as e:
        local_db.log_download("pharmgkb", "failed", error=str(e))
        progress_cb(f"[PharmGKB] Download failed: {e}")
        return 0

    import zipfile

    rows_to_insert = []
    count = 0

    # PharmGKB clinical_annotations.tsv columns:
    # Clinical Annotation ID, Variant/Haplotypes, Gene, Level of Evidence,
    # Level Override, Level Modifiers, Score, Phenotype Category,
    # PMID Count, Evidence Count, Drug(s), Phenotype(s), Latest History Date,
    # URL, Specialty Population

    LEVEL_PLAIN = {
        "1A": "Very high confidence — included in FDA drug labeling",
        "1B": "High confidence — validated in multiple studies",
        "2A": "Moderate confidence — replicated association",
        "2B": "Moderate confidence — single study association",
        "3":  "Preliminary evidence — limited studies",
        "4":  "Early evidence — case reports only",
    }

    CATEGORY_PLAIN = {
        "Metabolism/PK": "How fast your body processes this medication",
        "Efficacy":      "How well this medication works for you",
        "Toxicity":      "Your risk of side effects from this medication",
        "Dosage":        "Whether you may need a different dose",
    }

    try:
        with zipfile.ZipFile(path) as zf:
            tsv_name = next(
                (n for n in zf.namelist() if "clinical_annotations" in n.lower() and n.endswith(".tsv")),
                None,
            )
            if not tsv_name:
                # Try any TSV in the zip
                tsv_name = next((n for n in zf.namelist() if n.endswith(".tsv")), None)

            if not tsv_name:
                progress_cb("[PharmGKB] No TSV found in zip. Skipping.")
                local_db.log_download("pharmgkb", "failed", error="No TSV in zip")
                return 0

            with zf.open(tsv_name) as raw:
                reader = csv.DictReader(
                    io.TextIOWrapper(raw, encoding="utf-8", errors="replace"),
                    delimiter="\t",
                )
                for row in reader:
                    gene = row.get("Gene", "").strip()
                    variant = row.get("Variant/Haplotypes", "").strip()
                    level = row.get("Level of Evidence", "").strip()
                    drugs_raw = row.get("Drug(s)", "").strip()
                    phenotypes = row.get("Phenotype(s)", "").strip()
                    category = row.get("Phenotype Category", "").strip()

                    if not drugs_raw or not gene:
                        continue

                    # Extract rsID from variant field if present
                    rsid = None
                    for part in variant.split(";"):
                        part = part.strip()
                        if part.startswith("rs"):
                            rsid = part.split()[0]
                            break

                    level_plain = LEVEL_PLAIN.get(level, f"Evidence level {level}")
                    cat_plain = CATEGORY_PLAIN.get(category, category)

                    plain_english = (
                        f"{cat_plain}. "
                        f"Confidence: {level_plain}. "
                        f"Effect: {phenotypes}"
                    )

                    for drug in drugs_raw.split(";"):
                        drug = drug.strip()
                        if not drug:
                            continue
                        rows_to_insert.append((
                            rsid, gene, drug, None, phenotypes,
                            level, plain_english, cat_plain,
                        ))
                        count += 1

    except zipfile.BadZipFile as e:
        local_db.log_download("pharmgkb", "failed", error=str(e))
        progress_cb(f"[PharmGKB] Bad zip file: {e}")
        return 0

    progress_cb(f"[PharmGKB] Inserting {len(rows_to_insert):,} drug interactions...")
    local_db.executemany(
        """INSERT OR IGNORE INTO drug_interaction
           (rsid, gene, drug_name, drug_rxnorm, phenotype,
            significance, plain_english, category)
           VALUES (?,?,?,?,?,?,?,?)""",
        rows_to_insert,
    )

    local_db.log_download("pharmgkb", "completed", records=count)
    progress_cb(f"[PharmGKB] Done — {count:,} drug interactions stored.")
    return count


# ═══════════════════════════════════════════════════════════════════
#  Wellness traits seed (curated data, always loaded)
# ═══════════════════════════════════════════════════════════════════

def seed_wellness_traits(progress_cb: Callable[[str], None] = print) -> int:
    """Seed the wellness_trait table with curated SNP→trait mappings.

    These are hand-curated from published GWAS studies and population genetics
    research. Each entry explains not just what the SNP is, but what it means
    in everyday life and what a person can do about it.
    """
    progress_cb("[Wellness] Seeding curated trait database...")

    TRAITS = [
        # ── Nutrition ──────────────────────────────────────────────
        {
            "rsid": "rs1801133", "gene": "MTHFR", "category": "nutrition",
            "trait_name": "Folate & B-Vitamin Processing", "icon": "🥬",
            "entries": {
                "TT": ("elevated",
                       "Your body converts folate less efficiently than average",
                       "The MTHFR C677T variant reduces the activity of an enzyme that converts folic acid into its usable form (5-MTHF). This affects how well your body uses B vitamins for energy and cell repair.",
                       "Look for supplements labelled 'methylfolate' or '5-MTHF' rather than folic acid. Leafy greens, lentils, and avocado are rich natural sources.",
                       "PMID:9371533"),
                "CT": ("variant",
                       "You have one copy of the MTHFR variant — mild effect",
                       "One copy of this variant gives you slightly reduced folate conversion. Your body works a little harder to process B vitamins.",
                       "A diet rich in leafy greens and legumes helps. Methylated B-vitamins in supplements are also easily absorbed.",
                       "PMID:9371533"),
                "CC": ("typical",
                       "Your folate and B-vitamin processing is typical",
                       "You have the most common MTHFR genotype. Your body converts folate efficiently.",
                       "A balanced diet covering B-vitamin rich foods (leafy greens, eggs, legumes) supports good health.",
                       "PMID:9371533"),
            },
        },
        {
            "rsid": "rs4988235", "gene": "MCM6", "category": "nutrition",
            "trait_name": "Dairy & Lactose Tolerance", "icon": "🥛",
            "entries": {
                "TT": ("typical",
                       "You tend to digest dairy products comfortably",
                       "The rs4988235 T allele (near the LCT gene) is associated with continued lactase production into adulthood. This means your small intestine can break down lactose in milk.",
                       "You can enjoy dairy without issues for most people with this profile.",
                       "PMID:16251998"),
                "CT": ("variant",
                       "You may find large amounts of dairy cause some discomfort",
                       "One copy of the T allele means you retain some lactase activity, but it may diminish with age. Small portions of dairy are usually fine.",
                       "Fermented dairy (yogurt, aged cheese) is easier to digest as bacteria pre-digest much of the lactose.",
                       "PMID:16251998"),
                "CC": ("elevated",
                       "You are likely lactose intolerant — dairy may cause discomfort",
                       "The CC genotype is associated with reduced lactase production after childhood — the ancestral human pattern. Without lactase, lactose passes undigested into the colon where bacteria ferment it, causing bloating and discomfort.",
                       "Lactase enzyme supplements taken before meals help. Aged cheeses and yogurt are naturally lower in lactose. Plant-based milks are a comfortable alternative.",
                       "PMID:16251998"),
            },
        },
        {
            "rsid": "rs762551", "gene": "CYP1A2", "category": "nutrition",
            "trait_name": "Caffeine Metabolism Speed", "icon": "☕",
            "entries": {
                "AA": ("typical",
                       "You metabolise caffeine quickly — a 'fast metaboliser'",
                       "The CYP1A2 enzyme breaks down caffeine in your liver. The AA genotype produces a highly active form, clearing caffeine from your bloodstream about 4x faster than slow metabolisers.",
                       "Coffee effects wear off sooner for you. You're less likely to have sleep disruption from afternoon coffee, but may need more caffeine to feel the same effect.",
                       "PMID:16174292"),
                "AC": ("variant",
                       "You have mixed caffeine metabolism speed",
                       "One copy of the slow-metabolism allele means caffeine stays in your system longer than for fast metabolisers.",
                       "You may feel effects for 4-6 hours. An afternoon cutoff for caffeine (around 2pm) is a useful guideline.",
                       "PMID:16174292"),
                "CC": ("elevated",
                       "You metabolise caffeine slowly — it stays in your system longer",
                       "Slow CYP1A2 activity means caffeine lingers in your bloodstream much longer. Studies suggest slow metabolisers who drink 4+ coffees per day have a higher risk of cardiovascular effects.",
                       "Consider limiting to 1-2 cups per day and avoiding caffeine after midday. Green tea has a gentler caffeine profile that suits slow metabolisers better.",
                       "PMID:16174292"),
            },
        },
        {
            "rsid": "rs9939609", "gene": "FTO", "category": "nutrition",
            "trait_name": "Appetite & Hunger Signals", "icon": "🍽️",
            "entries": {
                "AA": ("elevated",
                       "Your genes may make you feel hungrier than average",
                       "The FTO gene influences levels of ghrelin (the hunger hormone) and affects brain circuits involved in food reward. The AA genotype is associated with stronger appetite signals and a preference for higher-calorie foods.",
                       "Mindful eating, high-protein meals, and regular mealtimes are especially effective strategies. Your appetite signals are louder — not a character flaw, just biology.",
                       "PMID:17434869"),
                "AT": ("variant",
                       "You have a moderate appetite signal variant",
                       "One copy of the FTO risk allele gives a modest increase in appetite signalling compared to the TT genotype.",
                       "Balanced meals with adequate protein and fibre help maintain satiety.",
                       "PMID:17434869"),
                "TT": ("typical",
                       "Your hunger and appetite signals are typical",
                       "The TT genotype is associated with typical appetite regulation and food reward signalling.",
                       "A balanced, varied diet supports good energy balance.",
                       "PMID:17434869"),
            },
        },
        # ── Fitness ────────────────────────────────────────────────
        {
            "rsid": "rs1815739", "gene": "ACTN3", "category": "fitness",
            "trait_name": "Muscle Fibre Type", "icon": "💪",
            "entries": {
                "CC": ("typical",
                       "Your muscles lean toward power and sprint performance",
                       "ACTN3 produces alpha-actinin-3, a protein found only in fast-twitch (Type II) muscle fibres. The CC genotype produces functional ACTN3 protein, giving these fibres better structure for explosive movements.",
                       "You may naturally excel at sports requiring speed and power. Strength training, HIIT, and sprint intervals tend to work very well for your muscle type.",
                       "PMID:12879365"),
                "CT": ("typical",
                       "You have a balanced mix of power and endurance muscle traits",
                       "One copy of each variant gives you functional ACTN3 in fast-twitch fibres alongside some characteristics that favour endurance. A versatile muscle profile.",
                       "Both strength training and aerobic exercise work well. You can adapt to a wide range of sports.",
                       "PMID:12879365"),
                "TT": ("variant",
                       "Your muscles lean toward endurance performance",
                       "The TT genotype produces no functional ACTN3 protein (about 18% of people worldwide). Rather than a deficiency, this seems to shift muscle metabolism toward more efficient oxygen use — an endurance advantage.",
                       "Endurance sports (running, cycling, swimming) may come more naturally. Long, steady-state cardio tends to suit you well. Olympic endurance athletes are disproportionately TT.",
                       "PMID:12879365"),
            },
        },
        {
            "rsid": "rs8192678", "gene": "PPARGC1A", "category": "fitness",
            "trait_name": "Aerobic Training Response", "icon": "🏃",
            "entries": {
                "GG": ("typical",
                       "Your aerobic fitness tends to improve quickly with training",
                       "PPARGC1A (PGC-1α) is the master regulator of mitochondrial biogenesis — it triggers your body to build more energy-producing mitochondria in response to exercise.",
                       "You're likely a good responder to aerobic training. Consistency pays off especially well for you.",
                       "PMID:12837945"),
                "GA": ("typical",
                       "Your aerobic training response is typical",
                       "One copy of each allele gives moderate PGC-1α activity.",
                       "Consistent aerobic exercise builds fitness well for your genotype.",
                       "PMID:12837945"),
                "AA": ("variant",
                       "Your aerobic fitness may build more gradually with training",
                       "The Ala variant of PGC-1α is associated with somewhat lower mitochondrial response to aerobic exercise. This doesn't mean you can't improve — it means consistency over a longer period matters more.",
                       "Longer warm-up periods, progressive overload, and adequate recovery time help maximise your aerobic gains.",
                       "PMID:12837945"),
            },
        },
        # ── Sleep ──────────────────────────────────────────────────
        {
            "rsid": "rs12736689", "gene": "CLOCK", "category": "sleep",
            "trait_name": "Morning Person vs Night Owl", "icon": "🌙",
            "entries": {
                "CC": ("typical",
                       "Your genes lean toward being a morning person",
                       "The CLOCK gene is the core driver of your circadian rhythm — your body's 24-hour internal clock. The CC genotype is associated with an earlier-shifting biological clock.",
                       "Your natural wake time is earlier. Working with this (early bedtime, morning exercise) tends to improve your energy and mood.",
                       "PMID:16002782"),
                "CT": ("typical",
                       "Your natural sleep timing is intermediate",
                       "Neither strongly morning nor evening-oriented. Your circadian rhythm is flexible.",
                       "You adapt reasonably well to different schedules, though consistency still helps sleep quality.",
                       "PMID:16002782"),
                "TT": ("variant",
                       "Your genes lean toward being a night owl",
                       "The TT genotype is associated with a later-shifting circadian clock. Your body naturally wants to sleep and wake later than average.",
                       "Fighting your chronotype causes chronic social jet-lag. Where possible, align your schedule with your natural rhythm. Evening exercise is fine for night owls — unlike morning types, it doesn't disrupt your sleep.",
                       "PMID:16002782"),
            },
        },
        {
            "rsid": "rs4680", "gene": "COMT", "category": "sleep",
            "trait_name": "Stress Response & Resilience", "icon": "🧘",
            "entries": {
                "GG": ("typical",
                       "You tend to handle stress calmly — 'Warrior' profile",
                       "COMT (catechol-O-methyltransferase) clears dopamine and adrenaline from the brain's prefrontal cortex. The GG (Val/Val) variant clears these faster, keeping the stress response shorter.",
                       "You stay calm under pressure but may feel less emotional reward from everyday pleasures. You thrive in high-pressure, competitive situations.",
                       "PMID:12161658"),
                "AG": ("typical",
                       "You have a balanced stress and reward profile",
                       "The intermediate COMT activity gives you moderate stress resilience alongside moderate emotional reward sensitivity.",
                       "You adapt well to varied situations. Mindfulness and regular exercise enhance both stress management and wellbeing.",
                       "PMID:12161658"),
                "AA": ("variant",
                       "You feel stress more acutely — 'Worrier' profile, but also more reward",
                       "Slower COMT activity means dopamine and adrenaline linger longer in the prefrontal cortex. This makes you more sensitive to stress but also more sensitive to pleasure and rewards — and often more creatively driven.",
                       "Stress management techniques (breathing, exercise, sleep) have an outsized positive effect for you. Avoid unnecessary stress exposure where you can. Many highly creative people share this profile.",
                       "PMID:12161658"),
            },
        },
        # ── Allergies & Immune ─────────────────────────────────────
        {
            "rsid": "rs1800629", "gene": "TNF", "category": "allergies",
            "trait_name": "Inflammatory Response Strength", "icon": "🤧",
            "entries": {
                "AA": ("elevated",
                       "Your immune system may respond more strongly to triggers",
                       "TNF-alpha is a key inflammatory signalling protein. The AA genotype is associated with higher TNF-alpha production, meaning your immune system can mount a stronger inflammatory response.",
                       "You may experience more pronounced allergy symptoms, longer recovery from infections, or more sensitivity to environmental triggers. Anti-inflammatory foods (oily fish, turmeric, berries) and avoiding smoke/pollution is especially worthwhile.",
                       "PMID:12490404"),
                "AG": ("variant",
                       "Your inflammatory response is moderately elevated",
                       "One copy of the TNF-alpha promoter variant gives moderately higher inflammatory signalling.",
                       "A diet rich in anti-inflammatory foods and regular moderate exercise supports a balanced immune response.",
                       "PMID:12490404"),
                "GG": ("typical",
                       "Your inflammatory response is typical",
                       "The GG genotype produces typical TNF-alpha levels, associated with a standard immune and inflammatory response.",
                       "A balanced diet and regular exercise maintain good immune health.",
                       "PMID:12490404"),
            },
        },
        {
            "rsid": "rs2243250", "gene": "IL4", "category": "allergies",
            "trait_name": "Hay Fever & Environmental Allergy Risk", "icon": "🌿",
            "entries": {
                "TT": ("elevated",
                       "You may be more prone to hay fever and environmental allergies",
                       "IL-4 drives the Th2 immune pathway responsible for allergic responses. Higher IL-4 production (associated with the T allele) shifts immunity toward allergy-prone responses, including higher IgE antibody production.",
                       "Air purifiers, regular nasal rinses, and tracking pollen counts help manage symptoms. Antihistamines tend to work well for this mechanism.",
                       "PMID:10790169"),
                "TC": ("variant",
                       "You have moderate environmental allergy sensitivity",
                       "One copy of the IL-4 variant gives moderate allergy tendency.",
                       "Monitoring seasonal triggers and reducing indoor allergens (dust, pet dander) is worthwhile.",
                       "PMID:10790169"),
                "CC": ("typical",
                       "Your environmental allergy sensitivity is typical",
                       "The CC genotype is associated with typical IL-4 signalling and standard allergy risk.",
                       "Standard sensible measures — keeping windows closed during high pollen, regular cleaning — are sufficient.",
                       "PMID:10790169"),
            },
        },
        # ── Supplements ────────────────────────────────────────────
        {
            "rsid": "rs2282679", "gene": "GC", "category": "supplements",
            "trait_name": "Vitamin D Transport Efficiency", "icon": "☀️",
            "entries": {
                "AA": ("elevated",
                       "Your body may struggle to transport Vitamin D efficiently",
                       "The GC gene produces Vitamin D Binding Protein (VDBP), which carries Vitamin D through the bloodstream. The AA genotype is associated with lower VDBP levels, meaning Vitamin D gets less efficiently delivered to tissues.",
                       "Vitamin D supplementation (1000-2000 IU daily) is widely recommended for this profile. Getting your blood Vitamin D level checked (25-OH-D test) gives a personal baseline.",
                       "PMID:20541252"),
                "AC": ("variant",
                       "Your Vitamin D transport is moderately reduced",
                       "One copy of the GC variant gives moderate reduction in VDBP levels.",
                       "Regular sunlight exposure and Vitamin D-rich foods (oily fish, eggs) support good levels.",
                       "PMID:20541252"),
                "CC": ("typical",
                       "Your Vitamin D transport is efficient",
                       "The CC genotype produces typical VDBP levels, efficiently transporting Vitamin D to where your body needs it.",
                       "Standard sun exposure and a balanced diet maintain good Vitamin D levels for most people with this profile.",
                       "PMID:20541252"),
            },
        },
        {
            "rsid": "rs1544410", "gene": "VDR", "category": "supplements",
            "trait_name": "Vitamin D Receptor Sensitivity", "icon": "🦴",
            "entries": {
                "CC": ("elevated",
                       "Your Vitamin D receptors may be less responsive",
                       "The Vitamin D Receptor (VDR) is how your cells actually respond to Vitamin D. The BsmI CC variant is associated with reduced receptor activity, meaning you may need higher circulating Vitamin D levels to get the same effect.",
                       "Combined with rs2282679, this doubles the case for monitoring your Vitamin D levels. Higher supplementation (2000-4000 IU, with medical guidance) may be appropriate.",
                       "PMID:17608567"),
                "CT": ("variant",
                       "Your Vitamin D receptor sensitivity is moderately reduced",
                       "One copy of this VDR variant gives moderate reduction in receptor response.",
                       "Adequate Vitamin D from sunlight, diet, and supplements covers most needs.",
                       "PMID:17608567"),
                "TT": ("typical",
                       "Your Vitamin D receptor sensitivity is typical",
                       "The TT genotype is associated with typical VDR responsiveness.",
                       "Standard sun exposure and diet maintain healthy Vitamin D activity.",
                       "PMID:17608567"),
            },
        },
        # ── Health history / awareness ─────────────────────────────
        {
            "rsid": "rs7903146", "gene": "TCF7L2", "category": "health_history",
            "trait_name": "Blood Sugar Regulation", "icon": "🩸",
            "entries": {
                "TT": ("elevated",
                       "Worth keeping an eye on your blood sugar balance",
                       "TCF7L2 is the strongest known genetic marker for blood sugar processing. The TT genotype is associated with reduced insulin secretion from the pancreas in response to carbohydrates. This is the most replicated type-2 diabetes genetic marker in research.",
                       "Reducing refined carbohydrates and sugar-sweetened drinks, increasing fibre intake, and regular walking after meals significantly reduce the effect of this variant. This is a lifestyle-responsive finding.",
                       "PMID:17463246"),
                "TC": ("variant",
                       "A mild blood sugar awareness marker",
                       "One copy of the TCF7L2 risk allele gives a moderate association with blood sugar processing.",
                       "A balanced diet moderate in refined carbohydrates and regular physical activity covers this well.",
                       "PMID:17463246"),
                "CC": ("typical",
                       "Your blood sugar regulation genetics are typical",
                       "The CC genotype is not associated with impaired insulin secretion.",
                       "A healthy diet and regular exercise maintain good blood sugar balance.",
                       "PMID:17463246"),
            },
        },
        {
            "rsid": "rs1333049", "gene": "CDKN2B-AS1", "category": "health_history",
            "trait_name": "Heart Health Awareness", "icon": "❤️",
            "entries": {
                "GG": ("elevated",
                       "A heart health awareness marker worth knowing about",
                       "This locus on chromosome 9p21 is one of the most replicated genetic markers for coronary artery disease in large population studies. The mechanism involves regulation of nearby genes that control arterial cell growth and senescence.",
                       "The good news: lifestyle factors override this genetic signal strongly. Not smoking, regular aerobic exercise, and a Mediterranean-style diet reduce risk substantially. Knowing this profile makes regular blood pressure and cholesterol checks extra worthwhile.",
                       "PMID:17478681"),
                "GC": ("variant",
                       "A mild heart health awareness marker",
                       "One copy of the 9p21 risk variant gives a moderate association in population studies.",
                       "Heart-healthy habits (exercise, diet, not smoking) are beneficial for everyone and especially offset this variant.",
                       "PMID:17478681"),
                "CC": ("typical",
                       "Your heart health genetics are typical for this marker",
                       "The CC genotype is not associated with elevated risk at this locus.",
                       "Heart-healthy habits benefit everyone regardless of genetics.",
                       "PMID:17478681"),
            },
        },
        # ── Personal traits ────────────────────────────────────────
        {
            "rsid": "rs1800497", "gene": "ANKK1", "category": "traits",
            "trait_name": "Motivation & Reward Drive", "icon": "🎯",
            "entries": {
                "GG": ("variant",
                       "You may be more driven to seek novelty and stimulation",
                       "This variant near the DRD2 dopamine receptor gene is associated with reduced D2 receptor density in the brain's reward system. Fewer receptors means the reward signal is less intense — motivating a drive to seek more stimulation, challenge, and novelty.",
                       "You may thrive in dynamic, varied environments. Boredom can be a real challenge. Channel this into challenging projects, social connection, and physical activity.",
                       "PMID:9795319"),
                "GA": ("typical",
                       "Your reward and motivation drive is moderate",
                       "One copy of each allele gives a balance between reward sensitivity and novelty-seeking.",
                       "A mix of routine and novelty works well for your motivation profile.",
                       "PMID:9795319"),
                "AA": ("typical",
                       "Your brain's reward system is typical in sensitivity",
                       "The AA genotype is associated with typical D2 receptor density. Rewards feel rewarding at normal intensity.",
                       "Standard goals and routines with regular variety maintain motivation well.",
                       "PMID:9795319"),
            },
        },
        {
            "rsid": "rs6265", "gene": "BDNF", "category": "traits",
            "trait_name": "Learning Style & Memory", "icon": "🧠",
            "entries": {
                "CC": ("typical",
                       "You tend to learn well from repeated practice",
                       "BDNF (Brain-Derived Neurotrophic Factor) supports the growth and maintenance of neurons. The Val66Val (CC) genotype is associated with higher activity-dependent BDNF release — meaning the more you practise a skill, the stronger the neural connections become.",
                       "Spaced repetition and deliberate practice are especially effective for you. Your memory for procedural and skill-based learning is a strength.",
                       "PMID:12546884"),
                "CT": ("variant",
                       "You have a mixed learning and memory profile",
                       "One Met allele reduces activity-dependent BDNF release slightly, but is also associated with stronger emotional memory consolidation.",
                       "You may remember emotionally significant events particularly vividly. Varied learning approaches (visual, verbal, hands-on) work well.",
                       "PMID:12546884"),
                "TT": ("variant",
                       "Your memory may be especially shaped by emotional context",
                       "The Met/Met (TT) genotype is associated with lower activity-dependent BDNF secretion, which reduces some forms of hippocampal memory — but strengthens amygdala-mediated emotional memory. Experiences with emotional weight are remembered very vividly.",
                       "Creating emotional hooks (stories, personal connections) improves learning and recall for you. Anxiety management is particularly beneficial since stress affects memory more for this profile.",
                       "PMID:12546884"),
            },
        },
    ]

    rows: list[tuple] = []
    for trait in TRAITS:
        for genotype, (risk_level, your_result, detail, what_to_do, source) in trait["entries"].items():
            rows.append((
                trait["rsid"], trait["gene"], trait["category"],
                trait["trait_name"], trait["icon"], genotype,
                risk_level, your_result, detail,
                f"This variant affects {trait['gene']} function.",
                what_to_do, source,
            ))

    local_db.executemany(
        """INSERT OR IGNORE INTO wellness_trait
           (rsid, gene, category, trait_name, icon, genotype,
            risk_level, your_result, detail, why_it_matters, what_to_do, source_study)
           VALUES (?,?,?,?,?,?,?,?,?,?,?,?)""",
        rows,
    )
    progress_cb(f"[Wellness] Seeded {len(rows)} trait entries.")
    return len(rows)


# ═══════════════════════════════════════════════════════════════════
#  CLI entry point
# ═══════════════════════════════════════════════════════════════════

# ═══════════════════════════════════════════════════════════════════
#  DISGENET
# ═══════════════════════════════════════════════════════════════════

def download_disgenet(progress_cb: Callable[[str], None] = print) -> int:
    """Download DisGeNET gene-disease associations.

    What it is:
        DisGeNET aggregates gene-disease associations from curated databases
        (OMIM, Orphanet, ClinVar, UniProt) and text-mined literature, scoring
        each association by evidence quality (0–1 score).

    Why it matters for users:
        It enriches trait cards with disease context — not just which gene is
        involved, but which diseases that gene has been scientifically linked to,
        giving users a broader picture of their genetic health landscape.

    Returns:
        Number of association records inserted.
    """
    local_db.log_download("disgenet", "started")
    progress_cb("[DisGeNET] Downloading gene-disease associations...")

    try:
        path = _fetch(DISGENET_ALL_URL, "DisGeNET", "disgenet_associations.tsv.gz")
    except Exception as e:
        local_db.log_download("disgenet", "failed", error=str(e))
        progress_cb(f"[DisGeNET] Download failed: {e}")
        return 0

    progress_cb("[DisGeNET] Parsing associations...")

    # Columns: geneId, geneSymbol, DSI, DPI, diseaseId, diseaseName,
    #          diseaseType, diseaseClass, diseaseSemanticType, score,
    #          EI, YearInitial, YearFinal, NofPmids, NofSnps, source

    rows_to_insert = []
    count = 0

    with _open_maybe_gz(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gene = row.get("geneSymbol", "").strip()
            disease = row.get("diseaseName", "").strip()
            if not gene or not disease:
                continue

            try:
                score = float(row.get("score", 0) or 0)
            except ValueError:
                score = 0.0

            # Only keep associations with meaningful evidence
            if score < 0.1:
                continue

            try:
                ei = float(row.get("EI", 0) or 0)
            except ValueError:
                ei = 0.0

            rows_to_insert.append((
                gene,
                row.get("geneId", "").strip(),
                disease,
                row.get("diseaseId", "").strip(),
                row.get("diseaseClass", "").strip(),
                score,
                ei,
            ))
            count += 1

            if count % 20000 == 0:
                progress_cb(f"[DisGeNET] Parsed {count:,} records...")

    progress_cb(f"[DisGeNET] Inserting {len(rows_to_insert):,} associations...")
    local_db.executemany(
        """INSERT OR IGNORE INTO disgenet_association
           (gene_symbol, gene_id, disease_name, disease_id, disease_class, score, ei)
           VALUES (?,?,?,?,?,?,?)""",
        rows_to_insert,
    )

    local_db.log_download("disgenet", "completed", records=count)
    progress_cb(f"[DisGeNET] Done — {count:,} gene-disease associations stored.")
    return count


# ═══════════════════════════════════════════════════════════════════
#  ENSEMBL REST (per-SNP variant consequences)
# ═══════════════════════════════════════════════════════════════════

def download_ensembl_consequences(
    rsids: list[str] = None,
    progress_cb: Callable[[str], None] = print,
) -> int:
    """Query Ensembl REST API for functional consequences of tracked SNPs.

    What it is:
        The Ensembl Variant Effect Predictor (VEP) annotates each variant with
        its predicted molecular impact: missense (amino acid change), synonymous,
        splice region, regulatory region, etc. SIFT and PolyPhen scores predict
        whether a missense change is damaging.

    Why it matters for users:
        It explains *how* a variant affects a gene — "this changes an amino acid
        in the protein" is more meaningful than "this is a missense variant".

    Args:
        rsids: List of rsIDs to annotate. Uses tracked wellness SNPs if None.

    Returns:
        Number of consequence records inserted.
    """
    # Default to the wellness SNPs we track
    if rsids is None:
        rows = local_db.query("SELECT DISTINCT rsid FROM wellness_trait")
        rsids = [r["rsid"] for r in rows]

    if not rsids:
        progress_cb("[Ensembl] No SNPs to annotate.")
        return 0

    local_db.log_download("ensembl", "started")
    progress_cb(f"[Ensembl] Fetching VEP consequences for {len(rsids)} SNPs...")

    rows_to_insert = []
    # Ensembl allows batches of up to 200 rsIDs per POST
    batch_size = 50

    for i in range(0, len(rsids), batch_size):
        batch = rsids[i: i + batch_size]
        url   = f"{ENSEMBL_REST_BASE}/vep/human/id"
        payload = json.dumps({"ids": batch}).encode()
        req = urllib.request.Request(
            url,
            data=payload,
            headers={
                "Content-Type": "application/json",
                "Accept":       "application/json",
                "User-Agent":   "HealthAgent/1.0",
            },
            method="POST",
        )
        try:
            with urllib.request.urlopen(req, timeout=30) as resp:
                results = json.loads(resp.read())
        except Exception as e:
            progress_cb(f"[Ensembl] Batch {i//batch_size + 1} failed: {e}")
            continue

        for variant in results:
            rsid = variant.get("id", "")
            for tc in variant.get("transcript_consequences", []):
                conseqs = ",".join(tc.get("consequence_terms", []))
                rows_to_insert.append((
                    rsid,
                    tc.get("gene_symbol", ""),
                    tc.get("transcript_id", ""),
                    conseqs,
                    tc.get("impact", ""),
                    tc.get("biotype", ""),
                    tc.get("sift_prediction", ""),
                    tc.get("polyphen_prediction", ""),
                ))

        time.sleep(0.34)   # Ensembl rate limit: ~3 req/s

    progress_cb(f"[Ensembl] Inserting {len(rows_to_insert):,} consequence records...")
    if rows_to_insert:
        local_db.executemany(
            """INSERT OR IGNORE INTO ensembl_consequence
               (rsid, gene_symbol, transcript_id, consequence, impact,
                biotype, sift_prediction, polyphen_pred)
               VALUES (?,?,?,?,?,?,?,?)""",
            rows_to_insert,
        )

    local_db.log_download("ensembl", "completed", records=len(rows_to_insert))
    progress_cb(f"[Ensembl] Done — {len(rows_to_insert):,} consequences stored.")
    return len(rows_to_insert)


# ═══════════════════════════════════════════════════════════════════
#  FINNGEN (public R10 phenome summary statistics)
# ═══════════════════════════════════════════════════════════════════

def download_finngen(progress_cb: Callable[[str], None] = print) -> int:
    """Download FinnGen phenotype manifest and top SNP associations.

    What it is:
        FinnGen is a Finnish biobank study of ~500,000 participants linking
        genomic data to health records. Their public R10 release includes
        summary statistics for ~2,600 disease endpoints — one of the most
        powerful biobank datasets in the world.

    Why it matters for users:
        FinnGen provides some of the most statistically robust associations
        between SNPs and disease endpoints, particularly for autoimmune,
        cardiovascular, and metabolic conditions common in Northern European
        ancestry populations.

    Returns:
        Number of phenotype-SNP associations inserted.
    """
    local_db.log_download("finngen", "started")
    progress_cb("[FinnGen] Fetching phenotype manifest from R10 release...")

    # Fetch the phenotype manifest (JSON list of all studied endpoints)
    try:
        req = urllib.request.Request(
            FINNGEN_MANIFEST_URL,
            headers={"User-Agent": "HealthAgent/1.0", "Accept": "application/json"},
        )
        with urllib.request.urlopen(req, timeout=30) as resp:
            phenos = json.loads(resp.read())
    except Exception as e:
        local_db.log_download("finngen", "failed", error=str(e))
        progress_cb(f"[FinnGen] Manifest fetch failed: {e}")
        return 0

    progress_cb(f"[FinnGen] Found {len(phenos)} phenotype endpoints.")

    # For tracked wellness SNPs, query the FinnGen association API
    snp_rows = local_db.query("SELECT DISTINCT rsid FROM wellness_trait")
    tracked_rsids = [r["rsid"] for r in snp_rows]

    rows_to_insert = []
    count = 0
    FINNGEN_ASSOC_BASE = "https://r10.finngen.fi/api/pheno"

    for rsid in tracked_rsids:
        try:
            url = f"https://r10.finngen.fi/api/variants/{rsid}"
            req = urllib.request.Request(
                url, headers={"User-Agent": "HealthAgent/1.0", "Accept": "application/json"}
            )
            with urllib.request.urlopen(req, timeout=15) as resp:
                data = json.loads(resp.read())
        except Exception:
            continue   # SNP not in FinnGen; skip silently

        for assoc in data.get("phewas", [])[:20]:  # top 20 phenotypes per SNP
            rows_to_insert.append((
                rsid,
                assoc.get("phenostring", ""),
                assoc.get("category", ""),
                "FinnGen R10",
                assoc.get("num_cases"),
                assoc.get("num_controls"),
                assoc.get("pval"),
                assoc.get("beta"),
            ))
            count += 1

        time.sleep(0.2)

    if rows_to_insert:
        local_db.executemany(
            """INSERT OR IGNORE INTO biobank_phenotype
               (rsid, phenotype, phenotype_category, study_name,
                n_cases, n_controls, p_value, beta, source)
               VALUES (?,?,?,?,?,?,?,?,'finngen')""",
            rows_to_insert,
        )

    local_db.log_download("finngen", "completed", records=count)
    progress_cb(f"[FinnGen] Done — {count:,} phenotype associations stored.")
    return count


# ═══════════════════════════════════════════════════════════════════
#  OPENTARGETS (gene-disease association scores)
# ═══════════════════════════════════════════════════════════════════

def download_opentargets(progress_cb: Callable[[str], None] = print) -> int:
    """Query OpenTargets Platform for gene-disease association scores.

    What it is:
        OpenTargets aggregates evidence from genetics, somatic mutations,
        drugs, literature, and animal models to score gene-disease associations.
        The platform integrates GWAS, ClinVar, ChEMBL, and 15 other sources.

    Why it matters for users:
        It gives a single composite 'confidence score' (0-1) for how strongly
        a gene has been linked to a disease across *all* evidence types —
        not just genetics, but also drugs, animal studies, and clinical data.

    Returns:
        Number of associations inserted.
    """
    local_db.log_download("opentargets", "started")
    progress_cb("[OpenTargets] Querying gene-disease associations via GraphQL API...")

    # Use the OpenTargets GraphQL API (no auth required, free tier)
    GRAPHQL_URL = "https://api.platform.opentargets.org/api/v4/graphql"

    # Get genes from tracked wellness SNPs
    snp_rows = local_db.query(
        "SELECT DISTINCT gene FROM wellness_trait WHERE gene IS NOT NULL AND gene != ''"
    )
    genes = [r["gene"] for r in snp_rows]

    if not genes:
        progress_cb("[OpenTargets] No genes to query.")
        local_db.log_download("opentargets", "completed", records=0)
        return 0

    rows_to_insert = []
    count = 0

    QUERY = """
    query GeneAssociations($symbol: String!) {
      target(ensemblId: $symbol) {
        approvedSymbol
        id
        associatedDiseases(page: {index: 0, size: 20}) {
          rows {
            disease { name id }
            score
            datatypeScores {
              id score
            }
          }
        }
      }
    }
    """

    # Also try by gene symbol using search
    SEARCH_QUERY = """
    query SearchGene($q: String!) {
      search(queryString: $q, entityNames: ["target"], page: {index: 0, size: 1}) {
        hits { id object { ... on Target { approvedSymbol } } }
      }
    }
    """

    for gene in genes:
        # First resolve gene symbol → Ensembl ID
        try:
            payload = json.dumps({"query": SEARCH_QUERY, "variables": {"q": gene}}).encode()
            req = urllib.request.Request(
                GRAPHQL_URL,
                data=payload,
                headers={"Content-Type": "application/json", "User-Agent": "HealthAgent/1.0"},
            )
            with urllib.request.urlopen(req, timeout=20) as resp:
                result = json.loads(resp.read())

            hits = result.get("data", {}).get("search", {}).get("hits", [])
            if not hits:
                continue
            ensembl_id = hits[0].get("id", "")
        except Exception:
            continue

        # Fetch associations using Ensembl ID
        try:
            payload = json.dumps({"query": QUERY, "variables": {"symbol": ensembl_id}}).encode()
            req = urllib.request.Request(
                GRAPHQL_URL,
                data=payload,
                headers={"Content-Type": "application/json", "User-Agent": "HealthAgent/1.0"},
            )
            with urllib.request.urlopen(req, timeout=20) as resp:
                result = json.loads(resp.read())

            target = result.get("data", {}).get("target") or {}
            assoc_diseases = target.get("associatedDiseases", {}).get("rows", [])

            for assoc in assoc_diseases:
                disease = assoc.get("disease", {})
                scores  = {s["id"]: s["score"] for s in assoc.get("datatypeScores", [])}
                rows_to_insert.append((
                    gene, ensembl_id,
                    disease.get("name", ""),
                    disease.get("id", ""),
                    assoc.get("score", 0),
                    scores.get("genetic_association", 0),
                    scores.get("somatic_mutation", 0),
                    scores.get("known_drug", 0),
                    scores.get("literature", 0),
                ))
                count += 1

        except Exception:
            pass

        time.sleep(0.5)   # be polite to the API

    if rows_to_insert:
        local_db.executemany(
            """INSERT OR IGNORE INTO opentargets_association
               (gene_symbol, ensembl_id, disease_name, disease_id,
                overall_score, genetic_score, somatic_score, drug_score, literature_score)
               VALUES (?,?,?,?,?,?,?,?,?)""",
            rows_to_insert,
        )

    local_db.log_download("opentargets", "completed", records=count)
    progress_cb(f"[OpenTargets] Done — {count:,} gene-disease associations stored.")
    return count


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    parser = argparse.ArgumentParser(description="HealthAgent database downloader")
    parser.add_argument("--all",      action="store_true", help="Download all sources")
    parser.add_argument("--gwas",     action="store_true", help="Download GWAS Catalog")
    parser.add_argument("--clinvar",  action="store_true", help="Download ClinVar")
    parser.add_argument("--pharmgkb", action="store_true", help="Download PharmGKB")
    parser.add_argument("--wellness",    action="store_true", help="Seed wellness traits")
    parser.add_argument("--disgenet",    action="store_true", help="Download DisGeNET")
    parser.add_argument("--ensembl",     action="store_true", help="Fetch Ensembl VEP consequences")
    parser.add_argument("--finngen",     action="store_true", help="Download FinnGen R10 associations")
    parser.add_argument("--opentargets", action="store_true", help="Download OpenTargets associations")
    parser.add_argument("--status",      action="store_true", help="Show DB stats")
    args = parser.parse_args()

    local_db.init_db()

    if args.status:
        stats = local_db.get_db_stats()
        print("\n── HealthAgent Local Database Status ──")
        for table, count in stats.items():
            print(f"  {table:<25} {count:>10,} rows")
        return

    if args.wellness or args.all:
        seed_wellness_traits()

    if args.gwas or args.all:
        download_gwas()

    if args.clinvar or args.all:
        download_clinvar()

    if args.pharmgkb or args.all:
        download_pharmgkb()

    if args.disgenet or args.all:
        download_disgenet()

    if args.ensembl or args.all:
        download_ensembl_consequences()

    if args.finngen or args.all:
        download_finngen()

    if args.opentargets or args.all:
        download_opentargets()


if __name__ == "__main__":
    main()
