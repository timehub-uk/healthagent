"""Health trait analysis — queries the local DB for a loaded DNA profile.

Combines data from:
  - wellness_trait  (curated SNP → trait mappings)
  - gwas_association (open GWAS Catalog data)
  - clinvar_variant  (clinical significance)
  - drug_interaction (PharmGKB pharmacogenomics)

All independent DB queries run in parallel via ThreadPoolExecutor.
Results are cached by profile hash to avoid redundant re-analysis.
"""

from __future__ import annotations

import hashlib
import threading
from concurrent.futures import ThreadPoolExecutor
from typing import TYPE_CHECKING

from healthagent.databases import local_db
from healthagent.databases.tcga_client import get_cached_tcga

if TYPE_CHECKING:
    from healthagent.dna_importer import DNAProfile


CATEGORY_META = {
    "nutrition":      {"label": "Food & Nutrition",        "icon": "🥗", "color": "#22c55e"},
    "fitness":        {"label": "Fitness & Energy",         "icon": "💪", "color": "#3b82f6"},
    "sleep":          {"label": "Sleep & Stress",           "icon": "🌙", "color": "#8b5cf6"},
    "allergies":      {"label": "Allergies & Immunity",     "icon": "🤧", "color": "#f59e0b"},
    "supplements":    {"label": "Vitamins & Supplements",   "icon": "💊", "color": "#06b6d4"},
    "health_history": {"label": "Health Awareness",         "icon": "❤️", "color": "#f43f5e"},
    "traits":         {"label": "Personal Traits",          "icon": "🧬", "color": "#a855f7"},
    "medications":    {"label": "Medication Interactions",  "icon": "💊", "color": "#ec4899"},
}

RISK_META = {
    "typical":  {"label": "Typical",      "color": "#22c55e", "bg": "#dcfce7"},
    "reduced":  {"label": "Reduced",      "color": "#3b82f6", "bg": "#dbeafe"},
    "variant":  {"label": "Variant",      "color": "#f59e0b", "bg": "#fef3c7"},
    "elevated": {"label": "Worth Noting", "color": "#f43f5e", "bg": "#fee2e2"},
}

# ── Profile analysis cache ─────────────────────────────────────────
_analysis_cache: dict = {}
_cache_lock = threading.Lock()


def _profile_hash(rsid_list: list[str]) -> str:
    """Deterministic hash of the rsid list for cache keying."""
    key = "|".join(sorted(rsid_list))
    return hashlib.md5(key.encode()).hexdigest()


def invalidate_cache() -> None:
    """Clear the analysis cache (call after a fresh DNA upload)."""
    with _cache_lock:
        _analysis_cache.clear()


def analyze_profile(profile: "DNAProfile") -> dict:
    """
    Run a full analysis of the DNA profile against the local database.

    Returns a dict with:
      - wellness:   list of matched wellness trait cards
      - gwas:       list of GWAS associations for profile SNPs
      - clinvar:    list of ClinVar clinical findings
      - drugs:      list of medication interaction findings
      - summary:    high-level counts and category breakdown

    All independent DB queries execute in parallel threads.
    Results are cached by profile hash; call invalidate_cache() after upload.
    """
    local_db.init_db()

    rsid_list = [s.rsid for s in profile.snps]
    geno_map  = {s.rsid: s.genotype.upper() for s in profile.snps}

    # ── Cache check ───────────────────────────────────────────────
    cache_key = _profile_hash(rsid_list)
    with _cache_lock:
        if cache_key in _analysis_cache:
            return _analysis_cache[cache_key]

    # ── Phase 1: all rsid-independent queries run in parallel ─────
    with ThreadPoolExecutor(max_workers=8, thread_name_prefix="ha-db") as pool:
        f_wellness  = pool.submit(local_db.get_wellness_traits, rsid_list, geno_map)
        f_gwas      = pool.submit(local_db.get_gwas_associations, rsid_list, 3)
        f_clinvar   = pool.submit(local_db.get_clinvar_variants, rsid_list)
        f_ensembl   = pool.submit(local_db.get_ensembl_consequences, rsid_list)
        f_biobank   = pool.submit(local_db.get_biobank_phenotypes, rsid_list)
        f_tcga      = pool.submit(get_cached_tcga, rsid_list)
        f_wt_genes  = pool.submit(
            local_db.query,
            "SELECT DISTINCT gene FROM wellness_trait WHERE gene IS NOT NULL AND gene != ''"
        )
        f_snp_genes = (
            pool.submit(local_db.query,
                        "SELECT DISTINCT gene FROM snp WHERE gene IS NOT NULL AND gene != ''")
            if rsid_list else None
        )

        wellness_rows = f_wellness.result()
        gwas_rows     = f_gwas.result()
        clinvar_rows  = f_clinvar.result()
        ensembl_rows  = f_ensembl.result()
        biobank_rows  = f_biobank.result()
        tcga          = f_tcga.result()
        wt_gene_rows  = f_wt_genes.result()
        snp_gene_rows = f_snp_genes.result() if f_snp_genes else []

    # Build gene set from phase-1 results
    genes: list[str] = list(
        {r["gene"] for r in wt_gene_rows}
        | {r["gene"] for r in snp_gene_rows if r["gene"]}
    )

    # ── Phase 2: gene-dependent queries run in parallel ───────────
    with ThreadPoolExecutor(max_workers=3, thread_name_prefix="ha-gene") as pool:
        f_drugs       = pool.submit(
            lambda: local_db.get_drug_interactions(rsids=rsid_list, genes=genes)
        )
        f_disgenet    = pool.submit(local_db.get_disgenet_diseases, genes, 0.3)
        f_opentargets = pool.submit(local_db.get_opentargets, genes)

        drug_rows        = f_drugs.result()
        disgenet_rows    = f_disgenet.result()
        opentargets_rows = f_opentargets.result()

    # ── Enrich wellness results ───────────────────────────────────
    wellness = []
    for row in wellness_rows:
        cat  = CATEGORY_META.get(row["category"], {})
        risk = RISK_META.get(row["risk_level"], RISK_META["typical"])
        wellness.append({
            **row,
            "category_label": cat.get("label", row["category"]),
            "category_color": cat.get("color", "#6b7280"),
            "risk_label":     risk["label"],
            "risk_color":     risk["color"],
            "risk_bg":        risk["bg"],
        })

    # ── Enrich GWAS results ───────────────────────────────────────
    gwas = []
    for row in gwas_rows:
        trait_name = row.get("trait") or row.get("trait_name") or ""
        if not trait_name:
            continue
        gwas.append({
            **row,
            "trait_name":  trait_name.title(),
            "category":    "gwas",
            "detail":      row.get("study_title") or "",
            "your_result": _gwas_plain_result(row.get("odds_ratio"), row.get("risk_allele")),
            "p_value_fmt": _fmt_pval(row.get("p_value")),
            "or_fmt":      _fmt_or(row.get("odds_ratio")),
        })

    # ── Enrich ClinVar results ────────────────────────────────────
    clinvar = []
    for row in clinvar_rows:
        sig       = row.get("clinical_sig", "")
        condition = row.get("condition") or ""
        gene      = row.get("gene") or ""
        trait_name = condition.split(";")[0].strip() if condition else (gene or "")
        if not trait_name:
            continue
        clinvar.append({
            **row,
            "trait_name":         trait_name,
            "category":           "clinvar",
            "your_result":        _clinvar_plain(sig),
            "detail":             condition,
            "plain_significance": _clinvar_plain(sig),
            "severity_color":     _clinvar_color(sig),
        })

    # ── Enrich drug interaction results ───────────────────────────
    drugs = []
    for row in drug_rows:
        drug_name = row.get("drug_name") or ""
        if not drug_name:
            continue
        phenotype = row.get("phenotype") or ""
        drugs.append({
            **row,
            "trait_name":  drug_name,
            "category":    "medications",
            "your_result": _drug_plain_result(phenotype, row.get("category")),
            "detail":      row.get("plain_english") or phenotype,
            "what_to_do":  (
                f"Speak to your doctor or pharmacist before taking {drug_name}. "
                "Your genetic profile may affect how this medication works for you."
            ),
        })

    # ── Enrich Ensembl results ────────────────────────────────────
    ensembl = []
    for row in ensembl_rows:
        ensembl.append({
            **row,
            "consequence_plain": _consequence_plain(row.get("consequence", "")),
            "impact_color":      _impact_color(row.get("impact", "")),
        })

    # ── Summary ───────────────────────────────────────────────────
    cat_counts: dict[str, int] = {}
    for w in wellness:
        cat_counts[w["category"]] = cat_counts.get(w["category"], 0) + 1

    worth_noting = sum(1 for w in wellness if w.get("risk_level") == "elevated")

    result = {
        "wellness":    wellness,
        "gwas":        gwas,
        "clinvar":     clinvar,
        "drugs":       drugs,
        "disgenet":    [dict(r) for r in disgenet_rows],
        "opentargets": [dict(r) for r in opentargets_rows],
        "ensembl":     ensembl,
        "biobank":     [dict(r) for r in biobank_rows],
        "tcga":        tcga,
        "summary": {
            "total_snps":      len(rsid_list),
            "traits_found":    len(wellness),
            "gwas_found":      len(gwas),
            "clinvar_found":   len(clinvar),
            "drugs_found":     len(drugs),
            "tcga_found":      len(tcga),
            "worth_noting":    worth_noting,
            "category_counts": cat_counts,
            "db_stats":        local_db.get_db_stats(),
        },
    }

    # Store in cache
    with _cache_lock:
        _analysis_cache[cache_key] = result

    return result


# ── Formatting helpers ─────────────────────────────────────────────

def _fmt_pval(p) -> str:
    if p is None:
        return "—"
    if p < 1e-50:
        return "< 10⁻⁵⁰"
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.4f}"


def _fmt_or(or_val) -> str:
    if or_val is None or or_val == 0:
        return "—"
    return f"{or_val:.2f}"


def _clinvar_plain(sig: str) -> str:
    mapping = {
        "Pathogenic":                   "Clinically significant — review with your doctor",
        "Likely pathogenic":            "Likely clinically significant",
        "Pathogenic/Likely pathogenic": "Clinically significant — review with your doctor",
        "Uncertain significance":       "Uncertain — science is still learning about this variant",
        "Likely benign":                "Likely no clinical concern",
        "Benign":                       "No clinical concern",
        "Benign/Likely benign":         "No clinical concern",
    }
    for key, val in mapping.items():
        if key in sig:
            return val
    return sig


def _consequence_plain(consequence: str) -> str:
    """Convert Ensembl consequence term to plain English."""
    mapping = {
        "missense_variant":           "Changes an amino acid in the protein",
        "synonymous_variant":         "Same amino acid — silent change",
        "stop_gained":                "Creates a premature stop codon",
        "splice_region_variant":      "May affect how the gene is read",
        "intron_variant":             "Located within a gene intron",
        "regulatory_region_variant":  "Affects gene regulatory region",
        "5_prime_UTR_variant":        "In the gene's start regulatory region",
        "3_prime_UTR_variant":        "In the gene's end regulatory region",
        "upstream_gene_variant":      "Near the start of a gene",
        "downstream_gene_variant":    "Near the end of a gene",
        "intergenic_variant":         "Between genes",
        "non_coding_transcript_exon_variant": "In a non-coding RNA region",
    }
    for key, val in mapping.items():
        if key in consequence:
            return val
    return consequence.replace("_", " ").capitalize() if consequence else "—"


def _impact_color(impact: str) -> str:
    return {
        "HIGH":     "#f43f5e",
        "MODERATE": "#f59e0b",
        "LOW":      "#3b82f6",
        "MODIFIER": "#6b7280",
    }.get(impact, "#6b7280")


def _clinvar_color(sig: str) -> str:
    if "Pathogenic" in sig:
        return "#f43f5e"
    if "Uncertain" in sig:
        return "#f59e0b"
    if "Benign" in sig:
        return "#22c55e"
    return "#6b7280"


def _gwas_plain_result(odds_ratio, risk_allele) -> str:
    """Turn a raw odds ratio into a plain-English sentence a non-scientist can understand."""
    try:
        or_val = float(odds_ratio)
    except (TypeError, ValueError):
        return "A DNA variant near this gene has been linked to this trait in research studies."

    allele_note = (
        f" (variant: {risk_allele})"
        if risk_allele and risk_allele not in ("—", "", "?")
        else ""
    )
    if or_val >= 2.0:
        return (
            f"Your DNA variant{allele_note} is associated with roughly {or_val:.1f}× higher "
            "likelihood of this trait in research studies."
        )
    elif or_val >= 1.2:
        pct = int((or_val - 1) * 100)
        return (
            f"Research links your DNA variant{allele_note} to about {pct}% higher likelihood "
            "of this trait compared to average."
        )
    elif or_val >= 0.95:
        return (
            f"Your DNA variant{allele_note} shows little difference from average "
            "for this trait in research studies."
        )
    elif or_val >= 0.5:
        pct = int((1 - or_val) * 100)
        return (
            f"Your DNA variant{allele_note} is associated with roughly {pct}% lower likelihood "
            "of this trait compared to average."
        )
    else:
        return (
            f"Your DNA variant{allele_note} is strongly associated with reduced likelihood "
            "of this trait in research studies."
        )


def _drug_plain_result(phenotype: str, category: str) -> str:
    """Convert a pharmacogenomics phenotype to a short plain-English finding."""
    if not phenotype:
        return "Your DNA may affect how your body responds to this medication."
    ph = phenotype.lower()
    if "poor" in ph or "slow" in ph or "decreased" in ph or "reduced" in ph:
        return "Your body may process this medication more slowly than average — it can build up to higher levels."
    if "rapid" in ph or "ultra" in ph or "fast" in ph or "increased" in ph:
        return "Your body may break down this medication faster than average — it may wear off sooner."
    if "intermediate" in ph:
        return "Your body processes this medication at a slightly below-average rate."
    if "toxicity" in ph or "adverse" in ph or "risk" in ph:
        return "Your DNA variant is linked to a higher risk of side effects with this medication."
    if "efficacy" in ph or "response" in ph or "effective" in ph:
        return "This medication may work differently for you based on your genetic profile."
    return phenotype
