"""Health trait analysis — queries the local DB for a loaded DNA profile.

Combines data from:
  - wellness_trait  (curated SNP → trait mappings)
  - gwas_association (open GWAS Catalog data)
  - clinvar_variant  (clinical significance)
  - drug_interaction (PharmGKB pharmacogenomics)
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from healthagent.databases import local_db
from healthagent.databases.tcga_client import get_cached_tcga, query_tcga_for_rsids

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


def analyze_profile(profile: "DNAProfile") -> dict:
    """
    Run a full analysis of the DNA profile against the local database.

    Returns a dict with:
      - wellness:   list of matched wellness trait cards
      - gwas:       list of GWAS associations for profile SNPs
      - clinvar:    list of ClinVar clinical findings
      - drugs:      list of medication interaction findings
      - summary:    high-level counts and category breakdown
    """
    local_db.init_db()

    # Build lookup maps from profile
    rsid_list  = [s.rsid for s in profile.snps]
    geno_map   = {s.rsid: s.genotype.upper() for s in profile.snps}
    gene_set   = {s.rsid: None for s in profile.snps}  # genes filled from DB

    # ── Wellness traits ───────────────────────────────────────────
    wellness_rows = local_db.get_wellness_traits(rsid_list, geno_map)

    # Enrich with category and risk metadata
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

    # ── GWAS associations ─────────────────────────────────────────
    gwas_rows = local_db.get_gwas_associations(rsid_list, limit=3)
    gwas = []
    for row in gwas_rows:
        trait_name = row.get("trait") or row.get("trait_name") or ""
        if not trait_name:
            continue  # skip unknown/empty GWAS records
        gwas.append({
            **row,
            "trait_name":   trait_name,
            "category":     "gwas",
            "detail":       f"{row.get('study_title', '')} (p={_fmt_pval(row.get('p_value'))})",
            "your_result":  f"OR: {_fmt_or(row.get('odds_ratio'))}  Risk allele: {row.get('risk_allele', '—')}",
            "p_value_fmt":  _fmt_pval(row.get("p_value")),
            "or_fmt":       _fmt_or(row.get("odds_ratio")),
        })

    # ── ClinVar ───────────────────────────────────────────────────
    clinvar_rows = local_db.get_clinvar_variants(rsid_list)
    clinvar = []
    for row in clinvar_rows:
        sig = row.get("clinical_sig", "")
        condition = row.get("condition") or ""
        gene = row.get("gene") or ""
        trait_name = condition.split(";")[0].strip() if condition else (gene or "")
        if not trait_name:
            continue  # skip empty ClinVar records
        clinvar.append({
            **row,
            "trait_name":        trait_name,
            "category":          "clinvar",
            "your_result":       _clinvar_plain(sig),
            "detail":            condition,
            "plain_significance": _clinvar_plain(sig),
            "severity_color":    _clinvar_color(sig),
        })

    # ── Gather gene names from all sources ────────────────────────
    genes: list[str] = []
    # From wellness_trait table (always populated)
    wt_genes = local_db.query(
        "SELECT DISTINCT gene FROM wellness_trait WHERE gene IS NOT NULL AND gene != ''"
    )
    genes = list({r["gene"] for r in wt_genes})
    # Also from snp table if populated (snp table is small — only tracked rsids)
    if rsid_list:
        snp_genes = local_db.query(
            "SELECT DISTINCT gene FROM snp WHERE gene IS NOT NULL AND gene != ''"
        )
        genes = list(set(genes) | {r["gene"] for r in snp_genes if r["gene"]})

    # ── Drug interactions ─────────────────────────────────────────
    drug_rows = local_db.get_drug_interactions(rsids=rsid_list, genes=genes)
    drugs = []
    for row in drug_rows:
        drug_name = row.get("drug_name") or ""
        if not drug_name:
            continue
        gene = row.get("gene") or ""
        drugs.append({
            **row,
            "trait_name": f"{drug_name}" + (f" ({gene})" if gene else ""),
            "category":   "medications",
            "your_result": row.get("phenotype") or "",
            "detail":      row.get("plain_english") or row.get("phenotype") or "",
        })

    # ── DisGeNET gene-disease associations ────────────────────────
    disgenet_rows = local_db.get_disgenet_diseases(genes, min_score=0.3)
    disgenet = [dict(r) for r in disgenet_rows]

    # ── OpenTargets gene-disease scores ──────────────────────────
    opentargets_rows = local_db.get_opentargets(genes)
    opentargets = [dict(r) for r in opentargets_rows]

    # ── Ensembl VEP consequences ──────────────────────────────────
    ensembl_rows = local_db.get_ensembl_consequences(rsid_list)
    ensembl = []
    for row in ensembl_rows:
        ensembl.append({
            **row,
            "consequence_plain": _consequence_plain(row.get("consequence", "")),
            "impact_color":      _impact_color(row.get("impact", "")),
        })

    # ── Biobank phenotypes (FinnGen etc.) ─────────────────────────
    biobank_rows = local_db.get_biobank_phenotypes(rsid_list)
    biobank = [dict(r) for r in biobank_rows]

    # ── TCGA cancer mutation context (cached, per-user SNPs only) ─
    tcga = get_cached_tcga(rsid_list)

    # ── Summary ───────────────────────────────────────────────────
    cat_counts: dict[str, int] = {}
    for w in wellness:
        cat_counts[w["category"]] = cat_counts.get(w["category"], 0) + 1

    worth_noting = sum(1 for w in wellness if w.get("risk_level") == "elevated")

    summary = {
        "total_snps":        len(rsid_list),
        "traits_found":      len(wellness),
        "gwas_found":        len(gwas),
        "clinvar_found":     len(clinvar),
        "drugs_found":       len(drugs),
        "tcga_found":        len(tcga),
        "worth_noting":      worth_noting,
        "category_counts":   cat_counts,
        "db_stats":          local_db.get_db_stats(),
    }

    return {
        "wellness":    wellness,
        "gwas":        gwas,
        "clinvar":     clinvar,
        "drugs":       drugs,
        "disgenet":    disgenet,
        "opentargets": opentargets,
        "ensembl":     ensembl,
        "biobank":     biobank,
        "tcga":        tcga,
        "summary":     summary,
    }


# ── Formatting helpers ────────────────────────────────────────────

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
