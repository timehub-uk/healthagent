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
        gwas.append({
            **row,
            "p_value_fmt": _fmt_pval(row.get("p_value")),
            "or_fmt":      _fmt_or(row.get("odds_ratio")),
        })

    # ── ClinVar ───────────────────────────────────────────────────
    clinvar_rows = local_db.get_clinvar_variants(rsid_list)
    clinvar = []
    for row in clinvar_rows:
        sig = row.get("clinical_sig", "")
        clinvar.append({
            **row,
            "plain_significance": _clinvar_plain(sig),
            "severity_color":     _clinvar_color(sig),
        })

    # ── Drug interactions ─────────────────────────────────────────
    # Gather gene names from DB for the profile's rsids
    genes = []
    if rsid_list:
        ph = ",".join("?" * len(rsid_list))
        rows = local_db.query(
            f"SELECT DISTINCT gene FROM snp WHERE rsid IN ({ph}) AND gene != ''",
            tuple(rsid_list),
        )
        genes = [r["gene"] for r in rows if r["gene"]]

    drug_rows = local_db.get_drug_interactions(rsids=rsid_list, genes=genes)
    drugs = [dict(r) for r in drug_rows]

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
        "worth_noting":      worth_noting,
        "category_counts":   cat_counts,
        "db_stats":          local_db.get_db_stats(),
    }

    return {
        "wellness": wellness,
        "gwas":     gwas,
        "clinvar":  clinvar,
        "drugs":    drugs,
        "summary":  summary,
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


def _clinvar_color(sig: str) -> str:
    if "Pathogenic" in sig:
        return "#f43f5e"
    if "Uncertain" in sig:
        return "#f59e0b"
    if "Benign" in sig:
        return "#22c55e"
    return "#6b7280"
