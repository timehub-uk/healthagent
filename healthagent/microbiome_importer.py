"""HUMAnN microbiome pathway abundance importer.

Accepts the *_pathabundance.tsv output from HUMAnN 3.x and stores
community-level pathway abundances in local SQLite.

What HUMAnN pathabundance.tsv looks like:
    # Pathway abundance
    UNMAPPED                                    1234567.0
    UNINTEGRATED                                987654.0
    UNINTEGRATED|g__Bacteroides.s__...          12345.0
    GLYCOLYSIS-PWY: Glycolysis I                45678.0
    GLYCOLYSIS-PWY: Glycolysis I|g__Bacte...    12345.0

We import only the community-level rows (no "|" species stratification).
Abundance is in RPK (reads per kilobase). We also compute a relative %
from the total mapped abundance for easier consumer display.
"""

import datetime
import io
import logging
from typing import Optional

from healthagent.databases import local_db

log = logging.getLogger(__name__)

# ── Pathway knowledge base ────────────────────────────────────────
# Maps MetaCyc pathway ID prefixes / name keywords → consumer context.
# Checked in order — first match wins.
# Format: (id_prefix_or_None, name_keyword, category, icon, plain_english, signal_logic)
# signal_logic: "high_good" | "high_bad" | "present_good" | "neutral"

PATHWAY_CATALOG: list[dict] = [
    # ── Short-chain fatty acids (SCFAs) — most important gut health signals
    {
        "id": "PWY-5022",
        "name_kw": "butanoate",
        "category": "scfa",
        "icon": "🛡️",
        "label": "Butyrate Production",
        "plain": (
            "Your gut bacteria are making butyrate — a powerful anti-inflammatory "
            "fat that feeds and protects your gut lining. Higher levels are linked "
            "to better gut health, reduced inflammation, and a healthier colon."
        ),
        "signal": "high_good",
    },
    {
        "id": "PWY-5971",
        "name_kw": "propionate",
        "category": "scfa",
        "icon": "🛡️",
        "label": "Propionate Production",
        "plain": (
            "Propionate is a short-chain fatty acid your gut microbes produce from "
            "fibre. It helps regulate blood sugar, reduces cholesterol synthesis, "
            "and signals fullness to your brain."
        ),
        "signal": "high_good",
    },
    {
        "id": "PWY-5177",
        "name_kw": "acetate",
        "category": "scfa",
        "icon": "🛡️",
        "label": "Acetate Production",
        "plain": (
            "Acetate is the most abundant short-chain fatty acid in your gut. "
            "It is the raw material other bacteria use to make butyrate, and "
            "helps keep your gut environment slightly acidic and healthy."
        ),
        "signal": "high_good",
    },
    {
        "id": "ANAEROFRUCAT-PWY",
        "name_kw": "homolactic",
        "category": "scfa",
        "icon": "🥛",
        "label": "Lactate Fermentation",
        "plain": (
            "Some of your gut bacteria ferment sugars into lactate. "
            "Moderate levels support a healthy gut pH and feed butyrate-producing bacteria."
        ),
        "signal": "neutral",
    },
    # ── Energy metabolism
    {
        "id": "GLYCOLYSIS-PWY",
        "name_kw": "glycolysis",
        "category": "gut_energy",
        "icon": "⚡",
        "label": "Glycolysis (Sugar Fermentation)",
        "plain": (
            "Your gut bacteria are actively breaking down sugars to generate energy — "
            "a sign of a metabolically busy microbiome."
        ),
        "signal": "neutral",
    },
    {
        "id": "OXIDATIVEPHOS-PWY",
        "name_kw": "oxidative phosphorylation",
        "category": "gut_energy",
        "icon": "⚡",
        "label": "Oxidative Energy Production",
        "plain": (
            "A small number of your gut microbes use oxygen-based energy pathways. "
            "This is normal in the outer gut layers."
        ),
        "signal": "neutral",
    },
    {
        "id": "TCA",
        "name_kw": "tca cycle",
        "category": "gut_energy",
        "icon": "⚡",
        "label": "TCA Cycle (Krebs Cycle)",
        "plain": (
            "The TCA cycle is a core energy pathway active in your gut bacteria. "
            "It produces building blocks for amino acids and vitamins."
        ),
        "signal": "neutral",
    },
    # ── Neurotransmitter / mood precursors
    {
        "id": "TRPSYN-PWY",
        "name_kw": "tryptophan biosynthesis",
        "category": "gut_brain",
        "icon": "🧠",
        "label": "Tryptophan Biosynthesis",
        "plain": (
            "Your gut microbes are making tryptophan — the amino acid your body "
            "uses to produce serotonin. About 90% of your serotonin is made in "
            "the gut, so microbiome tryptophan production can influence mood and sleep."
        ),
        "signal": "high_good",
    },
    {
        "id": "PWY-6731",
        "name_kw": "gaba biosynthesis",
        "category": "gut_brain",
        "icon": "🧠",
        "label": "GABA Biosynthesis",
        "plain": (
            "Some bacteria in your gut make GABA — your brain's main calming "
            "neurotransmitter. Gut-derived GABA may influence stress response "
            "and sleep quality via the gut-brain axis."
        ),
        "signal": "high_good",
    },
    {
        "id": "DOPASYN-PWY",
        "name_kw": "dopamine biosynthesis",
        "category": "gut_brain",
        "icon": "🧠",
        "label": "Dopamine Biosynthesis",
        "plain": (
            "Gut bacteria can produce dopamine precursors that travel to the brain. "
            "This pathway contributes to the gut-brain connection affecting motivation and mood."
        ),
        "signal": "high_good",
    },
    # ── Vitamins & cofactors
    {
        "id": "PWY66-400",
        "name_kw": "vitamin b12",
        "category": "vitamins",
        "icon": "💊",
        "label": "Vitamin B12 Biosynthesis",
        "plain": (
            "Some of your gut bacteria are synthesising vitamin B12 — essential for "
            "nerve function, red blood cell production, and DNA synthesis. "
            "Higher microbial B12 production may support your B12 levels."
        ),
        "signal": "high_good",
    },
    {
        "id": "FOLSYN-PWY",
        "name_kw": "folate biosynthesis",
        "category": "vitamins",
        "icon": "💊",
        "label": "Folate Biosynthesis",
        "plain": (
            "Your gut microbiome is making folate (vitamin B9) — critical for "
            "cell division, DNA repair, and (especially important) preventing "
            "neural tube defects during pregnancy."
        ),
        "signal": "high_good",
    },
    {
        "id": "PWY-6168",
        "name_kw": "vitamin k",
        "category": "vitamins",
        "icon": "💊",
        "label": "Vitamin K Biosynthesis",
        "plain": (
            "Gut bacteria produce vitamin K2, which is important for blood clotting "
            "and bone mineralisation. Gut-derived K2 is a meaningful contribution "
            "to your daily vitamin K needs."
        ),
        "signal": "high_good",
    },
    {
        "id": "RIBOSYN2-PWY",
        "name_kw": "riboflavin",
        "category": "vitamins",
        "icon": "💊",
        "label": "Riboflavin (B2) Biosynthesis",
        "plain": (
            "Your gut microbes are making riboflavin (vitamin B2), which supports "
            "energy metabolism and helps your body use other B vitamins."
        ),
        "signal": "high_good",
    },
    # ── Inflammation / immune
    {
        "id": "PWY-7111",
        "name_kw": "lipopolysaccharide",
        "category": "inflammation",
        "icon": "🔥",
        "label": "Lipopolysaccharide (LPS) Biosynthesis",
        "plain": (
            "Some bacteria make lipopolysaccharide (LPS), a molecule in their outer "
            "wall that can trigger inflammation if it enters your bloodstream. "
            "Moderate levels are normal — high levels may indicate a leaky-gut-related concern."
        ),
        "signal": "high_bad",
    },
    {
        "id": "PEPTIDOGLYCANSYN-PWY",
        "name_kw": "peptidoglycan",
        "category": "inflammation",
        "icon": "🔥",
        "label": "Peptidoglycan Biosynthesis",
        "plain": (
            "Peptidoglycan is a component of bacterial cell walls. Some is normal; "
            "very high levels can stimulate immune responses."
        ),
        "signal": "neutral",
    },
    # ── Protein / amino acid metabolism
    {
        "id": "ARGSYNBSUB-PWY",
        "name_kw": "arginine biosynthesis",
        "category": "amino_acids",
        "icon": "🧬",
        "label": "Arginine Biosynthesis",
        "plain": (
            "Your gut microbes are making arginine, an amino acid used to produce "
            "nitric oxide — which supports blood vessel health and immune function."
        ),
        "signal": "high_good",
    },
    {
        "id": "GLUTAMATE-SYN2-PWY",
        "name_kw": "glutamate biosynthesis",
        "category": "amino_acids",
        "icon": "🧬",
        "label": "Glutamate Biosynthesis",
        "plain": (
            "Glutamate is the most abundant amino acid in your body and the "
            "precursor to GABA. Gut microbial glutamate production feeds both "
            "protein synthesis and neurotransmitter pathways."
        ),
        "signal": "neutral",
    },
    # ── Carbohydrate / fibre fermentation
    {
        "id": "PWY-6588",
        "name_kw": "pyruvate fermentation",
        "category": "fibre_fermentation",
        "icon": "🌾",
        "label": "Pyruvate Fermentation",
        "plain": (
            "This pathway shows your gut bacteria fermenting carbohydrates — "
            "a healthy sign of active fibre digestion."
        ),
        "signal": "neutral",
    },
    {
        "id": "STARCH-DEGRADATION-PWY",
        "name_kw": "starch degradation",
        "category": "fibre_fermentation",
        "icon": "🌾",
        "label": "Starch Degradation",
        "plain": (
            "Your gut bacteria are breaking down starch — especially resistant "
            "starch that reaches your colon intact. This feeds SCFA production."
        ),
        "signal": "high_good",
    },
]

# Category metadata for UI rendering
CATEGORY_META = {
    "scfa":             {"label": "Short-Chain Fatty Acids", "color": "#22c55e",  "bg": "#dcfce7"},
    "gut_energy":       {"label": "Gut Energy Metabolism",   "color": "#3b82f6",  "bg": "#dbeafe"},
    "gut_brain":        {"label": "Gut-Brain Axis",          "color": "#8b5cf6",  "bg": "#ede9fe"},
    "vitamins":         {"label": "Vitamin Biosynthesis",    "color": "#06b6d4",  "bg": "#cffafe"},
    "inflammation":     {"label": "Inflammation Signals",    "color": "#f43f5e",  "bg": "#fee2e2"},
    "amino_acids":      {"label": "Amino Acid Metabolism",   "color": "#f59e0b",  "bg": "#fef3c7"},
    "fibre_fermentation": {"label": "Fibre Fermentation",   "color": "#84cc16",  "bg": "#ecfccb"},
    "other":            {"label": "Other Pathways",          "color": "#6b7280",  "bg": "#f3f4f6"},
}


def parse_pathabundance(content: str) -> list[dict]:
    """Parse a HUMAnN pathabundance.tsv string into a list of pathway dicts.

    Skips:
    - Comment lines (start with #)
    - UNMAPPED and UNINTEGRATED rows
    - Species-stratified rows (contain "|")

    Returns list of dicts with keys: pathway_id, pathway_name, abundance_rpk
    """
    pathways = []
    for line in content.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if "|" in line:
            continue  # skip species-stratified rows
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        raw_id = parts[0].strip()
        if raw_id in ("UNMAPPED", "UNINTEGRATED"):
            continue
        try:
            abundance = float(parts[1])
        except ValueError:
            continue
        # Split "PATHWAY-ID: Pathway Name" into id and name
        if ": " in raw_id:
            pathway_id, pathway_name = raw_id.split(": ", 1)
        else:
            pathway_id = raw_id
            pathway_name = raw_id
        pathways.append({
            "pathway_id":   pathway_id.strip(),
            "pathway_name": pathway_name.strip(),
            "abundance_rpk": abundance,
        })
    return pathways


def _annotate(pathway_id: str, pathway_name: str) -> dict:
    """Match a pathway to our catalog and return annotation fields."""
    name_lower = pathway_name.lower()
    for entry in PATHWAY_CATALOG:
        if entry.get("id") and pathway_id.startswith(entry["id"]):
            return entry
        if entry.get("name_kw") and entry["name_kw"] in name_lower:
            return entry
    return {
        "category": "other",
        "icon": "🧫",
        "label": pathway_name,
        "plain": f"This pathway ({pathway_id}) was detected in your microbiome sample.",
        "signal": "neutral",
    }


def _health_signal(entry: dict, relative_pct: float) -> str:
    """Compute a health signal label based on relative abundance."""
    signal_logic = entry.get("signal", "neutral")
    if relative_pct == 0:
        return "absent"
    if signal_logic == "high_good":
        if relative_pct >= 2.0:
            return "high"
        if relative_pct >= 0.5:
            return "normal"
        return "low"
    if signal_logic == "high_bad":
        if relative_pct >= 2.0:
            return "elevated"
        return "normal"
    return "normal"


def import_pathabundance(
    content: str,
    filename: str = "upload",
    session_id: Optional[str] = None,
) -> dict:
    """Parse and store a HUMAnN pathabundance.tsv.

    Args:
        content:    Raw file content as string.
        filename:   Original filename (for display).
        session_id: Optional session key; auto-generated if not provided.

    Returns:
        Summary dict with session_id, total_pathways, notable_pathways.
    """
    local_db.init_db()

    if session_id is None:
        session_id = datetime.datetime.utcnow().strftime("%Y%m%d_%H%M%S")

    pathways = parse_pathabundance(content)
    if not pathways:
        return {"error": "No pathway data found in file.", "session_id": session_id}

    # Compute total mapped RPK for relative %
    total_rpk = sum(p["abundance_rpk"] for p in pathways) or 1.0

    rows = []
    for p in pathways:
        ann = _annotate(p["pathway_id"], p["pathway_name"])
        rel_pct = round((p["abundance_rpk"] / total_rpk) * 100, 4)
        signal  = _health_signal(ann, rel_pct)
        rows.append((
            session_id,
            p["pathway_id"],
            p["pathway_name"],
            p["abundance_rpk"],
            rel_pct,
            ann.get("category", "other"),
            ann.get("plain", ""),
            signal,
            ann.get("icon", "🧫"),
        ))

    local_db.executemany(
        """INSERT OR REPLACE INTO microbiome_pathway
           (session_id, pathway_id, pathway_name, abundance_rpk,
            relative_pct, health_category, plain_english, health_signal, icon)
           VALUES (?,?,?,?,?,?,?,?,?)""",
        rows,
    )

    local_db.execute(
        """INSERT OR REPLACE INTO microbiome_session
           (session_id, filename, total_pathways, mapped_rpk)
           VALUES (?,?,?,?)""",
        (session_id, filename, len(pathways), total_rpk),
    )

    log.info(f"[Microbiome] Imported {len(pathways)} pathways for session {session_id}")

    return {
        "session_id":       session_id,
        "filename":         filename,
        "total_pathways":   len(pathways),
        "mapped_rpk":       round(total_rpk, 2),
    }


def get_microbiome_profile(session_id: str) -> dict:
    """Return annotated microbiome pathway profile for a session."""
    local_db.init_db()

    rows = local_db.query(
        """SELECT pathway_id, pathway_name, abundance_rpk, relative_pct,
                  health_category, plain_english, health_signal, icon
           FROM microbiome_pathway
           WHERE session_id = ?
           ORDER BY relative_pct DESC""",
        (session_id,),
    )

    if not rows:
        return {"pathways": [], "categories": {}, "session_id": session_id}

    pathways = []
    categories: dict[str, list] = {}
    for r in rows:
        d = dict(r)
        cat = d.get("health_category", "other")
        meta = CATEGORY_META.get(cat, CATEGORY_META["other"])
        d["category_label"] = meta["label"]
        d["category_color"] = meta["color"]
        d["category_bg"]    = meta["bg"]
        pathways.append(d)
        categories.setdefault(cat, []).append(d)

    # Notable = non-other category, non-neutral/absent, top 20
    notable = [
        p for p in pathways
        if p["health_category"] != "other"
        and p["health_signal"] not in ("absent",)
    ][:20]

    return {
        "session_id":  session_id,
        "pathways":    pathways,
        "notable":     notable,
        "categories":  {k: v for k, v in categories.items()},
        "category_meta": CATEGORY_META,
    }


def list_sessions() -> list[dict]:
    """Return all microbiome upload sessions."""
    local_db.init_db()
    rows = local_db.query(
        "SELECT * FROM microbiome_session ORDER BY uploaded_at DESC LIMIT 20"
    )
    return [dict(r) for r in rows]
