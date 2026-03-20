"""FastAPI web server for the HealthAgent DNA Wellness Dashboard."""

import asyncio
import base64
import logging
import tempfile
import time
from pathlib import Path
from typing import Optional

from fastapi import FastAPI, File, Form, UploadFile, HTTPException, BackgroundTasks
from pydantic import BaseModel
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from fastapi.requests import Request

from healthagent.dna_importer import import_dna_string, DNAProfile, SNP, DNAFormat
from healthagent.databases import local_db
from healthagent.databases.local_db import save_profile, load_profile
from healthagent.databases.downloader import (
    seed_wellness_traits,
    download_gwas,
    download_clinvar,
    download_pharmgkb,
    download_disgenet,
    download_ensembl_consequences,
    download_finngen,
    download_opentargets,
)
from healthagent.health_traits import analyze_profile
from healthagent.databases.tcga_client import query_tcga_for_rsids, get_cached_tcga
from healthagent.microbiome_importer import (
    import_pathabundance,
    get_microbiome_profile,
    list_sessions,
)

log = logging.getLogger(__name__)

app = FastAPI(title="HealthAgent DNA Wellness Dashboard")

BASE = Path(__file__).parent
templates = Jinja2Templates(directory=str(BASE / "templates"))
app.mount("/static", StaticFiles(directory=str(BASE / "static")), name="static")

# ── Session state ─────────────────────────────────────────────────
_profile: Optional[DNAProfile] = None

# ── Chunked upload state ───────────────────────────────────────────
# Maps upload_id → list of (chunk_index, encrypted_bytes) tuples
_chunk_store: dict[str, list] = {}

CHUNK_UPLOAD_DIR = Path(tempfile.gettempdir()) / "healthagent_chunks"
CHUNK_UPLOAD_DIR.mkdir(parents=True, exist_ok=True)


class CommitRequest(BaseModel):
    upload_id: str
    aes_key_b64: str
    filename: str = "upload.txt"

# ── DB update tracking ────────────────────────────────────────────
_DB_UPDATE_INTERVAL_H = 24          # hours between full database refreshes
_last_use_time: float = 0.0         # epoch seconds of last API call
_db_update_task: Optional[asyncio.Task] = None
_db_initialized: bool = False


# ═══════════════════════════════════════════════════════════════════
#  Startup — init DB and seed wellness traits
# ═══════════════════════════════════════════════════════════════════

@app.on_event("startup")
async def startup():
    global _db_initialized, _profile
    local_db.init_db()
    # Seed curated wellness traits (fast, always safe to re-run)
    seed_wellness_traits(progress_cb=log.info)
    _db_initialized = True

    # Restore last uploaded DNA profile from database (survives server restarts)
    restored = load_profile()
    if restored:
        _profile = restored
        log.info(f"HealthAgent: restored profile from DB — {restored.snp_count:,} SNPs ({restored.source_format.value})")
    else:
        log.info("HealthAgent DB initialised. No saved profile found.")

    # Schedule the 24-hour background update loop
    asyncio.create_task(_db_update_loop())


async def _db_update_loop():
    """Background task: refresh open-source databases every 24 hours,
    but only when the site is being actively used (last API call < 1 hour ago).
    If site is idle, defer until next use triggers a check."""
    while True:
        await asyncio.sleep(3600)   # check every hour
        now = time.time()
        idle_hours = (now - _last_use_time) / 3600

        if idle_hours > 1.0:
            log.info("[DB] Site idle — deferring database update.")
            continue

        # Check when DB was last updated
        rows = local_db.query(
            """SELECT MAX(finished_at) AS last_update FROM download_log
               WHERE status = 'completed'"""
        )
        last_update_str = rows[0]["last_update"] if rows else None

        if last_update_str:
            import datetime
            last_dt = datetime.datetime.fromisoformat(last_update_str)
            age_h = (datetime.datetime.utcnow() - last_dt).total_seconds() / 3600
            if age_h < _DB_UPDATE_INTERVAL_H:
                log.info(f"[DB] Last update {age_h:.1f}h ago — no refresh needed.")
                continue

        log.info("[DB] Starting scheduled 24-hour database refresh...")
        await asyncio.get_event_loop().run_in_executor(None, _run_full_update)


def _run_full_update():
    """Run all database downloads synchronously (called from executor).

    Order is deliberate:
      1. Wellness traits seed  — fast, always safe, curated data
      2. Ensembl VEP           — annotates our tracked SNPs first (small, fast)
      3. GWAS Catalog          — large bulk download
      4. ClinVar               — large bulk download
      5. PharmGKB              — medium, drug interactions
      6. DisGeNET              — gene-disease scored associations
      7. OpenTargets           — GraphQL per-gene queries (rate-limited)
      8. FinnGen               — biobank phenome per-SNP queries (rate-limited)
    """
    local_db.init_db()
    seed_wellness_traits(progress_cb=log.info)
    download_ensembl_consequences(progress_cb=log.info)
    download_gwas(progress_cb=log.info)
    download_clinvar(progress_cb=log.info)
    download_pharmgkb(progress_cb=log.info)
    download_disgenet(progress_cb=log.info)
    download_opentargets(progress_cb=log.info)
    download_finngen(progress_cb=log.info)
    log.info("[DB] Scheduled database refresh complete.")


def _touch_last_use():
    """Record that the site was just used."""
    global _last_use_time
    _last_use_time = time.time()


async def _maybe_trigger_update(background_tasks: BackgroundTasks):
    """Trigger a behind-the-scenes DB update on first use if DB is stale."""
    _touch_last_use()
    rows = local_db.query(
        """SELECT MAX(finished_at) AS last_update FROM download_log
           WHERE source != 'wellness' AND status = 'completed'"""
    )
    last_update = rows[0]["last_update"] if rows else None
    if last_update is None:
        # Never downloaded open-source data — kick off in background
        log.info("[DB] First use — triggering background database download.")
        background_tasks.add_task(_run_full_update)


# ═══════════════════════════════════════════════════════════════════
#  Routes
# ═══════════════════════════════════════════════════════════════════

@app.get("/", response_class=HTMLResponse)
async def index(request: Request):
    _touch_last_use()
    return templates.TemplateResponse("index.html", {"request": request})


@app.post("/api/upload")
async def upload_dna(
    file: UploadFile = File(...),
    background_tasks: BackgroundTasks = BackgroundTasks(),
):
    """Accept a raw DNA file and parse it into the session profile."""
    global _profile
    _touch_last_use()
    try:
        raw = (await file.read()).decode("utf-8", errors="replace")
        _profile = import_dna_string(raw, filename=file.filename)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc))

    # Persist to DB so profile survives server restarts
    background_tasks.add_task(save_profile, _profile)
    await _maybe_trigger_update(background_tasks)

    return {
        "format":    _profile.source_format.value,
        "snp_count": _profile.snp_count,
        "filename":  file.filename,
    }


@app.post("/api/upload/chunk")
async def upload_chunk(
    upload_id:    str        = Form(...),
    chunk_index:  int        = Form(...),
    total_chunks: int        = Form(...),
    filename:     str        = Form("upload.txt"),
    chunk:        UploadFile = File(...),
):
    """Receive one AES-GCM encrypted chunk of a large DNA file."""
    _touch_last_use()
    data = await chunk.read()
    # Store chunk on disk to avoid memory pressure for large files
    chunk_dir = CHUNK_UPLOAD_DIR / upload_id
    chunk_dir.mkdir(parents=True, exist_ok=True)
    (chunk_dir / f"{chunk_index:05d}.bin").write_bytes(data)
    return {"received": chunk_index, "total": total_chunks}


@app.post("/api/upload/commit")
async def commit_upload(
    req: CommitRequest,
    background_tasks: BackgroundTasks = BackgroundTasks(),
):
    """Reassemble and decrypt all chunks, then parse the DNA file.

    The client sends the AES-256-GCM key (base64) used to encrypt each chunk.
    Each chunk blob is: [12-byte IV][AES-GCM ciphertext].
    Chunks are sorted by index, decrypted, and concatenated before parsing.
    """
    global _profile
    _touch_last_use()

    chunk_dir = CHUNK_UPLOAD_DIR / req.upload_id
    if not chunk_dir.exists():
        raise HTTPException(status_code=404, detail="Upload session not found. Upload chunks first.")

    from cryptography.hazmat.primitives.ciphers.aead import AESGCM

    try:
        aes_key = base64.b64decode(req.aes_key_b64)
        aesgcm  = AESGCM(aes_key)
    except Exception:
        raise HTTPException(status_code=400, detail="Invalid AES key.")

    # Sort chunks by index and decrypt
    chunk_files = sorted(chunk_dir.glob("*.bin"), key=lambda p: int(p.stem))
    if not chunk_files:
        raise HTTPException(status_code=400, detail="No chunks found for this upload_id.")

    parts = []
    for cf in chunk_files:
        blob = cf.read_bytes()
        if len(blob) < 13:
            raise HTTPException(status_code=400, detail=f"Corrupt chunk: {cf.name}")
        iv         = blob[:12]
        ciphertext = blob[12:]
        try:
            plain = aesgcm.decrypt(iv, ciphertext, None)
        except Exception:
            raise HTTPException(status_code=400, detail=f"Decryption failed on chunk {cf.name}. Key mismatch?")
        parts.append(plain)

    # Clean up temp files
    import shutil
    shutil.rmtree(chunk_dir, ignore_errors=True)

    raw = b"".join(parts).decode("utf-8", errors="replace")

    try:
        _profile = import_dna_string(raw, filename=req.filename)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=f"DNA parse error: {exc}")

    # Persist to DB so profile survives server restarts
    background_tasks.add_task(save_profile, _profile)
    await _maybe_trigger_update(background_tasks)

    return {
        "format":    _profile.source_format.value,
        "snp_count": _profile.snp_count,
        "filename":  req.filename,
    }


@app.get("/api/demo")
async def load_demo(background_tasks: BackgroundTasks = BackgroundTasks()):
    """Load a synthetic demo profile with all wellness SNPs injected."""
    global _profile
    _touch_last_use()
    import random

    random.seed(42)
    bases  = "AGTC"
    chroms = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
    snps   = []

    # Random background SNPs
    for i in range(2000):
        chrom = random.choice(chroms)
        b1 = random.choice(bases)
        b2 = random.choice(bases)
        snps.append(SNP(
            rsid=f"rs{1000000 + i}",
            chromosome=chrom,
            position=random.randint(100_000, 50_000_000),
            genotype=b1 + b2,
        ))

    # Inject all known wellness SNPs with realistic demo genotypes
    DEMO_WELLNESS = {
        "rs1801133": ("1", 11796321,  "CT"),   # MTHFR — mild folate variant
        "rs4988235": ("2", 136608646, "TT"),   # LCT — lactose tolerant
        "rs762551":  ("15", 75041917, "AA"),   # CYP1A2 — fast caffeine
        "rs9939609":  ("16", 53820527, "AT"),  # FTO — moderate appetite
        "rs1815739":  ("11", 66328095, "CT"),  # ACTN3 — mixed muscle
        "rs8192678":  ("4", 23821950,  "GG"),  # PPARGC1A — good aerobic
        "rs12736689": ("4", 56095230,  "CT"),  # CLOCK — intermediate
        "rs4680":     ("22", 19951271, "AG"),  # COMT — balanced
        "rs1800629":  ("6", 31543031,  "AG"),  # TNF — moderate inflammation
        "rs2243250":  ("5", 132660696, "TC"),  # IL4 — moderate allergy
        "rs2282679":  ("4", 72618334,  "AC"),  # GC — moderate VitD transport
        "rs1544410":  ("12", 48272895, "CT"),  # VDR — typical receptor
        "rs7903146":  ("10", 114758349,"TC"),  # TCF7L2 — mild blood sugar
        "rs1333049":  ("9", 22115026,  "GC"),  # 9p21 — mild heart aware
        "rs1800497":  ("11", 113270828,"GA"),  # ANKK1 — moderate reward
        "rs6265":     ("11", 27679916, "CT"),  # BDNF — mixed memory
    }
    for rsid, (chrom, pos, geno) in DEMO_WELLNESS.items():
        snps.append(SNP(rsid=rsid, chromosome=chrom, position=pos, genotype=geno))

    _profile = DNAProfile(
        source_format=DNAFormat.TWENTYTHREE_AND_ME,
        snps=snps,
        metadata={"filename": "demo"},
    )
    await _maybe_trigger_update(background_tasks)

    return {"format": "23andme", "snp_count": len(snps), "filename": "demo"}


@app.get("/api/traits")
async def get_traits():
    """Return health trait cards derived from the loaded SNP profile."""
    _touch_last_use()
    if _profile is None:
        raise HTTPException(status_code=404, detail="No DNA profile loaded.")
    result = analyze_profile(_profile)
    return result


@app.get("/api/snps")
async def get_snps(
    chromosome: Optional[str] = None,
    limit: int = 500,
    offset: int = 0,
):
    """Return a page of SNP records, optionally filtered by chromosome."""
    _touch_last_use()
    if _profile is None:
        raise HTTPException(status_code=404, detail="No DNA profile loaded. Upload a file first.")

    snps = _profile.snps
    if chromosome:
        snps = [s for s in snps if s.chromosome == chromosome.upper().lstrip("CHR")]

    page = snps[offset: offset + limit]
    return {
        "total":  len(snps),
        "offset": offset,
        "limit":  limit,
        "snps": [
            {"rsid": s.rsid, "chromosome": s.chromosome,
             "position": s.position, "genotype": s.genotype}
            for s in page
        ],
    }


@app.get("/api/profile")
async def get_profile():
    """Return metadata about the currently loaded DNA profile."""
    _touch_last_use()
    if _profile is None:
        raise HTTPException(status_code=404, detail="No DNA profile loaded.")

    # Chromosome breakdown
    chrom_counts: dict[str, int] = {}
    for s in _profile.snps:
        chrom_counts[s.chromosome] = chrom_counts.get(s.chromosome, 0) + 1

    chroms_sorted = sorted(
        chrom_counts.keys(),
        key=lambda c: (not c.isdigit(), c.zfill(3)),
    )

    # Genotype stats
    homo_ref = sum(1 for s in _profile.snps if len(s.genotype) == 2 and s.genotype[0] == s.genotype[1])
    het      = sum(1 for s in _profile.snps if len(s.genotype) == 2 and s.genotype[0] != s.genotype[1])
    no_call  = sum(1 for s in _profile.snps if s.genotype in ("--", "00", ""))

    # DB persistence info
    meta = local_db.query(
        "SELECT status, records_added, started_at FROM download_log WHERE source='profile_meta' ORDER BY id DESC LIMIT 1"
    )
    saved_at = meta[0]["started_at"] if meta else None

    return {
        "format":          _profile.source_format.value,
        "snp_count":       _profile.snp_count,
        "chromosomes":     chroms_sorted,
        "chrom_counts":    {c: chrom_counts[c] for c in chroms_sorted},
        "homozygous":      homo_ref,
        "heterozygous":    het,
        "no_call":         no_call,
        "restored":        _profile.metadata.get("restored", False),
        "saved_at":        saved_at,
        "metadata":        {k: v for k, v in _profile.metadata.items() if k != "restored"},
    }


@app.get("/api/chromosomes")
async def get_chromosomes():
    """Return chromosomes present in the loaded profile."""
    _touch_last_use()
    if _profile is None:
        raise HTTPException(status_code=404, detail="No DNA profile loaded.")
    chroms = sorted(
        {s.chromosome for s in _profile.snps},
        key=lambda c: (not c.isdigit(), c.zfill(3)),
    )
    return {"chromosomes": chroms}


# Chromosome display names (AncestryDNA uses 23=X, 24=Y, 25=MT)
_CHROM_DISPLAY = {
    **{str(i): str(i) for i in range(1, 23)},
    "23": "X", "X": "X",
    "24": "Y", "Y": "Y",
    "25": "MT", "MT": "MT", "M": "MT",
}
_CHROM_ORDER = [str(i) for i in range(1, 26)] + ["X", "Y", "MT"]


@app.get("/api/scan")
async def scan_chromosomes():
    """Return per-chromosome breakdown of the loaded profile.

    Splits all SNPs into their 25 chromosomal groups and returns:
    - snp_count, wellness_hits, drug_hits, heterozygous, homozygous per strand
    - Ordered 1-22, X(23), Y(24), MT(25)
    """
    _touch_last_use()
    if _profile is None:
        raise HTTPException(status_code=404, detail="No DNA profile loaded.")

    # Build chromosome → SNP list map
    from collections import defaultdict
    chrom_snps: dict[str, list] = defaultdict(list)
    for s in _profile.snps:
        chrom_snps[s.chromosome].append(s)

    # Pre-load wellness rsid set and genotype map once
    wellness_rows = local_db.query("SELECT rsid, genotype FROM wellness_trait")
    _COMPLEMENT_TR = str.maketrans("ACGT", "TGCA")

    def _geno_variants(g: str):
        c = g.translate(_COMPLEMENT_TR)
        return {g, g[::-1], c, c[::-1]}

    # Build wellness lookup: rsid → set of matching genotypes
    wellness_genos: dict[str, set] = {}
    for r in wellness_rows:
        wellness_genos.setdefault(r["rsid"], set()).add(r["genotype"])

    # Drug rsids that match this profile (pre-computed set)
    drug_rsids = {
        r["rsid"] for r in local_db.query(
            "SELECT DISTINCT rsid FROM drug_interaction WHERE rsid IS NOT NULL AND rsid != ''"
        )
    }

    results = []
    all_chroms = sorted(chrom_snps.keys(), key=lambda c: (not c.isdigit(), c.zfill(3)))

    for chrom in all_chroms:
        snps = chrom_snps[chrom]
        display = _CHROM_DISPLAY.get(chrom, chrom)

        wellness_hits, drug_hits = 0, 0
        homo, het, nocall = 0, 0, 0

        for s in snps:
            g = s.genotype.upper()
            # Genotype stats
            if not g or g in ("--", "00"):
                nocall += 1
            elif len(g) == 2 and g[0] == g[1]:
                homo += 1
            elif len(g) == 2:
                het += 1

            # Wellness match
            if s.rsid in wellness_genos:
                expected = wellness_genos[s.rsid]
                variants = _geno_variants(g) if g else set()
                if expected & variants:
                    wellness_hits += 1

            # Drug match
            if s.rsid in drug_rsids:
                drug_hits += 1

        results.append({
            "chromosome":    chrom,
            "display":       display,
            "snp_count":     len(snps),
            "homozygous":    homo,
            "heterozygous":  het,
            "no_call":       nocall,
            "wellness_hits": wellness_hits,
            "drug_hits":     drug_hits,
            "has_findings":  (wellness_hits + drug_hits) > 0,
        })

    total_snps = sum(r["snp_count"] for r in results)
    return {
        "chromosomes":   results,
        "total_parts":   len(results),
        "total_snps":    total_snps,
        "source_format": _profile.source_format.value,
    }


@app.get("/api/db/status")
async def db_status():
    """Return database download status and row counts."""
    _touch_last_use()
    stats = local_db.get_db_stats()
    logs  = local_db.query(
        """SELECT source, status, records_added, finished_at
           FROM download_log ORDER BY id DESC LIMIT 20"""
    )
    return {
        "db_stats":  stats,
        "last_idle_h": round((time.time() - _last_use_time) / 3600, 2),
        "update_interval_h": _DB_UPDATE_INTERVAL_H,
        "recent_downloads": [dict(r) for r in logs],
    }


@app.post("/api/db/update")
async def trigger_db_update(background_tasks: BackgroundTasks):
    """Manually trigger a database refresh in the background."""
    _touch_last_use()
    background_tasks.add_task(_run_full_update)
    return {"status": "started", "message": "Database update started in background."}


@app.get("/api/tcga")
async def get_tcga(background_tasks: BackgroundTasks):
    """Return TCGA cancer mutation context for the loaded DNA profile.

    Returns cached data immediately if available.
    Queues a background GDC API fetch for any rsIDs not yet cached.
    """
    _touch_last_use()
    if _profile is None:
        raise HTTPException(status_code=404, detail="No DNA profile loaded.")

    rsid_list = [s.rsid for s in _profile.snps]

    # Return cached data immediately
    cached = get_cached_tcga(rsid_list)

    # Queue background fetch for any missing rsIDs
    cached_rsids = {r["rsid"] for r in cached}
    missing = [r for r in rsid_list if r not in cached_rsids]
    if missing:
        background_tasks.add_task(
            query_tcga_for_rsids, missing[:200], None, log.info
        )

    return {
        "tcga":          cached,
        "total":         len(cached),
        "pending_fetch": len(missing),
        "disclaimer": (
            "TCGA data shows somatic (tumour) mutations found in cancer research. "
            "Your germline DNA variants may differ — this is NOT a cancer risk prediction."
        ),
    }


# ═══════════════════════════════════════════════════════════════════
#  Microbiome routes
# ═══════════════════════════════════════════════════════════════════

_microbiome_session: Optional[str] = None   # most-recently uploaded session


@app.post("/api/microbiome/upload")
async def upload_microbiome(file: UploadFile = File(...)):
    """Accept a HUMAnN *_pathabundance.tsv file and store pathway data."""
    global _microbiome_session
    _touch_last_use()
    if not file.filename.endswith(".tsv"):
        raise HTTPException(
            status_code=400,
            detail="Please upload a HUMAnN *_pathabundance.tsv file.",
        )
    try:
        content = (await file.read()).decode("utf-8", errors="replace")
        result  = import_pathabundance(content, filename=file.filename)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc))

    if "error" in result:
        raise HTTPException(status_code=422, detail=result["error"])

    _microbiome_session = result["session_id"]
    return result


@app.get("/api/microbiome")
async def get_microbiome(session_id: Optional[str] = None):
    """Return annotated microbiome pathway profile."""
    _touch_last_use()
    sid = session_id or _microbiome_session
    if not sid:
        raise HTTPException(status_code=404, detail="No microbiome data uploaded yet.")
    profile = get_microbiome_profile(sid)
    if not profile["pathways"]:
        raise HTTPException(status_code=404, detail="No pathways found for this session.")
    return profile


@app.get("/api/microbiome/sessions")
async def get_microbiome_sessions():
    """List all microbiome upload sessions."""
    _touch_last_use()
    return {"sessions": list_sessions()}
