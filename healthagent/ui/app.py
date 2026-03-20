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

# ── Download progress tracker ─────────────────────────────────────
# Maps source name → {"pct": int, "msg": str, "done": bool}
_DOWNLOAD_TASKS = [
    "wellness", "ensembl", "gwas", "clinvar", "pharmgkb",
    "disgenet", "opentargets", "finngen",
]
_dl_progress: dict = {
    "active":       False,
    "overall_pct":  0,
    "tasks":        {t: {"pct": 0, "msg": "", "done": False} for t in _DOWNLOAD_TASKS},
    "current_task": "",
    "started_at":   None,
}


def _make_progress_cb(source: str):
    """Return a progress_cb that updates _dl_progress for the given source."""
    import re
    def cb(msg: str):
        log.info(msg)
        _dl_progress["tasks"][source]["msg"] = msg
        _dl_progress["current_task"] = source
        # Parse percentage from "[Label] Downloaded X/Y MB (N%)" messages
        m = re.search(r'\((\d+)%\)', msg)
        if m:
            _dl_progress["tasks"][source]["pct"] = int(m.group(1))
        elif "Done" in msg or "done" in msg or "complete" in msg.lower() or "stored" in msg.lower():
            _dl_progress["tasks"][source]["pct"] = 100
        _update_overall_pct()
    return cb


def _update_overall_pct():
    tasks = _dl_progress["tasks"]
    done_count = sum(1 for t in tasks.values() if t["done"])
    active_pct  = sum(t["pct"] for t in tasks.values() if not t["done"])
    # Each task is worth 1/N of total; add partial credit for in-progress task
    n = len(tasks)
    overall = int((done_count * 100 + active_pct) / n)
    _dl_progress["overall_pct"] = min(overall, 99 if _dl_progress["active"] else 100)


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
    import time as _time
    _dl_progress["active"] = True
    _dl_progress["overall_pct"] = 0
    _dl_progress["started_at"] = _time.time()
    for t in _dl_progress["tasks"].values():
        t["pct"] = 0
        t["msg"] = ""
        t["done"] = False

    def _step(source, fn, **kw):
        _dl_progress["current_task"] = source
        cb = _make_progress_cb(source)
        try:
            fn(progress_cb=cb, **kw)
        except Exception as e:
            log.error(f"[DB] {source} failed: {e}")
        _dl_progress["tasks"][source]["done"] = True
        _dl_progress["tasks"][source]["pct"] = 100
        _update_overall_pct()

    local_db.init_db()
    _step("wellness",    seed_wellness_traits)
    _step("ensembl",     download_ensembl_consequences)
    _step("gwas",        download_gwas)
    _step("clinvar",     download_clinvar)
    _step("pharmgkb",   download_pharmgkb)
    _step("disgenet",    download_disgenet)
    _step("opentargets", download_opentargets)
    _step("finngen",     download_finngen)

    _dl_progress["active"] = False
    _dl_progress["overall_pct"] = 100
    _dl_progress["current_task"] = ""
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


@app.get("/api/db/progress")
async def db_progress():
    """Return real-time download progress for all database sources."""
    _touch_last_use()
    import time as _time
    elapsed = None
    if _dl_progress["started_at"]:
        elapsed = round(_time.time() - _dl_progress["started_at"])
    return {
        "active":       _dl_progress["active"],
        "overall_pct":  _dl_progress["overall_pct"],
        "current_task": _dl_progress["current_task"],
        "elapsed_s":    elapsed,
        "tasks": {
            name: {
                "pct":  info["pct"],
                "done": info["done"],
                "msg":  info["msg"][-120:] if info["msg"] else "",
            }
            for name, info in _dl_progress["tasks"].items()
        },
    }


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


# ── Physical / appearance traits derived from known SNPs ──────────
_PHYSICAL_TRAIT_RSIDS: dict[str, dict] = {
    # Eye colour
    "rs12913832": {
        "trait": "Eye Colour",
        "icon": "👁️",
        "category": "appearance",
        "gene": "HERC2/OCA2",
        "body_part": "👁️ Eyes",
        "what_it_means": (
            "The HERC2 gene controls how much OCA2 protein is made in your iris. "
            "OCA2 produces the pigment eumelanin — more of it means darker (brown) eyes, "
            "less means lighter (blue or green) eyes. This SNP is the strongest single "
            "predictor of eye colour identified by genome-wide studies."
        ),
        "interpret": lambda g: (
            ("Likely blue or green eyes", "blue")   if "AA" in g else
            ("Likely brown or hazel eyes", "brown") if "GG" in g else
            ("Mixed — hazel or green likely", "hazel")
        ),
    },
    "rs1800407": {
        "trait": "Eye Colour (OCA2 modifier)",
        "icon": "👁️",
        "category": "appearance",
        "gene": "OCA2",
        "body_part": "👁️ Eyes",
        "what_it_means": (
            "A second OCA2 variant that can independently reduce iris pigmentation. "
            "Carrying the T allele tends to lighten eye colour, particularly shifting "
            "brown eyes toward hazel or green."
        ),
        "interpret": lambda g: (
            ("OCA2 variant — may lighten eye colour", "light") if "T" in g else
            ("Standard OCA2 — no lightening effect", "standard")
        ),
    },
    # Hair colour
    "rs1805007": {
        "trait": "Red Hair Tendency",
        "icon": "🦰",
        "category": "appearance",
        "gene": "MC1R",
        "body_part": "🧑 Hair & Skin",
        "what_it_means": (
            "The MC1R gene encodes the melanocortin-1 receptor, which switches melanocytes "
            "between making dark eumelanin (brown/black) and red pheomelanin. "
            "The rs1805007 T allele reduces receptor activity, dramatically increasing "
            "pheomelanin and producing red or auburn hair. It also reduces UV protection, "
            "increasing sunburn risk — use SPF daily."
        ),
        "interpret": lambda g: (
            ("Strong red-hair variant (MC1R)", "red") if "TT" in g else
            ("One copy of red-hair MC1R variant", "auburn") if "T" in g else
            ("No red-hair MC1R variant", "none")
        ),
    },
    "rs1805008": {
        "trait": "Red Hair (MC1R variant 2)",
        "icon": "🦰",
        "category": "appearance",
        "gene": "MC1R",
        "body_part": "🧑 Hair & Skin",
        "what_it_means": (
            "A second loss-of-function MC1R variant associated with red or auburn hair, "
            "fair skin, and increased sensitivity to sunlight. Combined with rs1805007, "
            "both copies significantly elevate melanoma risk compared to non-carriers."
        ),
        "interpret": lambda g: (
            ("MC1R variant — linked to red or auburn hair", "red") if "A" in g else
            ("No second MC1R red-hair variant", "none")
        ),
    },
    # Lactose tolerance
    "rs4988235": {
        "trait": "Lactose Tolerance",
        "icon": "🥛",
        "category": "digestion",
        "gene": "LCT",
        "body_part": "🫀 Digestive System",
        "what_it_means": (
            "Most mammals lose the ability to digest lactose (milk sugar) after weaning. "
            "In some human populations, a mutation near the LCT gene keeps the lactase "
            "enzyme active into adulthood. The G allele at rs4988235 is the European "
            "lactase-persistence variant. Without it, undigested lactose reaches the colon "
            "where bacteria ferment it, causing bloating, gas, and diarrhoea."
        ),
        "interpret": lambda g: (
            ("Lactose tolerant — LCT gene stays active", "tolerant") if "G" in g else
            ("Likely lactose intolerant in adulthood", "intolerant")
        ),
    },
    # Caffeine metabolism
    "rs762551": {
        "trait": "Caffeine Metabolism",
        "icon": "☕",
        "category": "metabolism",
        "gene": "CYP1A2",
        "body_part": "🫀 Liver",
        "what_it_means": (
            "CYP1A2 is the liver enzyme responsible for breaking down about 95% of caffeine. "
            "The rs762551 A allele produces a faster enzyme. Fast metabolisers can drink coffee "
            "later in the day with less sleep disruption, while slow metabolisers may feel "
            "jitteriness, anxiety, or insomnia even from moderate amounts. Studies also link "
            "slow CYP1A2 to higher heart disease risk with heavy coffee intake."
        ),
        "interpret": lambda g: (
            ("Fast caffeine metaboliser (CYP1A2)", "fast") if "AA" in g else
            ("Slow caffeine metaboliser — sensitivity likely", "slow")
        ),
    },
    # Earwax type
    "rs17822931": {
        "trait": "Earwax Type",
        "icon": "👂",
        "category": "traits",
        "gene": "ABCC11",
        "body_part": "👂 Ears & Sweat Glands",
        "what_it_means": (
            "ABCC11 encodes a transporter that secretes fatty compounds into ear-canal glands. "
            "The wet form (C or T allele) produces sticky, yellow earwax; the dry form (TT) "
            "produces grey, flaky wax common in East Asian populations. "
            "Interestingly, this same gene also affects body-odour intensity — dry-earwax carriers "
            "produce significantly less underarm odour."
        ),
        "interpret": lambda g: (
            ("Dry earwax (East Asian variant)", "dry") if "TT" in g else
            ("Wet earwax (typical)", "wet")
        ),
    },
    # Muscle composition
    "rs1815739": {
        "trait": "Muscle Fibre Type",
        "icon": "💪",
        "category": "fitness",
        "gene": "ACTN3",
        "body_part": "💪 Skeletal Muscle",
        "what_it_means": (
            "Alpha-actinin-3 (encoded by ACTN3) is a structural protein found only in fast-twitch "
            "(type II) muscle fibres used for explosive power. The CC genotype produces functional "
            "ACTN3 — associated with sprint and power performance. The TT genotype creates a "
            "non-functional protein; the body compensates with more efficient slow-twitch fibres, "
            "favouring endurance sports. About 18% of the global population is TT."
        ),
        "interpret": lambda g: (
            ("Power-oriented muscle fibres (ACTN3 CC)", "power")    if "CC" in g else
            ("Endurance-oriented muscle fibres (ACTN3 TT)", "endurance") if "TT" in g else
            ("Mixed muscle fibre type", "mixed")
        ),
    },
    # Bitter taste perception
    "rs713598": {
        "trait": "Bitter Taste Perception",
        "icon": "👅",
        "category": "traits",
        "gene": "TAS2R38",
        "body_part": "👅 Taste Buds & Tongue",
        "what_it_means": (
            "TAS2R38 encodes a bitter taste receptor that detects compounds like PTC and PROP, "
            "found in cruciferous vegetables (broccoli, Brussels sprouts, kale). "
            "Supertasters (GG) find these very bitter and often eat fewer vegetables, "
            "which may affect cancer-protective nutrient intake. Non-tasters (CC) may use "
            "more salt to compensate for blander food perception."
        ),
        "interpret": lambda g: (
            ("Supertaster — very sensitive to bitter foods", "supertaster") if "GG" in g else
            ("Non-taster — low bitter sensitivity", "nontaster")   if "CC" in g else
            ("Average bitter taste sensitivity", "average")
        ),
    },
    # Skin / UV sensitivity (SLC24A5)
    "rs1426654": {
        "trait": "Skin Pigmentation",
        "icon": "🌞",
        "category": "appearance",
        "gene": "SLC24A5",
        "body_part": "🧑 Skin",
        "what_it_means": (
            "SLC24A5 encodes a calcium transporter in melanosomes (pigment-producing organelles). "
            "The A allele (common in Europeans) reduces melanin production, resulting in lighter skin "
            "that is more efficient at synthesising vitamin D in low-sunlight environments but burns "
            "more easily. The G allele (common in Africans and South Asians) produces more melanin, "
            "offering better natural UV protection."
        ),
        "interpret": lambda g: (
            ("Lighter skin pigmentation (SLC24A5 AA)", "light") if "AA" in g else
            ("Darker skin pigmentation (SLC24A5 GG)", "dark")   if "GG" in g else
            ("Intermediate skin pigmentation", "medium")
        ),
    },
    # Alcohol flush
    "rs671": {
        "trait": "Alcohol Flush Reaction",
        "icon": "🍷",
        "category": "metabolism",
        "gene": "ALDH2",
        "body_part": "🫀 Liver & Blood Vessels",
        "what_it_means": (
            "ALDH2 breaks down acetaldehyde, the toxic by-product of alcohol metabolism. "
            "The A allele creates a near-non-functional ALDH2 enzyme common in East Asians (~35-40%). "
            "Acetaldehyde builds up, causing facial flushing, rapid heartbeat, nausea, and headache "
            "(the 'Asian flush'). ALDH2 deficiency also significantly increases oesophageal cancer "
            "risk with regular alcohol consumption — the WHO classifies acetaldehyde as a Group 1 carcinogen."
        ),
        "interpret": lambda g: (
            ("Strong flush reaction (ALDH2 deficiency)", "flush") if "AA" in g else
            ("Mild flush tendency (ALDH2 heterozygous)", "mild")  if "A" in g else
            ("No typical flush reaction", "none")
        ),
    },
    # Deep sleep / melatonin
    "rs57875989": {
        "trait": "Sleep Duration Tendency",
        "icon": "🌙",
        "category": "sleep",
        "gene": "ADRB1",
        "body_part": "🧠 Brain & Nervous System",
        "what_it_means": (
            "ADRB1 encodes the beta-1 adrenergic receptor involved in regulating arousal and sleep. "
            "A rare variant (A allele) identified in naturally short sleepers allows individuals to "
            "function well on 6 hours or less without apparent cognitive impairment. Most people with "
            "typical genetics need 7–9 hours; chronically sleeping less increases cardiovascular and "
            "metabolic risk regardless of how rested you feel."
        ),
        "interpret": lambda g: (
            ("Genetic short-sleep variant (ADRB1)", "short") if "A" in g else
            ("Typical sleep duration genetics", "typical")
        ),
    },
}


@app.get("/api/physical_traits")
async def get_physical_traits():
    """Return DNA-derived physical and appearance trait predictions."""
    _touch_last_use()
    if _profile is None:
        return {"traits": [], "available": False}

    geno_map = {s.rsid: s.genotype.upper() for s in _profile.snps}
    results = []
    for rsid, meta in _PHYSICAL_TRAIT_RSIDS.items():
        genotype = geno_map.get(rsid)
        if genotype is None:
            continue
        label, value_key = meta["interpret"](genotype)
        results.append({
            "rsid":          rsid,
            "trait":         meta["trait"],
            "icon":          meta["icon"],
            "category":      meta["category"],
            "gene":          meta.get("gene", ""),
            "genotype":      genotype,
            "result":        label,
            "value_key":     value_key,
            "what_it_means": meta.get("what_it_means", ""),
            "body_part":     meta.get("body_part", ""),
        })

    return {"traits": results, "available": True}
