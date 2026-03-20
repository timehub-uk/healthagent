"""FastAPI web server for the HealthAgent DNA Wellness Dashboard."""

import asyncio
import logging
import time
from pathlib import Path
from typing import Optional

from fastapi import FastAPI, File, UploadFile, HTTPException, BackgroundTasks
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from fastapi.requests import Request

from healthagent.dna_importer import import_dna_string, DNAProfile, SNP, DNAFormat
from healthagent.databases import local_db
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

log = logging.getLogger(__name__)

app = FastAPI(title="HealthAgent DNA Wellness Dashboard")

BASE = Path(__file__).parent
templates = Jinja2Templates(directory=str(BASE / "templates"))
app.mount("/static", StaticFiles(directory=str(BASE / "static")), name="static")

# ── Session state ─────────────────────────────────────────────────
_profile: Optional[DNAProfile] = None

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
    global _db_initialized
    local_db.init_db()
    # Seed curated wellness traits (fast, always safe to re-run)
    seed_wellness_traits(progress_cb=log.info)
    _db_initialized = True
    log.info("HealthAgent DB initialised.")

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

    await _maybe_trigger_update(background_tasks)

    return {
        "format":    _profile.source_format.value,
        "snp_count": _profile.snp_count,
        "filename":  file.filename,
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
