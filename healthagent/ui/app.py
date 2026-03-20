"""FastAPI web server for the HealthAgent DNA visualiser."""

import os
from pathlib import Path
from typing import Optional

from fastapi import FastAPI, File, UploadFile, HTTPException
from fastapi.responses import HTMLResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from fastapi.requests import Request

from healthagent.dna_importer import import_dna_string, DNAProfile

app = FastAPI(title="HealthAgent DNA Visualiser")

BASE = Path(__file__).parent
templates = Jinja2Templates(directory=str(BASE / "templates"))
app.mount("/static", StaticFiles(directory=str(BASE / "static")), name="static")

# In-memory store for the current session profile
_profile: Optional[DNAProfile] = None


@app.get("/", response_class=HTMLResponse)
async def index(request: Request):
    return templates.TemplateResponse("index.html", {"request": request})


@app.post("/api/upload")
async def upload_dna(file: UploadFile = File(...)):
    """Accept a raw DNA file upload and parse it into the session profile."""
    global _profile
    try:
        raw = (await file.read()).decode("utf-8", errors="replace")
        _profile = import_dna_string(raw, filename=file.filename)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc))

    return {
        "format": _profile.source_format.value,
        "snp_count": _profile.snp_count,
        "filename": file.filename,
    }


@app.get("/api/snps")
async def get_snps(
    chromosome: Optional[str] = None,
    limit: int = 500,
    offset: int = 0,
):
    """Return a page of SNP records, optionally filtered by chromosome."""
    if _profile is None:
        raise HTTPException(status_code=404, detail="No DNA profile loaded. Upload a file first.")

    snps = _profile.snps
    if chromosome:
        snps = [s for s in snps if s.chromosome == chromosome.upper().lstrip("CHR")]

    page = snps[offset : offset + limit]
    return {
        "total": len(snps),
        "offset": offset,
        "limit": limit,
        "snps": [
            {
                "rsid": s.rsid,
                "chromosome": s.chromosome,
                "position": s.position,
                "genotype": s.genotype,
            }
            for s in page
        ],
    }


@app.get("/api/chromosomes")
async def get_chromosomes():
    """Return the list of chromosomes present in the loaded profile."""
    if _profile is None:
        raise HTTPException(status_code=404, detail="No DNA profile loaded.")
    chroms = sorted({s.chromosome for s in _profile.snps},
                    key=lambda c: (c.isdigit() is False, c.zfill(3)))
    return {"chromosomes": chroms}


@app.get("/api/demo")
async def load_demo():
    """Load a synthetic demo profile so the UI works without a real file."""
    global _profile
    import random, string
    from healthagent.dna_importer import SNP, DNAProfile, DNAFormat

    random.seed(42)
    bases = "AGTC"
    chroms = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
    snps = []
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
    _profile = DNAProfile(source_format=DNAFormat.TWENTYTHREE_AND_ME, snps=snps,
                          metadata={"filename": "demo"})
    return {"format": "23andme", "snp_count": len(snps), "filename": "demo"}
