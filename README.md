# HealthAgent — Genomic Wellness Platform

A self-hosted, privacy-first DNA analysis platform. Upload a raw DNA file from 23andMe, AncestryDNA, or similar services and receive personalised genomic wellness insights powered by open-source scientific databases.

---

## Features

| Module | Description |
|---|---|
| 🧬 **DNA Upload** | Chunked, client-side validated upload for 23andMe, AncestryDNA, MyHeritage, FTDNA, VCF files |
| 🥗 **Health Insights** | Wellness traits across Nutrition, Fitness, Sleep, Supplements, Allergies, Health History |
| 💊 **Medications** | Pharmacogenomics lookup with US↔UK drug name aliasing and DNA interaction summary |
| 👤 **Profile** | Physical trait predictions (eye colour, hair, skin) from DNA |
| 🔬 **DNA Map** | Interactive chromosome scanner with streaming SNP data |
| 🦠 **Gut Health** | HUMAnN 3.x microbiome pathway analysis |
| 🔭 **Science View** | Raw GWAS, ClinVar, DisGeNET, OpenTargets, Ensembl VEP, FinnGen data |
| 🔒 **Privacy** | All data stays on your machine — no third-party cloud upload |

---

## Databases

All databases download automatically on first run and refresh every 24 hours in the background:

| Database | Data | Source |
|---|---|---|
| GWAS Catalog | Trait–SNP associations | EBI / NHGRI |
| ClinVar | Clinical significance | NCBI |
| PharmGKB | Pharmacogenomics / drug interactions | Stanford |
| DisGeNET | Gene–disease associations | disgenet.org |
| OpenTargets | Target–disease evidence scores | opentargets.org |
| Ensembl VEP | Variant functional consequences | Ensembl REST API |
| FinnGen / UK Biobank | Biobank phenome associations | FinnGen.fi |
| TCGA | Cancer mutation context | NCI GDC |
| Wellness Traits | Curated SNP wellness annotations | Built-in |

---

## Installation

### Requirements

- Python 3.10+
- pip

```bash
git clone https://github.com/timehub-uk/healthagent.git
cd healthagent
pip install -e ".[dev]"
```

### Run

```bash
python run_ui.py
```

Open **http://localhost:8000** in your browser.

On first start the server will:
1. Initialise the local SQLite database
2. Seed built-in wellness trait annotations
3. Restore any previously uploaded DNA profile
4. Trigger a background database refresh if data is stale (> 24 h)

---

## Usage

1. **Upload your DNA file** — click "⬆ Upload DNA" and select your raw file
2. **View Health Insights** — trait cards appear across Nutrition, Fitness, Sleep, and more
3. **Check Medications** — add a drug name and see how your DNA affects it (PharmGKB + US↔UK aliases)
4. **Explore Science** — switch to Science View for GWAS, ClinVar, DisGeNET, and Ensembl data
5. **Clear your data** — click "🗑 Clear Data" to permanently delete everything from this device

### Supported DNA formats

- 23andMe (v3, v4, v5)
- AncestryDNA
- MyHeritage
- Family Tree DNA (FTDNA)
- VCF (Variant Call Format)

---

## Architecture

```
healthagent/
├── ui/
│   ├── app.py              FastAPI server — all HTTP routes
│   └── templates/
│       └── index.html      Single-page application (all JS inline)
├── databases/
│   ├── local_db.py         Thread-safe SQLite connection pool (threading.local)
│   ├── downloader.py       Background DB download tasks
│   ├── schema.sql          SQLite schema (WAL mode)
│   └── tcga_client.py      TCGA/GDC API client
├── dna_importer.py         DNA file parser (all formats)
├── health_traits.py        Parallel profile analysis + LRU cache
└── microbiome_importer.py  HUMAnN 3.x pathway parser
```

### Backend

- **FastAPI** — async HTTP server with StreamingResponse for chromosome scanning
- **SQLite + WAL** — local database; `threading.local` gives each thread its own connection
- **ThreadPoolExecutor** — genomic DB queries run in parallel (Phase 1: 8 workers; Phase 2: 3 workers)
- **Analysis cache** — results keyed by profile MD5 hash; invalidated automatically on new upload
- **Background DB updates** — parallel downloader (4 workers), triggered at startup if stale

### Frontend

- **Vanilla JS + AJAX** — no framework; NDJSON streaming for real-time chromosome scanning
- **Three.js r128** — 3D rotating DNA double helix on the landing page
- **Web Crypto API** — per-device UUID + PBKDF2 → AES-GCM 256-bit key for local data encryption
- Medical clinical dark theme (navy/sky-blue)

### Security

- Client-side file validation — rejects non-DNA content and known abuse patterns
- Server-side: UUID path-traversal guard, chunk size/count limits, minimum SNP count check
- `Cache-Control: no-store` on the HTML route — prevents stale UI being served from browser cache
- AES-GCM 256-bit encryption for any locally stored sensitive data (per-device UUID key)

---

## Deployment on Plesk

See [`deploy/PLESK_INSTALL.md`](deploy/PLESK_INSTALL.md) for full Phusion Passenger WSGI instructions.

Quick steps:
1. Upload via Plesk File Manager or `git clone`
2. In Plesk → Python → create venv and `pip install -e .`
3. Set Passenger entry point to `deploy/passenger_wsgi.py`
4. Point document root to `public/`

---

## Privacy

- Your DNA file is **never sent to any external server**
- All analysis runs locally on the machine hosting HealthAgent
- Data is stored in `data/healthagent.db` (SQLite) on the server
- Use **Clear Data** to erase your profile at any time
- Databases download from public open-source APIs (GWAS Catalog, NCBI, etc.) — no user data is transmitted

---

## Licence

MIT — see [LICENSE](LICENSE)
