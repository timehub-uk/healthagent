# Deploying HealthAgent on Plesk (Phusion Passenger)

## Prerequisites

- Plesk Obsidian 18.x or later
- Python 3.10+ extension installed in Plesk
- A domain or subdomain set up in Plesk

---

## Step 1 — Upload the code

**Option A — Git**
```
cd /var/www/vhosts/yourdomain.com/httpdocs
git clone https://github.com/timehub-uk/healthagent.git .
```

**Option B — Plesk File Manager**
Upload the ZIP and extract to `/var/www/vhosts/yourdomain.com/httpdocs/`

---

## Step 2 — Create the Python virtual environment

In Plesk → **Domains** → your domain → **Python**:

| Setting | Value |
|---|---|
| Python version | 3.10 (or latest) |
| Application root | `/var/www/vhosts/yourdomain.com/httpdocs` |
| Application startup file | `deploy/passenger_wsgi.py` |
| Application URL | `/` |

Then click **Install requirements** (reads `requirements.txt` / `pyproject.toml`).

Or via SSH:
```bash
cd /var/www/vhosts/yourdomain.com/httpdocs
python3 -m venv venv
source venv/bin/activate
pip install -e ".[dev]"
```

---

## Step 3 — Set the document root

In Plesk → **Domains** → your domain → **Hosting Settings**:

- Document root: `public/`

The `public/` directory contains a placeholder `index.html` that Passenger replaces with the app.

---

## Step 4 — Create the data directory

```bash
mkdir -p /var/www/vhosts/yourdomain.com/httpdocs/data
chmod 755 /var/www/vhosts/yourdomain.com/httpdocs/data
```

The SQLite database (`healthagent.db`) will be created here automatically on first run.

---

## Step 5 — Restart the application

In Plesk → **Python** → click **Restart**.

Or via SSH:
```bash
touch /var/www/vhosts/yourdomain.com/httpdocs/tmp/restart.txt
```

---

## Step 6 — Verify

Open `https://yourdomain.com` — you should see the HealthAgent landing page.

Logs are in:
```
/var/www/vhosts/yourdomain.com/logs/error_log
```

---

## Environment variables

Set these in Plesk → **Python** → **Environment variables** (or in a `.env` file):

| Variable | Default | Description |
|---|---|---|
| `HEALTHAGENT_ENV` | `production` | Environment tag |
| `PORT` | `8000` | Port (Passenger manages this automatically) |

---

## Updating

```bash
cd /var/www/vhosts/yourdomain.com/httpdocs
git pull
source venv/bin/activate
pip install -e .
touch tmp/restart.txt
```

---

## Troubleshooting

**App shows "500 Internal Server Error"**
→ Check `/var/www/vhosts/yourdomain.com/logs/error_log`

**Database not found**
→ Ensure the `data/` directory exists and is writable by the Passenger user

**Passenger can't find the startup file**
→ Confirm the startup file path is exactly `deploy/passenger_wsgi.py` (relative to application root)

**Static files not serving**
→ Ensure Plesk's Apache/Nginx static file rule covers `/static/` → the app's `healthagent/ui/static/` directory
