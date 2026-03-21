"""Phusion Passenger WSGI entry point for Plesk hosting.

Plesk → Python App → Startup file: deploy/passenger_wsgi.py
Document root: public/
"""

import sys
import os

# ── Project root on the path ──────────────────────────────────────
# Adjust this if your Plesk home differs
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)

# ── Environment ───────────────────────────────────────────────────
os.environ.setdefault("HEALTHAGENT_ENV", "production")

# ── Import the FastAPI app and wrap it for WSGI ───────────────────
from healthagent.ui.app import app as _fastapi_app

# Passenger expects a WSGI callable named `application`
try:
    from asgiref.wsgi import WsgiToAsgi  # noqa: F401
    # Use uvicorn's ASGI adapter instead — Passenger supports ASGI via uvicorn
    raise ImportError("prefer uvicorn")
except ImportError:
    pass

# Recommended: use uvicorn as the ASGI server inside Passenger
# Add to requirements: uvicorn[standard]
from uvicorn.middleware.wsgi import WSGIMiddleware  # type: ignore
application = WSGIMiddleware(_fastapi_app)
