#!/usr/bin/env python3
"""Launch the HealthAgent DNA Visualiser web UI."""

import uvicorn

if __name__ == "__main__":
    uvicorn.run(
        "healthagent.ui.app:app",
        host="0.0.0.0",
        port=8000,
        reload=True,
    )
