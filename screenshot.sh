#!/usr/bin/env bash
# screenshot.sh — take a screenshot of the HealthAgent UI and print the file path
# Usage: ./screenshot.sh [url] [output.png]
#   url defaults to http://localhost:8000/
#   output defaults to /tmp/healthagent_$(date +%s).png

URL="${1:-http://localhost:8000/}"
OUT="${2:-/tmp/healthagent_$(date +%s).png}"
W="${3:-1440}"
H="${4:-900}"

google-chrome \
  --headless=new \
  --screenshot="$OUT" \
  --window-size="${W},${H}" \
  --no-sandbox \
  --disable-gpu \
  --hide-scrollbars \
  "$URL" 2>/dev/null

if [ -f "$OUT" ]; then
  echo "$OUT"
else
  echo "ERROR: screenshot failed" >&2
  exit 1
fi
