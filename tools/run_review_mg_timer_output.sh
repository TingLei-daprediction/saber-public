#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 3 ]]; then
  echo "Usage: run_review_mg_timer_output.sh <mg_timer_output_file> [output_dir] [top_n]"
  exit 1
fi

TIMER_FILE="$1"
OUTPUT_DIR="${2:-}"
TOP_N="${3:-25}"

CMD=(
  python training/docs/user-guide/examples/review_mg_timer_output.py
  "$TIMER_FILE"
  --top-n "$TOP_N"
)

if [[ -n "$OUTPUT_DIR" ]]; then
  CMD+=(--output-dir "$OUTPUT_DIR")
fi

"${CMD[@]}"
