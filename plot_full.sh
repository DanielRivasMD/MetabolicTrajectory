#!/usr/bin/env bash
set -euo pipefail

INPUT_DIR="dtw/full"
OUTPUT_DIR="png/full"
JULIA_SCRIPT="plot_full.jl"
META="sigma/meta.csv"

mkdir -p "$OUTPUT_DIR"

for CSV in "$INPUT_DIR"/RT*.csv; do
  [ -e "$CSV" ] || continue

  BASENAME=$(basename "$CSV" .csv)
  OUTPNG="$OUTPUT_DIR/${BASENAME}.png"

  echo "Processing $CSV → $OUTPNG"

  julia --project "$JULIA_SCRIPT" \
    --csv "$CSV" \
    --meta "$META" \
    --out "$OUTPNG" \
    --pad 10 \
    --title "$BASENAME"
done
