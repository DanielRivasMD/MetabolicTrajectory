#!/usr/bin/env bash

DTW_ROOT="dtw"
PNG_ROOT="png"
META="sigma/meta.csv"
JULIA_SCRIPT="src/plot_costmatrix.jl"

echo "Mirroring dtw → png and running Julia plots"
echo

# Ensure png root exists
mkdir -p "$PNG_ROOT"

# Process both raw/ and diff/
for SUB in raw diff; do
  echo "Processing $SUB/ ..."

  # Find all subdirectories like n100_len80.../
  find "$DTW_ROOT/$SUB" -mindepth 1 -maxdepth 1 -type d | while read -r DIR; do
    BASENAME=$(basename "$DIR")

    # Mirror directory under png/
    OUTDIR="$PNG_ROOT/$SUB/$BASENAME"
    mkdir -p "$OUTDIR"

    echo "  → Mirrored: $OUTDIR"

    # Loop over all cost matrices inside this directory
    for COST in "$DIR"/*_cost_matrix.csv; do
      [ -e "$COST" ] || continue

      IDS="${COST/_cost_matrix.csv/_order_ids.csv}"
      NAME=$(basename "$COST" _cost_matrix.csv)
      OUTPNG="$OUTDIR/${NAME}.png"

      echo "    Plotting $NAME"
      echo "      cost: $COST"
      echo "      ids:  $IDS"
      echo "      out:  $OUTPNG"

      julia --project "$JULIA_SCRIPT" \
        --cost "$COST" \
        --ids "$IDS" \
        --meta "$META" \
        --out "$OUTPNG" \
        --sex "" \
        --title "$NAME"
    done

    echo
  done
done

echo "All PNGs generated under: $PNG_ROOT/"
