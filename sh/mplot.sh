#!/usr/bin/env bash

# TODO: refactor & consolidate into one script to run sex differintiation & sex-agnostic analysis
BASE="dtw"
META="sigma/meta.csv"
OUTBASE="png"

# Create base output directory if missing
mkdir -p "$OUTBASE"

for DIR in "$BASE"/*/; do
    echo "Processing $DIR"

    # Compute relative subpath
    SUBDIR=$(basename "$DIR")

    # Create matching output directory
    OUTDIR="$OUTBASE/$SUBDIR"
    mkdir -p "$OUTDIR"

    # Loop over cost matrices
    for COST in "$DIR"/*_cost_matrix.csv; do
        IDS="${COST/_cost_matrix.csv/_order_ids.csv}"
        NAME=$(basename "$COST" _cost_matrix.csv)

        OUTPNG="$OUTDIR/${NAME}.png"

        echo "  → Plotting $NAME"
        julia --project src/plot_costmatrix.jl \
          --cost "$COST" \
          --ids "$IDS" \
          --meta "$META" \
          --out "$OUTPNG" \
          --sex "" \
          --title "$NAME"
            
            echo "$COST"
            echo "$IDS"
            echo "$META"
            echo "$OUTPNG"
            echo "$NAME"
            echo 
    done
done

