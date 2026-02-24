#!/usr/bin/env bash

SRC="$1"
DST="$2"

if [ -z "$SRC" ] || [ -z "$DST" ]; then
  echo "Usage: mirror.sh <source_dir> <destination_dir>"
  exit 1
fi

# Ensure destination exists
mkdir -p "$DST"

# Walk source and recreate directories
find "$SRC" -type d | while read -r dir; do
  # Compute relative path
  rel="${dir#$SRC/}"
  # Create corresponding directory in destination
  mkdir -p "$DST/$rel"
done
