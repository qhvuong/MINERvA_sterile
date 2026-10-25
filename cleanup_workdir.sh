#!/usr/bin/env bash

set -e

BASE="/exp/minerva/app/users/qvuong/MAT_AL9/CC-NuE-XSec"

cd "$BASE"

echo "Cleaning $BASE"

echo
echo "Clearing top-level grid_wrappers contents..."
find grid_wrappers -mindepth 1 -maxdepth 1 -exec rm -rf {} + 2>/dev/null || true

echo
echo "Clearing oscillations/grid_wrappers contents..."
find oscillations/grid_wrappers -mindepth 1 -maxdepth 1 -exec rm -rf {} + 2>/dev/null || true

echo
echo "Removing __pycache__ directories..."
find "$BASE" -type d -name "__pycache__" -prune -exec rm -rf {} +

echo
echo "Removing Python bytecode..."
find "$BASE" -type f \( -name "*.pyc" -o -name "*.pyo" \) -delete

echo
echo "Removing editor/temp files..."
find "$BASE" -type f \( \
    -name "*~" \
    -o -name "*.swp" \
    -o -name "*.swo" \
    -o -name ".DS_Store" \
\) -delete

echo
echo "Cleanup complete."

echo
echo "Current repository size:"
du -sh "$BASE"

echo
echo "Largest remaining top-level items:"
du -sh * .[^.]* 2>/dev/null | sort -h | tail -15