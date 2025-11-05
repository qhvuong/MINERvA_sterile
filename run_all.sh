#!/usr/bin/env bash
set -euo pipefail

# Helper to auto-press Enter at any prompt
run() { echo | python combine_file.py "$@"; }

run --playlist le1_p6 --i /pnfs/minerva/persistent/users/qvuong/CCNUE_selection_2025-10-22-100517_hists/ --cal_POT --ntuple_tag MAD --selection_tag elastic_noETheta2 --mc_only
run --playlist le13C_p6 --i /pnfs/minerva/persistent/users/qvuong/CCNUE_selection_2025-10-22-100751_hists/ --cal_POT --ntuple_tag MAD --selection_tag elastic_nodEdX --mc_only
run --playlist le13C_p6 --i /pnfs/minerva/persistent/users/qvuong/CCNUE_selection_2025-10-22-101840_hists/ --cal_POT --ntuple_tag MAD --selection_tag elastic_allCuts --mc_only
