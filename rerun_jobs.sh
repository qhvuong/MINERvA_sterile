#!/bin/bash
set -e
# Combining log: /exp/minerva/data/users/qvuong/runningNotes/Log_CCnumubar_FullCuts_v4p3p0_2026-11-02_012339.txt
# Submission log: /exp/minerva/data/users/qvuong/runningNotes/CCnumubar_FullCuts_v4p3p0_2026-11-02_012339.txt
# Exact original tarball: /pnfs/minerva/resilient/tarballs/qvuong-CCNUE_selection_CCnumubar_FullCuts_v4p3p0_2026-11-02_012339.tar.gz
# Commands must be run in CC-NuE-XSec/selection on GPVM.
# New grid outputs must be validated before replacing failed originals.
ls -lh /pnfs/minerva/resilient/tarballs/qvuong-CCNUE_selection_CCnumubar_FullCuts_v4p3p0_2026-11-02_012339.tar.gz

# le5_p8 mc: missing []; broken [16]; original count=1; original start=0
python selection/gridSelection.py --playlists le5_p8 --ntuple_tag MAD --truth --cal_POT --selection_tag CCnumubar_FullCuts_v4p3p0 --mc_only --start 16 --count 1 --njobs 1 --tarball /pnfs/minerva/resilient/tarballs/qvuong-CCNUE_selection_CCnumubar_FullCuts_v4p3p0_2026-11-02_012339.tar.gz

# le5_2p2h_p8 mc: missing []; broken []; original count=1; original start=0

# le5_p8 data: missing []; broken []; original count=20; original start=0

# Total: 1 rerun commands
