#!/bin/bash
set -e
# Combining log: /exp/minerva/data/users/qvuong/runningNotes/Log_CCnue_FullCuts_v4p3p0_2026-11-02_005952.txt
# Submission log: /exp/minerva/data/users/qvuong/runningNotes/CCnue_FullCuts_v4p3p0_2026-11-02_005952.txt
# Exact original tarball: /pnfs/minerva/resilient/tarballs/qvuong-CCNUE_selection_CCnue_FullCuts_v4p3p0_2026-11-02_005952.tar.gz
# Commands must be run in CC-NuE-XSec/selection on GPVM.
# New grid outputs must be validated before replacing failed originals.
ls -lh /pnfs/minerva/resilient/tarballs/qvuong-CCNUE_selection_CCnue_FullCuts_v4p3p0_2026-11-02_005952.tar.gz

# le1_p8 mc: missing []; broken [43]; original count=1; original start=0
python selection/gridSelection.py --playlists le1_p8 --ntuple_tag MAD --use-sideband dEdX --truth --cal_POT --selection_tag CCnue_FullCuts_v4p3p0 --mc_only --start 43 --count 1 --njobs 1 --tarball /pnfs/minerva/resilient/tarballs/qvuong-CCNUE_selection_CCnue_FullCuts_v4p3p0_2026-11-02_005952.tar.gz

# le7_p8 mc: missing []; broken []; original count=1; original start=0

# le9_p8 mc: missing []; broken [2]; original count=1; original start=0
python selection/gridSelection.py --playlists le9_p8 --ntuple_tag MAD --use-sideband dEdX --truth --cal_POT --selection_tag CCnue_FullCuts_v4p3p0 --mc_only --start 2 --count 1 --njobs 1 --tarball /pnfs/minerva/resilient/tarballs/qvuong-CCNUE_selection_CCnue_FullCuts_v4p3p0_2026-11-02_005952.tar.gz

# le13C_p8 mc: missing []; broken [13, 18]; original count=1; original start=0
python selection/gridSelection.py --playlists le13C_p8 --ntuple_tag MAD --use-sideband dEdX --truth --cal_POT --selection_tag CCnue_FullCuts_v4p3p0 --mc_only --start 13 --count 1 --njobs 1 --tarball /pnfs/minerva/resilient/tarballs/qvuong-CCNUE_selection_CCnue_FullCuts_v4p3p0_2026-11-02_005952.tar.gz
python selection/gridSelection.py --playlists le13C_p8 --ntuple_tag MAD --use-sideband dEdX --truth --cal_POT --selection_tag CCnue_FullCuts_v4p3p0 --mc_only --start 18 --count 1 --njobs 1 --tarball /pnfs/minerva/resilient/tarballs/qvuong-CCNUE_selection_CCnue_FullCuts_v4p3p0_2026-11-02_005952.tar.gz

# le13C_2p2h_p8 mc: missing []; broken []; original count=1; original start=0

# le1_p8 data: missing []; broken []; original count=20; original start=0

# le7_p8 data: missing []; broken []; original count=20; original start=0

# le9_p8 data: missing []; broken []; original count=20; original start=0

# le13A_p8 data: missing []; broken []; original count=20; original start=0

# le13B_p8 data: missing []; broken []; original count=20; original start=0

# le13C_p8 data: missing []; broken []; original count=20; original start=0

# le13D_p8 data: missing []; broken []; original count=20; original start=0

# le13E_p8 data: missing []; broken []; original count=20; original start=0

# Total: 4 rerun commands
