#!/bin/bash
set -e
# Combining log: /exp/minerva/data/users/qvuong/runningNotes/Log_CCnueswap_FullCuts_v4p3p0_fixedFluxNorm_2026-11-04_160121.txt
# Submission log: /exp/minerva/data/users/qvuong/runningNotes/CCnueswap_FullCuts_v4p3p0_fixedFluxNorm_2026-11-04_160121.txt
# Exact original tarball: /pnfs/minerva/resilient/tarballs/qvuong-CCNUE_selection_CCnueswap_FullCuts_v4p3p0_fixedFluxNorm_2026-11-04_160121.tar.gz
# Commands must be run in CC-NuE-XSec/selection on GPVM.
# New grid outputs must be validated before replacing failed originals.
ls -lh /pnfs/minerva/resilient/tarballs/qvuong-CCNUE_selection_CCnueswap_FullCuts_v4p3p0_fixedFluxNorm_2026-11-04_160121.tar.gz

# le1_p8 mc: missing []; broken []; original count=1; original start=0

# le7_p8 mc: missing []; broken []; original count=1; original start=0

# le9_p8 mc: missing [2]; broken []; original count=1; original start=0
python selection/gridSelection.py --playlists le9_p8 --ntuple_tag MAD --truth --cal_POT --selection_tag CCnueswap_FullCuts_v4p3p0_fixedFluxNorm --exclude-universes all --mc_only --start 2 --count 1 --njobs 1 --tarball /pnfs/minerva/resilient/tarballs/qvuong-CCNUE_selection_CCnueswap_FullCuts_v4p3p0_fixedFluxNorm_2026-11-04_160121.tar.gz

# le13C_p8 mc: missing []; broken [18]; original count=1; original start=0
python selection/gridSelection.py --playlists le13C_p8 --ntuple_tag MAD --truth --cal_POT --selection_tag CCnueswap_FullCuts_v4p3p0_fixedFluxNorm --exclude-universes all --mc_only --start 18 --count 1 --njobs 1 --tarball /pnfs/minerva/resilient/tarballs/qvuong-CCNUE_selection_CCnueswap_FullCuts_v4p3p0_fixedFluxNorm_2026-11-04_160121.tar.gz

# le13C_2p2h_p8 mc: missing []; broken []; original count=1; original start=0

# Total: 2 rerun commands
