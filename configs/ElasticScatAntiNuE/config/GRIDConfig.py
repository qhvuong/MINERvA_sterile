import argparse
import sys

PLAYLISTS = {
    "FHC": [
        "me1A_nx",
        "me1B_nx",
        "me1C_nx",
        "me1D_nx",
        "me1E_nx",
        "me1F_nx",
        "me1G_nx",
        "me1L_nx",
        "me1M_nx",
        "me1N_nx",
        "me1O_nx",
        "me1P_nx",
    ],

    "RHC": [
        "Minerva5",
    ],

    "FHC_NCDIF": [
        "me1A_NCDIF",
        "me1B_NCDIF",
        "me1C_NCDIF",
        "me1D_NCDIF",
    ],
    "FHC_Extended_2p2h":[
        "me1A_ext_2p2h"
    ]
}

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--RHC", "--rhc",
                    dest="helicity",
                    action="store_const",
                    const="RHC",
                    help="Use reverse horn current mode (RHC) instead of forward horn current (FHC).  FHC is the default.  If you supply both --FHC and --RHC, the trailing argument is the one that takes precedence.", 
                    default="FHC"
)

parser.add_argument("--FHC", "--fhc",
                    dest="helicity",
                    action="store_const",
                    const="FHC",
                    help="Use forward horn current (FHC) mode.  This is on by default; the option is supported in case you wish to be explicit.  If you supply both --FHC and --RHC, the trailing argument is the one that takes precedence."
)

parser.add_argument("-p", "--playlists",
                    dest="playlists",
                    nargs='*',
                    help="Only process given playlists. Defalut(None) is: for FHC: %s, and for RHC: %s" % (PLAYLISTS["FHC"],PLAYLISTS["RHC"])
)

parser.add_argument("--NCDIF", "--ncdif",
                    dest="special_sample",
                    action="store_const",
                    const="_NCDIF",
                    help="Use NCDIF special sample",
                    default=""
)

parser.add_argument("--EXT_2P2H", "--ext_2p2h",
                    dest="special_sample",
                    action="store_const",
                    const="_Extended_2p2h",
                    help="Use extended 2p2h special sample",
                    default=""
)

parser.add_argument("--data_only", "--data-only",
                    dest="data_only",
                    action="store_true",
                    help="Shortcut option to process only data from the 'data_types' option.  If you supply both '--data_only' and '--mc_only', '--mc-only' takes precedence.",
                    default=False,
)

parser.add_argument("--mc_only", "--mc-only",
                    dest="mc_only",
                    action="store_true",
                    help="Shortcut option to process only MC from the 'data_types' option.  If you supply both '--data_only' and '--mc_only', '--mc-only' takes precedence.",
                    default=False,
)
parser.add_argument("--ntuple_tag",
                    help="Use ntuple playlist tagged by given tag.",
                    default="LLR"
)

parser.add_argument("--count",
                    help="process n subruns per job",
                    type = int,
                    default = None)

parser.add_argument("--start",
                    help="starting from n-th entry of playlist file",
                    type =int,
                    default = 0)
parser.add_argument("--end",
                    help="ending at n-th entry of playlist file, Defalut is None, which means the end of playlist file, conflict with --njobs option",
                    type = int)

parser.add_argument("--njobs",
                    help="submit n jobs, conflict with --end option",
                    type = int)

parser.add_argument("--scratch",
                    help="save to scratch area",
                    dest="PNFS_switch",
                    action="store_const",
                    const="scratch",
                    default="persistent")

parser.add_argument("--cal_POT",
                    help = "count POT rather than submit jobs",
                    action="store_true",
                    default=False)

parser.add_argument("--memory","-m",
                    type = int,
                    help="amount of memory requested for this job, unit: MB",
                    dest="memory",
                    default=None)

parser.add_argument("--tarball",
                    help = "reuse tarball in given path",
                    default=None)

gridargs, anaargs = parser.parse_known_args()

if gridargs.playlists is None:
    gridargs.playlists = PLAYLISTS[gridargs.helicity+gridargs.special_sample]

if gridargs.njobs is not None and gridargs.end is not None:
    print("cann't fulfill --end and --njobs at the same time.")
    sys.exit(1)

if gridargs.data_only and gridargs.mc_only:
    print("You want data only AND mc only?")
    print("I'll process both.")
    gridargs.data_only=False
    gridargs.mc_only=False

if len(gridargs.special_sample) != 0 :
    gridargs.mc_only = True
