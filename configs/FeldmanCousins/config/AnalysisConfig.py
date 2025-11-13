"""
  AnalysisConfig.py:
   Centralized configuration for the nu_e CCQE analysis:
    signal definitions (nu vs. anti-nu),
    file locations,
    etc.
  
   Original author: J. Wolcott (jwolcott@fnal.gov)
                    January 2014
"""

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
#ROOT.gErrorIgnoreLevel = ROOT.kWarning
ROOT.gROOT.SetBatch()

import math
import argparse
import os
import sys
import pprint
import re

BLUEARC = "/exp/minerva/data/users/{}/surfaces/".format(os.environ["USER"])
GIRDOUTPUT ="/pnfs/minerva/persistent/"


SIDEBANDS=["dEdX"]
#SIDEBANDS=["Excess_High_Inline","Excess_Low_Inline","High_PsiEe_and_dEdX","Pi0","dEdX","ETheta2","True_Electron","True_Muon"]
#SIDEBANDS=["Excess_High_Inline","Excess_Low_Inline","High_Psi_and_dEdX"]

class _AnalysisConfig(object):
    Defaults = {
        "data_types": ["data", "mc",],
    }
    ALLOWED_HELICITIES = ["FHC", "RHC"]

    def __init__(self, **kwargs):
        params = _AnalysisConfig.Defaults.copy()
        params.update(kwargs)

        self._config_keys = set()
        for key, value in params.items():
            self._config_keys.add(key)
            setattr(self, key, value)

        retained_dts = []
        for dt in self.data_types:
            if (self.mc_only and "data" in dt.lower()) or (self.data_only and "mc" in dt.lower()):
                continue
            retained_dts.append(dt)

        self.data_types = retained_dts

    def __repr__(self):
        my_dict = dict([ (k, getattr(AnalysisConfig, k)) for k in AnalysisConfig._config_keys ])
        return pprint.pformat(my_dict)

    ####  properties ####
    #@property
    #def bknd_constraint_method(self):
    #    return BKND_CONSTRAINT_METHOD
    
    @property
    def helicity(self):
        return self._helicity

    @helicity.setter
    def helicity(self, value):
        value_upper = value.upper()
        if value_upper not in _AnalysisConfig.ALLOWED_HELICITIES:
            raise NameError("Allowed helicities are '%s', not '%s'" % (value, _AnalysisConfig.ALLOWED_HELICITIES))

        if "helicity" not in self._config_keys:
            self._config_keys.add("helicity")
        self._helicity = value_upper
    #@property
    #def POT(self):
    #    if self.processing_pass == "Resurrection":
    #        return POT
    #    else:
    #        raise Exception("Don't have POT data for processings other than Resurrection!")

    #@property
    #def POT_string(self, precision=2):
    #    exponent = int(math.log10(self.POT["data"]))
    #    mantissa = self.POT["data"] / 10**exponent
    #    fmt_string = "%%.%(precision)df #times 10^{%%d} P.O.T." % {"precision": precision}
    #    return fmt_string% (mantissa, exponent)

    # @property
    # def data_MC_POT_ratio(self):
    #     return self.POT["data"] / self.POT["MC"]

    # @property
    # def processing_pass(self):
    #     return self._processing_pass

    # @processing_pass.setter
    # def processing_pass(self, value):
    #     value_cap = value.lower().capitalize()
    #     if value_cap not in _AnalysisConfig.ALLOWED_PROCESSING_PASSES:
    #         raise NameError("Allowed processing passes are '%s', not '%s'" % (value, _AnalysisConfig.ALLOWED_PROCESSING_PASSES))
    #     if "processing_pass" not in self._config_keys:
    #         self._config_keys.add("processing_pass")
    #     self._processing_pass = value_cap

    @property
    def right_sign_electron_pdg(self):
        """ Sign of the PDG code for "right sign" in this configuration """
        return 1 if self.helicity == "FHC" else -1 # electron is +11; positron is -11



    #### public interface ####
    def DataTypeNtupleList(self):
        filename_config = {
            "proc_pass": self.processing_pass,
            "helicity": self.helicity,
        }
        for dt in self.data_types:
            filename_config.update( { "data_type": dt } )
            filename = NTUPLE_FILENAME_FORMAT % filename_config
#           print "  searching for filename:", filename
            fullpaths = []
#           print "looking here:", NTUPLE_PATH
            for root, dirs, files in os.walk(NTUPLE_PATH):
#               print "considering filenames:", files
                if filename in files:
                    fullpaths.append(os.path.join(root, filename))
            if len(fullpaths) == 1:
#               print "Found ntuple list for specification '%s':" % dt, fullpaths[-1]
                yield dt, fullpaths[-1]
            elif len(fullpaths) > 1:
                raise Exception("Multiple matches for data type specification '%s': %s" % (dt, fullpaths))
            else:
                continue

    def FilePath(self, top_dir,tier_tag, playlist, data_tag, type_tag):
        return "{}/{}_{}{}_{}".format(top_dir,tier_tag, data_tag,
                                       playlist, type_tag)
#        print signal_defn
#         signal_defn = signal_defn if signal_defn is not None else self.signal_defn
#         signal_defn = SIGNAL_DEFN_TEXT[signal_defn]
#        print signal_defn

#         path = os.path.join(TOP_DIR, signal_defn, dirname)
#         if subdir:
#             path = os.path.join(path, subdir)
#         if helicity:
#             path = os.path.join(path, self.helicity)
#         return path

    def SelectionHistoPath(self, playlist, is_data, is_output = True):
        return self.FilePath(self.output_dir if is_output else self.input_dir,
                             "kin_dist"+("test" if self.testing else ""), playlist,
                             "data" if is_data else "mc",
                             self.selection_tag+"_"+self.ntuple_tag+
                             ("_"+str(self.count[0]) if self.count[0] is not None else "")+
                             ".root")

    def CutStudyPath(self, playlist, is_data, is_output = True):
        return self.FilePath(self.output_dir if is_output else self.input_dir,
                             "cut_study"+("test" if self.testing else ""), playlist,
                             "data" if is_data else "mc",
                             self.selection_tag+"_"+self.ntuple_tag+
                             ("_"+str(self.count[0]) if self.count[0] is not None else "")+
                             ".root")


    def TruthHistoPath(self,playlist,is_output = True):
        return self.FilePath(self.output_dir if is_output else self.input_dir,
                             "truth_dist"+("test" if self.testing else ""),playlist,
                             "mc", self.selection_tag+"_"+self.ntuple_tag+("_"+str(self.count[0]) if self.count[0] is not None else "")+
                             ".root")

    def UnfoldedHistoPath(self,playlist,tag,is_output=True):
        return self.FilePath(self.output_dir if is_output else self.input_dir,
                             "unfolded"+("test" if self.testing else ""),playlist,
                             "", self.selection_tag+"_"+self.ntuple_tag+"_"+tag+
                             ".root")

    def XSecHistoPath(self,playlist,is_output=True):
        return self.FilePath(self.output_dir if is_output else self.input_dir,
                             "xsec"+("test" if self.testing else ""),playlist,
                             "", self.selection_tag+"_"+self.ntuple_tag+
                             ".root")

    
    def BackgroundFitPath(self, playlist, tag, is_output = True):
        return self.FilePath(self.output_dir if is_output else self.input_dir,
                             "bkgfit", playlist, "" , tag+"_"+self.selection_tag+"_"+self.ntuple_tag+".root")
    
    def PlotPath(self, plot_name, sideband,tag=""):
        return self.FilePath(self.output_dir,"plot/"+plot_name, sideband, "" , self.selection_tag+"_"+tag)

#### entry point ####

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-p", "--playlist",
                    dest="playlist",
                    help="Process given playlists."
)

parser.add_argument("--grid",
                    action="store_true",
                    default = False,
                    help = "Run macro on grid, Input/Output path must be updated to avoid direct access to BlueArc"
)

parser.add_argument("-d", "--data_types",
                    dest="data_types",
                    action="append",
                    help="Data types to process.  Supply this once for every data type you want to use.",
                    default=argparse.SUPPRESS,
)

parser.add_argument("--use-sideband", "--use_sideband",
                    dest="sidebands",
                    nargs='*',
                    help="Use this sideband rather than defalut.  (If you use this option at all, you must use it to specify ALL the sidebands you want.) ",
                    default=SIDEBANDS,
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

parser.add_argument("--pc",
                    dest="is_pc",
                    action="store_true",
                    help="running over particle cannon sample",
                    default=False,
)

parser.add_argument("--signal",
                    dest="signal",
                    action="store_true",
                    help="Use only the extracted true signal event samples.",
                    default=argparse.SUPPRESS,
)

parser.add_argument("-t", "--test", "--testing",
                    dest="testing",
                    action="store_true",
                    help="Use a smaller sample for testing purposes. ",
                    default=False
)

parser.add_argument("-o", "--output",
                    dest = "output_dir",
                    help="Use alternate location for output file.",
                    default=BLUEARC
)

parser.add_argument("-i", "--input",
                    dest = "input_dir",
                    help="Use alternate location for input files other than ntuple.",
                    default=BLUEARC
)

parser.add_argument("--ntuple_tag", "--ntuple-tag",
                    help="Use ntuple playlist tagged by given tag.",
                    default="LLR"
)

parser.add_argument("--selection_tag","--selection-tag",
                    help="Use event selection histograms tagged by given tag.",
                    default="collab1"
)

parser.add_argument("--bkgTune_tag","--bkgTune-tag",
                    help="Use event selection histograms tagged by given tag.",
                    default="Global"
)

parser.add_argument("--count",
                    help="process arg1 subruns starting from arg0 entry of playlist.",
                    type = int,
                    nargs = 2,
                    default = [None,None])

parser.add_argument("--cal_POT","--cal-POT",
                    help="recount POT even if POT info is available",
                    dest = "POT_cal",
                    action="store_true",
                    default=False)

parser.add_argument("--only_cal_POT","--only-cal-POT",
                    help="do not run selection but only count POT",
                    dest = "run_reco",
                    action="store_false",
                    default=True)

parser.add_argument("--no-reco","--no_reco",
                    help="do not run selection but only count POT",
                    dest = "run_reco",
                    action="store_false",
                    default=True)

parser.add_argument("--exclude_universes","--exclude-universes",
                    help="do not want some systematics universes, list their ShortName()",
                    nargs="*",
)

parser.add_argument("--skip_2p2h","--skip_2p2h",
                    help="do not want 2p2h events,(use this when you are going to run delicate 2p2h sample.)",
                    action="store_true",
                    default=False)

parser.add_argument("--truth",
                    help="run truth loop: more correct efficiency demominator",
                    action="store_true",
                    default=False)

parser.add_argument("--animate",
                    help="make contour plots as a video",
                    action="store_true",
                    default=False)

parser.add_argument("--compare_profiles",
                    help="compare different methods of flux profiling",
                    action="store_true",
                    default=False)

parser.add_argument("--pseudodata",
                    help="Use MC predicted pseudodata in signal region.",
                    action="store_true",
                    dest = 'pseudodata',
                    default=False)
parser.add_argument("--blinded",
                        dest = "blinded",
                        default = False,
                        action="store_true",
)
parser.add_argument("--extra_weighter",
                    help="Name of extra weighter you want to use",
                    default=None
)
parser.add_argument("--bin_width",
                    dest = "binwidth",
                    default = False,
                    action="store_true",
)
parser.add_argument("--thesis_plots",
                    dest = "thesis_plots",
                    default = False,
                    action="store_true",
)
parser.add_argument("--exclude",
                    dest ="exclude",
                    default = "None",
)
parser.add_argument("--ratio",
                    dest ="ratio",
                    default = False,
                    action = "store_true",
)
parser.add_argument("--lambda",
                    dest ="lambdaValue",
                    default = 1,
                    type=float
)
parser.add_argument("-m", "--delta_m",
                    dest = "delta_m",
                    help="Delta m^2 value to probe.",
                    type=int,
                    default=0
)
parser.add_argument("-Ue4", "--U_e4",
                    dest = "U_e4",
                    help="U_e4 parameter to probe.",
                    type=int,
                    default=0
)
parser.add_argument("-Umu4", "--U_mu4",
                    dest = "U_mu4",
                    help="U_mu4 parameter to probe.",
                    type=float,
                    default=0
)
parser.add_argument("-U", "--U_tau4",
                    dest = "U_tau4",
                    help="U_tau4 parameter to probe.",
                    type=float,
                    default=0
)
parser.add_argument("-excludeSys", "--exclude_systematic",
                    dest = "exclude_systematic",
                    help="Systematic category to exclude",
                    type=str,
                    default=""
)

options = parser.parse_args()

if options.data_only and options.mc_only:
    options.data_only = False

options.exclude = options.exclude.lower()
#if options.playlist is None:
    #print "Please specify a playlist."
    #sys.exit(1)

# If only cal POT, doesn't make sense to read from record.
if not options.run_reco:
    options.POT_cal=True

#if options.grid:
    #override start variable by $PROCESS variable
    #nth_job = int(os.environ["PROCESS"])
    #options.count[0]=nth_job*options.count[1]+options.count[0]

if options.testing:
    options.count = [0,1]

AnalysisConfig = _AnalysisConfig(**vars(options))

print("Using analysis configuration:")
print(AnalysisConfig)

