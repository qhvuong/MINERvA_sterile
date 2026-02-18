from collections import OrderedDict
from config import PlotConfig
import ROOT


USE_NLL = True
#HIST_TO_FIT= "Lepton Energy"
#HIST_TO_FIT= "Q3"
#HIST_TO_FIT= "Lepton Energy"
HIST_TO_FIT= "Biased Neutrino Energy"
HIST_OBSERVABLE= {"name":["Biased Neutrino Energy"]}
#HIST_OBSERVABLE= {"variables":["PsiEe","Lepton Energy"]}

REGULATION_PARAMETER = 0.01
SCALE_FACTORS = OrderedDict()
#SCALE_FACTORS["Signal"] =  [0,0.4,0.6,0.8,1.0,1.2]
SCALE_FACTORS["DIS"] =  [0,20]
#SCALE_FACTORS["Pi0"] =  [0,2.5,4,6,9,12,15,20]
#SCALE_FACTORS["Pi0"] =  [0,2.0]
#SCALE_FACTORS["Pi0"] =  [0,1.2,2.0]
#SCALE_FACTORS["Pi0"] =  [0,1.6]
#SCALE_FACTORS["Pi0"] =  [0,0.4,0.6,0.8,1.0,1.2,1.6,2.0]
#SCALE_FACTORS["Pi0"] =  [0,0.4,0.6,0.8,1.0,1.2,1.6]
#SCALE_FACTORS["Excess"] =  [0,2.5,4,6,9,12,15,20]
#SCALE_FACTORS["NCCoh"] =  [0,2.5,4,6,9,12,15,20]
SCALE_FACTORS["NCCoh"] =  [0,2.5,5,7.5,10,12.5,15,20]
SCALE_FACTORS["Excess"] =  [0,2.5,5,7.5,10,12.5,15,20]
SCALE_FACTORS["Signal"] =  [0,20]
#SCALE_FACTORS["NCCoh"] =  [0,2.0]
#SCALE_FACTORS["Excess"] =  [0,2.0]
Yaxis =True

CATEGORY_FACTORS= {
   "CCNuE":"Signal",
   "CCNuEQE":"Signal",
   "CCNuEDIS":"Signal",
   "CCNuEDelta":"Signal",
   "CCNuE2p2h":"Signal",
   #"NonPhaseSpace":"Signal",
   #"CCNuEAntiNu":"Signal",
    "NCDIS":"DIS",
#    "NCMultiMeson": "Pi0",
    "CCDIS":"DIS",
    "NCCOH":"NCCoh",
    "ExcessModel":"Excess"
}

ThreeFactorsFit = {
    "HIST_TO_FIT" : "Visible Energy",
    "REGULATION_PARAMETER" : 0.001,
    "SCALE_FACTORS" :OrderedDict(),
    "CATEGORY_FACTORS" : {
        "NCDIS":"Pi0",
        "NCCOH":"NCCoh",
        "ExcessModel":"Excess"
    }
}

ThreeFactorsFit["SCALE_FACTORS"]["Pi0"]=[0,2.0]
ThreeFactorsFit["SCALE_FACTORS"]["NCCoh"]=[0,2.0]
ThreeFactorsFit["SCALE_FACTORS"]["Excess"]=[0,2.0]

PT_TUNE = {
    "HIST_TO_FIT" : "Lepton Pt",
    "REGULATION_PARAMETER" : 0.01,
    "HIST_OBSERVABLE": "Visible Energy vs Lepton Pt",
    "Yaxis":True,
    "SCALE_FACTORS" :OrderedDict(),
    "CATEGORY_FACTORS" : {
        "NCDIS":"Pi0",
        "CCDIS":"Pi0",
        "NCCOH":"NCCoh",
        "ExcessModel":"Excess",
        "CCNuE":"Signal",
        "CCNuEQE":"Signal",
        "CCNuEDIS":"Signal",
        "CCNuEDelta":"Signal",
        "CCNuE2p2h":"Signal",
    }
}

PT_TUNE["SCALE_FACTORS"]["Excess"]=[0,1.6]
PT_TUNE["SCALE_FACTORS"]["NCCoh"]=[0,1.6]
PT_TUNE["SCALE_FACTORS"]["Pi0"]=[0,0.2,0.4,0.6,0.8,1.0,1.2,1.6]
PT_TUNE["SCALE_FACTORS"]["Signal"]=[0,1.6]

PI0PT_TUNE = {
    "HIST_TO_FIT" : "Lepton Pt",
    "HIST_OBSERVABLE": "Visible Energy vs Lepton Pt",
    "Yaxis":True,
    "REGULATION_PARAMETER" : 0.01,
    "SCALE_FACTORS" :OrderedDict(),
    "CATEGORY_FACTORS" : {
        "NCDIS":"Pi0",
        "CCDIS":"Pi0",
        "CCNuE":"Signal",
        "CCNuEQE":"Signal",
        "CCNuEDIS":"Signal",
        "CCNuEDelta":"Signal",
        "CCNuE2p2h":"Signal",
    }
}

PI0PT_TUNE["SCALE_FACTORS"]["Pi0"]=[0,0.2,0.4,0.6,0.8,1.0,1.2,1.6]
PI0PT_TUNE["SCALE_FACTORS"]["Signal"]=[0,1.6]

PI0PT_TUNE2 = {
    "HIST_TO_FIT" : "Lepton Pt",
    "HIST_OBSERVABLE":  {"variables":["Visible Energy","Lepton Pt"],
     "name":"Eavail_Lepton_Pt_AltTune",
     },
    "Yaxis":True,
    "REGULATION_PARAMETER" : 0.01,
    "SCALE_FACTORS" :OrderedDict(),
    "CATEGORY_FACTORS" : {
        "NCDIS":"Pi0",
        "CCDIS":"Pi0",
    }
}
PI0PT_TUNE2["SCALE_FACTORS"]["Pi0"]=[0,0.2,0.4,0.6,0.8,1.0,1.2,1.6]

REGPT_TUNE = {
    "HIST_TO_FIT" : "Lepton Pt",
    "REGULATION_PARAMETER" : 10,
    "HIST_OBSERVABLE": "Visible Energy vs Lepton Pt",
    "Yaxis":True,
    "SCALE_FACTORS" :OrderedDict(),
    "CATEGORY_FACTORS" : {
        "NCDIS":"Pi0",
        "NCCOH":"NCCoh",
        "ExcessModel":"Excess",
        "CCNuE":"Signal",
        "CCNuEQE":"Signal",
        "CCNuEDIS":"Signal",
        "CCNuEDelta":"Signal",
        "CCNuE2p2h":"Signal",
        "CCNuEAntiNu":"Signal",
    }
}

REGPT_TUNE["SCALE_FACTORS"]["Pi0"]=[0,0.4,0.6,0.8,1.0,1.2,1.6]
REGPT_TUNE["SCALE_FACTORS"]["NCCoh"]=[0,1.6]
REGPT_TUNE["SCALE_FACTORS"]["Excess"]=[0,1.6]
REGPT_TUNE["SCALE_FACTORS"]["Signal"]=[0,0.4,0.6,0.8,1.0,1.2,1.6]

Q3_TUNE = {
    "HIST_TO_FIT" : "Q3",
    "Yaxis":True,
    "HIST_OBSERVABLE": "Visible Energy vs q3",
    "REGULATION_PARAMETER" : 0.001,
    "SCALE_FACTORS" :OrderedDict(),
    "CATEGORY_FACTORS" : {
        "NCDIS":"Pi0",
    }
}

Q3_TUNE["SCALE_FACTORS"]["Pi0"]=[0,0.6,0.8,1.0,1.2,1.6,2.0]

NO_TUNE = {
    "HIST_TO_FIT" : "Lepton Pt",
    "Yaxis":False,
    "HIST_OBSERVABLE": "Visible Energy vs Lepton Pt",
    "REGULATION_PARAMETER" : 1000000000,
    "SCALE_FACTORS" :OrderedDict(),
    "CATEGORY_FACTORS" : {
        "NCDIS":"Pi0",
    }
}

NO_TUNE["SCALE_FACTORS"]["Pi0"]=[0,1.6]


VISE_TUNE = {
    "HIST_TO_FIT" : "Visible Energy",
    "HIST_OBSERVABLE": "Visible Energy vs Lepton Pt",
    "Yaxis":False,
    "REGULATION_PARAMETER" : 0.001,
    "SCALE_FACTORS" :OrderedDict(),
    "CATEGORY_FACTORS" : {
        "NCDIS":"Pi0",
        "NCCOH":"NCCoh",
        "ExcessModel":"Excess"
    }
}

VISE_TUNE["SCALE_FACTORS"]["Pi0"]=[0,0.2,0.4,0.6,0.8,1.0,1.2,1.6,2.0]
VISE_TUNE["SCALE_FACTORS"]["NCCoh"]=[0,2.0]
VISE_TUNE["SCALE_FACTORS"]["Excess"]=[0,2.0]

EEL_TUNE ={
    "HIST_OBSERVABLE": "Lepton Energy",
    "HIST_TO_FIT" : "Lepton Energy",
    "Yaxis":False,
    "REGULATION_PARAMETER" : 0.5,
    "SCALE_FACTORS" :OrderedDict(),
    "CATEGORY_FACTORS" : {
        "CCNuEQE":"Signal",
        "CCNuEDelta":"Signal",
        "CCNuEDIS":"Signal",
        "CCNuE2p2h":"Signal",
        "CCNuE":"Signal",
        "NCDIS":"Pi0",
        "CCDIS":"Pi0",
        #"NCCOH":"NCCoh",
        #"ExcessModel":"Excess"
    }
}
EEL_TUNE["SCALE_FACTORS"]["Pi0"]=[0,20]
EEL_TUNE["SCALE_FACTORS"]["Signal"]=[0,20]

N4_TUNE ={
    "HIST_OBSERVABLE": "Biased Neutrino Energy",
    "HIST_TO_FIT" : "Biased Neutrino Energy",
    "Yaxis":False,
    "REGULATION_PARAMETER" : 0.0,
    "SCALE_FACTORS" :OrderedDict(),
    "CATEGORY_FACTORS" : {
        "CCPi0Proton":"dEdX",
        "CCPi0":"dEdX",
        "NCCohPi0":"dEdX",
        "NCPi0Proton":"dEdX",
        "NCPi0":"dEdX",
        "CCPi":"dEdX",
        "NCPi":"dEdX",
        "NCDiff":"dEdX",
        #"NuEElastic":"dEdX",
        "Other":"dEdX",
    }
}

#N4_TUNE["SCALE_FACTORS"]["Pi0"]=[0,2.5,4,6,9,12,15,20]
#N4_TUNE["SCALE_FACTORS"]["dEdX"]=[0,20]
N4_TUNE["SCALE_FACTORS"]["dEdX"]=PlotConfig.NEUTRINO4_EE_BINNING

Tune_Strategy_map = {
    "pt_tune":PT_TUNE,
    "pi0pt_tune":PI0PT_TUNE,
    "q3_tune":Q3_TUNE,
    "Eel_tune":EEL_TUNE,
    "visE_tune":VISE_TUNE,
    "N4_tune":N4_TUNE,
    "no_tune":NO_TUNE,
    "3factors_tune":ThreeFactorsFit,
    "regpt_tune":REGPT_TUNE,
    "pi0pt_tune2":PI0PT_TUNE2
}

def SetGlobalParameter(iput):
    if iput in Tune_Strategy_map:
        print(("appling predefined strategy: {}".format(iput)))
        tune = Tune_Strategy_map[iput]
        global HIST_TO_FIT
        try:
            HIST_TO_FIT = tune["HIST_TO_FIT"]
        except KeyError:
            pass
        global REGULATION_PARAMETER
        try:
            REGULATION_PARAMETER = tune["REGULATION_PARAMETER"]
        except KeyError:
            pass
        global SCALE_FACTORS
        try:
            SCALE_FACTORS = tune["SCALE_FACTORS"]
        except KeyError:
            pass
        global CATEGORY_FACTORS
        try:
            CATEGORY_FACTORS = tune["CATEGORY_FACTORS"]
        except KeyError:
            pass
        global HIST_OBSERVABLE
        try:
            HIST_OBSERVABLE = tune["HIST_OBSERVABLE"]
        except KeyError:
            pass
        global Yaxis
        try:
            Yaxis = tune["Yaxis"]
        except KeyError:
            pass
        return True
    else:
        return False


# -------------------------
# NEW: Matrix-fit recipes
# -------------------------

class MatrixFitRecipe:
    def __init__(
        self, name, regions, components, cate_to_comp,
        fixed_cates=None,
        kreg=None,
        hist_to_fit=None,         # <-- NEW (1D scale-factor axis variable)
        hist_observable=None,     # <-- NEW (what you load in each region for the fit)
    ):
        self.name = name
        self.regions = regions
        self.components = components
        self.cate_to_comp = cate_to_comp
        self.fixed_cates = fixed_cates or set()
        self.kreg = kreg
        self.hist_to_fit = hist_to_fit
        self.hist_observable = hist_observable



# -------- CCnue recipe --------
# Safer than relying on PlotConfig.SIGNAL_DEFINITION unless you KNOW it matches HistHolder.hists keys.
# If you do have an authoritative list of category keys, replace CCNUE_SIGNAL_KEYS with that.
CCNUE_SIGNAL_KEYS = set([
    "CCNuE", "CCNuEQE", "CCNuEDelta", "CCNuEDIS", "CCNuE2p2h",
    # include other CCnue signal buckets you actually have as HistHolder keys
])

# If CCnue uses split elastic categories instead of "NuEElastic", add them here too.
CCNUE_FIXED_KEYS = set([
    "NuEElastic",  # keep if present
    # "NuEElasticE","NuEElasticMu","NuEElasticOther",  # uncomment if these are the actual keys
])

def cate_to_comp_CCNUE(cate):
    if cate in CCNUE_SIGNAL_KEYS:
        return "SignalLike"
    if cate in CCNUE_FIXED_KEYS:
        return None
    if cate == "Total":
        return None
    return "BkgLike"

CCNUE_MATRIX_RECIPE = MatrixFitRecipe(
    name="CCnue_matrix",
    regions=["Signal", "dEdX"],
    components=["BkgLike", "SignalLike"],
    cate_to_comp=cate_to_comp_CCNUE,
    fixed_cates=set(["NuEElastic"]),
    kreg=2,
    hist_to_fit="Biased Neutrino Energy",                 # EN4 axis
    hist_observable={"name":["Biased Neutrino Energy"]},  # EN4 observable
)


# -------- nu+e recipe --------
# FIT_SIGNAL_CATES_NUE = set(["NuEElasticE", "NuEElasticMu", "NuEElasticOther"])

# FIT_BKG_GROUPS_NUE = {
#     "BkgPhoton": set(["NCPi0", "NCCohPi0"]),
#     "BkgCC": set([
#         "CCNuEWrongSign", "CCNuEQE", "CCNuEDelta", "CCNuEDIS", "CCNuE2p2h", "CCNuE",
#         "CCPi0", "CCPi", "CCNuMu", "CCOther",
#     ]),
#     "BkgOther": set(["NCPi", "NCOther"]),
# }

def cate_to_comp_NUE(cate: str):
    # --- signal ---
    if cate.startswith("NuEElastic"):
        return "Signal"

    # --- CC backgrounds (everything CC) ---
    if cate.startswith("CC"):
        return "BkgCC"

    # --- photon-like backgrounds (pi0 / coherent pi0 / photon misID if you have it) ---
    if cate in ("NCNuECohPi0", "NCNuMuCohPi0", "NCPi0"):
        return "BkgPhoton"

    # --- everything else (NC non-pi0, unknown/diff, etc.) ---
    if cate in ("NCPi", "NCOther", "NCDiff", "Other"):
        return "BkgOther"

    # If something new appears, fail loudly (better than silently unmapped)
    return None

NUE_MATRIX_RECIPE = MatrixFitRecipe(
    name="nueElastic_matrix",
    regions=["Signal", "dEdX", "Eavail", "Etheta"],
    components=["BkgPhoton", "BkgCC", "BkgOther", "Signal"],
    cate_to_comp=cate_to_comp_NUE,
    fixed_cates=set(),
    kreg=3,
    hist_to_fit="Lepton Energy",                 # Eel axis
    hist_observable="Lepton Energy",             # or {"name":["Lepton Energy"]} to match your HistHolder API
)


# Helper map so the fitter can select by tag
MATRIX_RECIPE_MAP = {
    "CCnue_matrix": CCNUE_MATRIX_RECIPE,
    "nueElastic_matrix": NUE_MATRIX_RECIPE,
}

def GetMatrixRecipe(tag):
    return MATRIX_RECIPE_MAP.get(tag, None)



#make sure all groups have a color and defination
#assert set(BACKGROUND_GROUP.keys())==set(BACKGROUND_GROUP_COLOR.keys())
