import ROOT
from collections import OrderedDict
from config import PlotConfig
from tools.PlotLibrary import HistHolder
from tools import PlotTools
from config.SystematicsConfig import CONSOLIDATED_ERROR_GROUPS,DETAILED_ERROR_GROUPS



Default_Plot_Type="stacked"
Default_Scale = lambda histHolder:histHolder.POTScale(True)
DefaultSlicer = PlotTools.PrepareSlicer


# the order and color of each band in stacked plots defined here.
COLORS=ROOT.MnvColors.GetColors()
Categories = OrderedDict()
Categories["ExcessModel"] = {
    "title":"NC DFR",
    "color" : COLORS[2]
}
Categories["CCNuE"]={
    "title" : "CC #nu_{e}",
    "cate":{"CCNuEQE","CCNuEDelta","CCNuEDIS","CCNuE","CCNuE2p2h"},
    "color": COLORS[0]
}

# Categories["2p2h"]={
#     "title": "CC #nu_{e} 2p2h",
#     "cate" :{"CCNuE2p2h"},
#     "color": COLORS[1]
# }
Categories["CCDIS"]= {
    "title": "CC #nu_{#mu} #pi^{0}",
    "color" : COLORS[3]
}
Categories["NCCOH"]= {
    "title": "NC coh",
    "color" : COLORS[4]
}

Categories["NCDIS"]= {
    "title": "NC #pi^{0}",
    "color" : COLORS[5]
}

Categories["CCNuEAntiNu"] = {
    "title":"Wrong Sign",
    "color" : COLORS[6]
}
Categories["NuEElastic"] = {
    "title":"#nu+e",
    "color" : COLORS[1]
}

Categories["Other"] = {
    "title":"Others",
    "cate" : {"Other","NCOther","CCOther","NonPhaseSpace","NonFiducial"},
    "color" : COLORS[7]
}


SignalDecomposition = {
    "CCNuEQE" :
    {
        "title" : "CC #nu_{e}-QE",
        "color": COLORS[0]
    },
    "CCNuEDelta" : {
        "title" : "CC #nu_{e}-Delta",
        "color": COLORS[1]
    },
    "CCNuEDIS" : {
        "title" : "CC #nu_{e}-DIS",
        "color": COLORS[2]
    },
    "CCNuE2p2h" : {
        "title" : "CC #nu_{e}-2p2h",
        "color": COLORS[5]
    },
    "CCNuE": {
        "title" : "CC #nu_{e}-Other",
        "color": COLORS[3]
    },
    "Background" : {
        "title" : "Backgrounds",
        "cate": {"ExcessModel","NonFiducial","CCDIS","CCOther","NCCOH","NCRES","NCDIS","NCOther","NuEElastic","Other","CCNuEAntiNu","NonPhaseSpace"},
        "color": COLORS[4]
    }
    # "CCNuEAntiNu" : {
    #     "title" : "CC #bar{#nu_{e}}",
    #     "color": COLORS[4]
    # }
}

SignalBackground = {
    "Signal" :
    {
        "title" : "CC #nu_{e}",
        "cate": {
            "CCNuEQE","CCNuEDelta","CCNuEDIS","CCNuE2p2h","CCNuE"
        },
        "color": COLORS[0]
    },
    "Background" : {
        "title" : "Backgrounds",
        "cate": {"ExcessModel","NonFiducial","CCDIS","CCOther","NCCOH","NCRES","NCDIS","NCOther","NuEElastic","Other","CCNuEAntiNu","NonPhaseSpace"},
        "color": COLORS[4]
    }
}

NuEElasticCategory = {
    "Nu+e" :{
        "title": "nu+e",
        "cate" :{"NuEElastic"},
        "color":COLORS[1]
    },
    "CCNuE":
    {
        "title" : "CC #nu_{e}",
        "cate": {
            "CCNuEQE","CCNuEDelta","CCNuEDIS","CCNuE2p2h","CCNuE"
        },
        "color": COLORS[0]
    },
    "NCCOH":{
        "title" : "NCCoh",
        "color": COLORS[3]
    },
    "ExcessModel": {
        "title": "NC Diffractive",
        "color":COLORS[2]
    },
    "Other":{
        "title" : "Others",
        "cate": {"ExcessModel","NonFiducial","CCDIS","CCOther","NCRES","NCDIS","NCOther","Other","CCNuEAntiNu","NonPhaseSpace"},
        "color": COLORS[4]
    }
}

MODELS = {
    "MnvTune v1": {
        "errorband":(None,None),
        "color":COLORS[0]
    },
    "2p2h Tune (QE)": {
        "errorband":("Low_Recoil_2p2h_Tune",2),
        "color":COLORS[1]
    },
    "SuSA 2p2h" : {
        "errorband":("SuSA_Valencia_Weight",0),
        "color":COLORS[2]
    },
    "MK model": {
        "errorband": ("MK_model",0),
        "color":COLORS[7]
    },
    "Low Q2 Pion Joint": {
        "errorband" : ("LowQ2Pi_Joint",0),
        "color": COLORS[5]
    },
    "Low Q2 Pion NuPi0": {
        "errorband" : ("LowQ2Pi_NUPI0",0),
        "color": COLORS[6]
    }
}

DefaultPlotters={
    "comp":{"func": PlotTools.PrepareComp},
    "ratio":{"func": PlotTools.PrepareRatio},
    "err": {"func": PlotTools.PrepareErr,
            "args":(True,False,CONSOLIDATED_ERROR_GROUPS)},
    "stacked":{"func":PlotTools.PrepareStack,
               "args": (Categories,)},
    "diff":{"func":PlotTools.PrepareDiff},
    "migration":{"func":PlotTools.PrepareMigration},
    "sigdep":{"func":PlotTools.PrepareSignalDecompose,
              "args": (SignalDecomposition,True,False)},
    "errband":{"func":PlotTools.PrepareErrorBand},
    "model":{},
    "model_ratio" :{},
    "sigdep":{},
}

PLOTS_TO_MAKE = [
    # {"variables":["Visible Energy","Lepton Pt"],
    #   "plot_type":"errband",
    #   "args":("Leakage_Uncertainty",),
    #   },
    # {"name":"Lepton Energy",
    #  "plot_type":"stacked"},
    {"name":"Neutrino Energy",
      "plot_type":"stacked"},
    {"variables":["Neutrino Energy QE"],
      "plot_type":"stacked"},

    # {"name":"Neutrino Energy",
    #   "sideband_group": ("ex_and_pi0",("Excess","Pi0"))},
    # {"name":"Visible Energy",
    #   "sideband_group": ("ex_and_pi0",("Excess","Pi0"))},
    # {"name": "Lepton Theta"},
    # {"name": "Lepton Theta"},
    # {"name": "Psi"},
    {"name":"Q0"},
    {"name":"Q3"},
    {"name":"Visible Energy",
      "plot_type" : "comp",},
    # {"name":"Q3",
    #  "plot_type" : "sigdep",},
    {"name":"Visible Energy",
     "plot_type" : "err",},
    {"variables":["Leading Pi0 E","PsiEe"],
     "slicer":lambda hist: PlotTools.Make2DSlice(hist,True,interval=2),
     "sideband_group": ("Excess",("Excess_High_Inline","Excess_Low_Inline")),
     },

    {"variables":["Lepton Energy","Q3"],
      "slicer": lambda hist: [hist.ProjectionX("{}_px".format(hist.GetName()),0,6)]},
    {"variables":["Lepton Energy"]},
    {"variables":["Biased Neutrino Energy","Q3"],
     "slicer": lambda hist: [hist.ProjectionX("{}_px".format(hist.GetName()),0,6)]},
    {"variables":["Biased Neutrino Energy"]},
    {"variables":["Lepton Energy","Lepton Pt"]},
    {"variables":["Biased Neutrino Energy","Lepton Pt"]},
 
    {"variables":["Visible Energy","Lepton Energy"]},
    {"variables":["Visible Energy","Lepton Pt"],
     "plot_type":"stacked",
     "args":(Categories,("bkg_tune",0)),
     "tag":"bkg_tune1"},
    {"variables":["Visible Energy","Lepton Pt"],
     "plot_type":"stacked",
     "args":(Categories,("bkg_tune",1)),
     "tag":"bkg_tune2"},
    # {"variables":["PsiEe","Lepton Energy"],
    #  "slicer": PlotTools.ProjectionY},
    #  {"variables":["PsiEe","Lepton Energy"],
    #  },
    #  {"variables":["Visible Energy","Lepton Energy"],
    #  },
    # {"variables": ["Inline Upstream Energy","Lepton Energy"],
    #  },
    #  {"name":"Lepton Theta"},
    #   {"name":"Lepton Theta",
    #    "plot_type":"diff"},
    # {"variables": ["Vertex Difference","Lepton Energy"],
    #  },
    # {"variables": ["UIE Mean Pos","Lepton Energy"],
    #  },
    # {"variables": ["Shower Width","Lepton Energy"],
    #  },
    # {"variables": ["UIE Pos RMS","Lepton Energy"],
    #  },
    #  # {"variables": ["Rock Muon Removed","Lepton Energy"],
    #  #  "slicer": lambda hist: PlotTools.Make2DSlice(hist,interval=5),
    #  #  },
    #    {"variables": ["Vertex Energy","Lepton Energy"],
    #   "slicer": lambda hist: PlotTools.Make2DSlice(hist,interval=5),
    #   },
    # {"name":"Q0",
    #  "plot_type" : "err",},
    # {"name":"Q3",
    #  "plot_type" : "err",},
    # {"name":"Q3 Migration",
    #  "plot_type":"migration",
    # },
    # {"name":"Lepton Pt Migration",
    #  "plot_type":"migration"},
    #  {"name":"Visible Energy Migration",
    #   "plot_type":"migration"},
    # {"name":"Visible Energy vs q3 Migration",
    #  "plot_type":"migration"},
    # {"name":"Lepton Energy Migration",
    #  "plot_type":"migration"},
    # {"name":"Lepton Energy Resolution",
    #  "plot_type":"comp"},
    # {"name":"Lepton Theta Resolution",
    #  "plot_type":"comp"},
    # {"name":"Lepton Theta Migration",
    #  "plot_type":"migration"},
    # {"name":"Visible Energy vs Lepton Pt Migration",
    #  "plot_type":"migration"},
    {"name":"Visible Energy vs q3"},
    {"name":"Visible Energy vs q3",
     "plot_type":"stacked",
     "args":(Categories,("bkg_tune",0)),
     "tag":"bkg_tune1"},
       {"name":"Visible Energy vs q3",
     "plot_type":"stacked",
     "args":(Categories,("bkg_tune",1)),
     "tag":"bkg_tune2"},
    {"name":"Visible Energy vs q3",
     "plot_type":"comp"},
     {"name":"Visible Energy vs q3",
     "plot_type":"err"},
    {"variables":["Q0","True Visible Energy"],
     "scale":lambda hist:hist.AreaScale(True),
     "slicer":lambda hist: PlotTools.Make2DSlice(hist,True),
     "plot_type":"comp",
     "args":(True,False,True),
     },


    # {"name":"True Signal Visible Energy vs q3"
    #  "plot_type":"err"},
    # {"name":"Visible Energy vs q3",
    #  "plot_type" : "comp"},
    # {"name":"Visible Energy vs q3",
    #  "plot_type" : "ratio"},
    # {"name":"Visible Energy vs q3",
    #  "plot_type" : "err"},
    # {"name":"Visible Energy vs q3",
    #  "plot_type" : "sigdep",},
    # {"name":"q0 vs q3",
    #  "plot_type" : "sigdep",},
    # {"name":"Visible Energy"},
    {"variables":["Visible Zoomin","Lepton Energy"],
     "slicer": lambda hist: PlotTools.Make2DSlice(hist,interval=1),
     "args":(NuEElasticCategory,),
     "scale":lambda histHolder:histHolder.POTScale(False)},
    # {"name":"FrontdEdX_position" },
    {"name":"Front dEdX",
      "sideband_group": ("All",("Excess_High_Inline","Excess_Low_Inline","Pi0","Signal")),
     },
    {"variables":["PsiEe"],
      "sideband_group": ("Sidebands",("Excess_High_Inline","Excess_Low_Inline","Pi0")),
     },
    # {"name":"Inline Upstream Energy Weighted Position",
    #  "canvasconfig":PlotTools.Logy},
    # {"name":"Inline Upstream Energy",
    #  "canvasconfig":PlotTools.Logx,
    #  "scale" : lambda histHolder:histHolder.POTScale(False)},
    {"variables":["Visible Energy","Lepton Pt"]},
    {"variables":["Visible Energy","Lepton Pt"],
     "plot_type":"comp"},
     {"variables":["Visible Energy","Lepton Pt"],
     "plot_type":"err"},
    # {"variables":["Visible Energy","Lepton Pt"],
    #  "slicer": PlotTools.IdentitySlicer,
    #  "scale":lambda histHolder:histHolder.POTScale(False),
    #  "plot_type":"diff"},
    # {"variables":["Visible Energy","Lepton Pt"],
    #  "sideband_group": ("Excess",("Excess_High_Inline","Excess_Low_Inline")),
    #  "slicer": PlotTools.IdentitySlicer,
    #  "scale":lambda histHolder:histHolder.POTScale(False),
    #  "plot_type":"diff"},
    {"variables":["Visible Energy","Lepton Pt"],
     "sideband_group": ("Excess",("Excess_High_Inline","Excess_Low_Inline"))},
    {"name":"Visible Energy vs q3",
     "sideband_group": ("Excess",("Excess_High_Inline","Excess_Low_Inline"))},
    #  {"variables":["Lepton Energy"],
    #  "sideband_group": ("Excess",("Excess_High_Inline","Excess_Low_Inline"))},
    # {"variables":["Sum Visible Energy","Lepton Pt"],
    #  "sideband_group": ("Excess",("Excess_High_Inline","Excess_Low_Inline")),
    #  "slicer":PlotTools.ProjectionX},
    # {"variables":["Sum Visible Energy","Lepton Pt"],
    #  "slicer":PlotTools.ProjectionX},
    {"variables":["Visible Energy","Lepton Pt"],
     "plot_type":"err"},
    {"name":"Lepton Pt"},
    # {"variables":["Lepton Pt"],
    #  "plot_type":"err"},
    # {"variables":["Visible Energy","Lepton Energy"]},
    # {"variables":["N_EMLike_prongs"]},
    # {"variables":["Proton Angle"]},
    # {"variables":["Proton Angle Residual"],
    #  "plot_type":"sigdep"},
    # {"variables":["Proton Angle Residual"]},
    # {"name":"Lepton Energy Low Inline",
    #  "plot_type" : "comp"},
    # {"name":"Lepton Energy High Inline",
    #  "plot_type" : "comp"},
    #  {"name":"Lepton Energy Low Inline",
    #  "plot_type" : "err"},
    # {"name":"Lepton Energy High Inline",
    #  "plot_type" : "err"},
    # {"name":"Pseudo q3"},
    # {"name":"Psi"},
    # {"variables":["W"]},
    # {"variables":["Qsquare"]},
    #{"variables":["Vertex Multiplicity"]},
    #{"variables":["Start Multiplicity"]},
    #{"variables":["Lepton Visible Energy"],
    # },
    #{"variables":["Sum Visible Energy"],
    # },

    {"variables":["Exuv"]},
    {"variables":["Euv"]},
    # {"name":"True Signal Visible Energy vs q3",
    #  "plot_type":"comp"},
    #{"variables":["Visible Energy","Lepton Visible Energy"],
    # "plot_type":"comp"},
    # {"variables":["Sum Visible Energy"],
    #  "plot_type":"stacked",
    #  },
    # {"name":"Sum Visible Energy low inline"},
    # {"name":"Sum Visible Energy high inline"},
    # {"variables":["Front dEdX", "Visible Energy"],
    #  "sideband_group": ("Merged",("Signal","Excess","Pi0"))},
    #  {"variables":["Front dEdX", "Q3"],
    #  "sideband_group": ("Merged",("Signal","Excess","Pi0"))},
    # {"variables":["Front dEdX", "Lepton Pt"],
    #  "sideband_group": ("Merged",("Signal","Excess","Pi0"))},
    # {"variables":["Front dEdX", "Sum Visible Energy"],
    #  "sideband_group": ("Merged",("Signal","Excess","Pi0"))},

    # {"variables":["Lepton Visible Energy","Visible Energy"]},
    # {"variables":["Sum Visible Energy","Visible Energy"]},
    # #  {"variables":["Front dEdX", "Q3"],
    # # "sideband_group": ("Merged",("Excess","Pi0","Signal"))},
    # {"variables":["Vertex Apothem","Visible Energy"]},
    #  #"sideband_group": ("Merged",("Signal","Excess","Pi0"))},
    # {"variables":["Vertex Y","Visible Energy"]},
    # {"variables":["Vertex Z","Visible Energy"]}
]

