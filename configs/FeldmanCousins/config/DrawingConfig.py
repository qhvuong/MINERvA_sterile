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
COLORS=ROOT.MnvColors.GetColors(ROOT.MnvColors.k36Palette)
Categories = OrderedDict()

Categories["CCNuE"]= {
    "title": "CC Anti#nu_{e}",
    "cate":{"CCNu","CCNuEQE","CCNuEDelta","CCNuE2p2h","CCNuEDIS","CCNuE"},
    "color" : COLORS[0]
}

#Categories["CCQElike"]={
#    "title" : "CCQElike #bar#nu_{e}",
#   #"cate":{"CCNuEQE","CCNuEDelta","CCNuEDIS","CCNuE2p2h","CCNuE","NonPhaseSpace"},
#    "color": COLORS[0]
#}
#Categories["notCCQElike"]= {
#    "title": "notCCQElike #bar{#nu}_{e}",
#    #"cate" : {"CCPi0","CCOther"},
#    "color" : COLORS[1]
#}
Categories["NCDiff"] = {
     "title":"NC Diffractive #pi^{0}",
     "color" : COLORS[2]
}
Categories["CCPi0"]= {
    "title": "Not #nu_{e} CC #pi^{0}",
    "color" : COLORS[3]
}
Categories["NCCohPi0"]= {
    "title": "NC Coh #pi^{0}",
    "color" : COLORS[4]
}
Categories["NCPi0"]= {
    "title": "NC #pi^{0}",
    "color" : COLORS[5]
}
Categories["NCPi"]= {
    "title": "NC #pi+-",
    "color" : COLORS[8]
}
Categories["CCPi"]= {
    "title": "CC #pi+-",
    "color" : COLORS[9]
}
Categories["NuEElastic"] = {
    "title":"#nu + e elastic",
    "color" : COLORS[6]
}
Categories["NCOther"] = {
    "title":"NCOthers",
    "color" : COLORS[7]
}
Categories["CCOther"] = {
    "title":"CCOthers",
    "color" : COLORS[10]
}

BackgroundDecomposition = {
    "NCDiff" :
    {
        "title" : "NC Diffractive",
        "color": COLORS[0]
    },
    "NuEElastic" : {
        "title" : "#nu + e",
        "color": COLORS[1]
    },
    "CCPi0" : {
        "title" : "CC #pi^{0}",
        "color": COLORS[3]
    },
    "NCPi0" : {
        "title" : "NC #pi^{0}",
        "color": COLORS[5]
    },
    "CCPi" : {
        "title" : "CC #pi",
        "color": COLORS[6]
    },
    "NCPi" : {
        "title" : "NC #pi",
        "color": COLORS[7]
    },
    "NCCohPi0": {
        "title" : "NC Coh #pi^{0}", 
        "color": COLORS[8]
    },
    "NCOther" : {
        "title" : "NCOther",
        "color": COLORS[9]
    },
    "CCOther" : {
        "title" : "CCOther",
        "color": COLORS[10]
    }
    }

BackgroundDecomposition = Categories.copy()
del BackgroundDecomposition["CCNuE"]

#Categories = OrderedDict()
#Categories["NCnumu"] = {
#     "title":"NC #nu_{#mu}",
#     "color" : COLORS[2]
#}
#Categories["NCnue"] = {
#     "title":"NC #nu_{e}",
#     "color" : COLORS[3]
#}
#Categories["CCnue"] = {
#     "title":"CC #nu_{e}",
#     "color" : COLORS[4]
#}
#Categories["CCantinue"] = {
#     "title":"CC anti #nu_{e}",
#     "color" : COLORS[5]
#}
#Categories["NueOther"] = {
#    "title":"NueOthers",
#    "color" : COLORS[0]
#}

#Categories["ExcessModel"] = {
#    "title":"Ad hoc excess model",
#    "color" : COLORS[2]
#}



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
        "color": COLORS[6]
    },
    "CCNuE": {
        "title" : "CC #nu_{e}-Other", 
        "color": COLORS[4]
    },
    "CCNuEWrongSign": {
        "title" : "CC #nu_{e} wrong sign", 
        "color": COLORS[7]
    },
     "Background" : {
        "title" : "Backgrounds",
        "cate" : {"NonPhaseSpace","NCDiff", "CCPi0", "NCPi0", "NCCohPi0","NuEElastic","CCPi","NCPi","CCOther","NCOther","Other"}, 
        "color": COLORS[5]
    }
    }

#SignalDecomposition = {
#    "NueBar" :
#    {
#        "title" : "Anti#nu_{e}", 
#        "color": COLORS[0]
#    },
#    "Nue" :
#    {
#        "title" : "#nu_{e}", 
#        "color": COLORS[1]
#    },
#    "Background" : {
#        "title" : "Backgrounds",
#        "cate": {"Other"},
#        "color": COLORS[2]
#    }
#}
#}


#SignalDecomposition = {
#    "CCNuE" :
#    {
#        "title" : "Signal",
#        "color": COLORS[0]
#    },
#    "Background" : {
#        "title" : "Backgrounds",
#        "cate": {"NCDiff", "NuEElastic", "NonPhaseSpace", "NonFiducial", "CCPi0", "NCCohPi0", "NCPi0", "Other"},
#        "color": COLORS[2]
#    }
#}

SignalBackground = {
    "Signal" :
    {
        "title" : "CC #nu_{e}",
        "cate":{"CCNu","CCNuEQE","CCNuEDelta","CCNuE2p2h","CCNuEDIS","CCNuE"},
        "color": COLORS[0]
    },
    "Background" : {
        "title" : "Backgrounds",
        #"cate": {"CCPi0","Other","NCDiff","NuEElastic","NonPhaseSpace","NonFiducial","NCCohPi0","NCPi0","CCNu"},
        #"cate": {"CCPi0","Other","NCDiff","NuEElastic","NCCohPi0","NCPi0","CCNu"},
        "cate": {"CCPi0","CCPi","NCPi","NCOther","CCOther","NCDiff","NuEElastic","NCCohPi0","NCPi0"},
        "color": COLORS[4]
    }
}

SignalChargedBackground = {
    "Signal" :
    {
        "title" : "CC #nu_{e}",
        "cate":{"CCNu","CCNuEQE","CCNuEDelta","CCNuE2p2h","CCNuEDIS","CCNuE"},
        "color": COLORS[0]
    },
    "CCBackground" : {
        "title" : "Charged Current Backgrounds",
        "cate": {"CCPi0","CCPi0Proton","CCOther","CCPi"},
        "color": COLORS[1]
    },
    "NCBackground" : {
        "title" : "Neutral Current Backgrounds",
        "cate": {"Other","NCOther","NCDiff","NCCohPi0","NCPi0","NCPi0Proton","NCPi"},
        "color": COLORS[2]
    },
    "NuEEl" :{
        "title" : "#nu + e elastic",
        "cate": {"NuEElastic"},
        "color": COLORS[3]
    }
}

ChargedBackground = {
    "CCBackground" : {
        "title" : "Charged Current Backgrounds",
        "cate": {"CCPi0","CCPi0Proton","CCOther","CCPi"},
        #"cate": {"CCPi0","CCPi"},
        "color": COLORS[1]
    },
    "NCBackground" : {
        "title" : "Neutral Current Backgrounds",
        "cate": {"Other","NCOther","NCDiff","NCCohPi0","NCPi0","NCPi0Proton","NCPi"},
        "color": COLORS[2]
    },
    "NuEEl" :{
        "title" : "#nu + e elastic",
        "cate": {"NuEElastic"},
        "color": COLORS[3]
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
        #"cate": {"ExcessModel","NonFiducial","CCDIS","CCOther","NCRES","NCDIS","NCOther","Other","CCNuEAntiNu"},
        "cate": {"ExcessModel","CCDIS","CCOther","NCRES","NCDIS","NCOther","Other","CCNuEAntiNu"},
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
    "bkgratio":{"func": PlotTools.PrepareBkgRatio},
    "err": {"func": PlotTools.PrepareErr,
            "args":(SignalDecomposition,True,False,CONSOLIDATED_ERROR_GROUPS)},
    "stacked":{"func":PlotTools.PrepareStack,
        "args": (SignalDecomposition,)},
    "diff":{"func":PlotTools.PrepareDiff},
    "migration":{"func":PlotTools.PrepareMigration},
    "category_hist":{"func":PlotTools.CategoryHist,
                     "args": (SignalChargedBackground,)},
    "sigdep":{"func":PlotTools.PrepareSignalDecompose,
              "args": (ChargedBackground,True,False)},
    "sigdepratio":{"func":PlotTools.PrepareSignalDecomposeRatio,
              "args": (SignalDecomposition,True)},
    "errband":{"func":PlotTools.PrepareErrorBand},
    "model":{},
    "model_ratio" :{},
    #"sigdep":{},
}

PROB_PLOTS = [
    {"name":"True Energy vs Neutrino Length Travelled"},
    {"name":["Biased Neutrino Energy"]},
    {"name":"True Energy Inverse vs Biased Neutrino Energy",
      "plot_type" : "migration"},
]

PLOTS_TO_MAKE = [
####Taken from old script##########
    #{"name":"E Theta SquaredVisible Energy vs q3"}, 
    #{"name":"Q3 Migration",
    #    "plot_type" : "migration"},
    #{"name":"Visible Energy"},
    #{"name":"Visible Energy vs q3 Migration",
    # "plot_type":"migration"},
    #{"name":"Visible Energy vs Lepton Pt Migration",
    # "plot_type":"migration"},
    #{"variables":["Biased Neutrino Energy"]},
    #{"name":"Visible Energy vs Lepton Pt",
    #    "plot_type" : "migration",
    # "slicer": lambda x: [x.Clone()]},
    #{"name":"Biased Neutrino 4 Energy vs q3"},
    #{"name":"Biased Neutrino 4 Energy vs Lepton Pt"},
      #"scale" : lambda histHolder:histHolder.POTScale(False)}, 
    #{"name":"Biased Neutrino 4 Energy vs Available Energy"},
      #"scale" : lambda histHolder:histHolder.POTScale(False)}, 
    #{"name":"E Theta Squared"},
      #"scale" : lambda histHolder:histHolder.POTScale(False)}, 
   # {"name":"Biased Neutrino 4 Energy vs E Theta Squared"},
  #  {"name":"Biased Neutrino 4 Energy vs E Theta Squared"},
#    {"name":"Visible Energy vs q3",
#        "plot_type" : "migration",
#     "slicer": lambda x: [x.Clone()]},
#    {"name":"True Neutrino Energy"},
#    {"name":"Electron Energy vs q3"},
#    {"name":"Electron Theta vs q3"},
    #{"name":"True Signal Visible Energy vs q3",
    # "plot_type" : "sigdep"},
    #{"name":"True Signal Visible Energy vs Lepton Pt",
    # "plot_type" : "sigdep"},
#    {"name":"True Neutrino Energy",
#     "plot_type" : "sigdep",
#     "args": (PlotConfig.SignalDecomposition,True,False)}, 

####Hang Plots#####################
    # {"name":"Neutrino Energy",
    #  "plot_type":"stacked"},
    # {"name":"Neutrino Energy",
    #   "sideband_group": ("ex_and_pi0",("Excess","Pi0"))},
    # {"name":"Visible Energy",
    #   "sideband_group": ("ex_and_pi0",("Excess","Pi0"))},
    # {"name": "Lepton Theta"},
    # {"name": "Lepton Theta"},
    #{"name": "Psi"},
    #{"name":"Visible Energy vs Lepton Pt"},
    #{"variables":["Visible Energy","Lepton Pt"]},
    #{"variables":["Visible Energy","Lepton Pt"],
    # "sideband_group": ("Excess",("Excess_High_Inline","Excess_Low_Inline"))}, 
    #{"name":"Q0"},
    #{"name":"Q3"},
    #{"name":"Lepton Pt"},
    #{"name":"Lepton Pt Migration",
    #   "plot_type":"migration",
    #  "slicer": lambda x: [x.Clone()]}, 
    #{"name":"Lepton Energy"},  
    #{"name":"Vertex Z coordinate"},
    #{"variables":["Vertex Z coordinate"],
    #"sideband_group": ("Excess",("Excess_High_Inline","Excess_Low_Inline"))},
    #{"name":"Vertex R"},
    #{"name":"Transverse Shower Asymmetry"},
    #{"name":"Extra Energy Long"},
    #{"name":"Extra Energy Trans"},
    #{"name":"EeTheta2"}, 
    #{"name":"Vertex Energy"},
    #{"name":"Inline Upstream Energy"},
    #{"name":"Inline Upstream Energy",
    #   "sideband_group": ("ex_and_pi0",("Excess","Pi0"))},   
    # {"name":"Inline Upstream Energy",
    #  "canvasconfig":PlotTools.Logx,
    #  "scale" : lambda histHolder:histHolder.POTScale(False)}, 
    #{"name":"Extra Energy"},
   # {"name":"Lepton Energy",
   #   "plot_type" : "sigdep",},
    #{"name":"Lepton Theta"},
   #{"name":"True Visible Energy vs Visible Difference"},
    #{"name":"Visible Energy",
    #  "plot_type" : "comp",},
    #{"name":"Visible Energy Migration",
    #  "plot_type":"migration",
    #  "slicer": lambda x: [x.Clone()]}, 
    # {"name":"Q3",
    #  "plot_type" : "sigdep",},
    #{"name":"Visible Energy",
    # "plot_type" : "sigdep",},
    #{"variables":["Leading Pi0 E","PsiEe"],
    # "slicer":lambda hist: PlotTools.Make2DSlice(hist,True,interval=2),
    # "sideband_group": ("Excess",("Excess_High_Inline","Excess_Low_Inline")),
    # },
    #{"variables":["PsiEe","Lepton Energy"],
    # "slicer": lambda hist: PlotTools.Make2DSlice(hist,interval=5),
    # "sideband_group": ("Excess",("Excess_High_Inline","Excess_Low_Inline")),
    # },
    # {"variables":["PsiEe","Lepton Energy"],
    # "slicer": lambda hist: PlotTools.Make2DSlice(hist,interval=5),
    # },
    #{"variables":["Lepton Energy","Q3"]},
    #{"variables":["Lepton Energy"]},
  #    "scale" : lambda histHolder:histHolder.POTScale(False)}, 


  
    #{"name":"Signal Reco Energy vs L/E",
    #  "plot_type" : "migration"},
    #{"name":"Signal Reco Energy vs sin",
    #  "plot_type" : "migration"},
    #E Theta Squared",


    #{"name":"Available Energy vs Lepton Energy"},
    #{"name":"Available Energy vs Lepton Pt"},
    #{"name":"Lepton Pt vs Lepton Energy"},
    #{"name":"Lepton Pt vs Available Energy"},

    
    #{"name":"Lepton Energy vs Available Energy",
    #   "plot_type" : "migration"},
    #{"name":"True Lepton Energy vs Available Energy"},
    #{"name":"True Lepton Energy vs Available Energy",
    #   "plot_type" : "migration"},

    #"Biased Neutrino 4 Energy vs E Theta Squared",
    #{"name":"Biased Neutrino Energy"}, 

    ###### NEUTRINO4 PLOTS ########
    #{"name":"Biased Neutrino Energy",
    #  "plot_type" : "stacked"},
    #{"name":"Reco Energy vs L/E"},
    #{"name":"Signal Biased Neutrino Energy",
    #  "plot_type" : "err"},
    #{"name":"Signal Biased Neutrino Energy",
    #  "plot_type" : "sigdep"},
    #{"name":"Biased Neutrino Energy",
    #  "plot_type" : "ratio"},
    #{"name":"Biased Neutrino Energy",
    #  "plot_type" : "err"},
    #{"name":"Biased Neutrino Energy",
    #  "plot_type" : "sigdep"},
    #{"variables":["Neutrino Length Travelled"]},
    #{"variables":["Visible Energy"]},
    #{"variables":["Lepton Pt"]},
    #{"name":"Available Energy vs Lepton Energy",
    #    "plot_type" : "stacked",},
    #{"name":"Available Energy vs Lepton Energy",
    #    "plot_type" : "bkgratio",},

    #{"name":"Biased Neutrino Energy"},
    #{"name":"Biased Neutrino Energy",
    #  "plot_type" : "ratio"},
    #{"name":"Biased Neutrino Energy",
    #  "plot_type" : "err"},
    #{"name":"Biased Neutrino Energy",
    #  "plot_type" : "sigdep"},
    
    #{"variables":["Neutrino Length Travelled"]},
    #{"variables":["Visible Energy"]},
    #{"variables":["Lepton Pt"]},

    #{"name":"Available Energy vs Lepton Energy",
    #    "plot_type" : "stacked",
    #    "canvasconfig":PlotTools.Logx,},
    #{"name":"Available Energy vs Lepton Pt",
    #    "plot_type" : "stacked",
    #    "canvasconfig":PlotTools.Logx,},
    #{"name":"Lepton Pt vs Lepton Energy",
    #    "plot_type" : "stacked",
    #    "canvasconfig":PlotTools.Logx,},
    # {"name":"Lepton Pt vs Available Energy",
    #    "plot_type" : "stacked",
    #    "canvasconfig":PlotTools.Logx,},

    #{"name":"True Signal Biased Neutrino Energy"},
    #{"name":"Signal Biased Neutrino Energy"},

    ###### NEUTRINO4 TEMPLATES ########
    #{"name":"template 0",
      #"plot_type" : "migration"},
    #{"name":"template 1",
      #"plot_type" : "migration"},
    #{"name":"template 2",
      #"plot_type" : "migration"},
    #{"name":"template 3",
      #"plot_type" : "migration"},
    #{"name":"template 4",
      #"plot_type" : "migration"},
    #{"name":"template 5",
      #"plot_type" : "migration"},
    #{"name":"template 6",
      #"plot_type" : "migration"},
    #{"name":"template 7",
      #"plot_type" : "migration"},
    #{"name":"template 8",
      #"plot_type" : "migration"},
    #{"name":"template 9",
      #"plot_type" : "migration"},
    #{"name":"template 10",
      #"plot_type" : "migration"},
    #{"name":"template 11",
      #"plot_type" : "migration"},
    #{"name":"template 12",
      #"plot_type" : "migration"},
    #{"name":"template 13",
      #"plot_type" : "migration"},
    #{"name":"template 14",
      #"plot_type" : "migration"},
    ##{"name":"template 15",
      #"plot_type" : "migration"},
    #{"name":"template 20",
      #"plot_type" : "migration"},
    ##{"name":"template 30",
      #"plot_type" : "migration"},
    #{"name":"template 50",
      #"plot_type" : "migration"},
    #{"name":"nue_EL",
    #    "plot_type" : "migration"},
    #{"name":"electron_energy",
    #    "plot_type" : "sigdep"},
    #{"name":"Reco Energy vs L/E",
    #    "plot_type" : "migration"},
    #{"name":"E Theta Squared"},
    #{"name":"True Energy vs Biased Neutrino Energy",
    #        "plot_type" : "category_hist"},
    {"name":"Biased Neutrino Energy"},
    #{"name":"Biased Neutrino Energy",
    #     "plot_type" : "sigdepratio"},
    #{"name":"Biased Neutrino Energy",
    #     "plot_type" : "err"},
    #{"name":"Delta Phi"},
    #{"name":"Q2"},
    #{"name":"Front dEdX"},
    #{"name":"Front dEdX",
    #    "plot_type" : "sigdepratio"},

    #{"name":"Biased Neutrino Energy",
    #     "plot_type" : "bkgratio"},
    #{"name":"Biased Neutrino Energy",
    #     "plot_type" : "sigdepratio"},
    #{"name":"Front dEdX"},
    #{"name":"Front dEdX",
    #     "plot_type" : "bkgratio"},
    #{"name":"Front dEdX",
    #     "plot_type" : "sigdepratio"},
    #{"name":"E Theta Squared",
    #     "plot_type" : "sigdepratio"},
    #{"name":"Biased Neutrino 4 Energy vs E Theta Squared",
    #     "plot_type" : "category_hist"},
    #{"name":"flux",
    #     "plot_type":"err"},
    ]


#for i in PlotConfig.LOW_RECOIL_BIN_Q0:
    #PLOTS_TO_MAKE.append({"name":"Eavail bin "+str(i),"plot_type":"migration"})
    #PLOTS_TO_MAKE.append({"name":"E Estimator Eavail bin "+str(i),"plot_type":"stacked"})
    #PLOTS_TO_MAKE.append({"name":"E Estimator Eavail bin "+str(i),"plot_type":"ratio"})
    #PLOTS_TO_MAKE.append({"name":"E Estimator Eavail bin "+str(i),"plot_type":"err"})
