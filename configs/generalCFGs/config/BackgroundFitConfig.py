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
#        "CCNuEQE":"Signal",
#        "CCNuEDelta":"Signal",
#        "CCNuEDIS":"Signal",
#        "CCNuE2p2h":"Signal",
#        "CCNuE":"Signal",
        "CCPi0":"dEdX",
        "NCCohPi0":"dEdX",
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




#make sure all groups have a color and defination
#assert set(BACKGROUND_GROUP.keys())==set(BACKGROUND_GROUP_COLOR.keys())
