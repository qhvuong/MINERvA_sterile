"""
  PlotLibrary.py:
   Define plots to be handled by MakePlotProcessors.
   Configure by PlotConfig

   Original author: J. Wolcott (jwolcott@fnal.gov)
                    May 2014
"""
import ROOT
import math
from config import PlotConfig,SignalDef,SystematicsConfig
from tools import TruthTools,Utilities


#tags are treated in MyHistograms.MakePlotProcessors
# sideband: make the plot in sideband sample as well. Default is only signal sample
# truth_class : make the plot for different truth_categroy seperately, only works for mc of course
# mc_only: only make this plot for mc sample
# signal_only: only make this plot for true signal, only works for mc of course.

#tags for histograms that will be used to produce cross section.
reco_tags={"sideband","truth_class"}
migration_tags={"mc_only","signal_only"} 
signal_tags={"sideband","mc_only","signal_only"} 
truth_signal_tags={"mc_only","signal_only","ignore_selection","truth_class"}

#tags for histograms for various study
resolution_tags={"sideband","mc_only"}
truth_tags={"sideband","truth_class","mc_only"}

skip_sys = lambda universe: universe.ShortName() == "cv"

"""
How to define a new plot:
1) add a entry to PLOT_SETTINGS. format:
"name of plot":
{
"name" : the name of mnvHXD,
"title": the title of mnvHXD, "title;x-axis-lable;y-axis-lable",
"binning": the binning(s) of mnvHXD. [ binning for x-axis, binning for y-axis if 2D histogram ]
           binning is a list of bin edges for example: [0,1,2,] means two bins, 0-1 and 1-2.
"value_getter" : the function to get variable(s) to be filled to mnvHXD. [lambda function for x-variable, lambda function for y-variable if 2D,]
"tags" : giving tags to plot such that it will be treated specially in MakePlotProcessors.
         available tags at the time of writing are:
         "sideband" : make this plot for sideband region ( default is for signal region only)
         "mc_only" : only make this plot for mc samples.
         "truth_class" : make variants of this plot for each truth_class defined in SignalDef.py
"cuts" : list of lambda functions decide if a entry will be fill per universe/event. can be used to make additional cut for individual plot rather all plots.
}

2) add the name of this plot to HIST_TO_MAKE in PlotConfig.py

"""



DecMap = {
    "low" : lambda func, binning: Utilities.decorator_ReLU(func,binning[0]),
    "high" : lambda func, binning: Utilities.decorator_Cap(func,binning[-2]),
}

def VariantPlotsNamingScheme(*args):
    return "_".join(args)

def TranslateSettings(key): 
    if isinstance(key,str):
        settings = PLOT_SETTINGS[key].copy()  
    elif isinstance(key,dict):
        #make a new settings out of the dict 
        settings = key.copy() 
        for _ in settings["variables"]:
            var = VARIABLE_DICT[_] 
            settings.setdefault("binning",[]).append(var["binning"])
            tmp = var["value_getter"]
            if "dec" in var:
                for key in var["dec"]:
                    tmp = DecMap[key](tmp,var["binning"])

            settings.setdefault("value_getter",[]).append(tmp)
        settings["title"] = ";"+";".join(VARIABLE_DICT[_]["title"] for _ in settings["variables"])
        settings["name"] = "_".join(VARIABLE_DICT[_]["name"] for _ in settings["variables"])
        del settings["variables"]

    return settings

def passHybridProtonNodeCut(event):
    means=[]
    means.append(31.302)
    means.append(11.418)
    means.append(9.769)
    means.append(8.675)
    means.append(7.949)
    sigmas=[]
    sigmas.append(8.997)
    sigmas.append(3.075)
    sigmas.append(2.554)
    sigmas.append(2.484)
    sigmas.append(2.232)

    nodeEnergyVal=0
    chi2=0
    n_nodes=event.MasterAnaDev_proton_nodes_nodesNormE_sz
    if n_nodes>5:
        for i in range(n_nodes):
            if i==6:
                break;
            if i==0:
                nodeEnergyVal+=event.MasterAnaDev_proton_nodes_nodesNormE[0]
            elif i==1:
                nodeEnergyVal+=event.MasterAnaDev_proton_nodes_nodesNormE[1]
            else:
                nodeEnergyVal=event.MasterAnaDev_proton_nodes_nodesNormE[i]
            if i>=1:
                chi2+=(nodeEnergyVal-means[i-1])*(nodeEnergyVal-means[i-1])/(sigmas[i-1]*sigmas[i-1])
    else:
        passes=passPrimaryProtonNodeCut(event)
        if passes:
            chi2=0;
        else:
            chi2=75
    return chi2

def passPrimaryProtonNodeCut(event):
    #CutValuesbasedon22302
    #Node0-1
    #Node2
    #Node3
    #Node4
    #Node5
    #Node6

    cutval1=19
    cutval2=10
    cutval3=9
    cutval4=8
    cutval5=5


    #Primaryproton
    vect_size=event.MasterAnaDev_proton_nodes_nodesNormE_sz

    if vect_size==0:
        return False#/nonodes
    if event.MasterAnaDev_proton_nodes_nodesNormE[0]+event.MasterAnaDev_proton_nodes_nodesNormE[1]<cutval1:
        return False
    if event.MasterAnaDev_proton_nodes_nodesNormE[2]<cutval2:
        return False
    if event.MasterAnaDev_proton_nodes_nodesNormE[3]<cutval3:
        return False
    if event.MasterAnaDev_proton_nodes_nodesNormE[4]<cutval4:
        return False
    if vect_size>5:
        if event.MasterAnaDev_proton_nodes_nodesNormE[5]<cutval5:
            return False

    #survivedloops?
    return True

def vertexDistance(event,n_prong=0):
    protonX = event.MasterAnaDev_proton_startPointX
    protonY = event.MasterAnaDev_proton_startPointY
    protonZ = event.MasterAnaDev_proton_startPointZ
    electronX = event.prong_axis_vertex[n_prong][0]
    electronY = event.prong_axis_vertex[n_prong][1]
    electronZ = event.prong_axis_vertex[n_prong][2]
    return(math.sqrt((protonX-electronX)**2 + (protonY-electronY)**2 + (protonZ-electronZ)**2))

def vertexDifference(event,n_prong=0):
    protonZ = event.MasterAnaDev_proton_startPointZ
    electronZ = event.prong_axis_vertex[n_prong][2]
    return(electronZ - protonZ)

def CalcApothem(x,y):
    x=abs(x)
    y=abs(y)
    if ( x == 0 or y/x > 1/math.sqrt(3)):
        return (y+x/math.sqrt(3))/2*math.sqrt(3)
    else:
        return x

VARIABLE_DICT = {
    "Biased Neutrino Energy":
    {
        "name" : "EN4",
        "title" : "E_{e}+E_{avail} (GeV)",
        #"binning" : PlotConfig.NEUTRINO_ENERGY_BINNING,
        "binning" : PlotConfig.NEUTRINO4_EE_BINNING,
        "value_getter" : lambda event: event.kin_cal.reco_E_lep+event.kin_cal.reco_visE,
    },
    "Neutrino Length Travelled":
    {
        "name" : "nu_length",
        "title" : "Length_{#nu} (km)",
        "binning" : PlotConfig.NEUTRINO4_LENGTH_BINNING,
        "value_getter" : lambda event: .9825+event.mc_vtx[2]/1e6 - event.mc_fr_nuParentDecVtx[2]/1e6,
        "tags": "mc_only"
    },
    "Visible Energy":
    {
        "name" : "Eavail",
        "title" : "E_{avail} (GeV)",
        "binning" : PlotConfig.LOW_RECOIL_BIN_Q0,
        "value_getter" : lambda event: event.kin_cal.reco_visE,
    },
    "Lepton Energy":
    {
        "name" : "Eel",
        "title" : "Reconstructed E_lep (GeV)",
        "binning" : PlotConfig.NEUTRINO4_EE_BINNING,
        "value_getter" : lambda event: event.kin_cal.reco_E_lep,
    },
     "Q3" :
    {
        "name" : "Q3",
        "title": "q3 (GeV)",
        "binning" : PlotConfig.LOW_RECOIL_BIN_Q3,
        "value_getter" : lambda event: event.kin_cal.reco_q3,
    },
    "Lepton Pt":
    {
        "name" : "Lepton_Pt",
        "title" : "Pt_lepton (GeV); NEvents",
        "binning" : PlotConfig.PT_BINNING,
        "value_getter" : lambda event: event.kin_cal.reco_P_lep*math.sin(event.kin_cal.reco_theta_lep_rad),
    },
    "Nu Parent Energy":
    {
        "name" : "nuParent_energy",
        "title" : "Parent Energy (GeV); NEvents",
        "binning" : PlotConfig.NEUTRINO4_EE_BINNING,
        "value_getter" : lambda event: event.mc_fr_nuParentProdP[3]/1000,
    },

}

PLOT_SETTINGS= {
    "Biased Neutrino Energy":
    {
        "name" : "EN4",
        "title" : "E_{e}+E_{avail} (GeV)",
        #"binning" : PlotConfig.NEUTRINO_ENERGY_BINNING,
        "binning" : [PlotConfig.NEUTRINO4_EE_BINNING],
        "value_getter" : [lambda event: event.kin_cal.reco_E_lep+event.kin_cal.reco_visE],
        "tags":reco_tags
    },
    "True Lepton Energy":
    {
        "name" : "true_Eel",
        "title" : "True E_lep (GeV)",
        "binning" : [PlotConfig.NEUTRINO4_EE_BINNING],
        "value_getter" : [lambda event: event.kin_cal.true_E_lep],
        "tags":truth_tags
    },
    "Lepton Pt":
    {
        "name" : "Lepton_Pt",
        "title" : "Pt_lepton (GeV); NEvents",
        "binning" : [PlotConfig.PT_BINNING],
        "value_getter" : [lambda event: event.kin_cal.reco_P_lep*math.sin(event.kin_cal.reco_theta_lep_rad)],
        "tags":reco_tags
    },
    "True Signal Biased Neutrino Energy":
    {
        "name" : "true_EN4_true_signal",
        "title" : "E_{e}+E_{avail} (GeV)",
       #"binning" : PlotConfig.NEUTRINO_ENERGY_BINNING,
        "binning" : [PlotConfig.NEUTRINO4_EE_BINNING],
        "value_getter" : [lambda event: event.mc_incomingE/1000],
        "tags" : truth_signal_tags
    },
    "Signal Biased Neutrino Energy":
    {
        "name" : "EN4_true_signal",
        "title" : "E_{e}+E_{avail} (GeV)",
        #"binning" : PlotConfig.NEUTRINO_ENERGY_BINNING,
        "binning" : [PlotConfig.NEUTRINO4_EE_BINNING],
        "value_getter" : [lambda event: event.kin_cal.reco_E_lep+event.kin_cal.reco_visE],
        "tags" : signal_tags
    },
    "Lepton Energy":
    {
        "name" : "Eel",
        "title" : "Reconstructed E_lep (GeV)",
        "binning" : PlotConfig.NEUTRINO4_EE_BINNING,
        "value_getter" : lambda event: event.kin_cal.reco_E_lep,
    },
     "Front dEdX":
    {
        "name": "frontdedx",
        "title": "Mean front dE/dx; dEdX (MeV/cm); NEvents",
        "binning": [PlotConfig.DEDX_BINNING],
        "value_getter" : [lambda event: event.prong_dEdXMeanFrontTracker[0]],
        "tags":reco_tags
    },
    "True Signal Lepton Pt CCQE":
    {
        "name" : "Lepton_Pt_CCQE",

        "title": "Lepton Pt; Lepton Pt (GeV); NEvents",
        "binning" : [PlotConfig.PT_BINNING],
        "value_getter" : [lambda event: TruthTools.TransverseMomentum(event,13,True)],
        "tags": truth_signal_tags
    },
    "True Signal Lepton Pt Aaron":
    {
        "name" : "Lepton_Pt_Aaron",

        "title": "Lepton Pt; Lepton Pt (GeV); NEvents",
        "binning" : [PlotConfig.PT_BINNING_AARON],
        "value_getter" : [lambda event: TruthTools.TransverseMomentum(event,13,True)],
        "tags": truth_signal_tags
    },
    "Delta Phi":
    {
        "name" : "dphi",
        "title" : "Muon Phi - Proton Phi",
        "binning" : [[i * 0.1 for i in range(-60,60)]],
        "value_getter" : [lambda event: event.muon_phi - event.MasterAnaDev_proton_phi],
        "tags": truth_tags
    },
    "Q2" :
    {
        "name" : "Q2",
        "title": "q2 ; Q2 (GeV^{2}); dNEvents/dq2",
        "binning" : [[i/10 for i in range(40)]],
        "value_getter" : [lambda event: event.mc_Q2/1e6],
        "tags": truth_tags
    },
    "Q3" :
    {
        "name" : "Q3",
        "title": "q3 ; q3 (GeV); dNEvents/dq3",
        "binning" : [PlotConfig.LOW_RECOIL_BIN_Q3],
        "value_getter" : [lambda event: event.kin_cal.reco_q3],
        "tags": reco_tags
    },
    "E Theta Squared":
    {
        "name" : "E_Theta_Squared",
        "title" : "E_{lepton} #theta^{2} ; E_{lepton} #theta^{2} (GeV) ; NEvents",
        "binning" : [PlotConfig.NEUTRINO4_EE_THETA_BINNING],
        "value_getter" : [lambda event: event.kin_cal.reco_E_lep*(event.kin_cal.reco_theta_lep_rad)**2],
        "tags":reco_tags
    },
    "Lepton Angle":
    {
        "name" : "Lep_Angle",
        "title" : "Lepton Angle ; Lepton Angle (deg) ; NEvents",
        "binning" : [PlotConfig.LEPTON_ANGLE_BINNING],
        "value_getter" : [lambda event: event.kin_cal.reco_theta_lep],
        "tags":reco_tags
    },
     "Estimator vs Front dEdX":
    {
        "name" : "estimator_vs_frontdedx",
        "title" : "Energy Estimator vs dE/dx; Energy_{Estimator} (GeV); Mean Front dE/dx (MeV/cm); NEvents",
        "binning" : [PlotConfig.NEUTRINO4_EE_BINNING,PlotConfig.DEDX_BINNING],
        "value_getter" : [lambda event: event.kin_cal.reco_E_lep+event.kin_cal.reco_visE,lambda event: event.prong_dEdXMeanFrontTracker[0]],
        "tags":reco_tags
    },
     "Estimator vs Available Energy":
    {
        "name" : "estimator_vs_eavail",
        "title" : "Energy Estimator vs Available Energy; Energy_{Estimator} (GeV); Energy_{avail} (GeV); NEvents",
        "binning" : [PlotConfig.NEUTRINO4_EE_BINNING,PlotConfig.LOW_RECOIL_BIN_Q0],
        "value_getter" : [lambda event: event.kin_cal.reco_E_lep+event.kin_cal.reco_visE,lambda event: event.kin_cal.reco_visE],
        "tags":reco_tags
    },
    "Biased Neutrino 4 Energy vs E Theta Squared":
    {
        "name" : "EN4_Ethetasquared",
        "title" : "Neutrino 4 Energy v.s. E #theta^2; E_{lep}+E_{avail} (GeV); E #theta^2 (GeV); d^{2}NEvents/dq3dE_{e}",
        "binning" : [PlotConfig.NEUTRINO4_EE_BINNING,
                     PlotConfig.NEUTRINO4_EE_THETA_BINNING],
        "value_getter" : [lambda event: event.kin_cal.reco_E_lep+event.kin_cal.reco_visE,
                          lambda event: event.kin_cal.reco_E_lep*(event.kin_cal.reco_theta_lep_rad)**2],
        "tags":reco_tags
    },
    "Biased Neutrino 4 Energy vs q3":
    {
        "name" : "EN4_q3",
        "title" : "Neutrino 4 Energy v.s. q3; E_{lep}+E_{avail} (GeV); q3 (GeV); d^{2}NEvents/dq3dE_{e}",
        "binning" : [PlotConfig.NEUTRINO4_EE_BINNING,
                     PlotConfig.BACKGROUND_FIT_Q3_BIN],
        "value_getter" : [lambda event: event.kin_cal.reco_E_lep+event.kin_cal.reco_visE,
                          lambda event: event.kin_cal.reco_q3],
        "tags":reco_tags
    },
    "Estimator vs Lepton Pt":
    {
        "name" : "estimator_Lepton_Pt",

        "title": "Neutrino 4 Energy v.s. Lepton p_{t}; E_{lep}+E_{avail} (GeV); Lepton p_{t} (GeV); NEvents",
        "binning" : [PlotConfig.NEUTRINO4_EE_BINNING,
                     PlotConfig.PT_BINNING],
        "value_getter" : [lambda event: event.kin_cal.reco_E_lep+event.kin_cal.reco_visE,
                          lambda event: event.kin_cal.reco_P_lep*math.sin(event.kin_cal.reco_theta_lep_rad)],
        "tags": reco_tags   
    },
    "EN4 vs Lepton Pt":
    {
        "name" : "EN4_Lepton_Pt",

        "title": "Neutrino 4 Energy v.s. Lepton p_{t}; E_{lep}+E_{avail} (GeV); Lepton p_{t} (GeV); NEvents",
        "binning" : [PlotConfig.NEUTRINO4_EE_BINNING,
                     PlotConfig.PT_BINNING],
        "value_getter" : [lambda event: event.kin_cal.reco_E_lep+event.kin_cal.reco_visE,
                          lambda event: event.kin_cal.reco_P_lep*math.sin(event.kin_cal.reco_theta_lep_rad)],
        "tags": reco_tags   
    },
    "Lepton Pt vs Lepton Energy":
    {
        "name" : "Leton_Pt_Elepton",

        "title": "Lepton Pt v.s. Lepton Energy; Lepton Pt (GeV); E_{lepton} (GeV); NEvents",
        "binning" : [PlotConfig.PT_BINNING,
                     PlotConfig.NEUTRINO4_EE_BINNING],
        "value_getter" : [lambda event: event.kin_cal.reco_P_lep*math.sin(event.kin_cal.reco_theta_lep_rad),
                          lambda event: event.kin_cal.reco_E_lep],
        "tags": reco_tags   
    },
    "Lepton Pt vs Available Energy":
    {
        "name" : "Leton_Pt_Eavail",

        "title": "Lepton Pt v.s. Available Energy; Lepton Pt (GeV); E_{avail} (GeV); NEvents",
        "binning" : [PlotConfig.PT_BINNING,
                     PlotConfig.LOW_RECOIL_BIN_Q0],
        "value_getter" : [lambda event: event.kin_cal.reco_P_lep*math.sin(event.kin_cal.reco_theta_lep_rad),
                          lambda event: event.kin_cal.reco_visE],
        "tags": reco_tags   
    },
    "Lepton Energy vs Available Energy":
    {
        "name" : "Elepton_Eavail",

        "title": "Lepton Energy v.s Available Energy; E_{lepton} (GeV); E_{avail} (GeV); NEvents",
        "binning" : [PlotConfig.NEUTRINO4_EE_BINNING,
                     PlotConfig.LOW_RECOIL_BIN_Q0],
        "value_getter" : [lambda event: event.kin_cal.reco_E_lep,
                          lambda event: event.kin_cal.reco_visE],
        "tags": reco_tags   
    },
    "Available Energy vs Lepton Energy":
    {
        "name" : "Eavail_Elepton",

        "title": "Available Energy v.s. Lepton Energy; E_{avail} (GeV); E_{lepton} (GeV); NEvents",
        "binning" : [PlotConfig.LOW_RECOIL_BIN_Q0,
                     PlotConfig.NEUTRINO4_EE_BINNING],
        "value_getter" : [lambda event: event.kin_cal.reco_visE,
                          lambda event: event.kin_cal.reco_E_lep],
        "tags": reco_tags   
    },
    "Available Energy vs Lepton Pt":
    {
        "name" : "Eavail_Lepton_Pt",

        "title": "Available Energy v.s. Lepton Pt; E_{avail} (GeV); Lepton Pt (GeV); NEvents",
        "binning" : [PlotConfig.LOW_RECOIL_BIN_Q0,
                     PlotConfig.PT_BINNING],
        "value_getter" : [lambda event: event.kin_cal.reco_visE,
                          lambda event: event.kin_cal.reco_P_lep*math.sin(event.kin_cal.reco_theta_lep_rad)],
        "tags": reco_tags   
    },
    "True Lepton Energy vs Available Energy":
    {
        "name" : "TrueELepton_Eavail",

        "title": "True Lepton Energy v.s. Available Energy; E_{lep} (GeV); E_{avail} (GeV); NEvents",
        "binning" : [PlotConfig.NEUTRINO4_EE_BINNING,
                     PlotConfig.LOW_RECOIL_BIN_Q0],
        "value_getter" : [lambda event: event.kin_cal.true_E_lep,
                          lambda event: event.kin_cal.reco_visE],
        "tags": reco_tags   
    },
    "True Energy vs Biased Neutrino Energy":
    {
        "name" : "ETrue_EReco",

        "title": "E_{l} + E_{avail} vs True E_{#nu}; E_{l} + E_{avail} (GeV); True E_{#nu} (GeV); NEvents",
        "binning" : [PlotConfig.NEUTRINO4_EE_BINNING,
                     PlotConfig.NEUTRINO4_EE_BINNING],
        "value_getter" : [lambda event: event.kin_cal.reco_E_lep+event.kin_cal.reco_visE,
                          lambda event: event.mc_incomingE/1000],
        "tags": truth_tags   
    },
    "Nu Parent Energy vs Length Travelled":
    {
        "name" : "nuParent_energy_vs_Length",
        "title": "True E v.s. True Length; Length (km); True E (GeV); NEvents",
        "binning" : [PlotConfig.NEUTRINO4_LENGTH_BINNING,
                     PlotConfig.NEUTRINO4_EE_BINNING],
        "value_getter" : [lambda event: .9825+event.mc_vtx[2]/1e6 - event.mc_fr_nuParentDecVtx[2]/1e6,
                          lambda event: event.mc_fr_nuParentProdP[3]/1000],
        "tags": truth_tags   
    },
    "Nu Parent Energy vs Nu Energy":
    {
        "name" : "nuParent_energy_vs_nuEnergy",
        "title": "True Parent Energy vs True Neutrino Energy; Nu Parent Energy (GeV); Nu Energy (GeV); NEvents",
        "binning" : [PlotConfig.NEUTRINO4_EE_BINNING,
                     PlotConfig.NEUTRINO4_EE_BINNING],
        "value_getter" : [lambda event: event.mc_incomingE/1000,
                          lambda event: event.mc_fr_nuParentProdP[3]/1000],
        "tags": truth_tags   
    },
    "Reco Energy vs Longitudinal Distance":
    {
        "name" : "recoE_vs_longDist",
        "title": "Energy Estimator vs Longitudinal Distance; E_{avail}+E_{lep} (GeV); Distance (cm); NEvents",
        "binning" : [[i for i in range(0,500,500)],
                     PlotConfig.NEUTRINO4_EE_BINNING],
        "value_getter" : [lambda event: event.vtx[2]/1e1,
                          lambda event: event.kin_cal.reco_E_lep+event.kin_cal.reco_visE],
        "tags": truth_tags   
    },
    "Reco Energy vs Transverse Distance":
    {
        "name" : "recoE_vs_transDist",
        "title": "Energy Estimator vs Transverse Distance; E_{avail}+E_{lep} (GeV); Distance (cm); NEvents",
        "binning" : [[i for i in range(0,400,400)],
                     PlotConfig.NEUTRINO4_EE_BINNING],
        "value_getter" : [lambda event: math.sqrt((event.vtx[0]/1e1)**2+(event.vtx[1]/1e1)**2),
                          lambda event: event.kin_cal.reco_E_lep+event.kin_cal.reco_visE],
        "tags": truth_tags   
    },
    "True Energy vs Neutrino Length Travelled":
    {
        "name" : "ETrue_Length",

        "title": "True E v.s. True Length; Length (km); True E (GeV); NEvents",
        "binning" : [PlotConfig.NEUTRINO4_LENGTH_BINNING,
                     PlotConfig.NEUTRINO4_EE_BINNING],
        "value_getter" : [lambda event: .9825+event.mc_vtx[2]/1e6 - event.mc_fr_nuParentDecVtx[2]/1e6,
                          lambda event: event.mc_incomingE/1000],
        "tags": truth_tags   
    },
    "Reco Energy vs L/E":
    {
        "name" : "EReco_LE",

        "title": "Reco E v.s. True L/E; True L/E (km/GeV); E_{e} + E_{avail} GeV; NEvents",
        "binning" : [PlotConfig.NEUTRINO4_LE_BINNING,
                     PlotConfig.NEUTRINO4_EE_BINNING],
        "value_getter" : [lambda event: (.9825+event.mc_vtx[2]/1e6 - event.mc_fr_nuParentDecVtx[2]/1e6)/(event.mc_incomingE/1000),
                          lambda event: event.kin_cal.reco_E_lep+event.kin_cal.reco_visE], 
        "tags": truth_tags   
    },
    "True Energy vs L/E":
    {
        "name" : "E_LE",

        "title": "True E v.s. True L/E; True L/E (km/GeV); E_{e} + E_{avail} GeV; NEvents",
        "binning" : [PlotConfig.NEUTRINO4_LE_BINNING,
                     PlotConfig.NEUTRINO4_EE_BINNING],
        "value_getter" : [lambda event: (.9825+event.mc_vtx[2]/1e6 - event.mc_fr_nuParentDecVtx[2]/1e6)/(event.mc_incomingE/1000),
                          lambda event: event.mc_incomingE/1000], 
        "tags": truth_tags   
    },
    "Signal Reco Energy vs L/E":
    {
        "name" : "true_EReco_LE",

        "title": "Reco E v.s. True L/E; True L/E (km/GeV); E_{e} + E_{avail} GeV; NEvents",
        "binning" : [PlotConfig.NEUTRINO4_LE_BINNING,
                     PlotConfig.NEUTRINO4_EE_BINNING],
        "value_getter" : [lambda event: (.9825+event.mc_vtx[2]/1e6 - event.mc_fr_nuParentDecVtx[2]/1e6)/(event.mc_incomingE/1000),
                          lambda event: event.kin_cal.reco_E_lep+event.kin_cal.reco_visE], 
        "tags": signal_tags
    },
    "Signal Reco Energy vs sin":
    {
        "name" : "true_EReco_SIN",

        "title": "Reco E v.s. sin^{2}(1.27 \Delta m^{2} L/E); True sin^{2}(1.27 \Delta m^{2} L/E); E_{e} + E_{avail} GeV; NEvents",
        "binning" : [PlotConfig.NEUTRINO4_P_BINNING,
                     PlotConfig.NEUTRINO4_EE_BINNING],
        "value_getter" : [lambda event: math.sin(1.27 * 7.34 * (.9825+event.mc_vtx[2]/1e6 - event.mc_fr_nuParentDecVtx[2]/1e6)/(event.mc_incomingE/1000))**2,
                          lambda event: event.kin_cal.reco_E_lep+event.kin_cal.reco_visE], 
        "tags": signal_tags
    },
    "Predicted MC":
    {
        "name" : "EN4_predicted_Signal",
        "title" : "E_{e}+E_{avail} (GeV)",
        #"binning" : PlotConfig.NEUTRINO_ENERGY_BINNING,
        "binning" : [PlotConfig.NEUTRINO4_EE_BINNING],
        "value_getter" : [lambda event: event.kin_cal.reco_E_lep+event.kin_cal.reco_visE],
        "tags":reco_tags
    },
    "Background Subbed Data":
    {
        "name" : "EN4_ERROR_data_bkgSubbed",
        "title" : "E_{e}+E_{avail} (GeV)",
        #"binning" : PlotConfig.NEUTRINO_ENERGY_BINNING,
        "binning" : [PlotConfig.NEUTRINO4_EE_BINNING],
        "value_getter" : [lambda event: event.kin_cal.reco_E_lep+event.kin_cal.reco_visE],
        "tags":reco_tags
    },
    "Proton Chi2":
    {
        "name" : "proton_chi2",
        "title" : "TKI Proton Chi2; chi2; NEvents",
        #"binning" : PlotConfig.NEUTRINO_ENERGY_BINNING,
        "binning" : [[i for i in range(50)]],
        "value_getter" : [lambda event: passHybridProtonNodeCut(event)],
        "tags":reco_tags
    },
    "Proton Chi2 vs Lepton Energy":
    {
        "name" : "protonchi2_ELepton",
        "title" : "chi2",
        #"binning" : PlotConfig.NEUTRINO_ENERGY_BINNING,
        "binning" : [[i for i in range(50)],
                     PlotConfig.NEUTRINO4_EE_BINNING],
        "value_getter" : [lambda event: passHybridProtonNodeCut(event),
                          lambda event: event.kin_cal.reco_E_lep],
        "tags":reco_tags
    },
    "Electron Proton Distance vs Electron Vertex Distance":
    {
        "name" : "electron_proton_vs_electron_vertex",
        "title" : "Prong Vertex Distances; Z_{electron} - Z_{proton} [ 5 mm ]; Z_{electron} - Z_{reco vertex} [ 5 mm ]; NEvents",
        #"binning" : PlotConfig.NEUTRINO_ENERGY_BINNING,
        "binning" : [[i/10 for i in range(-500,305,5)],
                     [i/10 for i in range(-200,205,5)]],
        "value_getter" : [lambda event: vertexDifference(event),
                          lambda event: event.prong_axis_vertex[0][2] - event.vtx[2]],
        "tags":reco_tags
    },
    "Prong Distance":
    {
        "name" : "proton_electron_distance",
        "title" : "distance between proton and electron candidate; distance [ mm ]; NEvents",
        "binning" : [[i for i in range(300)]],
        "value_getter" : [lambda event: vertexDistance(event)],
        "tags":reco_tags
    },
    "Prong Z Difference":
    {
        "name" : "proton_electron_z_difference",
        "title" : "z-difference between proton and electron candidate; Z_{proton} - Z_{electron} [ mm ]; NEvents",
        "value_getter" : [lambda event: vertexDifference(event)],
        "tags":reco_tags
    },
    "Electron Vertex Z Difference":
    {
        "name" : "electron_vertex_z_difference",
        "title" : "z-difference between electron and vertex; Z_{electron} - Z_{vertex} [ mm ]; NEvents",
        "binning" : [[i for i in range(-300,300,10)]],
        "value_getter" : [lambda event: event.prong_axis_vertex[0][2] - event.vtx[2]],
        "tags":reco_tags
    },
    "True Electron Vertex Z Difference":
    {
        "name" : "true_electron_vertex_z_difference",
        "title" : "z-difference between electron and true vertex; Z_{electron} - Z_{true #nu vertex} [ mm ]; NEvents",
        "binning" : [[i for i in range(-300,300,10)]],
        "value_getter" : [lambda event: event.prong_axis_vertex[0][2] - event.mc_vtx[2]],
        "tags":truth_tags
    },
    "Proton Vertex Z Difference":
    {
        "name" : "proton_vertex_z_difference",
        "title" : "z-difference between proton and vertex; Z_{proton} - Z_{vertex} [ mm ]; NEvents",
        "binning" : [[i for i in range(-300,300,10)]],
        "value_getter" : [lambda event: event.MasterAnaDev_proton_startPointZ - event.vtx[2]],
        "tags":reco_tags
    },
    "Proton End Z":
    {
        "name" : "proton_end_z",
        "title" : "z-distribution of proton end; Z_{proton} [ mm ]; NEvents",
        "binning" : [[i for i in range(4000,10000,10)]],
        "value_getter" : [lambda event: event.MasterAnaDev_proton_endPointZ],
        "tags":reco_tags
     },
    "Proton Start Z":
    {
        "name" : "proton_start_z",
        "title" : "z-distribution of proton start; Z_{proton} [ mm ]; NEvents",
        "binning" : [[i for i in range(4000,10000,10)]],
        "value_getter" : [lambda event: event.MasterAnaDev_proton_startPointZ],
        "tags":reco_tags
    },
    "Estimator vs Proton Length":
    {
        "name" : "estimator_vs_proton_length",
        "title" : "Energy Estimator vs Proton Length; Energy_{Estimator} (GeV); L_{proton} [ 100 mm ]; NEvents",
        "binning" : [PlotConfig.NEUTRINO4_EE_BINNING,[i for i in range(0,4000,100)]],
        "value_getter" : [lambda event: event.kin_cal.reco_E_lep+event.kin_cal.reco_visE,lambda event: event.proton_track_length],
        "tags":reco_tags
    },
    "True Proton Vertex Z Difference":
    {
        "name" : "true_proton_vertex_z_difference",
        "title" : "z-difference between proton and true vertex; Z_{proton} - Z_{true #nu vertex} [ mm ]; NEvents",
        "binning" : [[i for i in range(-300,300,10)]],
        "value_getter" : [lambda event: event.MasterAnaDev_proton_startPointZ - event.mc_vtx[2]],
        "tags":truth_tags
    },
    "Proton Electron Angle":
    {
        "name" : "proton_electron_angle",
        "title" : "angle between proton and electron; #theta_{e,p} [ degrees ]; NEvents",
        "binning" : [[i for i in range(0,180,10)]],
        "value_getter" : [lambda event: event.ElectronProtonAngle()],
        "tags":reco_tags
    },
    "Neutrino Z Vertex":
    {
        "name" : "neutrino_vertex_z",
        "title" : "Reco Neutrino Z Vertex; Z [ 10 mm ]; NEvents",
        "binning" : [[i for i in range(4000,8500,10)]],
        "value_getter" : [lambda event: event.vtx[2]],
        "tags":reco_tags
    },
    "nue_EL":
    {
        "name" : "nue_scattering_template",
        "title" : "#nu + e Scattering; True L/E; Electron Energy; NEvents",
        "binning" : [PlotConfig.NEUTRINO4_LE_BINNING,PlotConfig.NUE_SCATTERING_BINNING],
        "value_getter" : [lambda event: (.9825+event.mc_vtx[2]/1e6 - event.mc_fr_nuParentDecVtx[2]/1e6)/(event.mc_incomingE/1000),
            lambda event: event.kin_cal.true_E_lep],
        "tags":reco_tags
    },
    "Available Energy vs True W":
    {
        "name" : "Eavail_trueW",
        "title": "Reco Available Energy v.s. True W; E_{avail} (GeV); True W (GeV); NEvents",
        "binning" : [PlotConfig.LOW_RECOIL_BIN_Q0,
                     PlotConfig.W_BINNING],
        "value_getter" : [lambda event: event.kin_cal.reco_visE,
                          lambda event: event.mc_w/1e3],
        "tags": truth_tags,
    },
    "electron_energy":
    {
        "name" : "electron_energy",
        "title" : "#nu + e Scattering; Electron Energy; NEvents",
        "binning" : [PlotConfig.NUE_SCATTERING_BINNING],
        "value_getter" : [lambda event: event.kin_cal.true_E_lep],
        "tags":reco_tags
    },
    "electron_energy_nue":
    {
        "name" : "electron_energy_nue",
        "title" : "#nu + e Scattering; Electron Energy; NEvents",
        "binning" : [PlotConfig.NUE_SCATTERING_BINNING],
        "value_getter" : [lambda event: event.kin_cal.true_E_lep],
        "tags":reco_tags
    },
    "electron_energy_numu":
    {
        "name" : "electron_energy_numu",
        "title" : "#nu + e Scattering; Electron Energy; NEvents",
        "binning" : [PlotConfig.NUE_SCATTERING_BINNING],
        "value_getter" : [lambda event: event.kin_cal.true_E_lep],
        "tags":reco_tags
    },
    "electron_energy_anue":
    {
        "name" : "electron_energy_anue",
        "title" : "#nu + e Scattering; Electron Energy; NEvents",
        "binning" : [PlotConfig.NUE_SCATTERING_BINNING],
        "value_getter" : [lambda event: event.kin_cal.true_E_lep],
        "tags":reco_tags
    },
    "electron_energy_anumu":
    {
        "name" : "electron_energy_anumu",
        "title" : "#nu + e Scattering; Electron Energy; NEvents",
        "binning" : [PlotConfig.NUE_SCATTERING_BINNING],
        "value_getter" : [lambda event: event.kin_cal.true_E_lep],
        "tags":reco_tags
    },
    "flux":
    {
        "name" : "Flux",
        "title" : "#nu Flux; True Neutrino Energy; NEvents",
        "binning" : [PlotConfig.NEUTRINO4_EE_BINNING],
        "value_getter" : [lambda event: event.mc_incomingE/1000],
        "tags":reco_tags
    },
    "Visible Energy":
    {
        "name" : "E_avail",
        "title" : "Available Energy; E_{avail} (GeV); NEvents",
        "binning" : [PlotConfig.LOW_RECOIL_BIN_Q0],
        "value_getter" : [lambda event: event.kin_cal.reco_visE],
        "tags":reco_tags
    },
    "Lepton Energy":
    {
        "name" : "Eel",
        "title" : "Reconstructed E_lep (GeV)",
        "binning" : [PlotConfig.NEUTRINO4_EE_BINNING],
        "value_getter" : [lambda event: event.kin_cal.reco_E_lep],
        "tags":reco_tags
    },
    "True vs Reconstructed Lepton Energy":
    {
        "name" : "TrueVsRecoLeptonEnergy",

        "title": "True vs Reco Lepton Energy; True E_{lep} (GeV); Reco E_{lep} (GeV); NEvents",
        "binning" : [PlotConfig.NEUTRINO4_EE_BINNING,
                     PlotConfig.NEUTRINO4_EE_BINNING],
        "value_getter" : [lambda event: event.kin_cal.true_E_lep,
                          lambda event: event.kin_cal.reco_E_lep],
        "tags": reco_tags   
    },
    "True vs Reconstructed Lepton Theta":
    {
        "name" : "TrueVsRecoLeptonTheta",

        "title": "True vs Reco Lepton Theta; True theta_{lep} (rad); Reco theta_{lep} (rad); NEvents",
        "binning" : [PlotConfig.ELECTRON_ANGLE_BINNING,
                     PlotConfig.ELECTRON_ANGLE_BINNING],
        "value_getter" : [lambda event: event.kin_cal.true_theta_lep_rad,
                          lambda event: event.kin_cal.reco_theta_lep_rad],
        "tags": reco_tags   
    },
    "True vs Reconstructed Lepton Theta X":
    {
        "name" : "TrueVsRecoLeptonThetaX",

        "title": "True vs Reco Lepton Theta X; True thetaX (rad); Reco thetaX (rad); NEvents",
        "binning" : [PlotConfig.ELECTRON_ANGLE_BINNING,
                     PlotConfig.ELECTRON_ANGLE_BINNING],
        "value_getter" : [lambda event: event.kin_cal.true_thetaX_lep_rad,
                          lambda event: event.kin_cal.reco_thetaX_lep_rad],
        "tags": reco_tags   
    },
    "True vs Reconstructed Lepton Theta Y":
    {
        "name" : "TrueVsRecoLeptonThetaY",

        "title": "True vs Reco Lepton Theta Y; True thetaY (rad); Reco thetaY (rad); NEvents",
        "binning" : [PlotConfig.ELECTRON_ANGLE_BINNING,
                     PlotConfig.ELECTRON_ANGLE_BINNING],
        "value_getter" : [lambda event: event.kin_cal.true_thetaY_lep_rad,
                          lambda event: event.kin_cal.reco_thetaY_lep_rad],
        "tags": reco_tags   
    },
    "True vs Reconstructed Lepton Theta 2D":
    {
        "name" : "TrueVsRecoLeptonTheta2D",

        "title": "True vs Reco Lepton Theta 2D; True theta2D (rad); Reco theta2D (rad); NEvents",
        "binning" : [PlotConfig.ELECTRON_ANGLE_BINNING,
                     PlotConfig.ELECTRON_ANGLE_BINNING],
        "value_getter" : [lambda event: event.kin_cal.true_theta2D_lep_rad,
                          lambda event: event.kin_cal.reco_theta2D_lep_rad],
        "tags": reco_tags   
    },
    "True vs Reconstructed Lepton Phi":
    {
        "name" : "TrueVsRecoLeptonPhi",

        "title": "True vs Reco Lepton Phi; True phi (rad); Reco phi (rad); NEvents",
        "binning" : [PlotConfig.ELECTRON_PHI_BINNING,
                     PlotConfig.ELECTRON_PHI_BINNING],
        "value_getter" : [lambda event: event.kin_cal.true_phi_lep_rad,
                          lambda event: event.kin_cal.reco_phi_lep_rad],
        "tags": reco_tags   
    },
    "Reco ThetaX vs Vertex X":
    {
        "name" : "RecoThetaXVsVertexX",

        "title": "Reco ThetaX vs Vertex X; Vertex X (mm); Reco thetaX (rad); NEvents",
        "binning" : [PlotConfig.ELECTRON_VERTEXX_BINNING,
                     PlotConfig.ELECTRON_ANGLE_BINNING],
        "value_getter" : [lambda event: event.vtx[0],
                          lambda event: event.kin_cal.reco_thetaX_lep_rad],
        "tags": reco_tags   
    },
    "Reco ThetaY vs Vertex Y":
    {
        "name" : "RecoThetaYVsVertexY",

        "title": "Reco ThetaY vs Vertex Y; Vertex Y (mm); Reco thetaY (rad); NEvents",
        "binning" : [PlotConfig.ELECTRON_VERTEXY_BINNING,
                     PlotConfig.ELECTRON_ANGLE_BINNING],
        "value_getter" : [lambda event: event.vtx[1],
                          lambda event: event.kin_cal.reco_thetaY_lep_rad],
        "tags": reco_tags   
    },
    "True ThetaX vs Vertex X":
    {
        "name" : "TrueThetaXVsVertexX",

        "title": "True ThetaX vs Vertex X; Vertex X (mm); True thetaX (rad); NEvents",
        "binning" : [PlotConfig.ELECTRON_VERTEXX_BINNING,
                     PlotConfig.ELECTRON_ANGLE_BINNING],
        "value_getter" : [lambda event: event.mc_vtx[0],
                          lambda event: event.kin_cal.true_thetaX_lep_rad],
        "tags": truth_tags   
    },
    "True ThetaY vs Vertex Y":
    {
        "name" : "TrueThetaYVsVertexY",

        "title": "True ThetaY vs Vertex Y; Vertex Y (mm); True thetaY (rad); NEvents",
        "binning" : [PlotConfig.ELECTRON_VERTEXY_BINNING,
                     PlotConfig.ELECTRON_ANGLE_BINNING],
        "value_getter" : [lambda event: event.mc_vtx[1],
                          lambda event: event.kin_cal.true_thetaY_lep_rad],
        "tags": truth_tags   
    },
    "Reco Q2" :
    {
        "name" : "Q2_reco",
        "title": "Reco q2 ; Q2 (GeV^{2}); dNEvents/dq2",
        "binning" : [[i/100 for i in range(20)]],
        "value_getter" : [lambda event: event.kin_cal.reco_q2_cal],
        "tags": reco_tags
    },
    "True E Theta Squared":
    {
        "name" : "true_E_Theta_Squared",
        "title" : "E_{lepton} #theta^{2} ; E_{lepton} #theta^{2} (GeV) ; NEvents",
        "binning" : [PlotConfig.NEUTRINO4_EE_THETA_BINNING],
        "value_getter" : [lambda event: event.kin_cal.true_E_lep*(event.kin_cal.true_theta_lep_rad)**2],
        "tags":truth_tags
    },
    "Neutrino Vertex Apothem":
    {
        "name" : "neutrino_apothem",
        "title" : "Reco Neutrino Vertex Apothem; Vertex Apothem; NEvents",
        "binning" : [[30*i for i in range(40)]],
        "value_getter" : [lambda event: CalcApothem(event.vtx[0],event.vtx[1])],
        "tags":reco_tags
    },
    "Neutrino Vertex Z":
    {
        "name" : "neutrino_vertex_z",
        "title" : "Reco Neutrino Z Vertex; Z [ 10 mm ]; NEvents",
        "binning" : [[i for i in range(5500,9000,70)]],
        "value_getter" : [lambda event: event.vtx[2]],
        "tags":reco_tags
    },
    "EMLikeTrackScore":
    {
        "name" : "EMLikeScore",
        "title" : "EM-Like Score; EM-Like Score; NEvents",
        "binning" : [PlotConfig.EMLIKETRACKSCORE_BINNING],
        "value_getter" : [lambda event: event.prong_part_score],
        "tags":reco_tags
    },
    "TransverseGapScore":
    {
        "name" : "TransGapScore",
        "title" : "Transverse Gap Score; Transverse Gap Score; NEvents",
        "binning" : [PlotConfig.TRANSVERSEGAPSCORE_BINNING],
        "value_getter" : [lambda event: event.prong_TransverseGapScore],
        "tags":reco_tags
    },
    "NonMIPClusFrac":
    {
        "name" : "NonMIPClusFracion",
        "title" : "Non-MIP Cluster Energy Fraction; NonMIPClusFrac; NEvents",
        "binning" : [PlotConfig.NONMIPCLUSFRAC_BINNING],
        "value_getter" : [lambda event: event.prong_NonMIPClusFrac],
        "tags":reco_tags
    },
    "ODCalVisE": 
    {
        "name" : "ODCalVisE",
        "title" : "ODCalVisE; ODCalVisE; NEvents",
        "binning" : [PlotConfig.ODCALVISE_BINNING],
        "value_getter": [
            lambda event: [
                (od / side if side > 0.1 else (1e10 if od > 0 else 0))
                for od, side in zip(event.prong_ODVisE, event.prong_SideECALVisE)
            ]
        ],
        "tags":reco_tags
    },
    "DSCalVisE": 
    {
        "name" : "DSCalVisE",
        "title" : "DSCalVisE; DSCalVisE; NEvents",
        "binning" : [PlotConfig.DSCALVISE_BINNING],
        "value_getter": [
            lambda event: [
                (h / ecal if ecal > 0.1 else (1e10 if h > 0 else 0))
                for h, ecal in zip(event.prong_HCALVisE, event.prong_ECALVisE)
            ]
        ],
        "tags":reco_tags
    },
    "Afterpulsing": 
    {
        "name" : "Afterpulsing",
        "title" : "Afterpulsing; Afterpulsing; NEvents",
        "binning" : [PlotConfig.AFTERPULSING_BINNING],
        "value_getter": [lambda event: event.prong_FirstFireFraction],
        "tags":reco_tags
    },
    "DeadTime": 
    {
        "name" : "DeadTime",
        "title" : "DeadTime; DeadTime; NEvents",
        "binning" : [[-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5]],
        "value_getter": [lambda event: event.phys_n_dead_discr_pair_upstream_prim_track_proj],
        "tags":reco_tags
    },
    "VertexTrackMultiplicity": 
    {
        "name" : "VertexTrackMultiplicity",
        "title" : "VertexTrackMultiplicity; VertexTrackMultiplicity; NEvents",
        "binning" : [[-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5]],
        "value_getter": [lambda event: event.VertexTrackMultiplicity],
        "tags":reco_tags
    },
    "StartPointVertexMultiplicity": 
    {
        "name" : "StartPointVertexMultiplicity",
        "title" : "StartPointVertexMultiplicity; StartPointVertexMultiplicity; NEvents",
        "binning" : [[-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5]],
        "value_getter": [lambda event: event.StartPointVertexMultiplicity],
        "tags":reco_tags
    },
    "HasNoVertexMismatch": 
    {
        "name" : "HasNoVertexMismatch",
        "title" : "HasNoVertexMismatch; HasNoVertexMismatch; NEvents",
        "binning" : [[-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5]],
        "value_getter": [lambda event: event.HasNoVertexMismatch],
        "tags":reco_tags
    },
    "HasTracks": 
    {
        "name" : "HasTracks",
        "title" : "HasTracks; HasTracks; NEvents",
        "binning" : [[-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5]],
        "value_getter": [lambda event: event.HasTracks],
        "tags":reco_tags
    },
    "HasNoBackExitingTracks": 
    {
        "name" : "HasNoBackExitingTracks",
        "title" : "HasNoBackExitingTracks; HasNoBackExitingTracks; NEvents",
        "binning" : [[-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5]],
        "value_getter": [lambda event: event.HasNoBackExitingTracks],
        "tags":reco_tags
    },
}

for i in PlotConfig.LOW_RECOIL_BIN_Q0:
    key = "Eavail bin "+str(i)
    name = "Eavail_bin_"+str(i).replace(".","p")
    title = "Reco E vs L/E"
    binning = [PlotConfig.NEUTRINO4_EE_BINNING,PlotConfig.NEUTRINO4_L_OVER_E_BINNING]
    value_getter = [lambda event: event.kin_cal.reco_E_lep+event.kin_cal.reco_visE, 
                    lambda event: (.9825+event.mc_vtx[2]/1e6 - event.mc_fr_nuParentDecVtx[2]/1e6)/(event.mc_incomingE/1000)]
    tags = signal_tags
    cuts = [(lambda i=i: lambda event: event.kin_cal.reco_visE < i)()]
    entry = {"name" : name,
             "title": title,
             "binning": binning,
             "value_getter":value_getter,
             "tags": tags,
             "cuts": cuts}

    PLOT_SETTINGS[key] = entry

    key = "E Estimator Eavail bin "+str(i)
    name = "E_Estimator_Eavail_bin_"+str(i).replace(".","p")
    title = "Reco E"
    binning = [PlotConfig.NEUTRINO4_EE_BINNING]
    value_getter = [lambda event: event.kin_cal.reco_E_lep+event.kin_cal.reco_visE]
    tags = reco_tags
    cuts = [(lambda i=i: lambda event: event.kin_cal.reco_visE < i)()]
    entry = {"name" : name,
             "title": title,
             "binning": binning,
             "value_getter":value_getter,
             "tags": tags,
             "cuts": cuts}

    PLOT_SETTINGS[key] = entry

    key = "E Estimator Eavail bin "+str(i) + " Background Subbed MC"
    name = "E_Estimator_Eavail_bin_"+str(i).replace(".","p") + "_mc_bkgSubbed"
    title = "Reco E"
    binning = [PlotConfig.NEUTRINO4_EE_BINNING]
    value_getter = [lambda event: event.kin_cal.reco_E_lep+event.kin_cal.reco_visE]
    tags = reco_tags
    cuts = [(lambda i=i: lambda event: event.kin_cal.reco_visE < i)()]
    entry = {"name" : name,
             "title": title,
             "binning": binning,
             "value_getter":value_getter,
             "tags": tags,
             "cuts": cuts}

    PLOT_SETTINGS[key] = entry

    key = "E Estimator Eavail bin "+str(i) + " Background Subbed Data"
    name = "E_Estimator_Eavail_bin_"+str(i).replace(".","p") + "_data_bkgSubbed"
    title = "Reco E"
    binning = [PlotConfig.NEUTRINO4_EE_BINNING]
    value_getter = [lambda event: event.kin_cal.reco_E_lep+event.kin_cal.reco_visE]
    tags = reco_tags
    cuts = [(lambda i=i: lambda event: event.kin_cal.reco_visE < i)()]
    entry = {"name" : name,
             "title": title,
             "binning": binning,
             "value_getter":value_getter,
             "tags": tags,
             "cuts": cuts}

    PLOT_SETTINGS[key] = entry
#like histfolio, connect related histograms together
class HistHolder:
    def __init__(self, name, f, sideband, is_mc, pot=1.0,data_pot=None): 
        self.plot_name = TranslateSettings(name)["name"]
        self.dimension = len(TranslateSettings(name)["value_getter"])
        self.is_mc = is_mc
        self.hists= {}
        self.sideband=sideband
        self.pot_scale = data_pot/pot if (data_pot is not None and pot is not None) else pot
        self.valid = False
        self.bin_width_scaled = False
        self.POT_scaled =False 
        
        if f is not None: 
            self.valid = self.ReadHistograms(f)

    def ReadHistograms(self,f):
        variant_arg = [self.plot_name]
        if self.sideband != "Signal":
            variant_arg.append(self.sideband)
        namestring = VariantPlotsNamingScheme(*variant_arg)
        self.hists["Total"] = Utilities.GetHistogram(f, namestring)
        if self.is_mc:
            for cate in list(SignalDef.TRUTH_CATEGORIES.keys())+["Other"]:
                cate_namestring = VariantPlotsNamingScheme(*(variant_arg+[cate]))
                self.hists[cate]=Utilities.GetHistogram(f, cate_namestring) 
        return self.hists["Total"] is not None


    def Scale(self, scale, bin_width_normalize):
        if not self.valid:
            return None
        if self.bin_width_scaled:
            raise ValueError("hist holder alreay scaled by bin_width")

        for k, v in self.hists.items():
            if v is not None:
                v.Scale(scale,"width" if bin_width_normalize else "")
        self.bin_width_scaled = bin_width_normalize

    def AreaScale(self,bin_width_normalize = False):
        if not self.valid:
            return None
        scale = 1.0/self.hists["Total"].Integral()
        self.Scale(scale,bin_width_normalize)


    def POTScale(self,bin_width_normalize = False):
        if not self.valid:
            return None
        if self.POT_scaled and self.pot_scale != 1.0:
            raise ValueError("hist holder alreay scaled by POT")

        for k, v in self.hists.items():
            if v is not None:                
                v.Scale(self.pot_scale if self.is_mc else 1.0,"width" if bin_width_normalize else "")

        self.POT_scaled = True

    def GetCateList(self, grouping=None, with_raw=False):
        if not self.valid:
            return None

        def _raw_integral(h):
            if not h:
                return 0.0
            # MnvH1D/MnvH2D -> use CV hist for a stable integral
            cv = h.GetCVHistoWithStatError() if hasattr(h, "GetCVHistoWithStatError") else h

            # 2D vs 1D
            if cv.InheritsFrom("TH2"):
                return float(cv.Integral(0, cv.GetNbinsX() + 1, 0, cv.GetNbinsY() + 1))
            return float(cv.Integral(0, cv.GetNbinsX() + 1))

        _mc_ints, _colors, _titles, _raws = [], [], [], []
        local_grouping = grouping or {}

        for cate in list(local_grouping.keys())[::-1]:
            config = local_grouping[cate]
            hist = None
            raw = 0.0

            if "cate" in config:
                for fine_cate in config["cate"]:
                    h = self.hists.get(fine_cate, None)
                    if not h:
                        continue

                    # build summed hist
                    if hist is None:
                        hist = h.Clone()
                    else:
                        hist.Add(h)

                    # raw: prefer cached raw_counts if available, otherwise integrate
                    if getattr(self, "raw_counts", None):
                        raw += float(self.raw_counts.get(fine_cate, 0.0))
                    else:
                        raw += _raw_integral(h)

            else:
                h = self.hists.get(cate, None)
                if not h:
                    continue

                # IMPORTANT: clone here to avoid mutating originals downstream
                hist = h.Clone()

                if getattr(self, "raw_counts", None):
                    raw = float(self.raw_counts.get(cate, 0.0))
                else:
                    raw = _raw_integral(h)

            if hist:
                hist.SetTitle(config["title"])
                _mc_ints.append(hist)
                _titles.append(config["title"])
                _colors.append(config["color"])
                _raws.append(raw)

        if with_raw:
            return _mc_ints, _colors, _titles, _raws
        return _mc_ints, _colors, _titles

    # def GetCateList(self,grouping = None):
    #     if not self.valid:
    #         return None
    #     _mc_ints = []
    #     _colors= []
    #     _titles= []
    #     local_grouping = grouping
    #     for cate in list(local_grouping.keys())[::-1]: 
    #         config = local_grouping[cate]
    #         hist = None
    #         if "cate" in config:
    #             for fine_cate in config["cate"]:
    #                 if fine_cate not in self.hists.keys(): 
    #                     continue
    #                 if hist is None and self.hists[fine_cate]:
    #                     hist = self.hists[fine_cate].Clone()
    #                 elif self.hists[fine_cate]:
    #                     hist.Add(self.hists[fine_cate])
    #                 else:
    #                     continue
    #         else:
    #             if not self.hists[cate]:
    #                 continue
    #             else:
    #                 hist = self.hists[cate] if self.hists[cate] else None
                
    #         if hist:
    #             hist.SetTitle(config["title"])
    #             _mc_ints.append(hist)
    #             _titles.append(config["title"])
    #             _colors.append(config["color"])
    #         else:
    #             continue
    #         del hist
    #     return _mc_ints,_colors,_titles

    def GetHist(self): 
        if not self.valid:  
            return None 
        return self.hists["Total"]

    def GetTrueSignalHist(self): 
        if not self.valid:  
            return None 
        return self.hists["CCNuE"]

    def ResumTotal(self):
        self.hists["Total"].Reset("ICESM")
        for i in self.hists:
            if i=="Total":
                continue
            else:
                self.hists["Total"].Add(self.hists[i])

    def GetSumCate(self,cate_list):
        if len(cate_list) == 0:
            return None
        hnew = self.hists[cate_list[0]].Clone()
        for i in range(1,len(cate_list)):
            htmp = self.hists[cate_list[i]]
            if htmp:
                hnew.Add(htmp)
            else:
                continue
        return hnew

    def Add(self,hist_holder):
        if not isinstance(hist_holder,HistHolder):
            raise ValueError("Can only add histholder to histholder")
        for i in self.hists:
            if self.hists[i] and hist_holder.hists[i]:
                self.hists[i].Add(hist_holder.hists[i])
            else:
                print(i,self.hists[i],hist_holder.hists[i])

