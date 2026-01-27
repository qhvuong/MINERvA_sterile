"""
    Analysis Cut Library that provide functions/class to be used in EventClassification.
    Cut value is configured in config/CutConfig.py

    Adapted from Jeremy's CutConfig.

    author: H. Su
    date: Jul 2019
"""


import array
import functools
import math
import ROOT

from config import CutConfig


class SelectionCut(object):
    def __init__(self, **kwargs):
        self.name = kwargs["name"]
        try:
            self._value_getter = kwargs["value_getter"]
        except KeyError:
            self._value_getter = functools.partial(lambda name, event, nprong: getattr(event, name), self.name)

        self._cut_fn = kwargs["cut_fn"]

        try:
            self._variable_range = kwargs["variable_range"]
        except KeyError:
            self._variable_range = None

    def Values(self, event, nprong=0):
        return self._value_getter(event, nprong)

    def DoesEventPass(self, event, nprong=0):
        return self.DoValuesPass( self.Values(event, nprong) )

    def DoValuesPass(self, vals):
        return self._cut_fn(vals)

# these are used often.  just consolidate them here.
REQUIRE_POSITIVE_INT = lambda val: val >= 1
REQUIRE_UNITY_INT = lambda val: int(val) == 1  #val >= 1 and val < 2

VARIABLE_RANGE_01 = [-0.1,0.9,1.9] # for true/false varialbe, avoiding ambiguity


"""
  How to define a new cut:
1) add an entry to CUT_CONFIGS, format is :
  "name of cut" : {
   "value_getter" : lambada function for getting the cut variable,
   "cut_fn" : lambda function for determing if the cut varialbe passes.
  }
2) add the cut to the list of sample cuts in CutConfig.py 
"""

def CalcApothem(x,y):
    x=abs(x)
    y=abs(y)
    if ( x == 0 or y/x > 1/math.sqrt(3)):
        return (y+x/math.sqrt(3))/2*math.sqrt(3)
    else:
        return x

def CalTheta2Hybrid(x,y):
    return x*x+y*y

def FHC_Cut(available_energy,lepton_energy):
    if lepton_energy < 1.00:
        return(0 < available_energy < 2.00)
    elif 1.00 <= lepton_energy <= 1.75:
        return(0 < available_energy < 0.9)
    elif 1.75 < lepton_energy <= 3.00:
        return(0 < available_energy < 1.3)
    elif 3.00 < lepton_energy:
        return(0 < available_energy < 2.00)
    else:
        return(False)

def RHC_Cut(available_energy,lepton_energy):
    if lepton_energy < 0.50:
        return(0 < available_energy < 0.1)
    elif 0.50 <= lepton_energy <= 0.75:
        return(0 < available_energy < 1.3)
    elif 0.75 < lepton_energy <= 1.00:
        return(0 < available_energy < 0.35)
    elif 1.00 < lepton_energy <= 1.25:
        return(0 < available_energy < 0.90)
    elif 1.25 < lepton_energy <= 1.50:
        return(0 < available_energy < 2.00)
    elif 1.50 < lepton_energy <= 2.25:
        return(0 < available_energy < 0.90)
    elif 2.25 < lepton_energy <= 3.50:
        return(0 < available_energy < 1.30)
    elif 3.50 < lepton_energy:
        return(0 < available_energy < 2.00)
    else:
        return(False)

def passSingleProtonCut(event, scoreShift = 0): 

    # Keep around 81% efficiency in each Q2 region
    pass_proton = False

    #q2_reco = event.MasterAnaDev_Q2 + q2Shift 
    #q2_reco /= 10**6 
    q2_reco = event.kin_cal.reco_q2_cal
    protonScore1_reco = event.MasterAnaDev_proton_score1 + scoreShift 
    if q2_reco < 0.2: 
        if protonScore1_reco > 0.2:
            pass_proton = True
        elif q2_reco >= 0.2 and q2_reco < 0.6:
            if protonScore1_reco > 0.1:
                pass_proton = True
        elif q2_reco >= 0.6: 
            if protonScore1_reco > 0.0:
                pass_proton = True 

    return(pass_proton)

def passHybridProtonNodeCut(event, cutval):
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
    return chi2<cutval

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
        return False #/nonodes
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

CUT_CONFIGS = {
    "NoCut": {
        "value_getter": lambda event, nprong: None,
        "cut_fn": lambda val: True,
        "variable_range": VARIABLE_RANGE_01,
    },
    "True NuMu": {
        "value_getter": lambda event, nprong: abs(event.mc_incoming),
        "cut_fn": lambda val: val==14,
    },
    "True NuE": {
        "value_getter": lambda event, nprong: abs(event.mc_incoming),
        "cut_fn": lambda val: val==12,
    },
    "True Scattering": {
        "value_getter": lambda event, nprong: abs(event.mc_intType),
        "cut_fn": lambda val: val==7,
    },
    "HasFiducialVertex": {
        "cut_fn": REQUIRE_POSITIVE_INT,
        "variable_range": VARIABLE_RANGE_01,
    },
    "Low UIE": {
        "value_getter": lambda event,nprong: event.UpstreamInlineEnergy,
        "cut_fn":lambda val : val<10,
    },
    "High UIE": {
        "value_getter": lambda event,nprong: event.UpstreamInlineEnergy,
        "cut_fn":lambda val : val>10,
    },
    "PsiEe": {
        "value_getter": lambda event, nprong: (event.prong_TotalVisE[nprong]/1e3, event.Psi*event.kin_cal.reco_E_lep) ,
        "cut_fn": lambda vals: vals[1] < CutConfig.PSIEE_FLAT_CUT,
    },
     "LowPsiEe": {
        "value_getter": lambda event, nprong: event.kin_cal.reco_E_lep*event.Psi,
        "cut_fn": lambda vals: vals < CutConfig.PsiEe_CUT,
        "variable_range": [0.2* i for i in range(0,21)]
    },
    "InversePsiEe": {
        "value_getter": lambda event, nprong: (event.prong_TotalVisE[nprong]/1e3, event.Psi*event.kin_cal.reco_E_lep),
        "cut_fn": lambda vals: vals[1] > CutConfig.PSIEE_FLAT_CUT,
    },
    "Vertex_Z": {
        "value_getter": lambda event, nprong: event.vtx[2],
        "cut_fn": lambda val: CutConfig.FIDUCIAL_Z_RANGE[0] <= val <= CutConfig.FIDUCIAL_Z_RANGE[1],
        "variable_range": [100 * i for i in range(100)]
    },
    "Vertex_Z_Centralized" : {
        "value_getter": lambda event,nprong: event.classifier.ZRange.passesCut(event,ROOT.PlotUtils.detail.empty()),
        "cut_fn" : lambda val:val
    },
    "Vertex_Apothem_Centralized" : {
        "value_getter": lambda event,nprong: event.classifier.Apothem.passesCut(event,ROOT.PlotUtils.detail.empty()),
        "cut_fn" : lambda val:val
    },
    "Vertex_Apothem": {
        "value_getter":lambda event,nprong: CalcApothem(event.vtx[0],event.vtx[1]),
        "cut_fn": lambda val: val <= CutConfig.FIDUCIAL_APOTHEM,
        "variable_range": [10*i for i in range(120)]
    },

    "Truth_Vertex_Z": {
        "value_getter": lambda event, nprong: event.mc_vtx[2],
        "cut_fn": lambda val: CutConfig.FIDUCIAL_Z_RANGE[0] <= val <= CutConfig.FIDUCIAL_Z_RANGE[1],
        "variable_range": [100 * i for i in range(100)]
    },
    "Truth_Vertex_Apothem": {
        "value_getter":lambda event,nprong: CalcApothem(event.mc_vtx[0],event.mc_vtx[1]),
        "cut_fn": lambda val: val <= CutConfig.FIDUCIAL_APOTHEM,
        "variable_range": [10*i for i in range(120)]
    },

    "HasTracks": {
        "value_getter": lambda event, nprong: event.n_prongs>0,
        "cut_fn": REQUIRE_POSITIVE_INT,
        "variable_range": VARIABLE_RANGE_01,
    },
    
    "HasNoBackExitingTracks": {
        "cut_fn": REQUIRE_POSITIVE_INT,
        "variable_range": VARIABLE_RANGE_01,
    },

    "EMLikeTrackScore": {
        "value_getter": lambda event,n_prongs: event.prong_part_score,
        "cut_fn": lambda val: [_>=CutConfig.PID_SCORE_CUT for _ in val].count(True)>0,
        "variable_range": [0.1*i for i in range(5,11)]
    },

    "UniqueEMLikeTrackScore": {
        "value_getter": lambda event,n_prongs: event.prong_part_score,
        "cut_fn": lambda val: [_>=CutConfig.PID_SCORE_CUT for _ in val].count(True)==1,
        "variable_range": [0.1*i for i in range(5,11)]
    },
    
    "DSCalVisE" : {
        "value_getter": lambda event, nprong:
            # return None if there's HCAL energy but not ECAL energy...
            event.prong_HCALVisE[nprong]/event.prong_ECALVisE[nprong] if event.prong_ECALVisE[nprong] > 0.1 else (1e10 if event.prong_HCALVisE[nprong] > 0 else 0),
        "cut_fn": lambda val: val <= CutConfig.DS_CAL_VISE_CUT,
        "variable_range": [0.1* i for i in range(0,11)]
    },
    
    "ODCalVisE": {
        "value_getter": lambda event, nprong:
            # return None if there's OD energy but not side ECAL energy...
            event.prong_ODVisE[nprong]/event.prong_SideECALVisE[nprong] if event.prong_SideECALVisE[nprong] > 0.1 else (1e10 if event.prong_ODVisE[nprong] > 0 else 0),
        "cut_fn": lambda val: val <= CutConfig.OD_CAL_VISE_CUT,
        "variable_range": [0.1* i for i in range(0,11)]
    },

    "StartPointVertexMultiplicity": {
        "cut_fn": REQUIRE_UNITY_INT,
        "variable_range": list(range(0,6))
    },


    "HasNoVertexMismatch": {
        "cut_fn": REQUIRE_UNITY_INT,
        "variable_range": VARIABLE_RANGE_01,
    },


    "VertexTrackMultiplicity": {
        "cut_fn": lambda val: CutConfig.MIN_VERTEX_TRACK_MULTIPLICITY <=  val <= CutConfig.MAX_VERTEX_TRACK_MULTIPLICITY,
        "variable_range": list(range(0,6))
    },

    # "NuEVertexTrackMultiplicity": {
    #     "value_getter": lambda event: event.VertexTrackMultiplicity,
    #     "cut_fn": lambda val: val < CutConfig.MIN_VERTEX_TRACK_MULTIPLICITY,
    #     "variable_range": list(range(0,2))
    # },

    "HasNoNonEMExitingTracks": {
        "cut_fn": REQUIRE_POSITIVE_INT,
        "variable_range": VARIABLE_RANGE_01,
    },

    "Psi": {
        "value_getter": lambda event, nprong: (event.prong_TotalVisE[nprong]/1e3, event.Psi) ,
        "cut_fn": lambda vals: vals[1] < CutConfig.PSI_FLAT_CUT,
        "variable_range": [0.1* i for i in range(0,11)]
    },

    "Wexp": {
        "value_getter": lambda event, nprong: event.kin_cal.reco_W,
        "cut_fn": lambda vals: vals>0 and vals <= CutConfig.WEXP_CUT,
        "variable_range": [0.5* i for i in range(0,11)]
    },
    "FHC_proton": {
        "value_getter": lambda event, nprong: event,
        "cut_fn": lambda val: passHybridProtonNodeCut(val,10),
        #"cut_fn": lambda val: passSingleProtonCut(val,0),
        "variable_range": [0.1* i for i in range(0,11)]
    },
    "Pt": {
        "value_getter": lambda event, nprong: event.kin_cal.reco_Pt_lep,
        "cut_fn": lambda val: CutConfig.RECO_PT_RANGE[0] <= val < CutConfig.RECO_PT_RANGE[1],
        "variable_range": [0.1*i for i in range(0,21)]
    },
    "Etheta": {
        # "value_getter": lambda event, nprong: event.kin_cal.reco_E_lep * (event.kin_cal.reco_theta_lep_rad)**2,
        "value_getter": lambda event, nprong: event.kin_cal.reco_E_lep * CalTheta2Hybrid(event.kin_cal.reco_thetaX_lep_rad, event.kin_cal.true_thetaY_lep_rad),
        "cut_fn": lambda vals: vals >= CutConfig.Ethetasquared_CUT,
    },
    "EthetaSideband": {
        "value_getter": lambda event, nprong: event.kin_cal.reco_E_lep * (event.kin_cal.reco_theta_lep_rad)**2,
        "cut_fn": lambda vals: CutConfig.Ethetasquared_SIDEBAND_LOW <= vals <= CutConfig.Ethetasquared_SIDEBAND_HIGH,
    },
    "Eavail": {
        "value_getter": lambda event, nprong: event.kin_cal.reco_visE,
        "cut_fn": lambda val: CutConfig.visE_RANGE[0] <= val < CutConfig.visE_RANGE[1]
    },
    "EavailSideband": {
        "value_getter": lambda event, nprong: event.kin_cal.reco_visE,
        "cut_fn": lambda val: val >= CutConfig.visE_RANGE[1]
    },
    "LeptonEnergy": {
        "value_getter": lambda event,nprong: event.kin_cal.reco_E_lep,
        "cut_fn": lambda val: val > CutConfig.ELECTRON_ENERGY_CUT
    },
    "Q2": {
        "value_getter": lambda event,nprong: event.kin_cal.reco_q2_cal,
        "cut_fn": lambda val: val < CutConfig.Q2_CUT
    },
    "OpeningAngle": {
        "value_getter": lambda event, nprong: event.ElectronProtonAngle(),
        "cut_fn": lambda val: CutConfig.RECO_ANGLE[0] <= val <= CutConfig.RECO_ANGLE[1],
    },
    "InversePsi": {
         "value_getter": lambda event, nprong: (event.prong_TotalVisE[nprong]/1e3, event.Psi) ,
         "cut_fn": lambda vals: vals[1] > CutConfig.PSI_FLAT_CUT,
         "variable_range": [0.1* i for i in range(0,11)]
     },

    "MeanFrontdEdX": {
        "value_getter": lambda event, nprong: event.prong_dEdXMeanFrontTracker[nprong],
        "cut_fn": lambda val: 0 < val <= CutConfig.FRONT_DEDX_CUT,
        "variable_range": [0.1* i for i in range(0,51)]
    },

    "MidMeanFrontdEdX" : {
        "value_getter": lambda event, nprong: event.prong_dEdXMeanFrontTracker[nprong],
        "cut_fn": lambda val: val > CutConfig.FRONT_DEDX_CUT and val<CutConfig.FRONT_DEDX_PI0_UPPERBOUND,
        "variable_range": [0.1* i for i in range(0,51)]
    },

    "HighMeanFrontdEdX" : {
        "value_getter": lambda event, nprong: event.prong_dEdXMeanFrontTracker[nprong],
        "cut_fn": lambda val: val > CutConfig.FRONT_DEDX_PI0_UPPERBOUND,
        "variable_range": [0.1* i for i in range(0,51)]
    },

    "NonMIPClusFrac": {
        "value_getter": lambda event, nprong: event.prong_NonMIPClusFrac[nprong],
        "cut_fn": lambda val: val > CutConfig.NONMIP_CLUS_FRAC_CUT,
        "variable_range": [0.1* i for i in range(0,11)]
    },

    "TransverseGapScore": {
        "value_getter": lambda event, nprong: event.prong_TransverseGapScore[nprong],
        "cut_fn": lambda val: val > CutConfig.TRANSVERSE_GAP_SCORE_CUT,
        "variable_range": [1.5* i for i in range(0,21)]
    },
    
    "ProtonScore": {
        "value_getter": lambda event, nprong: event.MasterAnaDev_proton_score1,
        "cut_fn": lambda val: 0 <= val <= 1,
    },

    "ZDifference": {
        "value_getter": lambda event, nprong: (event.MasterAnaDev_proton_startPointZ,event.prong_axis_vertex[nprong][2]),
        "cut_fn": lambda vals: CutConfig.RECO_ZDIFF_RANGE[0]<= vals[0]-vals[1] <= CutConfig.RECO_ZDIFF_RANGE[1],
    },

    "ZDifference_vtx": {
        "value_getter": lambda event, nprong: (event.prong_axis_vertex[nprong][2],event.vtx[2]),
        "cut_fn": lambda vals: CutConfig.RECO_VTX_ZDIFF_RANGE[0]<= vals[0]-vals[1] <= CutConfig.RECO_VTX_ZDIFF_RANGE[1],
    },
    
    "ProtonEnd": {
        "value_getter": lambda event, nprong: event.MasterAnaDev_proton_endPointZ,
        "cut_fn": lambda val: val <= CutConfig.RECO_PROTON_END,
    },
    
    #"HasNoMichelElectrons": {
    #"value_getter": lambda event, nprong: (event.michel_digits, event.michel_energy, event.michel_slice_energy),
    #   "cut_fn": lambda vals: all([not(vals[0][i] < 35 and vals[1][i] < 55 and vals[2][i] < 100 and vals[1][i]/vals[0][i] > 0.8) for i in range(len(vals[0]))])
    #},#Better way for this?
    
    "DeadTime": {
        "value_getter": lambda event, nprong: event.phys_n_dead_discr_pair_upstream_prim_track_proj,
        "cut_fn": lambda val: val >= 0 and val <= 1,
        "variable_range": [0.1* i for i in range(0,11)]
    },
    
    "Afterpulsing": {
        "value_getter": lambda event, nprong: event.prong_FirstFireFraction[nprong],
        "cut_fn": lambda val: val >= CutConfig.FIRST_FIRE_FRACTION_CUT,
        "variable_range": [0.1* i for i in range(0,11)]
    },

    "Exuv" : {
        "value_getter": lambda event, nprong: event.kin_cal.Exuv,
        "cut_fn": lambda val:abs(val) <= CutConfig.EXUV_CUT,
        "variable_range": [0.1* i -0.5 for i in range(0,11)]
    },

    "Euv" : {
        "value_getter": lambda event, nprong: event.kin_cal.Euv,
        "cut_fn": lambda val: abs(val) <= CutConfig.EUV_CUT,
        "variable_range": [0.1* i -0.5 for i in range(0,11)]
    },

    "LLR" : {
        "value_getter": lambda event, nprong: max(-20,min(19.9, event.kin_cal.LLR)) if event.kin_cal.LLR is not None else None,
        "cut_fn": lambda val: val is not None and val <= 0,
        "variable_range": [i for i in range(-20,20)]
    },

    "Neutrino Helicity" : {
        "value_getter": lambda event, nprong: -1 if event.MasterAnaDev_nuHelicity ==1 else 1,
        "cut_fn": lambda val: CutConfig.HELICITY is None or val * CutConfig.HELICITY>0,
        "variable_range": [-1.1,0,1.1]
    },
    "HasMINOSMatch" : {
        "value_getter": lambda event,nprong: event.IsMinosMatchMuon(),
        "cut_fn": lambda val: val,
        "variable_range": [-1.1,0,1.1]
    }

}  # CUT_CONFIGS


KINEMATICS_CUT_CONFIGS = {
    "RecoEavail": {
        "value_getter": lambda event, nprong: event.kin_cal.reco_visE,
        "cut_fn": lambda val: CutConfig.visE_RANGE[0] <= val < CutConfig.visE_RANGE[1],
        "variable_range": [0.4* i for i in range(0,11)]
    },
    "RecoLeptonEnergy": {
        "value_getter": lambda event,nprong: event.kin_cal.reco_E_lep,
        "cut_fn": lambda val: CutConfig.ELECTRON_ENERGY_RANGE[0] <= val < CutConfig.ELECTRON_ENERGY_RANGE[1],
        "variable_range": [0.5*i for i in range(0,21)]
    },
    
    "RecoLeptonAngle": {
        "value_getter": lambda event,nprong: event.kin_cal.reco_theta_lep,
        "cut_fn": lambda val: CutConfig.LEPTON_ANGLE_RANGE[0] <= val < CutConfig.LEPTON_ANGLE_RANGE[1],
    },

    "RecoNeutrinoEnergy": {
        "value_getter": lambda event,nprong: event.kin_cal.reco_E_nu_cal,
        "cut_fn": lambda val: CutConfig.NEUTRINO_ENERGY_RANGE[0] <= val < CutConfig.NEUTRINO_ENERGY_RANGE[1],
    },
    "RecoQ3": {
        "value_getter": lambda event,nprong: event.kin_cal.reco_q3,
        "cut_fn": lambda val: CutConfig.RECO_Q3_RANGE[0] <= val < CutConfig.RECO_Q3_RANGE[1],
        "variable_range": [0.1*i for i in range(0,21)]
    },
    "RecoPt": {
        "value_getter": lambda event,nprong: event.kin_cal.reco_Pt_lep,
        "cut_fn": lambda val: CutConfig.RECO_PT_RANGE[0] <= val < CutConfig.RECO_PT_RANGE[1],
        "variable_range": [0.1*i for i in range(0,21)]
    },
    "TrueEavail": {
        "value_getter": lambda event, nprong: event.kin_cal.true_visE,
        "cut_fn": lambda val: CutConfig.visE_RANGE[0] <= val < CutConfig.visE_RANGE[1],
        "variable_range": [0.4* i for i in range(0,11)]
    },
    "TrueLeptonEnergy": {
        "value_getter": lambda event,nprong: event.kin_cal.true_E_lep,
        "cut_fn": lambda val: CutConfig.ELECTRON_ENERGY_RANGE[0] <= val < CutConfig.ELECTRON_ENERGY_RANGE[1],
    },
    
    "TrueLeptonAngle": {
        "value_getter": lambda event,nprong: event.kin_cal.true_theta_lep,
        "cut_fn": lambda val: CutConfig.LEPTON_ANGLE_RANGE[0] <= val < CutConfig.LEPTON_ANGLE_RANGE[1],
    },

    "TrueNeutrinoEnergy": {
        "value_getter": lambda event,nprong: event.kin_cal.true_enu_genie,
        "cut_fn": lambda val: CutConfig.NEUTRINO_ENERGY_RANGE[0] <= val < CutConfig.NEUTRINO_ENERGY_RANGE[1],
    },
    "TrueQ3": {
        "value_getter": lambda event,nprong: event.kin_cal.true_q3,
        "cut_fn": lambda val: CutConfig.TRUE_Q3_RANGE[0] <= val < CutConfig.TRUE_Q3_RANGE[1],
        "variable_range": [0.1*i for i in range(0,21)]
    },
    "TruePt": {
        "value_getter": lambda event,nprong: event.kin_cal.true_Pt_lep,
        "cut_fn": lambda val: CutConfig.RECO_PT_RANGE[0] <= val < CutConfig.RECO_PT_RANGE[1],
        "variable_range": [0.1*i for i in range(0,21)]
    },
}

if CutConfig.HACK_R2:
    print("\n\n=====================================================================")
    print("WARNING: using R < 850mm and z cuts instead of ntuple fiducial volume cut...")
    print("=====================================================================\n\n")
    
    CUT_CONFIGS["HasFiducialVertex"] = {
        "value_getter": lambda event, nprong: (math.sqrt(event.vtx[0]**2 + event.vtx[1]**2), event.vtx[2]),
        "cut_fn": lambda vals: vals[0] <= 850 and vals[1] > 5990 and vals[1] < 8340,
    }
    CUT_CONFIGS["ContainedProton"] = {
        "value_getter": lambda event, nprong: (math.sqrt(event.MasterAnaDev_proton_endPointX**2 + event.MasterAnaDev_proton_endPointY**2), event.MasterAnaDev_proton_endPointZ),
        "cut_fn": lambda vals: vals[0] <= 850 and vals[1] > 5990 and vals[1] < 8500,
    }

# make SelectionCut objects now

class Cuts(object):
    def __init__(self):
        self.cutmap = {}
    def __getitem__(self,key):
        
        if key in self.cutmap:
            return self.cutmap[key]
        elif key in CUT_CONFIGS:
            config = CUT_CONFIGS[key]
            config["name"]=key
            self.cutmap[key]=SelectionCut(**config)
            return self.cutmap[key]
        elif key in KINEMATICS_CUT_CONFIGS:
            config = KINEMATICS_CUT_CONFIGS[key]
            config["name"]=key
            self.cutmap[key]=SelectionCut(**config)
            return self.cutmap[key]
        else:
            raise KeyError("Cut {} not defined".format(key))


CUTS = Cuts()
