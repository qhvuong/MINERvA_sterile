
"""
  SystematicsConfig.py:
   Centralization of common systematics configurations.
  
   Original author: J. Wolcott (jwolcott@fnal.gov)
                    May 2014
"""


import itertools
import math

############################ Config Error Bands  ########################
#number of random shifted universes
NUM_UNIVERSE = 100
#number of flux uniberses
USE_NUE_CONSTRAINT = False
AnaNuPDG = 14
USE_SWAPPED=False
NUM_FLUX_UNIVERSE = 100
# detector mass uncertainty
MASS_UNCERTAINTY = 0.014  # = 1.4% (it's fractional).  Laura (Doc7615) says she got this from Ron (Doc6016).

# data EM scale shift in ECAL:
EM_ENERGY_SCALE_SHIFT_ECAL = -0.058 # downward 5.8%

# EM scale uncertainty in ECAL,HCAL, quoted from nu+e paper
EM_ENERGY_SCALE_UNCERTAINTY = {
    "ECAL"   : 0.015,
    "HCAL"   : 0.05,
    "Tracker": 0.02
}

BEAM_ANGLE = math.radians(-3.3)
BEAM_XANGLE_UNCERTAINTY = 1*1e-3 #radians
BEAM_YANGLE_UNCERTAINTY = 0.9*1e-3

LEAKAGE_CORRECTION = lambda E: 0.008*E
AVAILABLE_E_CORRECTION = 1.17

LEAKAGE_SYSTEMATICS = 2 # MeV
LEAKAGE_BIAS = 5#MeV

# electron angle uncertainty
LEPTON_ANGLE_UNCERTAINTY = 1e-3 # this is muon angular resolution. I am worry about this.

# Genie knobs
GENIE_UNIVERSES = [
    "AGKYxF1pi",
    "AhtBY",
    "BhtBY",
    "CCQEPauliSupViaKF",
    "CV1uBY",
    "CV2uBY",
    "EtaNCEL",
    "FrAbs_N",
    "FrAbs_pi",
    "FrCEx_N",
    "FrCEx_pi",
    "FrElas_N",
    "FrElas_pi",
    "FrInel_N",
    "FrPiProd_N",
    "FrPiProd_pi",
    "MFP_N",
    "MFP_pi",
    "MaNCEL",
    "NormDISCC",
    "NormNCRES",
    "RDecBR1gamma",
    "Rvn2pi",
    "Rvp2pi",
    "Theta_Delta2Npi",
    "VecFFCCQEshape"
]

# minerva tune errorbands
UNIVERSES_2P2H = [1,2,3] #2p2h universe variations
RPA_UNIVERSES = {
    "HighQ2":[1,2],
    "LowQ2" :[3,4]
}

NonResPi=True
LowQ2PiWeightChannel = "MENU1PI"
LowQ2PiWeightSysChannel = [None]
NumZExpansionUniverses = 100 #Means Don't use Zexpansion. 100 is default Z expansion

RESPONSE_BRANCHES = [
    "p",
    "meson",
    "em",
    "other",
    "xtalk",
]

NEUTRON_RESPONSE = False

if NEUTRON_RESPONSE:
    RESPONSE_BRANCHES.extend([
        "low_neutron",
        "mid_neutron",
        "high_neutron"
    ])

GEANT_PARTICLES = [
    2212,2112,211
]

################################# Error summary plot config ############################

DETECTOR_RESPONSE_ERROR_GROUPS = {
    "Angular resolution": ["eltheta",],
    "Beam Angle": ["beam_angle",],
    "EM energy scale": ["elE_ECAL","elE_HCAL"],
    "Birk's Constant" : ["birks"],
    "Particle Response":["response_"+i for i in RESPONSE_BRANCHES],
    "Leakage Estimation" : ["Leakage_Uncertainty"],
    "Target Mass" : ["Target_Mass_CH"]
}

MINERVA_TUNNING_ERROR_GROUPS = {
    "RPA" : ["RPA_"+i for i in RPA_UNIVERSES],
    "Low Recoil 2p2h Tune" : ["Low_Recoil_2p2h_Tune"],
    "Low Q2 Pion": ["LowQ2Pi"],
    "FSI bugfix" : ["fsi_weight"],
    "SuSA 2p2h" : ["SuSA_Valencia_Weight"],
    "MK model" : ["MK_model"],
}

MINERVA_TUNNING_ERROR_GROUPS2 = {
    "MK model" : ["MK_model"],
    "FSI bugfix" : ["fsi_weight"],
    "SuSA 2p2h" : ["SuSA_Valencia_Weight"],
}


GENIE_ERROR_GROUPS = {
    "GENIE" : ["GENIE_"+ i for i in (GENIE_UNIVERSES+["EP_MvRES","MaRES","NormCCRES","D2_MaRES","D2_NormCCRES","MaZExpCCQE","MaCCQE", "Rvn1pi", "Rvp1pi"] ) if not (i.startswith("Fr") or i.startswith("MFP")) ]
}

FSI_ERROR_GROUPS = {
    "GENIE-FSI" : ["GENIE_"+ i for i in GENIE_UNIVERSES  if (i.startswith("Fr") or i.startswith("MFP")) ]
}

GEANT_ERROR_GROUPS = {
    "GEANT" : ["GEANT_" +i for i in ("Neutron","Pion","Proton")]
}

BKG_TUNNING_ERROR_GROUPS = {
    "BKG_TUNNING" : ["bkg_tune"]
}


CONSOLIDATED_ERROR_GROUPS_CONFIG = {
    "Detector model": [DETECTOR_RESPONSE_ERROR_GROUPS,GEANT_ERROR_GROUPS],
    "Interaction model": [GENIE_ERROR_GROUPS,FSI_ERROR_GROUPS],
    "MnvTunes" :[MINERVA_TUNNING_ERROR_GROUPS],
    "Muon Reconstruction" :[{"Muon Energy":["MuonAngleXResolution","MINOS_Reconstruction_Efficiency","Muon_Energy_Resolution","MuonAngleYResolution","Muon_Energy_MINERvA","Muon_Energy_MINOS"]}],
    "Alternative Tunning methods" : [BKG_TUNNING_ERROR_GROUPS]
}



CONSOLIDATED_ERROR_GROUPS = {
    key:[e for group in CONSOLIDATED_ERROR_GROUPS_CONFIG[key] for e in itertools.chain.from_iterable(iter(group.values()))] for key in CONSOLIDATED_ERROR_GROUPS_CONFIG
}


DETAILED_ERROR_GROUPS = DETECTOR_RESPONSE_ERROR_GROUPS.copy()
DETAILED_ERROR_GROUPS.update(GENIE_ERROR_GROUPS)
DETAILED_ERROR_GROUPS.update(FSI_ERROR_GROUPS)
DETAILED_ERROR_GROUPS.update(GEANT_ERROR_GROUPS)
