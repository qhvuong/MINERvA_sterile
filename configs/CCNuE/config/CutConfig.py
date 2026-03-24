"""
  Analysis cut values and the functions that represent them.
  Tune them using studies/CutTuning.py
  
   author: J. Wolcott <jwolcott@fnal.gov>
   date: December 2013

"""

############################################################################

# hack the fiducial vertex cut,
HACK_R2 = False

# tuned cut values
# DS_CAL_VISE_CUT = 0.2
DS_CAL_VISE_CUT = 0.1
OD_CAL_VISE_CUT = 0.005
PSI_FLAT_CUT = 0.1
FRONT_DEDX_CUT = 2.64  # in MeV/cm
PID_SCORE_CUT = 0.7
MIN_VERTEX_TRACK_MULTIPLICITY = 1
MAX_VERTEX_TRACK_MULTIPLICITY = 5

NONMIP_CLUS_FRAC_CUT = 0.4
TRANSVERSE_GAP_SCORE_CUT = 15
# FIRST_FIRE_FRACTION_CUT = 0.25
FIRST_FIRE_FRACTION_CUT = 0.7
UPSTREAM_OD_ENERGY_CUT = 5000
# EXUV_CUT = 0.2
# EUV_CUT = 0.3

# DEDX_PLANES = [5,9]
# HELICITY= 1

#OLD
FIDUCIAL_APOTHEM = 850
FIDUCIAL_Z_RANGE = [5980,8422]
# FIDUCIAL_APOTHEM = 881.25
# FIDUCIAL_Z_RANGE = [5980,8422]


# Kinematics cutoffs
ELECTRON_ENERGY_RANGE = [1.5, float("inf")] # in GeV
NEUTRINO_ENERGY_RANGE = [0, 100] # in GeV.
LEPTON_ANGLE_RANGE = [0, 20] # in deg
RECO_Q3_RANGE = [0,4]
RECO_PT_RANGE= [.2,1.6]
TRUE_PT_RANGE= [.2,1.6]
TRUE_Q3_RANGE = [0,4]

PSIEE_FLAT_CUT = 0.5
WEXP_CUT = 2
visE_RANGE = [0.0,0.7]
Ethetasquared_CUT = .0032
FRONT_DEDX_PI0_UPPERBOUND = 5

EAVAIL_LOW = [0.0,0.05]
EAVAIL_HIGH = [0.0,0.2]
ELECTRON_ENERGY_CUTOFF = 2.5

############################################################################
# choose the cuts you want from cut library

SAMPLE_CUTS = {
    "Signal" : [ 
        "NoCut",
        "HasNoBackExitingTracks",
        "HasNoVertexMismatch", 
        "StartPointVertexMultiplicity",
        "VertexTrackMultiplicity",
        "HasTracks",
        "Vertex_Z",
        "Vertex_Apothem",
        "DeadTime",
        "Afterpulsing",
        "DSCalVisE",
        "ODCalVisE",
        "EMLikeTrackScore",
        "MeanFrontdEdX",
        "Etheta",
        "Eavail",
        # "Pt",
        # "LeptonAngle",

        ## REMOVED CUTS
        # "NonMIPClusFrac",
        # "TransverseGapScore",
    ],
    "dEdX" : [
        "NoCut",
        "HasNoBackExitingTracks",
        "HasTracks",
        "Vertex_Z",
        "Vertex_Apothem",
        "EMLikeTrackScore",
        "DSCalVisE",
        "ODCalVisE",
        "DeadTime",
        "Afterpulsing",
        "NonMIPClusFrac",
        "TransverseGapScore",
        "HasNoVertexMismatch",
        "StartPointVertexMultiplicity",
        "VertexTrackMultiplicity",
        "Etheta",
        "InverseMeanFrontdEdX",
        "LeptonAngle",
        "Eavail",
        "Pt",
    ]
}

KINEMATICS_CUTS = [
    # "LeptonAngle",
    # "Eavail",
    # "Pt",
]
#######################################
