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
DS_CAL_VISE_CUT = 0.2
OD_CAL_VISE_CUT = 0.05
#PSI_CUT = { 2: 0.1, 3: 0.06, None: 0.03 } # energy threshold in GeV, cut value
PSI_FLAT_CUT = 0.1
#FRONT_DEDX_CUT = 2.2  # in MeV/cm
FRONT_DEDX_CUT = 2.4  # in MeV/cm
#FRONT_DEDX_EXCESS_REGION = (2.2, 3.4)  # in MeV/cm
# determined by studying electron/photon particle cannon
# where the search window was varied
#FRONT_DEDX_CUT_BY_POS = [
#    (0, 2.2),
#    (400, 2.2),
#   (50, 2.2),
#   (100, 2.6),
#   (150, 2.8),
#   (200, 3.0),
#   (250, 5.2),
#   (300, 6.4),
#   (350, 7.8),
#   (400, 8.6),
#]
#FRONTDEDX_CUT_GRAPH = ROOT.TGraph( len(FRONT_DEDX_CUT_BY_POS), array.array('d', [p[0] for p in FRONT_DEDX_CUT_BY_POS]), array.array('d', [p[1] for p in FRONT_DEDX_CUT_BY_POS]) )
PID_SCORE_CUT = 0.7
#PID_SCORE_CUT = 0.0
MIN_VERTEX_TRACK_MULTIPLICITY = 1
MAX_VERTEX_TRACK_MULTIPLICITY = 6

NONMIP_CLUS_FRAC_CUT = 0.4
TRANSVERSE_GAP_SCORE_CUT = 15
FIRST_FIRE_FRACTION_CUT = 0.25
UPSTREAM_OD_ENERGY_CUT = 5000
EXUV_CUT = 0.2
EUV_CUT = 0.3

DEDX_PLANES = [5,9]
HELICITY= -1

FIDUCIAL_APOTHEM = 850
FIDUCIAL_Z_RANGE = [5980,8422]

# Kinematics cutoffs
ELECTRON_ENERGY_RANGE = [2.5, float('inf')] # in GeV
NEUTRINO_ENERGY_RANGE = [0, 100] # in GeV.
ELECTRON_ANGLE_RANGE = [0, 30] # in deg
RECO_Q3_RANGE = [0,4]
TRUE_Q3_RANGE = [0,4]
RECO_PT_RANGE= [0,1.0]
RECO_VISE_RANGE= [0,1.2]

WEXP_CUT = 2
Reco_visEcut = 2
FRONT_DEDX_PI0_UPPERBOUND = 5
PsiEe_CUT = 0.5


############################################################################
# choose the cuts you want from cut library

SAMPLE_CUTS = {
    "Signal" : [
        #precuts
        "NoCut",
        "HasTracks",
        "HasNoBackExitingTracks",
        "EMLikeTrackScore",
        #quality cuts
        "DSCalVisE",
        "ODCalVisE",
        "VertexTrackMultiplicity",
        "Afterpulsing",
        "DeadTime",
        #"HasNoNonEMExitingTracks",
        #EM shower quality
        "StartPointVertexMultiplicity",
        "MeanFrontdEdX",
        "NonMIPClusFrac",
        "TransverseGapScore",
        #phase space
        "Vertex_Z",
        "Vertex_Apothem",
        "Eavail",
        # "HasNoVertexMismatch",
        #"RecoQ3",
        #"Exuv",
        #"Euv",

        #"LLR",
        #"Psi",
    ],
    # "High_Psi" : [
    #     "NoCut",
    #     "HasFiducialVertex",
    #     "HasTracks",
    #     "HasNoBackExitingTracks",
    #     "EMLikeTrackScore",
    #     "DSCalVisE",
    #     "ODCalVisE",
    #     "StartPointVertexMultiplicity",
    #     "HasNoVertexMismatch",
    #     "VertexTrackMultiplicity",
    #     "HasNoNonEMExitingTracks",
    #     "InversePsi",  # different from signal
    #     "MeanFrontdEdX",
    #     "NonMIPClusFrac",
    #     "TransverseGapScore",
    #     "DeadTime",
    #     "Afterpulsing",
    # ],
    "Excess_Low_Inline" : [
        #"HasNoNonEMExitingTracks",
        "NoCut",
        "HasTracks",
        "HasNoBackExitingTracks",
        "Vertex_Z",
        "Vertex_Apothem",
        "EMLikeTrackScore",
        "DSCalVisE",
        "ODCalVisE",
        #"HasNoVertexMismatch",
        "VertexTrackMultiplicity",
        "StartPointVertexMultiplicity",
        #"Psi",         #differnt from signal
        "NonMIPClusFrac",
        "TransverseGapScore",
        "DeadTime",
        "Afterpulsing",
        "Eavail",
        "InverseMeanFrontdEdX",
        "Low UIE",
        "LowPsiEe"
        #"Exuv",
        #"Euv"
    ],
     "Excess_High_Inline" : [
        #"HasNoNonEMExitingTracks",
        "NoCut",
        "HasTracks",
        "HasNoBackExitingTracks",
        "Vertex_Z",
        "Vertex_Apothem",
        "EMLikeTrackScore",
        "DSCalVisE",
        "ODCalVisE",
        #"HasNoVertexMismatch",
        "VertexTrackMultiplicity",
         "StartPointVertexMultiplicity",
        #"Psi",         #differnt from signal
        "NonMIPClusFrac",
        "TransverseGapScore",
        "DeadTime",
        "Afterpulsing",
        "Eavail",
        "InverseMeanFrontdEdX",
        "High UIE",
        "LowPsiEe"
        #"Exuv",
        #"Euv"
    ],
    "High_PsiEe_and_dEdX" : [
        #"HasNoNonEMExitingTracks",
        "NoCut",
        "Vertex_Z",
        "Vertex_Apothem",
        "HasTracks",
        "HasNoBackExitingTracks",
        "EMLikeTrackScore",
        "DSCalVisE",
        "ODCalVisE",
        #"HasNoVertexMismatch",
        "VertexTrackMultiplicity",
        "StartPointVertexMultiplicity",
        "NonMIPClusFrac",
        "TransverseGapScore",
        "DeadTime",
        "Afterpulsing",
        "Eavail",
        "InverseMeanFrontdEdX",
        "InverseLowPsiEe"
        #"Exuv",
        #"Euv"
    ],
    "Muon" : [
        "HasMINOSMatch",
        "Vertex_Z",
        "Vertex_Apothem",
        "Neutrino Helicity",
     ]
}

KINEMATICS_CUTS = [
    "LeptonEnergy",
    #"ElectronAngle",
    #"NeutrinoEnergy",
    #"Q3"
]
#######################################



# def _dedx_cut_getter(tree):
#     if tree.n_prongs != 1:
#         return None

#     first_bin = min(tree.prong_binned_energy_bin_indices[0])
#     bin_width = tree.prong_projection_bin_width[0]
#     dE = 0.
#     for i, bin_idx in enumerate(tree.prong_binned_energy_bin_indices[0]):
#         if DEDX_PLANES[0] <= bin_idx - first_bin <= DEDX_PLANES[1]:
#             dE += tree.prong_binned_energy_bin_contents[0][i]
#             return dE / ( (DEDX_PLANES[1] - DEDX_PLANES[0] + 1) * bin_width) * 10   # x10 because we want MeV/cm, not MeV/mm



# def _asym_cut_getter(event, np):
#     if event.n_prongs != 1:
#         return None

#     asymx_num[np] = event.prong_TransverseShowerAsymmetryNumeratorX[np]

#     asymu_num[np] = event.prong_TransverseShowerAsymmetryNumeratorU[np]
#     asymv_num[np] = event.prong_TransverseShowerAsymmetryNumeratorV[np]
#     asymx_den[np] = event.prong_TransverseShowerAsymmetryDenominatorX[np]
#     asymu_den[np] = event.prong_TransverseShowerAsymmetryDenominatorU[np]
#     asymv_den[np] = event.prong_TransverseShowerAsymmetryDenominatorV[np]
    
#     pasymx = asymx_num[np]/asymx_den[np] if asymx_den[np]!=0.0 else 0.0
#     pasymu = asymu_num[np]/asymu_den[np] if asymu_den[np]!=0.0 else 0.0
#     pasymv = asymv_num[np]/asymv_den[np] if asymv_den[np]!=0.0 else 0.0
    
#     pasymx_coord = ((2.0*asymx_num[np] + asymu_num[np] + asymv_num[np])/(asymx_den[np] + asymu_den[np] + asymv_den[np]))/3.0 if (asymx_den[np] + asymu_den[np] + asymv_den[np]) != 0.0 else 0.0
#     pasymy_coord = ((asymu_num[np] - asymv_num[np])/(asymu_den[np] + asymv_den[np]))/math.sqrt(3.0) if (asymu_den[np] + asymv_den[np])!=0.0 else 0.0
#     pasym = math.sqrt(math.pow(pasymx_coord,2) + math.pow(pasymy_coord,2))

#     return pasym

