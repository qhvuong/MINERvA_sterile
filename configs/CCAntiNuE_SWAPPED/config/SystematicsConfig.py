
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
AnaNuPDG=-12
USE_SWAPPED=False
NUM_FLUX_UNIVERSE = 100
# detector mass uncertainty
MASS_UNCERTAINTY = 0.014  # = 1.4% (it's fractional).  Laura (Doc7615) says she got this from Ron (Doc6016).

# data EM scale shift in ECAL:
EM_ENERGY_SCALE_SHIFT_ECAL = -0.058 # downward 5.8%

# EM scale uncertainty in ECAL,HCAL, quoted from nu+e paper
EM_ENERGY_SCALE_UNCERTAINTY = {
    "ECAL": 0.015,
    "HCAL": 0.05
}

BEAM_ANGLE = math.radians(-3.3)
BEAM_XANGLE_UNCERTAINTY = 1*1e-3 #radians
BEAM_YANGLE_UNCERTAINTY = 0.9*1e-3

LEAKAGE_CORRECTION = lambda E: 0.008*E
#LEAKAGE_CORRECTION = lambda E: 0
AVAILABLE_E_CORRECTION = 1.17

LEAKAGE_SYSTEMATICS = 2 # MeV
LEAKAGE_BIAS = 5#MeV

# electron angle uncertainty
ELECTRON_ANGLE_UNCERTAINTY = 1e-3 # this is muon angular resolution. I am worry about this.

#DEDX_ENERGY_SCALE_UNCERTAINTY = 0.015 # = 1.5% (still fractional).  Based on eyeballing the dE/dx distribution for this analysis in the extra energy sideband.

#COHERENT_UNCERTAINTY = 0.2  # = 20% (fractional).  Covers the measurement errors and residual difference between GENIE COH and MINERvA's measurement after the reweighting done in tools/Corrections.py
#COHERENT_EB_NAME = "CoherentModel"

#MYSTERY_ALT_EB_NAME = "MysteryProcessModel"
#SIDEBAND_MODEL_EB_NAME = "SidebandModel"  # systematic addressing data-MC disagreement in EE sideband at low energy.  used in BackgroundFitting.py

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
    #"FrInel_pi",
    "FrPiProd_N",
    "FrPiProd_pi",
    "MFP_N",
    "MFP_pi",
    #"MaCCQE",
    #"MaCCQEshape",
    #"NormCCQE",  these three are taken care by seperate class
    "MaNCEL",
    #"MaRES",
    #"MvRES",
    #"NormCCRES", these three are taken care of by seperate class for deuterium and EP updates
    "NormDISCC",
    "NormNCRES",
    "RDecBR1gamma",
    #"Rvn1pi",
    #"Rvp1pi",  these two are taken care by seperate class
    "Rvn2pi",
    "Rvp2pi",
    "Theta_Delta2Npi",
    "VecFFCCQEshape"
]


# Two non-standard genie knobs? should revisit in the future
  #"efnucr_flat", "EFNUCR",

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
# not all the variations are relevant for the PC excess samples
#PC_EXCESS_VARIATIONS = VARIATIONS[:]
#PC_EXCESS_VARIATIONS.remove("BirksUp")
#PC_EXCESS_VARIATIONS.remove("BirksDown")
#PC_EXCESS_VARIATIONS.remove("MichelEScaleUp")
#PC_EXCESS_VARIATIONS.remove("MichelEScaleDown")
#PC_EXCESS_VARIATIONS.remove("GENIE_EFNUCR-1Sigma")
#PC_EXCESS_VARIATIONS.remove("GENIE_EFNUCR+1Sigma")


# EXTERNAL_ERROR_BANDS = {
# 	"BirksConstant": ["BirksDown", "BirksUp",],
# 	"CrossTalk": [ "XtalkSmear", ],
# 	"EnergyScale": [ "EMEScaleUp", "EMEScaleDown", ],
# 	"MichelEScale": [ "MichelEScaleUp", "MichelEScaleDown",],
# 	"dEdXEScale": ["dEdXEScaleUp", "dEdXEScaleDown",],
# 	"EFNUCR": ["GENIE_EFNUCR-1Sigma", "GENIE_EFNUCR+1Sigma"],
# 	MYSTERY_ALT_EB_NAME: [MYSTERY_ALT_EB_NAME, ]
# }


# for err_name, err_group in BKND_UNIV_SYSTEMATICS.iteritems():
# 	# these are constraint err bands, not spectators
# 	if "flux" in err_name.lower():
# 		continue
		
	# for err_name in err_group:
	# 	ERROR_CONSTRAINT_STRATEGIES[err_name] = PlotUtils.MnvHistoConstrainer.PRESERVE_FRACTIONAL_ERR

# if USE_GENIE_INDIV_KNOBS:
# 	BKND_UNIV_SYSTEMATICS["GENIE"] += [
# 		"truth_genie_wgt_AGKYxF1pi",
# 		"truth_genie_wgt_AhtBY",
# 		"truth_genie_wgt_BhtBY",
# 		"truth_genie_wgt_CCQEPauliSupViaKF",
# 		"truth_genie_wgt_CV1uBY",
# 		"truth_genie_wgt_CV2uBY",
# 		"truth_genie_wgt_EtaNCEL",
# 		"truth_genie_wgt_FrAbs_N",
# 		"truth_genie_wgt_FrAbs_pi",
# 		"truth_genie_wgt_FrCEx_N",
# 		"truth_genie_wgt_FrCEx_pi",
# 		"truth_genie_wgt_FrElas_N",
# 		"truth_genie_wgt_FrElas_pi",
# 		"truth_genie_wgt_FrInel_N",
# 		"truth_genie_wgt_FrInel_pi",
# 		"truth_genie_wgt_FrPiProd_N",
# 		"truth_genie_wgt_FrPiProd_pi",
# 		"truth_genie_wgt_MFP_N",
# 		"truth_genie_wgt_MFP_pi",
# 		"truth_genie_wgt_MaCCQE",
# 		"truth_genie_wgt_MaNCEL",
# 		"truth_genie_wgt_MaRES",
# 		"truth_genie_wgt_MvRES",
# 		"truth_genie_wgt_NormDISCC",
# 		"truth_genie_wgt_RDecBR1gamma",
# 		"truth_genie_wgt_Rvn1pi",
# 		"truth_genie_wgt_Rvn2pi",
# 		"truth_genie_wgt_Rvp1pi",
# 		"truth_genie_wgt_Rvp2pi",
# 		"truth_genie_wgt_Theta_Delta2Npi",
# 		"truth_genie_wgt_VecFFCCQEshape",
# 	]

# GENIE_ERROR_GROUPS = {
# 	#"Elastic": ["NormCCQE", "MaCCQEshape", "MaCCQE", "VecFFCCQEshape", "CCQEPauliSupViaKF", "MaNCEL", "EtaNCEL"],
# 	#"Resonance": ["NormCCRES", "MaRES", "MvRES", "NormNCRES", "RDecBR1gamma", "Theta_Delta2Npi"],
# 	#"DIS": ["NormDISCC", "Rvp1pi", "Rvn1pi", "Rvp2pi", "Rvn2pi", "AGKYxF1pi", "AhtBY", "BhtBY", "CV1uBY", "CV2uBY",],
#     "Primary interaction" :[
#         # Elastic
#         "NormCCQE", "MaCCQEshape", "MaCCQE", "VecFFCCQEshape", "CCQEPauliSupViaKF", "MaNCEL", "EtaNCEL",
#         # Resonance
#         "NormCCRES", "MaRES", "MvRES", "NormNCRES", "RDecBR1gamma", "Theta_Delta2Npi",
#         # DIS
#         "NormDISCC", "Rvp1pi", "Rvn1pi", "Rvp2pi", "Rvn2pi", "AGKYxF1pi", "AhtBY", "BhtBY", "CV1uBY", "CV2uBY",
#     ],
# 	"FSI model": ["FrAbs_N", "FrAbs_pi", "FrCEx_N", "FrCEx_pi", "FrElas_N", "FrElas_pi", "FrInel_N", "FrInel_pi", "FrPiProd_N", "FrPiProd_pi", "MFP_N", "MFP_pi",],
#     #"efnucr_flat", "EFNUCR"],
# }

# OTHER_INTERACTION_ERROR_GROUPS = {
# 	"Sideband model": ["BackgroundFit", "SidebandModel" ],
# 	"Coherent": ["CoherentModel",], 
# 	"Excess process model": ["MysteryProcessModel",],
# }


# DETECTOR_RESPONSE_ERROR_GROUPS = {
# 	"Angular resolution": ["ElectronAngle", "theta_bias", "theta_smear"],
# 	"EM energy scale": ["EnergyScale", "dEdXEScale"],
# 	"Birks's constant": ["BirksConstant", "birks_constant_flat"],
# 	"Cross-talk": ["CrossTalk", "xtalk_alt"],
# 	"Global energy scale": ["meu",],
# 	"Hadron response": ["pion_response", "proton_response", "other_response", "neutron_pathlength"],
# 	"Michel electron energy scale": ["MichelEScale",],
# 	"MINOS": ["minos_overlap",],
# 	"Muon energy scale": ["muon_energy_scale", ],
# 	"Target mass": ["Mass", "target_mass"],
# 	"Timing": ["WideTimeWindow",],
# 	"Reconstruction": ["no_isoblobs_cut", "rock_subtraction", "tracking_efficiency", "unfolding", "vertex_smear"],
# }

################################# Error summary plot config ############################

DETECTOR_RESPONSE_ERROR_GROUPS = {
    "Angular resolution": ["eltheta",],
    "Beam Angle": ["beam_angle","BeamAngleX","BeamAngleY","BeamAngleY"],
    "EM energy scale": ["elE_ECAL","elE_HCAL",""],
    "Birk's Constant" : ["birks"],
    "Particle Response":["response_"+i for i in RESPONSE_BRANCHES],
    "Tracking Efficiency":["Proton_TrackEff"],
    "Leakage Estimation" : ["Leakage_Uncertainty"],
    "Target Mass" : ["Target_Mass","Target_Mass_CH"]
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
    "Muon Reconstruction" :[{"Muon Energy":["MINOS_Reconstruction_Efficiency","Michel_Efficiency","MuonReconstruction","MuonAngleXResolution","MuonAngleYResolution","Muon_Energy_MINERvA","Muon_Energy_MINOS","Muon_Energy_Resolution"]}],
    "Alternative Tunning methods" : [BKG_TUNNING_ERROR_GROUPS]
   # "Others":[DETECTOR_RESPONSE_ERROR_GROUPS,GEANT_ERROR_GROUPS,GENIE_ERROR_GROUPS,FSI_ERROR_GROUPS,MINERVA_TUNNING_ERROR_GROUPS],
}



CONSOLIDATED_ERROR_GROUPS = {
    key:[e for group in CONSOLIDATED_ERROR_GROUPS_CONFIG[key] for e in itertools.chain.from_iterable(iter(group.values()))] for key in CONSOLIDATED_ERROR_GROUPS_CONFIG
}


DETAILED_ERROR_GROUPS = DETECTOR_RESPONSE_ERROR_GROUPS.copy()
DETAILED_ERROR_GROUPS.update(GENIE_ERROR_GROUPS)
DETAILED_ERROR_GROUPS.update(FSI_ERROR_GROUPS)
DETAILED_ERROR_GROUPS.update(GEANT_ERROR_GROUPS)
#DETAILED_ERROR_GROUPS["MnvTunes"]=CONSOLIDATED_ERROR_GROUPS["MnvTunes"]


# ERROR_SUBGROUPS = {
# 	"Detector model": DETECTOR_RESPONSE_ERROR_GROUPS,
# 	"Flux model": FLUX_ERROR_GROUPS,
# 	"Interaction model": dict(GENIE_ERROR_GROUPS.items() + OTHER_INTERACTION_ERROR_GROUPS.items()),
# }

# systematics from universes generated from external calculations
# (e.g., dedicated variation samples).
# these should be in the form of histogram(s) in the same binning
# as the CV histogram, stored in a file somewhere

# EXTERNAL_UNIV_CONFIG = {
# 	"BkndShape": {
# 		"err_band_type": PlotUtils.MnvVertErrorBand,
# 		"histo_file": os.path.join( AnalysisConfig.PlotPath("BackgroundFitting"), AnalysisConfig.bknd_constraint_method, "background_prediction_histos.root" ),
# 		"histo_names": ["mc_bknd_pred",]
# 	},
# }

# DEFAULT_VERT_SYSTEMATICS = BKND_UNIV_SYSTEMATICS.copy()


# ########################

# def GroupErrors(plotter, grouping=CONSOLIDATED_ERROR_GROUPS):
# 	Vector_string_type = getattr(ROOT, "vector<string>")
	
# 	# so.......... I can't seem to get genreflex
# 	# to create the proper dictionaries to allow me
# 	# to interrogate the error_summary_group_map for its list of keys.
# 	# thus I will store them separately.
# 	if not hasattr(plotter, "error_group_list"):
# 		plotter.error_group_list = []

# #	for prefix, grouping in zip( ("truth_genie_wgt", "mc_wgt_Flux"), (GENIE_ERROR_GROUPS, FLUX_ERROR_GROUP) ):

# 	for group_name, error_list in grouping.iteritems():
# 		vector = Vector_string_type()
# 		for error_name in error_list:
# 			vector.push_back(error_name)
# 		plotter.error_summary_group_map[group_name] = vector
# 		plotter.error_group_list.append(group_name)

# # make an MnvPlotter that has the right error groups
# plotter = PlotUtils.MnvPlotter()
# GroupErrors(plotter)
# plotter.text_size = 0.03
# #plotter.legend_n_columns = 2
# plotter.legend_offset_x = 0.02
# plotter.legend_text_size = 0.03
# plotter.height_nspaces_per_hist = 1
# plotter.width_xspace_per_letter = 0.3
# plotter.legend_offset_y = -0.1

# # also make an MnvHistoConstrainer with the right error strategies
# constrainer = PlotUtils.MnvHistoConstrainer()
# constrainer.LoadConstraint("nu+e", os.path.expandvars("$MPARAMFILESROOT/data/FluxConstraints/nu+e.txt"))
# for err_name, strategy in ERROR_CONSTRAINT_STRATEGIES.iteritems():
# 	constrainer.SetSpectatorCorrectionStrategy(err_name, strategy)

