"""
   Define your signal and background here.

   author: H. Su
"""
from collections import OrderedDict
from config.CutConfig import KINEMATICS_CUTS
from tools.CutLibrary import CUTS
from tools import TruthTools

#helper functions for signal defination

def IsFiducial(event):
    return all(map((lambda cutname:CUTS[cutname].DoesEventPass(event)),["Truth_Vertex_Z","Truth_Vertex_Apothem"] ))

IsCC = lambda event: event.mc_current == 1
IsNC = lambda event: event.mc_current == 2

IsQE =  lambda event : event.mc_intType==1 and event.mc_charm!=1
IsDelta = lambda event : event.mc_intType == 2 and event.mc_resID==0
IsDIS = lambda event : event.mc_intType == 3
IsCoherent = lambda event: event.mc_intType == 4
IsElastic = lambda event: event.mc_intType == 7
Is2p2h = lambda event: event.mc_intType == 8

IsPC = lambda event: event.mc_processType ==5
IsUnknown  = lambda event : event.mc_intType == 10

IsNuE = lambda event: abs(event.mc_incoming) == 12
IsNuMu = lambda event: abs(event.mc_incoming) == 14
IsAntiNu = lambda event: event.mc_incoming < 0 # assuming only neutrino incomes
IsNue = lambda event: event.mc_incoming == 12
IsNueBar = lambda event: event.mc_incoming == -12
IsNumu = lambda event: event.mc_incoming == 14
IsNumuBar = lambda event: event.mc_incoming == -14
IsPi0InFinalState = lambda event: 111 in event.mc_FSPartPDG
IsChargedPionInFinalState = lambda event: 211 in event.mc_FSPartPDG or -211 in event.mc_FSPartPDG
IsProtonInFinalState = lambda event: 2212 in event.mc_FSPartPDG
IsMultiMeson = lambda event: countPDG(event.mc_FSPartPDG, [211,-211, 321,-321,323,-323,111,130,310,311])>1
IsSinglePion = lambda event: countPDG(event.mc_FSPartPDG, [211]) == 1 and countPDG(event.mc_FSPartPDG, [-211, 321,-321,323,-323,111,130,310,311]) == 0

IsHeavyBaryon = lambda event: 3112 in event.mc_FSPartPDG or 3122 in event.mc_FSPartPDG or 3212 in event.mc_FSPartPDG or 3222 in event.mc_FSPartPDG or 4112 in event.mc_FSPartPDG or 4122 in event.mc_FSPartPDG or 4212 in event.mc_FSPartPDG or 4222 in event.mc_FSPartPDG
IsMeson = lambda event: 211 in event.mc_FSPartPDG or -211 in event.mc_FSPartPDG or 321 in event.mc_FSPartPDG or -321 in event.mc_FSPartPDG or 323 in event.mc_FSPartPDG or -323 in event.mc_FSPartPDG  or 111 in event.mc_FSPartPDG or 130 in event.mc_FSPartPDG or 310 in event.mc_FSPartPDG or 311 in event.mc_FSPartPDG
IsNucleon = lambda event: 2212 in event.mc_FSPartPDG or 2112 in event.mc_FSPartPDG
IsDeexcitationPhoton =  lambda event: event.mc_FSPartPDG[0] == 22 and event.mc_FSPartE[0] < 10
IsPhoton = lambda event: 22 in event.mc_FSPartPDG and event.mc_FSPartPDG[0] != 22

def IsInKinematicPhaseSpace(event):
    return all(CUTS["True{}".format(cut)].DoesEventPass(event) for cut in KINEMATICS_CUTS)


# In case a event satisfy multiple definations, the first takes priority.
TRUTH_CATEGORIES = OrderedDict()
TRUTH_CATEGORIES["NCDiff"] = lambda event: IsUnknown(event)
TRUTH_CATEGORIES["NuEElastic"] = lambda event: IsElastic(event)
TRUTH_CATEGORIES["NonPhaseSpace"] = lambda event: IsCC(event) and IsNuMu(event) and not IsInKinematicPhaseSpace(event)

TRUTH_CATEGORIES["CCNuMuWrongSign"] = lambda event: IsCC(event) and IsNuMu(event) and not IsAntiNu(event)
TRUTH_CATEGORIES["CCNuMuQE"] = lambda event: IsCC(event) and IsNuMu(event) and IsQE(event)
TRUTH_CATEGORIES["CCNuMuDelta"] = lambda event: IsCC(event) and IsNuMu(event) and IsDelta(event)
TRUTH_CATEGORIES["CCNuMuDIS"] = lambda event: IsCC(event) and IsNuMu(event) and IsDIS(event)
TRUTH_CATEGORIES["CCNuMu2p2h"] = lambda event: IsCC(event) and IsNuMu(event) and Is2p2h(event)
TRUTH_CATEGORIES["CCNuMu"] = lambda event: IsCC(event) and IsNuMu(event)

# My signal is one or more of the listed categories.
SIGNAL_DEFINATION = [
    "CCNuMuQE",
    "CCNuMuDelta",
    "CCNuMuDIS",
    "CCNuMu",
    "CCNuMu2p2h",
    "CCNuMuWrongSign"
]

EXTRA_OTHER = [
]