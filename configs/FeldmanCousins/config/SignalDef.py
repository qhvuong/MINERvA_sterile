"""
   Define your signal and background here.

   author: H. Su
"""
from collections import OrderedDict
from config.CutConfig import KINEMATICS_CUTS
from tools.CutLibrary import CUTS
from tools import TruthTools
import math

#helper functions for signal defination
def countPDG(fsparticles,pdgs):
    return sum([1 if x in pdgs else 0 for x in fsparticles])

def IsFiducial(event):
    return all(map((lambda cutname:CUTS[cutname].DoesEventPass(event)),["Truth_Vertex_Z","Truth_Vertex_Apothem"] ))

def IsLowPhoton(event):
    for i in range(len(event.mc_FSPartPDG)):
        if event.mc_FSPartPDG[i] == 22:
            if event.mc_FSPartE[i] > 10:
                return False
    return True

def CCQESignal(event):
    if 13 in event.mc_FSPartPDG:
        muonloc = event.mc_FSPartPDG.index(13)
    elif -13 in event.mc_FSPartPDG:
        muonloc = event.mc_FSPartPDG.index(-13)
    else:
        return False
    muonphi = event.kin_cal.true_theta_lep_rad
    return (muonphi < 0.349066 and muonphi > -0.349066)


def SinglePionSignal(event):
    if 13 in event.mc_FSPartPDG:
        muonloc = event.mc_FSPartPDG.index(13)
    elif -13 in event.mc_FSPartPDG:
        muonloc = event.mc_FSPartPDG.index(-13)
    else:
        return False

    if countPDG(event.mc_FSPartPDG, [211]) == 1:
        piloc = event.mc_FSPartPDG.index(211)
    elif countPDG(event.mc_FSPartPDG, [-211]) == 1:
        piloc = event.mc_FSPartPDG.index(-211)
    else:
        return False

    muonphi = event.kin_cal.true_theta_lep_rad
    muon_p = math.sqrt(event.mc_FSPartPx[muonloc]**2 + event.mc_FSPartPy[muonloc]**2 + event.mc_FSPartPz[muonloc]**2)
    pionKE = event.mc_FSPartE[piloc] - 139.6
    Wexp = event.mc_w
    return (muonphi < 0.226893 and muonphi > -0.226893 and muon_p < 20e3 and muon_p > 1.5e3 and pionKE < 350 and pionKE > 35 and Wexp < 1.4e3)

def IsAaronSignal(event):
    passed = IsCC(event) and IsNumu(event) and IsSinglePion(event) and IsLowPhoton(event) and SinglePionSignal(event)
    return passed

def Is2DSignal(event):
    passed = IsCC(event) and IsNumu(event) and not IsHeavyBaryon(event) and not IsMeson(event) and IsNucleon(event) and CCQESignal(event) and IsLowPhoton(event)
    return passed

IsCC = lambda event: event.mc_current == 1
IsNC = lambda event: event.mc_current == 2

IsQE =  lambda event : event.mc_intType==1 and event.mc_charm!=1
IsDelta = lambda event : event.mc_intType == 2
IsDIS = lambda event : event.mc_intType == 3
IsCoherent = lambda event: event.mc_intType == 4
IsElastic = lambda event: event.mc_intType == 7
Is2p2h = lambda event: event.mc_intType == 8
IsNotNue = lambda event: abs(event.mc_incoming) != 12

IsPC = lambda event: event.mc_processType ==5
IsUnknown  = lambda event : event.mc_intType == 10


IsNuE = lambda event: event.mc_incoming == 12
IsNuMu = lambda event: abs(event.mc_incoming) == 14
IsAntiNu = lambda event: event.mc_incoming < 0 # assuming only neutrino incomes
IsNu = lambda event: event.mc_incoming > 0
IsNue = lambda event: event.mc_incoming == 12
IsNueBar = lambda event: event.mc_incoming == -12
IsAntiNuE = lambda event: event.mc_incoming == -12
IsNumu = lambda event: event.mc_incoming == 14
IsNumuBar = lambda event: event.mc_incoming == -14

IsPi0InFinalState = lambda event: 111 in event.mc_FSPartPDG
IsChargedPionInFinalState = lambda event: 211 in event.mc_FSPartPDG or -211 in event.mc_FSPartPDG
IsProtonInFinalState = lambda event: 2212 in event.mc_FSPartPDG
IsMultiMeson = lambda event: countPDG(event.mc_FSPartPDG, [211,-211, 321,-321,323,-323,111,130,310,311])>1
IsSinglePion = lambda event: countPDG(event.mc_FSPartPDG, [211]) == 1 and countPDG(event.mc_FSPartPDG, [-211, 321,-321,323,-323,111,130,310,311]) == 0
IsNucleon = lambda event: 2212 in event.mc_FSPartPDG or 2112 in event.mc_FSPartPDG


IsHeavyBaryon = lambda event: 3112 in event.mc_FSPartPDG or 3122 in event.mc_FSPartPDG or 3212 in event.mc_FSPartPDG or 3222 in event.mc_FSPartPDG or 4112 in event.mc_FSPartPDG or 4122 in event.mc_FSPartPDG or 4212 in event.mc_FSPartPDG or 4222 in event.mc_FSPartPDG
IsMeson = lambda event: 211 in event.mc_FSPartPDG or -211 in event.mc_FSPartPDG or 321 in event.mc_FSPartPDG or -321 in event.mc_FSPartPDG or 323 in event.mc_FSPartPDG or -323 in event.mc_FSPartPDG  or 111 in event.mc_FSPartPDG or 130 in event.mc_FSPartPDG or 310 in event.mc_FSPartPDG or 311 in event.mc_FSPartPDG
IsDeexcitationPhoton =  lambda event: event.mc_FSPartPDG[0] == 22 and event.mc_FSPartE[0] < 10
IsPhoton = lambda event: 22 in event.mc_FSPartPDG and event.mc_FSPartPDG[0] != 22

def IsInKinematicPhaseSpace(event):
    return all(CUTS["True{}".format(cut)].DoesEventPass(event) for cut in KINEMATICS_CUTS)


# In case a event satisfy multiple definations, the first takes priority.
TRUTH_CATEGORIES = OrderedDict()
TRUTH_CATEGORIES["NCDiff"] = lambda event: IsUnknown(event)
TRUTH_CATEGORIES["NuEElastic"] = lambda event: IsElastic(event)
TRUTH_CATEGORIES["NonPhaseSpace"] = lambda event: IsCC(event) and (IsNuE(event) or IsNue(event)) and not IsInKinematicPhaseSpace(event)

TRUTH_CATEGORIES["CCNuEWrongSign"] = lambda event: IsCC(event) and IsNue(event)
TRUTH_CATEGORIES["CCNuEQE"] = lambda event: IsCC(event) and IsNuE(event) and IsQE(event)
TRUTH_CATEGORIES["CCNuEDelta"] = lambda event: IsCC(event) and IsNuE(event) and IsDelta(event)
TRUTH_CATEGORIES["CCNuEDIS"] = lambda event: IsCC(event) and IsNuE(event) and IsDIS(event)
TRUTH_CATEGORIES["CCNuE2p2h"] = lambda event: IsCC(event) and IsNuE(event) and Is2p2h(event)
TRUTH_CATEGORIES["CCNuE"] = lambda event: IsCC(event) and IsNuE(event)

TRUTH_CATEGORIES["CCPi0"] = lambda event: IsCC(event) and IsPi0InFinalState(event) and IsNotNue(event)
TRUTH_CATEGORIES["NCCohPi0"] = lambda event: IsCoherent(event) and IsNC(event) and IsPi0InFinalState(event)
TRUTH_CATEGORIES["NCPi0"] = lambda event: IsNC(event) and IsPi0InFinalState(event)
TRUTH_CATEGORIES["CCPi"] = lambda event: IsCC(event) and IsChargedPionInFinalState(event)
TRUTH_CATEGORIES["NCPi"] = lambda event: IsNC(event) and IsChargedPionInFinalState(event)
TRUTH_CATEGORIES["CCOther"] = lambda event: IsCC(event)
TRUTH_CATEGORIES["NCOther"] = lambda event: IsNC(event)

TRUTH_CATEGORIES["NC"] = lambda event: IsNC(event)
TRUTH_CATEGORIES["CCWrongSign"] = lambda event: IsCC(event) and IsNu(event)
#TRUTH_CATEGORIES["NonFiducial"] = lambda event: not IsFiducial(event)
TRUTH_CATEGORIES["CCQE"] = lambda event: IsQE(event)
TRUTH_CATEGORIES["CCDelta"] = lambda event: IsDelta(event)
TRUTH_CATEGORIES["CC2p2h"] = lambda event: Is2p2h(event)
TRUTH_CATEGORIES["CCDIS"] = lambda event: IsDIS(event)

SIGNAL_DEFINITION = [
    "CCNuEQE",
    "CCNuEDelta",
    "CCNuEDIS",
    "CCNuE",
    "CCNuE2p2h",
    "CCNuEWrongSign"
]
SWAP_SIGNAL_DEFINITION = [
    "CCNu",
    "CCNuEQE",
    "CCNuEDelta",
    "CCNuEDIS",
    "CCNuE",
    "CCNuE2p2h",
    "CCNuEAntiNu",
]
SWAP_SIGNAL_DEFINITION = [
    "CCNu",
    "CCNuE",
    "CCNuEQE",
    "CCNuEDelta",
    "CCNuEDIS",
    "CCNuE",
    "CCNuE2p2h",
    "CCNuEAntiNu",
    "CCNuEWrongSign",
]

#SIGNAL_DEFINITION = [
#    "CCQE",
#    "SinglePion"
#]

EXTRA_OTHER = [
]
