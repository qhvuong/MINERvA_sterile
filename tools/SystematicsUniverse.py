"""
    SystematicsUniverse.py:
    Define systematics shift universe. 
"""

import math
import os,sys
import random
import itertools
import operator
import random
import re
import numpy as np
import math
import ROOT
import PlotUtils
import ctypes

from config import SystematicsConfig
from tools import Utilities
from collections import OrderedDict
from tools import pcweight
#from .pcweight import GetModelWeight,MyWeighter
from tools import AvailableECorrection


OneSigmaShift = [-1.0,1.0]
M_e = 0.511 # MeV
M_mu =  105.658 #MeV
M_e_sqr = M_e**2
M_mu_sqr = M_mu**2

USE_MNVTUNE_V1 = False
USE_MNVTUNE_V4P3P0 = False
USE_MNVTUNE_V4P4P1 = True


V1_ADD_DEUTERIUM = False
V1_ADD_ZEXP = False

V1_ADD_LOWQ2PI = False
V1_LOWQ2PI_CHANNEL = None

V1_ADD_COH = False          # ME coherent tune implementation
V1_ADD_LE_COH = False        # LE coherent tune implementation
V1_ADD_DIFFCOH = False

# MnvTune v4.3.0
V4P3P0_ADD_LOWQ2PI = True
V4P3P0_LOWQ2PI_CHANNEL = "MENU1PI"

V4P3P0_ADD_ME_COH = True
V4P3P0_ADD_DIFFRACTIVE = True

# Detector weights, separately configurable
V4P3P0_ADD_GEANT_HADRON = True
V4P3P0_ADD_MINOS_EFFICIENCY = True


# v4.4.1 generator corrections
V4P4P1_ADD_ZEXP = True
V4P4P1_ADD_LOWQ2PI = False
V4P4P1_LOWQ2PI_CHANNEL = None
V4P4P1_ADD_LE_COH = False
V4P4P1_ADD_ME_COH = False
V4P4P1_ADD_DIFFRACTIVE = False
V4P4P1_ADD_ELASTIC_FSI = True


TEST_TAG = "V1D2_LECOH_DiffCOH_noLowQ2"



#Alias:
#ROOT.PythonMinervaUniverse = PlotUtils.DefaultCVPythonUniverse

#MethodName: UpperCaseCamel
#universeShortName: lower_case_underscore

def _clamp(x, lo=-1.0, hi=1.0):
    return lo if x < lo else hi if x > hi else x

def _hex_edge_distance_mm(x, y, apothem_mm=881.25):
    """
    Distance (mm) from point (x,y) to the nearest side of a regular hexagon centered at (0,0) with given apothem.
    Positive inside, 0 on boundary, negative outside.
    """
    sqrt3 = math.sqrt(3.0)
    u1 = x
    u2 = 0.5*x + 0.5*sqrt3*y
    u3 = -0.5*x + 0.5*sqrt3*y
    umax = max(abs(u1), abs(u2), abs(u3))
    return apothem_mm - umax



def apply_theta_bias_correction_det(
    p3_det,
    vtxX_mm: float,
    vtxY_mm: float,
    t0x: float,
    t1x: float,
    t0y: float,
    t1y: float,
    *,
    clamp_vtx: bool = True,
    vtx_min_mm: float = -1200.0,
    vtx_max_mm: float = 1200.0,
    keep_pmag: bool = True,
    max_abs_angle_rad = None,
    make_vec=None,
):
    """
    Correct direction using angular biases (projected angles):
      <dThetaX>(vtxX) = t0x + t1x*vtxX_mm
      <dThetaY>(vtxY) = t0y + t1y*vtxY_mm

    Definitions:
      thetaX = atan2(px, pz)
      thetaY = atan2(py, pz)

    After correcting thetaX/thetaY, reconstruct a UNIT direction using:
      tx = tan(thetaX_corr), ty = tan(thetaY_corr)
      u  ∝ (tx, ty, 1), then normalize.

    Returns:
      (p3_corr, thetaX_corr, thetaY_corr, theta_corr)
    where theta_corr is the 3D polar angle: atan2(sqrt(px^2+py^2), pz).
    """
    if p3_det is None:
        return None, None, None, None

    px, py, pz = float(p3_det.X()), float(p3_det.Y()), float(p3_det.Z())
    pmag = float(p3_det.R())
    if pmag <= 0:
        return p3_det, None, None, None

    if clamp_vtx:
        vtxX_mm = max(vtx_min_mm, min(vtx_max_mm, vtxX_mm))
        vtxY_mm = max(vtx_min_mm, min(vtx_max_mm, vtxY_mm))

    # Current projected angles (radians)
    thetaX = math.atan2(px, pz)
    thetaY = math.atan2(py, pz)

    # Bias models (radians)
    bx = t0x + t1x * vtxX_mm
    by = t0y + t1y * vtxY_mm

    # Apply correction (subtract mean bias)
    thetaX_corr = thetaX - bx
    thetaY_corr = thetaY - by

    # Optional safety clamp on corrected angles
    if max_abs_angle_rad is not None:
        thetaX_corr = max(-max_abs_angle_rad, min(max_abs_angle_rad, thetaX_corr))
        thetaY_corr = max(-max_abs_angle_rad, min(max_abs_angle_rad, thetaY_corr))

    # Rebuild a UNIT direction from corrected projected angles
    tx = math.tan(thetaX_corr)
    ty = math.tan(thetaY_corr)

    inv_norm = 1.0 / math.sqrt(1.0 + tx*tx + ty*ty)
    ux = tx * inv_norm
    uy = ty * inv_norm
    uz = 1.0 * inv_norm

    # 3D polar angle from corrected direction
    theta_corr = math.atan2(math.sqrt(ux*ux + uy*uy), uz)

    # Scale to momentum magnitude (or not)
    scale = pmag if keep_pmag else 1.0
    px_c, py_c, pz_c = scale * ux, scale * uy, scale * uz

    if make_vec is None:
        try:
            p3_corr = p3_det.__class__(px_c, py_c, pz_c)
        except Exception:
            raise RuntimeError("apply_theta_bias_correction_det: please pass make_vec (e.g. ROOT.Math.XYZVector)")
    else:
        p3_corr = make_vec(px_c, py_c, pz_c)

    return p3_corr, thetaX_corr, thetaY_corr, theta_corr



def apply_pmag_frac_correction_det(
    p3_det,
    vtxY_mm: float,
    a0: float,
    a1: float,
    *,
    clamp_vtx: bool = True,
    vtx_min_mm: float = -1200.0,
    vtx_max_mm: float = 1200.0,
    max_abs_frac: float = 0.2,   # safety clamp: 20% by default
    make_vec=None,
):
    """
    Correct |p| bias using a linear model in vtxY:
      <dPmagFrac>(vtxY) = a0 + a1*vtxY_mm
    where dPmagFrac = (|p|reco - |p|true)/|p|true

    We correct reco magnitude:
      |p|corr = |p|reco / (1 + bias)

    Direction is preserved.
    """
    if p3_det is None:
        return None

    pmag = p3_det.R()
    if pmag <= 0:
        return p3_det

    if clamp_vtx:
        vtxY_mm = max(vtx_min_mm, min(vtx_max_mm, vtxY_mm))

    bias = a0 + a1 * vtxY_mm  # dimensionless
    # optional clamp to avoid crazy corrections
    if max_abs_frac is not None:
        bias = max(-max_abs_frac, min(max_abs_frac, bias))

    denom = 1.0 + bias
    if denom <= 0:
        return p3_det

    scale = 1.0 / denom
    px_c = p3_det.X() * scale
    py_c = p3_det.Y() * scale
    pz_c = p3_det.Z() * scale

    if make_vec is None:
        try:
            return p3_det.__class__(px_c, py_c, pz_c)
        except Exception:
            raise RuntimeError("apply_pmag_frac_correction_det: pass make_vec (e.g. ROOT.Math.XYZVector)")
    return make_vec(px_c, py_c, pz_c)


_flux_cache = {}

def _get_flux_hist(path, hname="flux_E_cvweighted"):
    key = (path, hname)
    if key not in _flux_cache:
        f = ROOT.TFile.Open(path)
        if not f or f.IsZombie():
            raise RuntimeError(f"Could not open flux file: {path}")
        h = f.Get(hname)
        if not h:
            raise RuntimeError(f"Could not find histogram {hname} in {path}")
        h = h.Clone(f"{h.GetName()}_{len(_flux_cache)}")
        h.SetDirectory(0)
        _flux_cache[key] = (f, h)
    return _flux_cache[key][1]

def get_flux_ratio_me_to_le(enu_gev, pdg, me_flux_tag,
                            flux_dir="/exp/minerva/app/users/qvuong/MAT_AL9/CC-NuE-XSec/custom_plotutils/data/flux",
                            min_den=1e-12,
                            max_ratio=10.0):
    pdg = int(pdg)

    h_le = _get_flux_hist(
        f"{flux_dir}/flux-gen2thin-pdg{pdg}-minerva1.root",
        "flux_E_cvweighted"
    )
    h_me = _get_flux_hist(
        f"{flux_dir}/flux-gen2thin-pdg{pdg}-minervame{me_flux_tag}.root",
        "flux_E_cvweighted"
    )

    x_le = h_le.GetXaxis()
    x_me = h_me.GetXaxis()

    # Reject energies outside histogram range instead of clamping
    if enu_gev < x_le.GetXmin() or enu_gev > x_le.GetXmax():
        return 0.0
    if enu_gev < x_me.GetXmin() or enu_gev > x_me.GetXmax():
        return 0.0

    b_le = x_le.FindBin(enu_gev)
    b_me = x_me.FindBin(enu_gev)

    if b_le < 1 or b_le > h_le.GetNbinsX():
        return 0.0
    if b_me < 1 or b_me > h_me.GetNbinsX():
        return 0.0

    le_val = h_le.GetBinContent(b_le)
    me_val = h_me.GetBinContent(b_me)

    if me_val <= min_den:
        return 0.0

    ratio = le_val / me_val

    if max_ratio is not None:
        ratio = min(ratio, max_ratio)

    return ratio

_trueEL_cache = {}
_trueEL_pdf_cache = {}

def _get_trueEL_hist(path, hname="ETrue_Length_CCNuE"):
    key = (path, hname)
    if key not in _trueEL_cache:
        f = ROOT.TFile.Open(path)
        if not f or f.IsZombie():
            raise RuntimeError(f"Could not open trueE-trueL file: {path}")
        h = f.Get(hname)
        if not h:
            raise RuntimeError(f"Could not find histogram {hname} in {path}")
        h = h.Clone(f"{h.GetName()}_{len(_trueEL_cache)}")
        h.SetDirectory(0)
        _trueEL_cache[key] = (f, h)
    return _trueEL_cache[key][1]

def GetNeutrinoTravelledLength(enu_gev, pdg,
                          nue_file="/exp/minerva/data/users/qvuong/nu_e/kin_dist_mcleFHC_CCnue_trueEL_fullStatsFluxes_MAD.root",
                          numu_file="/exp/minerva/data/users/qvuong/nu_mu/kin_dist_mcleFHC_CCnumu_trueEL_fullStatsFluxes_MAD.root",
                          nue_hname="ETrue_Length_CCNuE",
                          numu_hname="ETrue_Length_CCNuMu"):
    """
    Return a sampled true L for a neutrino with true energy enu_gev and species pdg.

    Uses:
      - nue template file/hist for |pdg| == 12
      - numu template file/hist for |pdg| == 14

    Sampling is done from the true-L distribution in the corresponding true-E slice.
    """

    apdg = abs(int(pdg))
    if apdg == 12:
        h2 = _get_trueEL_hist(nue_file, nue_hname)
    elif apdg == 14:
        h2 = _get_trueEL_hist(numu_file, numu_hname)
    else:
        raise ValueError(f"Unsupported PDG for GetNeutrinoTravelledLength: {pdg}")

    by = h2.GetYaxis().FindBin(enu_gev)
    by = max(1, min(by, h2.GetNbinsY()))

    hL = h2.ProjectionX(f"hL_tmp_{apdg}_{by}_{ROOT.TUUID().AsString()}", by, by)
    hL.SetDirectory(0)

    if hL.Integral() <= 0:
        by1 = max(1, by - 1)
        by2 = min(h2.GetNbinsY(), by + 1)
        hL = h2.ProjectionX(f"hL_tmp_wide_{apdg}_{by1}_{by2}_{ROOT.TUUID().AsString()}", by1, by2)
        hL.SetDirectory(0)

    if hL.Integral() <= 0:
        hL = h2.ProjectionX(f"hL_tmp_all_{apdg}_{ROOT.TUUID().AsString()}", 1, h2.GetNbinsY())
        hL.SetDirectory(0)

    if hL.Integral() <= 0:
        raise RuntimeError(f"No valid L distribution found for pdg={pdg}, Enu={enu_gev}")

    return hL.GetRandom()

def GetNeutrinoTravelledLengthPDF(
    enu_gev,
    pdg,
    nue_file=f"/exp/minerva/data/users/qvuong/nu_e/kin_dist_mcleFHC_CCnue_trueEL_{TEST_TAG}_MAD.root",
    numu_file=f"/exp/minerva/data/users/qvuong/nu_mu/kin_dist_mcleFHC_CCnumu_trueEL_{TEST_TAG}_MAD.root",
    nuebar_file=f"/exp/minerva/data/users/qvuong/antinu_e/kin_dist_mcle5_CCnuebar_trueEL_{TEST_TAG}_MAD.root",
    numubar_file=f"/exp/minerva/data/users/qvuong/antinu_mu/kin_dist_mcle5_CCnumubar_trueEL_{TEST_TAG}_MAD.root",
    nue_hname="ETrue_Length_CCNuE",
    numu_hname="ETrue_Length_CCNuMu",
):
    """
    Return a normalized TH1 of P(L | Enu, pdg).

    Uses separate files for particles and antiparticles:
      pdg =  12  -> nue_file,     nue_hname
      pdg = -12  -> nuebar_file,  nue_hname
      pdg =  14  -> numu_file,    numu_hname
      pdg = -14  -> numubar_file, numu_hname

    The returned histogram is cached by file, hist name, signed PDG, and E bin.
    Treat the returned histogram as read-only.
    """

    pdg = int(pdg)

    if pdg == 12:
        this_file = nue_file
        this_hname = nue_hname
    elif pdg == -12:
        this_file = nuebar_file
        this_hname = nue_hname
    elif pdg == 14:
        this_file = numu_file
        this_hname = numu_hname
    elif pdg == -14:
        this_file = numubar_file
        this_hname = numu_hname
    else:
        raise ValueError(
            "Unsupported PDG for GetNeutrinoTravelledLengthPDF: {}".format(pdg)
        )

    h2 = _get_trueEL_hist(this_file, this_hname)

    by = h2.GetYaxis().FindBin(enu_gev)
    by = max(1, min(by, h2.GetNbinsY()))

    cache_key = (this_file, this_hname, pdg, by)

    if cache_key in _trueEL_pdf_cache:
        return _trueEL_pdf_cache[cache_key]

    hL = h2.ProjectionX(
        "hL_pdf_{}_{}_{}".format(pdg, by, ROOT.TUUID().AsString()),
        by,
        by,
    )
    hL.SetDirectory(0)

    if hL.Integral() <= 0:
        hL.Delete()

        by1 = max(1, by - 1)
        by2 = min(h2.GetNbinsY(), by + 1)

        hL = h2.ProjectionX(
            "hL_pdf_wide_{}_{}_{}_{}".format(
                pdg,
                by1,
                by2,
                ROOT.TUUID().AsString(),
            ),
            by1,
            by2,
        )
        hL.SetDirectory(0)

    if hL.Integral() <= 0:
        hL.Delete()

        hL = h2.ProjectionX(
            "hL_pdf_all_{}_{}".format(pdg, ROOT.TUUID().AsString()),
            1,
            h2.GetNbinsY(),
        )
        hL.SetDirectory(0)

    integ = hL.Integral()
    if integ <= 0:
        hL.Delete()
        raise RuntimeError(
            "No valid L PDF found for pdg={}, Enu={}, file={}, hist={}".format(
                pdg,
                enu_gev,
                this_file,
                this_hname,
            )
        )

    hL.Scale(1.0 / integ)
    _trueEL_pdf_cache[cache_key] = hL

    return hL



# The base universe, define functions to be used for all universes.
class CVPythonUniverse():
    is_pc = False
    def __init__(self, chain, nsigma = None): #nsigma is None for data because we don't shift data
        self.weight = None
        self.tuning_weight = None
        self.InitWithoutSuper(chain,nsigma)

    def InitWithoutSuper(self,chain,nsigma):
        self.weight = None
        self.tuning_weight = None
        self.nEntries = chain.GetEntries()
        self.mc = nsigma is not None
        self.chain=chain
        self.pcweighter = pcweight.MyWeighter

    def ResetWeight(self):
        self.weight = None

    def _is_truth_mode(self):
        """Truth ntuples lack reco branches (e.g. HasNoBackExitingTracks)."""
        try:
            if self.IsTruth():
                return True
        except Exception:
            pass
        try:
            tree = self.chain.GetTree()
            return tree is not None and tree.GetName() == "Truth"
        except Exception:
            return False

    def HasInclusiveMINOSMatch(self):
        """
        Inclusive CC MINOS matching definition.
        Matches MINERvA-101-Cross-Section:
            IsMinosMatchMuon()
            -> GetIsMinosMatchTrack()
            -> isMinosMatchTrack == 1
        """
        return self.GetInt("isMinosMatchTrack") == 1

    # magic function that allow direct access to TBranch by Universe.TBranch, but prohibit any python magic variable being generated.
    def __getattr__(self,attrName):
        def makelist(x):
            return list(makelist(_) for _ in x) if hasattr(x, "__getitem__") else x

        if attrName == "kin_cal" or attrName == "classifier":
            #the two attrs can't be in chainwrapper 
            raise AttributeError()

        if attrName ==  "pcweighter":
            print (self.ShortName())

        # determine the dimension of leaf.
        try:
            return makelist(self.LeafGetters[attrName](attrName))
        except KeyError:
            #havn't see this function. Got to figure out how to get it.
            i = ctypes.c_int(1)
            leaf = self.chain.GetLeaf(attrName,-1)
            type = leaf.GetTypeName()
            if leaf.GetLeafCounter(i) or i.value > 1:
                dimension=1
            else:
                dimension=type.count("vector")
            # Assuming type names are Double_t,Int_t,Bool_t,vector<vector<Double_t>>, use regex to match the type
            regex="(<|^)([A-Za-z]+)(_t|>|$)"
            type = re.search(regex,type).group(2).capitalize()

            self.LeafGetters[attrName] = getattr(self, "Get{}{}".format("" if dimension==0 else ("Vec" if dimension==1 else "VecOfVec"),type))

        return makelist(self.LeafGetters[attrName](attrName))

    # return text when you print this class object. must reimplement or __getattr__ will hijack the method resolution
    def __repr__(self):
        return self.LatexName()+" Universe at sigma: "+ str(self.GetSigma())

    @property
    def n_prongs(self):
        return self.GetInt("n_prongs")

    # keep a reference to kinematicsCalulator and eventClassifier such that other function can pull info from universe
    def LoadTools(self,kinematicCalculator,eventClassifier):
        self.kin_cal = kinematicCalculator
        self.classifier = eventClassifier

    def SetLeptonType(self):
        if abs(SystematicsConfig.AnaNuPDG) == 12: 
            self.LeptonTheta       = self.ElectronTheta ## <===== fixed beam 3D angle
            self.LeptonEnergy      = self.ElectronEnergy
            self.LeptonP3D_det     = self.ElectronP3D_det
            self.LeptonP3D         = self.ElectronP3D
            # self.LeptonThetaX  = self.ElectronThetaX
            # self.LeptonThetaY  = self.ElectronThetaY
            # self.LeptonTheta2D = self.ElectronTheta2D ## <===== fixed beam 2D angle
            # self.LeptonPhi     = self.ElectronPhi
            # self.LeptonTheta   = self.ElectronTheta_vtxcorr ## <===== vertex corrected 3D angle
            # self.LeptonThetaX  = self.ElectronThetaX_vtxcorr
            # self.LeptonThetaY  = self.ElectronThetaY_vtxcorr
            # self.LeptonTheta2D = self.ElectronTheta2D_vtxcorr ## <===== vertex corrected 2D angle
            # self.LeptonPhi     = self.ElectronPhi_vtxcorr

            # self.Vertex3D_beam     = self.Vertex_beam
            # self.LeptonThetaY_det  = self.ElectronThetaY_det
            # self.LeptonThetaY_beam = self.ElectronThetaY_beam
            self.M_lep_sqr         =  M_e_sqr
            self.GetCorrection     = self.GetLeakageCorrection
            return True
        else:
            self.LeptonTheta       = self.GetThetamu        # angle in BEAM coord
            self.LeptonEnergy      = self.GetEmu
            self.LeptonP3D_det     = self.MuonP3D_det
            self.LeptonP3D         = self.MuonP3D_beam
            # # self.LeptonP3D         = self.GetPmu
            # # self.Vertex3D_beam     = self.Vertex_beam
            self.M_lep_sqr         =  M_mu_sqr
            self.GetCorrection     = self.GetNuEFuzz
            return False

    def New_True_L_Over_E(self):
        if "_cached_new_true_l_over_e" in self.__dict__:
            return self.__dict__["_cached_new_true_l_over_e"]

        enu_gev = self.mc_incomingE * 1e-3
        L = GetNeutrinoTravelledLength(enu_gev, self.mc_incoming)
        val = L / enu_gev if enu_gev > 0 else 0.0
        self.__dict__["_cached_new_true_l_over_e"] = val
        return val

    @property
    def nsigma(self):
        return self.GetSigma() if self.mc else None

    def GetWeight(self,bkgtuning=True):
        if self.weight is None:
            if self.nsigma is None:
                self.weight = 1
            else:
                self.weight = self.GetModelWeight() if self.is_pc else self.GetStandardWeight()
        if bkgtuning and self.nsigma is None:
            if self.tuning_weight is None:
                self.tuning_weight = self.pcweighter.GetWeight(self)
            return self.weight *self.tuning_weight
        else:
            return self.weight




    def DebugGenieCV(self):
        import math

        if self.ShortName() != "cv":
            return

        if not hasattr(CVPythonUniverse, "_debug_genie_counts"):
            CVPythonUniverse._debug_genie_counts = {
                "total": 0,
                "nonrespi": 0,
                "ccres": 0,
                "nonrespi_and_ccres": 0,
                "nonfinite_mat": 0,
                "nonfinite_expected": 0,
                "mismatch": 0,
                "sum_v1_nonrespi": 0.0,
                "sum_v1_deuterium": 0.0,
                "sum_mat": 0.0,
            }

        c = CVPythonUniverse._debug_genie_counts
        c["total"] += 1

        int_type = int(self.mc_intType)
        current = int(self.mc_current)

        rvn = self.GetVecElem("truth_genie_wgt_Rvn1pi", 2)
        rvp = self.GetVecElem("truth_genie_wgt_Rvp1pi", 2)

        is_nonrespi = rvn < 1.0 or rvp < 1.0
        is_ccres = int_type == 2 and current == 1

        if is_nonrespi:
            c["nonrespi"] += 1

        if is_ccres:
            c["ccres"] += 1

        if is_nonrespi and is_ccres:
            c["nonrespi_and_ccres"] += 1

        # V1 baseline: NonResPi only.
        w_nonrespi = 0.43 if is_nonrespi else 1.0

        # Deuterium resonance correction:
        # MAT: GENIE_MaRES=1.12, DEUTERIUM_MaRES=0.94,
        #      GENIE_MaRES_1Sig=0.22, RES_NORM=1.15
        w_deuterium = 1.0
        ma_res = None

        if is_ccres:
            ma_res = self.GetVecElem("truth_genie_wgt_MaRES", 2)
            w_deuterium = 1.15 * (
                1.0 + (0.18 / 0.22) * (ma_res - 1.0)
            )

        w_v1 = w_nonrespi
        w_v1_d2 = w_nonrespi * w_deuterium

        # Actual MAT implementation, with the active C++ switches.
        w_mat = self.GetGenieWeight()

        if not math.isfinite(w_mat):
            c["nonfinite_mat"] += 1
            if c["nonfinite_mat"] <= 10:
                print(
                    "[GENIE NONFINITE MAT]",
                    f"intType={int_type}",
                    f"current={current}",
                    f"mat={w_mat}",
                    flush=True
                )
            return

        if not (math.isfinite(w_v1) and math.isfinite(w_v1_d2)):
            c["nonfinite_expected"] += 1
            return

        c["sum_v1_nonrespi"] += w_v1
        c["sum_v1_deuterium"] += w_v1_d2
        c["sum_mat"] += w_mat

        # Match the currently configured C++ GENIE switches.
        if USE_MNVTUNE_V1:
            use_deuterium = V1_ADD_DEUTERIUM
            use_zexp = V1_ADD_ZEXP

        elif USE_MNVTUNE_V4P4P1:
            use_deuterium = True
            use_zexp = V4P4P1_ADD_ZEXP

        else:
            use_deuterium = True
            use_zexp = V1_ADD_ZEXP

        expected = w_v1_d2 if use_deuterium else w_v1

        if use_zexp:
            expected *= self.GetZExpWeight()

        if not math.isclose(
            w_mat, expected, rel_tol=1e-6, abs_tol=1e-8
        ):
            c["mismatch"] += 1

            if c["mismatch"] <= 20:
                print(
                    "[GENIE CV MISMATCH]",
                    f"intType={int_type}",
                    f"current={current}",
                    f"rvn={rvn}",
                    f"rvp={rvp}",
                    f"MaRES={ma_res}",
                    f"NonResPi={w_nonrespi}",
                    f"D2={w_deuterium}",
                    f"expected={expected}",
                    f"MAT={w_mat}",
                    flush=True
                )

        if c["total"] % 10000 == 0:
            print("[GENIE AUDIT PROGRESS]", c, flush=True)



    def Debug2p2hCV(self):
        import math

        if self.ShortName() != "cv":
            return

        # CV correction only affects MEC on nuclear targets.
        if self.mc_intType != 8 or self.mc_targetZ < 2:
            return

        if not hasattr(CVPythonUniverse, "_debug_2p2h_counts"):
            CVPythonUniverse._debug_2p2h_counts = {
                "total_mec": 0,
                "old_guard_pass": 0,
                "old_guard_fail": 0,
                "invalid_Q2": 0,
                "invalid_W": 0,
                "invalid_q0q3": 0,
                "mat_exceptions": 0,
                "mat_nonfinite": 0,
                "mat_nonunity": 0,
                "changed_weight": 0,
                "sum_mat": 0.0,
                "sum_current": 0.0,
            }

        counts = CVPythonUniverse._debug_2p2h_counts
        counts["total_mec"] += 1

        q2 = self.mc_Q2 * 1e-6
        W = self.mc_w * 1e-3

        q0 = self.Getq0True() * 1e-3
        q3 = self.Getq3True() * 1e-3

        # Exactly reproduce the current Python safeguard.
        valid_Q2 = math.isfinite(q2) and q2 >= 0.0
        valid_W = math.isfinite(W) and abs(W) < 90.0

        old_guard = valid_Q2 and valid_W

        if old_guard:
            counts["old_guard_pass"] += 1
        else:
            counts["old_guard_fail"] += 1

        if not valid_Q2:
            counts["invalid_Q2"] += 1

        if not valid_W:
            counts["invalid_W"] += 1

        # MAT evaluates its Gaussian using true q0 and q3.
        if not (math.isfinite(q0) and math.isfinite(q3)):
            counts["invalid_q0q3"] += 1

            if counts["invalid_q0q3"] <= 10:
                print(
                    "[2P2H INVALID Q0/Q3]",
                    "Q2 =", q2,
                    "W =", W,
                    "q0 =", q0,
                    "q3 =", q3,
                    flush=True
                )
            return

        # Compare against MAT directly, without changing CV weights.
        try:
            w_mat = self.GetLowRecoil2p2hWeight()
        except Exception as exc:
            counts["mat_exceptions"] += 1

            if counts["mat_exceptions"] <= 10:
                print("[2P2H MAT EXCEPTION]", exc, flush=True)
            return

        if not math.isfinite(w_mat):
            counts["mat_nonfinite"] += 1

            if counts["mat_nonfinite"] <= 10:
                print(
                    "[2P2H NONFINITE MAT]",
                    "Q2 =", q2,
                    "W =", W,
                    "q0 =", q0,
                    "q3 =", q3,
                    "w_mat =", w_mat,
                    flush=True
                )
            return

        # Reproduce what our existing Python code applies.
        w_current = w_mat if old_guard else 1.0

        counts["sum_mat"] += w_mat
        counts["sum_current"] += w_current

        if abs(w_mat - 1.0) > 1e-10:
            counts["mat_nonunity"] += 1

        if abs(w_mat - w_current) > 1e-10:
            counts["changed_weight"] += 1

            if counts["changed_weight"] <= 20:
                print(
                    "[2P2H GUARD DIFFERENCE]",
                    "Q2 =", q2,
                    "W =", W,
                    "q0 =", q0,
                    "q3 =", q3,
                    "MAT =", w_mat,
                    "CURRENT =", w_current,
                    flush=True
                )

        # Frequent progress reports for small diagnostic runs.
        if counts["total_mec"] % 100 == 0:
            print(
                "[2P2H AUDIT PROGRESS]",
                counts,
                flush=True
            )


    def GetStandardWeight(self):

        active_tunes = [
            USE_MNVTUNE_V1,
            USE_MNVTUNE_V4P3P0,
            USE_MNVTUNE_V4P4P1,
        ]

        if sum(active_tunes) > 1:
            raise RuntimeError(
                "Only one MINERvA tune can be active."
            )

        if USE_MNVTUNE_V1:
            return self.GetStandardWeightV1()

        if USE_MNVTUNE_V4P3P0:
            return self.GetStandardWeightV4p3p0()

        if USE_MNVTUNE_V4P4P1:
            return self.GetStandardWeightV4p4p1()

        # Leave existing legacy v4 implementation below.
        weight = 1.0
        weight *= self.GetGenieWeight()

        ### take care of flavor swapped sample below ###
        pdg = self.mc_incoming
        newpdg = pdg
        if SystematicsConfig.USE_SWAPPED:
            if pdg == 12:
                newpdg = 14
            elif pdg == -12:
                newpdg = -14
            elif pdg == 14:
                newpdg = 12
            elif pdg == -14:
                newpdg = -12
            else:
                print("Not neutrino PDG? {}".format(pdg))
        pdg = newpdg


        # enu = self.mc_incomingE * 1e-3
        # flux_cv = self.GetFluxAndCVWeight(enu, pdg)

        # # ---- DEBUG: true-E binned flux CV summary ----
        # if self.ShortName() == "cv":
        #     if not hasattr(CVPythonUniverse, "_trueE_flux_summary"):
        #         CVPythonUniverse._trueE_flux_summary = {}

        #     ebin = math.floor(enu / 0.5) * 0.5
        #     key = (pdg, ebin)

        #     if key not in CVPythonUniverse._trueE_flux_summary:
        #         CVPythonUniverse._trueE_flux_summary[key] = [0, 0.0, 0.0]

        #     # [N, sum fluxCV, sum weight before flux]
        #     CVPythonUniverse._trueE_flux_summary[key][0] += 1
        #     CVPythonUniverse._trueE_flux_summary[key][1] += flux_cv
        #     CVPythonUniverse._trueE_flux_summary[key][2] += weight

        #     if not hasattr(CVPythonUniverse, "_trueE_flux_print_counter"):
        #         CVPythonUniverse._trueE_flux_print_counter = 0

        #     if CVPythonUniverse._trueE_flux_print_counter < 10:
        #         CVPythonUniverse._trueE_flux_print_counter += 1
        #         print("[TRUEE_FLUX_BIN_SUMMARY]")
        #         for k, v in sorted(CVPythonUniverse._trueE_flux_summary.items()):
        #             pdg_bin, ebin_print = k
        #             n, sum_flux, sum_w_before = v
        #             if n < 20:
        #                 continue
        #             print(
        #                 "pdg=", pdg_bin,
        #                 "Ebin=", ebin_print, "-", ebin_print + 0.5,
        #                 "n=", n,
        #                 "avg_fluxCV=", sum_flux / n,
        #                 "sum_fluxCV=", sum_flux,
        #                 "avg_weight_before_flux=", sum_w_before / n,
        #             )

        # weight *= flux_cv
        # # ---- END DEBUG ----

        weight *= self.GetFluxAndCVWeight(self.mc_incomingE*1e-3,pdg)
        # weight *= get_flux_ratio_me_to_le(self.mc_incomingE * 1e-3, self.mc_incoming, "1N")

        # weight *= self.GetLowRecoil2p2hWeight()
        # Low-recoil 2p2h CV weight.
        #
        # Guard pathological truth kinematics. Some events have placeholder
        # values such as mc_Q2=-99999 and mc_w=99999; the 2p2h reweighter can
        # return NaN for those inputs. For those cases, use neutral weight 1.
        # w_2p2h = 1.0

        # q2_gev2 = self.mc_Q2 * 1e-6
        # w_gev = self.mc_w * 1e-3

        # valid_2p2h_kinematics = (
        #     math.isfinite(q2_gev2)
        #     and math.isfinite(w_gev)
        #     and q2_gev2 >= 0.0
        #     and abs(w_gev) < 90.0
        # )

        # if valid_2p2h_kinematics:
        #     w_2p2h = self.GetLowRecoil2p2hWeight()

        #     if not math.isfinite(w_2p2h):
        #         # Keep neutral weight for pathological reweighter output.
        #         w_2p2h = 1.0

        # weight *= w_2p2h

        # ------------------------------------------------------------
        # Low-recoil 2p2h CV weight
        # ------------------------------------------------------------
        # Compare the direct MAT weight with our safeguarded weight.
        # Apply the safeguarded version to the prediction.

        w_2p2h_safe = 1.0

        q2_gev2 = self.mc_Q2 * 1e-6
        w_gev = self.mc_w * 1e-3

        valid_2p2h_kinematics = (
            math.isfinite(q2_gev2)
            and math.isfinite(w_gev)
            and q2_gev2 >= 0.0
            and abs(w_gev) < 90.0
        )

        # Only potentially nonunity for MEC on nuclear targets.
        is_mec = self.mc_intType == 8 and self.mc_targetZ >= 2

        # Evaluate MAT directly for diagnostics.
        # For events rejected by our old guard, this lets us see
        # whether MAT actually returns a valid correction.
        w_2p2h_raw = None

        if is_mec:
            try:
                w_2p2h_raw = self.GetLowRecoil2p2hWeight()
            except Exception as exc:
                print(
                    "[2P2H MAT EXCEPTION]",
                    f"Q2={q2_gev2}",
                    f"W={w_gev}",
                    f"error={exc}",
                    flush=True
                )

        # Reproduce our original safeguarded implementation.
        if valid_2p2h_kinematics:
            if w_2p2h_raw is not None:
                if math.isfinite(w_2p2h_raw):
                    w_2p2h_safe = w_2p2h_raw
            elif not is_mec:
                # MAT gives unity for ineligible events.
                w_2p2h_safe = 1.0

        # Print cases where the raw MAT weight is nonfinite
        # or our guard changes a finite MAT weight.
        if is_mec and (
            w_2p2h_raw is None
            or not math.isfinite(w_2p2h_raw)
            or (
                math.isfinite(w_2p2h_raw)
                and abs(w_2p2h_safe - w_2p2h_raw) > 1e-10
            )
        ):
            if not hasattr(CVPythonUniverse, "_2p2h_debug_prints"):
                CVPythonUniverse._2p2h_debug_prints = 0

            if (
                self.ShortName() == "cv"
                and CVPythonUniverse._2p2h_debug_prints < 20
            ):
                CVPythonUniverse._2p2h_debug_prints += 1

                print(
                    "[2P2H SAFE VS RAW]",
                    f"intType={self.mc_intType}",
                    f"Q2={q2_gev2}",
                    f"W={w_gev}",
                    f"q0={self.Getq0True() * 1e-3}",
                    f"q3={self.Getq3True() * 1e-3}",
                    f"old_guard={valid_2p2h_kinematics}",
                    f"raw={w_2p2h_raw}",
                    f"safe={w_2p2h_safe}",
                    flush=True
                )

        weight *= w_2p2h_safe



        weight *= self.GetRPAWeight()
        weight *= self.GetMyLowQ2PiWeight() # using MENU1PI for Aaron's result
        weight *= self.GetGeantHadronWeight()
        weight *= self.GetMyMinosEfficiencyWeight()

        ### MnvTune v4.3.1 block below ###
        #weight *= self.GetTargetMassWeight() # 
        # weight *= self.GetCOHPionWeight() # Alejandro Ramirez's Result
        if V1_ADD_COH and V1_ADD_LE_COH:
            raise RuntimeError(
                "Choose either ME coherent weighting or LE coherent weighting, not both."
            )

        if V1_ADD_COH:
            # ME theta_pi x E_pi coherent correction
            weight *= self.GetCOHPionWeight()

        if V1_ADD_LE_COH:
            # LE prescription
            weight *= self.GetLECOHPionWeight()



        if self.mc_intType == 4:
            weight *= 1.4368 # diffractive coherent reweight from Aaron Mislivec

        return weight

    def GetStandardWeightV4p3p0(self):
        """
        MINERvA Tune v4.3.0 central-value event weight.

        Generator corrections:
        - NonResPi reduction
        - Full deuterium bubble-chamber fit
        - Valencia QE RPA
        - LE low-recoil Valencia 2p2h fit
        - MENU1PI low-Q2 pion correction
        - ME coherent pion correction
        - Diffractive correction (hydrogen only)

        No QE z-expansion or elastic-FSI bug-fix weight.

        Detector corrections are controlled separately.
        """

        weight = 1.0

        # --------------------------------------------------
        # GENIE CV: NonResPi + full deuterium fit
        # --------------------------------------------------
        # Configured through SetNonResPiReweight(True),
        # SetDeuteriumGeniePiTune(True), and
        # SetZExpansionFaReweight(False).
        weight *= self.GetGenieWeight()

        # --------------------------------------------------
        # Flux CV
        # --------------------------------------------------
        pdg = self.mc_incoming

        if SystematicsConfig.USE_SWAPPED:
            if pdg == 12:
                pdg = 14
            elif pdg == -12:
                pdg = -14
            elif pdg == 14:
                pdg = 12
            elif pdg == -14:
                pdg = -12
            else:
                raise ValueError(
                    "Unexpected incoming neutrino PDG: {}".format(pdg)
                )

        weight *= self.GetFluxAndCVWeight(
            self.mc_incomingE * 1e-3,
            pdg
        )

        # --------------------------------------------------
        # LE low-recoil Valencia 2p2h
        # --------------------------------------------------
        weight *= self.GetLowRecoil2p2hWeight()

        # --------------------------------------------------
        # Valencia QE RPA
        # --------------------------------------------------
        weight *= self.GetRPAWeight()

        # --------------------------------------------------
        # Low-Q2 pion correction: Aaron / MENU1PI
        # --------------------------------------------------
        if V4P3P0_ADD_LOWQ2PI:
            weight *= self.GetMyLowQ2PiWeight(
                V4P3P0_LOWQ2PI_CHANNEL
            )

        # --------------------------------------------------
        # ME coherent-pion reweight
        # --------------------------------------------------
        if V4P3P0_ADD_ME_COH:
            weight *= self.GetCOHPionWeight()

        # --------------------------------------------------
        # Diffractive coherent correction
        # Hydrogen only
        # --------------------------------------------------
        if (
            V4P3P0_ADD_DIFFRACTIVE
            and self.mc_intType == 4
            and self.mc_targetZ == 1
        ):
            weight *= 1.4368

        # --------------------------------------------------
        # Detector / reconstruction corrections
        # --------------------------------------------------
        if V4P3P0_ADD_GEANT_HADRON:
            weight *= self.GetGeantHadronWeight()

        if V4P3P0_ADD_MINOS_EFFICIENCY:
            weight *= self.GetMyMinosEfficiencyWeight()

        return weight

    def GetStandardWeightV4p4p1(self):
        """
        MINERvA Tune v4.4.1 central-value weight.

        Core corrections:
        - NonResPi reduction
        - Full deuterium bubble-chamber tune
        - QE axial form factor z-expansion
        - GENIE elastic-FSI bug fix
        - Valencia RPA
        - LE low-recoil Valencia 2p2h fit

        Optional corrections are controlled by the
        V4P4P1_ADD_* configuration flags.
        """

        weight = 1.0

        # --------------------------------------------------
        # GENIE CV
        # NonResPi + deuterium + z-expansion
        # --------------------------------------------------
        # Switches configured in GetAllSystematicsUniverses().
        weight *= self.GetGenieWeight()

        # --------------------------------------------------
        # GENIE elastic-FSI bug fix
        # iWeight=1: elastic correction only
        # --------------------------------------------------
        if V4P4P1_ADD_ELASTIC_FSI:
            weight *= self.GetFSIWeight(1)

        # --------------------------------------------------
        # Flux CV
        # --------------------------------------------------
        pdg = self.mc_incoming

        if SystematicsConfig.USE_SWAPPED:
            if pdg == 12:
                pdg = 14
            elif pdg == -12:
                pdg = -14
            elif pdg == 14:
                pdg = 12
            elif pdg == -14:
                pdg = -12
            else:
                raise ValueError(
                    "Unexpected incoming neutrino PDG: {}".format(pdg)
                )

        weight *= self.GetFluxAndCVWeight(
            self.mc_incomingE * 1e-3, pdg
        )

        # --------------------------------------------------
        # LE low-recoil Valencia 2p2h
        # Original MAT implementation, without safeguards
        # --------------------------------------------------
        weight *= self.GetLowRecoil2p2hWeight()

        # --------------------------------------------------
        # Valencia RPA
        # --------------------------------------------------
        weight *= self.GetRPAWeight()

        # --------------------------------------------------
        # Optional low-Q2 pion correction
        # --------------------------------------------------
        if V4P4P1_ADD_LOWQ2PI:
            weight *= self.GetMyLowQ2PiWeight(
                V4P4P1_LOWQ2PI_CHANNEL
            )

        # --------------------------------------------------
        # Optional coherent pion corrections
        # --------------------------------------------------
        if V4P4P1_ADD_LE_COH and V4P4P1_ADD_ME_COH:
            raise RuntimeError(
                "Choose either LE or ME coherent weighting, not both."
            )

        if V4P4P1_ADD_LE_COH:
            weight *= self.GetLECOHPionWeight()

        if V4P4P1_ADD_ME_COH:
            weight *= self.GetCOHPionWeight()

        # --------------------------------------------------
        # Optional diffractive correction
        # Hydrogen only
        # --------------------------------------------------
        if (
            V4P4P1_ADD_DIFFRACTIVE
            and self.mc_intType == 4
            and self.mc_targetZ == 1
        ):
            weight *= 1.4368

        # --------------------------------------------------
        # Detector / reconstruction corrections
        # --------------------------------------------------
        # weight *= self.GetGeantHadronWeight()
        weight *= self.GetMyMinosEfficiencyWeight()

        return weight



    def GetStandardWeightV1(self):

        weight = 1.0

        # ------------------------------------------------------------
        # GENIE / V1 pion tune
        # ------------------------------------------------------------

        if V1_ADD_DEUTERIUM:
            # With:
            #   NonResPi = ON
            #   Deuterium = ON
            #   z-exp = OFF
            #
            # GetGenieWeight() is exactly:
            #   NonResPi * Deuterium
            weight *= self.GetGenieWeight()

        else:
            # Original validated V1 baseline:
            # Non-resonant pion suppression only.
            rvn1pi = self.GetVecElem("truth_genie_wgt_Rvn1pi", 2)
            rvp1pi = self.GetVecElem("truth_genie_wgt_Rvp1pi", 2)

            is_nonrespi = (rvn1pi < 1.0) or (rvp1pi < 1.0)
            weight *= 0.43 if is_nonrespi else 1.0

        # ------------------------------------------------------------
        # Flux
        # ------------------------------------------------------------
        pdg = self.mc_incoming

        if SystematicsConfig.USE_SWAPPED:
            if pdg == 12:
                pdg = 14
            elif pdg == -12:
                pdg = -14
            elif pdg == 14:
                pdg = 12
            elif pdg == -14:
                pdg = -12
            else:
                print("Not neutrino PDG? {}".format(pdg))

        weight *= self.GetFluxAndCVWeight(self.mc_incomingE * 1e-3, pdg)
        # # weight *= get_flux_ratio_me_to_le(self.mc_incomingE * 1e-3, self.mc_incoming, "1M")
        # Flux CV: reproduce Maria's treatment
        # pdg = self.mc_incoming

        # if abs(pdg) == 12:
        #     flux_weight = 1.0
        # else:
        #     flux_weight = self.GetFluxAndCVWeight(
        #         self.mc_incomingE * 1e-3,
        #         pdg
        #     )

        # weight *= flux_weight


        # ------------------------------------------------------------
        # Low-recoil 2p2h
        # ------------------------------------------------------------
        # w_2p2h = 1.0

        # q2_gev2 = self.mc_Q2 * 1e-6
        # w_gev = self.mc_w * 1e-3

        # valid_2p2h_kinematics = (
        #     math.isfinite(q2_gev2)
        #     and math.isfinite(w_gev)
        #     and q2_gev2 >= 0.0
        #     and abs(w_gev) < 90.0
        # )

        # if valid_2p2h_kinematics:
        #     w_2p2h = self.GetLowRecoil2p2hWeight()

        #     if not math.isfinite(w_2p2h):
        #         w_2p2h = 1.0

        # weight *= w_2p2h

        weight *= self.GetLowRecoil2p2hWeight()

        # ------------------------------------------------------------
        # Low-recoil 2p2h CV weight
        # ------------------------------------------------------------
        # Compare the direct MAT weight with our safeguarded weight.
        # Apply the safeguarded version to the prediction.

        # w_2p2h_safe = 1.0

        # q2_gev2 = self.mc_Q2 * 1e-6
        # w_gev = self.mc_w * 1e-3

        # valid_2p2h_kinematics = (
        #     math.isfinite(q2_gev2)
        #     and math.isfinite(w_gev)
        #     and q2_gev2 >= 0.0
        #     and abs(w_gev) < 90.0
        # )

        # # Only potentially nonunity for MEC on nuclear targets.
        # is_mec = self.mc_intType == 8 and self.mc_targetZ >= 2

        # # Evaluate MAT directly for diagnostics.
        # # For events rejected by our old guard, this lets us see
        # # whether MAT actually returns a valid correction.
        # w_2p2h_raw = None

        # if is_mec:
        #     try:
        #         w_2p2h_raw = self.GetLowRecoil2p2hWeight()
        #     except Exception as exc:
        #         print(
        #             "[2P2H MAT EXCEPTION]",
        #             f"Q2={q2_gev2}",
        #             f"W={w_gev}",
        #             f"error={exc}",
        #             flush=True
        #         )

        # # Reproduce our original safeguarded implementation.
        # if valid_2p2h_kinematics:
        #     if w_2p2h_raw is not None:
        #         if math.isfinite(w_2p2h_raw):
        #             w_2p2h_safe = w_2p2h_raw
        #     elif not is_mec:
        #         # MAT gives unity for ineligible events.
        #         w_2p2h_safe = 1.0

        # # Print cases where the raw MAT weight is nonfinite
        # # or our guard changes a finite MAT weight.
        # if is_mec and (
        #     w_2p2h_raw is None
        #     or not math.isfinite(w_2p2h_raw)
        #     or (
        #         math.isfinite(w_2p2h_raw)
        #         and abs(w_2p2h_safe - w_2p2h_raw) > 1e-10
        #     )
        # ):
        #     if not hasattr(CVPythonUniverse, "_2p2h_debug_prints"):
        #         CVPythonUniverse._2p2h_debug_prints = 0

        #     if (
        #         self.ShortName() == "cv"
        #         and CVPythonUniverse._2p2h_debug_prints < 20
        #     ):
        #         CVPythonUniverse._2p2h_debug_prints += 1

        #         print(
        #             "[2P2H SAFE VS RAW]",
        #             f"intType={self.mc_intType}",
        #             f"Q2={q2_gev2}",
        #             f"W={w_gev}",
        #             f"q0={self.Getq0True() * 1e-3}",
        #             f"q3={self.Getq3True() * 1e-3}",
        #             f"old_guard={valid_2p2h_kinematics}",
        #             f"raw={w_2p2h_raw}",
        #             f"safe={w_2p2h_safe}",
        #             flush=True
        #         )

        # weight *= w_2p2h_safe




        # ------------------------------------------------------------
        # RPA
        # ------------------------------------------------------------
        weight *= self.GetRPAWeight()

        # ------------------------------------------------------------
        # Optional later tune components
        # ------------------------------------------------------------
        if V1_ADD_LOWQ2PI:
            weight *= self.GetMyLowQ2PiWeight(
                V1_LOWQ2PI_CHANNEL
            )

        if V1_ADD_COH and V1_ADD_LE_COH:
            raise RuntimeError(
                "Choose either ME coherent weighting or LE coherent weighting, not both."
            )

        if V1_ADD_COH:
            # ME theta_pi x E_pi coherent correction
            weight *= self.GetCOHPionWeight()

        if V1_ADD_LE_COH:
            # LE prescription
            weight *= self.GetLECOHPionWeight()

        if V1_ADD_DIFFCOH and self.mc_intType == 4:
            weight *= 1.4368

        # ------------------------------------------------------------
        # Detector / reconstruction corrections
        # ------------------------------------------------------------
        weight *= self.GetGeantHadronWeight()
        weight *= self.GetMyMinosEfficiencyWeight()

        return weight


    def GetMyMinosEfficiencyWeight(self):
        if self._is_truth_mode():
            return 1.0

        w_raw = self.GetMinosEfficiencyWeight()

        if self.HasNoBackExitingTracks:
            return 1.0
        else:
            return w_raw


    def GetMyLowQ2PiWeight(self, channel=None):
        if channel is None:
            channel = SystematicsConfig.LowQ2PiWeightChannel

        if channel is None:
            return 1.0

        return self.GetLowQ2PiWeight(channel.upper())

    def GetOldBroadMENU1PIWeight(self):
        if self.mc_current != 1 or self.mc_intType != 2:
            return 1.0

        return ROOT.PlotUtils.weight_lowq2pi().getWeight(
            self.mc_Q2 * 1e-6,
            "MENU1PI",
            0
        )


    def GetModelWeight(self):
        w = GetModelWeight(self, 'Eel',"Pi0") # (Ee, Theta, ETh), (PCElectron, PCPhoton, PCPi0) 
        return w

    def GetLECOHPionWeight(self):
        """
        LE coherent-pion prescription from the MINERvA tune documentation:

        coherent event with true E_pi < 0.5 GeV -> weight 0.5
        otherwise                              -> weight 1.0

        E_pi is the true final-state pion total energy.
        """

        # Only coherent interactions
        if self.mc_intType != 4:
            return 1.0

        # Find the charged final-state pion
        for i in range(self.mc_nFSPart):
            try:
                pdg = self.mc_FSPartPDG[i]
            except Exception:
                continue

            if abs(pdg) != 211:
                continue

            true_epi_gev = self.mc_FSPartE[i] / 1000.0

            if true_epi_gev < 0.5:
                return 0.5

            return 1.0

        # No charged pion found
        return 1.0

    def GetCOHPionWeight(self):
        def GetTrueHighEpi():
            nFSpi = self.mc_nFSPart
            pionE = -1.0
            pionKE = -1.0
            for i in range(nFSpi):
                try:
                    pdg = self.mc_FSPartPDG[i]
                except:
                    continue

                if pdg != 211:
                    continue
                energy = self.mc_FSPartE[i]
                mass = 139.569
                tpi = energy - mass
                if tpi >= pionKE:
                    pionKE = tpi
                    pionE = energy
            return pionE

        def GetTrueAngleHighTpi():
            nFSpi = self.mc_nFSPart
            angle = -9999 #WRTbeam and in degrees
            pionKE = 0.0
            idk = -9999
            for i in range(nFSpi):
                try:
                    pdg = self.mc_FSPartPDG[i]
                except:
                    continue

                if pdg != 211:
                    continue
                energy = self.mc_FSPartE[i]
                mass = 139.569
                tpi = energy - mass
                if tpi >= pionKE:
                    pionKE = tpi
                    pimomentumvec = ROOT.TVector3(self.mc_FSPartPx[i], self.mc_FSPartPy[i],self.mc_FSPartPz[i])
                    r = ROOT.Math.RotationX(SystematicsConfig.BEAM_ANGLE)
                    s=r(pimomentumvec)
                    deg_wrtb = s.Theta()
                    angle = deg_wrtb

            #Making sure angle is only between 0 and pi
            if angle < 0.0:
                angle = -1.0*angle
            if angle > math.pi:
                angle = 2.0*math.pi - angle
            return angle*180./math.pi # Degrees
       
        weight = 1.0
        if self.mc_intType != 4:
           return 1.0
        elif self.mc_intType == 4:
           npi = self.mc_FSPartPDG.count(211)
           if npi == 0:
               return 1.0
           angle = GetTrueAngleHighTpi()
           epi = GetTrueHighEpi()/1000. # This is supposed to be the energy of the pion!!!! 
           if epi < 0:
               return 1.0
           else:
               weight *= self.GetCoherentPiWeight(angle, epi) #Inputs are in Degrees and GeVi
               return weight
        else:
            return 1.0



    def IsCCNucleonPion(self):
        npion = 0

        for i in range(self.mc_nFSPart):
            pdg = self.mc_FSPartPDG[i]

            if self.mc_incoming == 14 and pdg == 211:
                npion += 1
            elif self.mc_incoming == -14 and pdg == -211:
                npion += 1

        return (
            self.mc_w < 1400.0
            and npion == 1
            and self.mc_current == 1
            and self.mc_intType in (2, 3)
        )

    def GetDiagnosticLowQ2PiWeight(self, channel):
        channel = channel.upper()

        if channel == "MENU1PI":
            if not self.IsCCNucleonPion():
                return 1.0
        else:
            if not (self.mc_current == 1 and self.mc_intType == 2):
                return 1.0

        return ROOT.PlotUtils.weight_lowq2pi().getWeight(
            self.mc_Q2 * 1e-6,
            channel,
            0,
            self.mc_targetNucleus
        )

    #################### functions to be overide in systematics shift universes #################
    #### No correction
    def ElectronEnergyRaw(self):
        return self.GetVecElem("prong_part_E",0,3)

    @Utilities.decorator_ReLU
    def ElectronEnergy(self):
        return self.ElectronEnergyRaw() + self.GetEMEnergyShift() #+ self.GetLeakageCorrection()

    ## ORIGINAL electron momentum calculation in DET coord, which assumes the momentum to be in DET coord
    def ElectronP3D_det(self):
        electronp = self.GetVecOfVecDouble("prong_part_E")
        scale = self.GetEMEnergyShift()/electronp[0][3] if (electronp[0][3]>0) else 0
        p = ROOT.Math.XYZVector(*tuple(list(electronp[0])[:3]))*(1+scale)
        return p

    # ## ORIGINAL electron momentum calculation, which assumes the momentum to be in DET coord and needs to be rotated
    # ## to be in beam coord
    # def ElectronP3D(self):              # momentum in DET coord already, need to rotate
    #     electronp = self.GetVecOfVecDouble("prong_part_E")
    #     scale = self.GetEMEnergyShift()/electronp[0][3] if (electronp[0][3]>0) else 0
    #     p = ROOT.Math.XYZVector(*tuple(list(electronp[0])[:3]))*(1+scale)
    #     r = ROOT.Math.RotationX(SystematicsConfig.BEAM_ANGLE)
    #     s=r(p)
    #     if s is None:
    #         print (p,r,tuple(list(electronp[0])[:3]))
    #     return s
    # ## ORIGINAL electron momentum calculation, which assumes the momentum to be in DET coord and needs to be rotated


    #### Yes correction: angle only, no momentum-magnitude correction
    def ElectronP3D_det(self):
        # ---- 0) build uncorrected p(det) from branches + existing EM energy shift ----
        electronp = self.GetVecOfVecDouble("prong_part_E")
        scale = self.GetEMEnergyShift()/electronp[0][3] if (electronp[0][3] > 0) else 0.0
        p_det = ROOT.Math.XYZVector(*tuple(list(electronp[0])[:3])) * (1.0 + scale)

        # ---- 1) get vertex (mm) ----
        vx = float(self.vtx[0])
        vy = float(self.vtx[1])

        make_vec = ROOT.Math.XYZVector

        # ---- 2) apply direction correction only; keep |p| unchanged ----
        p_det_corr, thetaX_corr, thetaY_corr, theta_corr = apply_theta_bias_correction_det(
            p_det,
            vx, vy,
            # This is nu+e fit parameters
            # We use nu+e angle correction coefficients for both nu+e and CCnue
            # # OLD p6 coefficients
            # t0x=0.0001289, t1x=1.402e-05,
            # t0y=0.001007,  t1y=1.87e-05,
            # # p8 coefficients before hadron reweights
            # t0x=0.000486, t1x=1.392e-05,
            # t0y=0.0008169,  t1y=1.969e-05,
            # # NEW p8 coefficients AFTER hadron reweight (ME FV)
            # t0x=0.0005425, t1x=1.389e-05,
            # t0y=0.0007356,  t1y=1.953e-05,
            # NEW p8 coefficients AFTER hadron reweight (inclusive FV)
            t0x=0.0005755, t1x=1.429e-05,
            t0y=0.0007915,  t1y=1.951e-05,
            clamp_vtx=True,
            vtx_min_mm=-1200.0,
            vtx_max_mm=1200.0,
            keep_pmag=True,
            make_vec=make_vec,
        )

        return p_det_corr


    def ElectronP3D_det_noEMshift(self):
        # Build unshifted p(det) from branches, with NO EM energy shift scaling
        electronp = self.GetVecOfVecDouble("prong_part_E")
        p_det = ROOT.Math.XYZVector(*tuple(list(electronp[0])[:3]))

        vx = float(self.vtx[0])
        vy = float(self.vtx[1])

        make_vec = ROOT.Math.XYZVector

        # Direction correction only; keep |p| unchanged
        p_det_corr, _, _, _ = apply_theta_bias_correction_det(
            p_det,
            vx, vy,
            # This is nu+e fit parameters
            # We use nu+e angle correction coefficients for both nu+e and CCnue
            # # OLD p6 coefficients
            # t0x=0.0001289, t1x=1.402e-05,
            # t0y=0.001007,  t1y=1.87e-05,
            # # p8 coefficients before hadron reweights
            # t0x=0.000486, t1x=1.392e-05,
            # t0y=0.0008169,  t1y=1.969e-05,
            # # NEW p8 coefficients AFTER hadron reweight (ME FV)
            # t0x=0.0005425, t1x=1.389e-05,
            # t0y=0.0007356,  t1y=1.953e-05,
            # NEW p8 coefficients AFTER hadron reweight (inclusive FV)
            t0x=0.0005755, t1x=1.429e-05,
            t0y=0.0007915,  t1y=1.951e-05,
            clamp_vtx=True,
            vtx_min_mm=-1200.0,
            vtx_max_mm=1200.0,
            keep_pmag=True,
            make_vec=make_vec,
        )

        return p_det_corr


    ## Electron momentum in BEAM coord:
    ## - get angle-corrected momentum in DET coord
    ## - rotate into beam coord
    def ElectronP3D(self):
        p_det_corr = self.ElectronP3D_det()

        r = ROOT.Math.RotationX(SystematicsConfig.BEAM_ANGLE)
        p_beam = r(p_det_corr)

        if p_beam is None:
            print(p_det_corr, r)

        return p_beam


    # Keep electron energy branch-based, since we are NOT applying
    # a momentum-magnitude correction anymore.
    def ElectronEnergyRaw(self):
        return self.GetVecElem("prong_part_E", 0, 3)


    @Utilities.decorator_ReLU
    def ElectronEnergy(self):
        return self.ElectronEnergyRaw() + self.GetEMEnergyShift()
    #### Yes correction: angle only, no momentum-magnitude correction





    # #### Yes correction both angle and magnitude
    # def ElectronP3D_det(self):
    #     # ---- 0) build uncorrected p(det) from branches + existing EM energy shift ----
    #     electronp = self.GetVecOfVecDouble("prong_part_E")
    #     scale = self.GetEMEnergyShift()/electronp[0][3] if (electronp[0][3] > 0) else 0.0
    #     p_det = ROOT.Math.XYZVector(*tuple(list(electronp[0])[:3])) * (1.0 + scale)

    #     # ---- 1) get vertex (mm) ----
    #     # If you already have self.vtx like event.vtx, use it:
    #     vx = float(self.vtx[0])
    #     vy = float(self.vtx[1])

    #     # If you *don’t* have self.vtx, replace the two lines above with however you access vtx:
    #     # vtx = self.GetVecDouble("vtx")  # example
    #     # vx, vy = float(vtx[0]), float(vtx[1])

    #     make_vec = ROOT.Math.XYZVector

    #     # ---- 2) apply direction correction (keep |p|) ----
    #     p_det_corr, thetaX_corr, thetaY_corr, theta_corr = apply_theta_bias_correction_det(
    #         p_det,
    #         vx, vy,
    #         # This is nu+e fit parameters
    #         # We use nu+e angle correction coefficients for both nu+e and CCnue
    #         t0x=0.0001289, t1x=1.402e-05,
    #         t0y=0.001007,  t1y=1.87e-05,
    #         clamp_vtx=True,
    #         vtx_min_mm=-1200.0,
    #         vtx_max_mm=1200.0,
    #         keep_pmag=True,
    #         make_vec=make_vec,
    #     )
        
    #     # print("SystematicsConfig.ELASTIC_MODE: ", SystematicsConfig.ELASTIC_MODE)
    #     # ---- 3) choose magnitude correction by analysis mode ----
    #     if SystematicsConfig.ELASTIC_MODE:
    #         # nu+e coefficients
    #         a0 = -0.01189
    #         a1 = 1.06e-04
    #     else:
    #         # CCnue coefficients
    #         a0 = 0.01914
    #         a1 = 6.03e-05


    #     # a0 = -0.01189
    #     # a1 = 1.06e-04
    #     # ---- 4) apply magnitude correction ----
    #     p_det_corr = apply_pmag_frac_correction_det(
    #         p_det_corr,
    #         vy,
    #         a0=a0,
    #         a1=a1,
    #         make_vec=make_vec
    #     )

    #     return p_det_corr

    # def ElectronP3D_det_noEMshift(self):
    #     # Build unshifted p(det) from branches (NO EM energy shift scaling here)
    #     electronp = self.GetVecOfVecDouble("prong_part_E")
    #     p_det = ROOT.Math.XYZVector(*tuple(list(electronp[0])[:3]))

    #     vx = float(self.vtx[0])
    #     vy = float(self.vtx[1])

    #     make_vec = ROOT.Math.XYZVector

    #     # direction correction (keep |p|)
    #     p_det_corr, _, _, _ = apply_theta_bias_correction_det(
    #         p_det,
    #         vx, vy,
    #         # This is nu+e fit parameters
    #         # We use nu+e angle correction coefficients for both nu+e and CCnue
    #         t0x=0.0001289, t1x=1.402e-05,
    #         t0y=0.001007,  t1y=1.87e-05,
    #         clamp_vtx=True,
    #         vtx_min_mm=-1200.0,
    #         vtx_max_mm=1200.0,
    #         keep_pmag=True,
    #         make_vec=make_vec,
    #     )

    #     # choose magnitude correction by analysis mode ----
    #     if SystematicsConfig.ELASTIC_MODE:
    #         # nu+e coefficients
    #         a0 = -0.01189
    #         a1 = 1.06e-04
    #     else:
    #         # CCnue coefficients
    #         a0 = 0.01914
    #         a1 = 6.03e-05

    #     # magnitude correction
    #     p_det_corr = apply_pmag_frac_correction_det(
    #         p_det_corr,
    #         vy,
    #         a0=a0,
    #         a1=a1,
    #         make_vec=make_vec
    #     )

    #     return p_det_corr

    # ## Electron momentum in BEAM coord:
    # ## - get corrected momentum in DET coord
    # ## - rotate into beam coord
    # def ElectronP3D(self):
    #     # corrected det-coordinate momentum (XYZqVector)
    #     p_det_corr = self.ElectronP3D_det()

    #     r = ROOT.Math.RotationX(SystematicsConfig.BEAM_ANGLE)
    #     p_beam = r(p_det_corr)

    #     if p_beam is None:
    #         print(p_det_corr, r)

    #     return p_beam

    # # Recalculate electron energy using corrected momentum
    # def ElectronEnergy(self):
    #     # Use corrected det momentum (or beam — magnitude is the same either way)
    #     p = self.ElectronP3D_det().R()   # |p| from corrected momentum
    #     return math.sqrt(p*p + M_e_sqr)

    # def ElectronEnergyRaw(self):
    #     p = self.ElectronP3D_det_noEMshift().R()   # |p| from corrected momentum
    #     return math.sqrt(p*p + M_e_sqr)
    # #### Yes correction both angle and magnitude







    def ElectronTheta(self):
        return self.ElectronP3D().Theta()



    def ElectronProtonAngle(self):
        electronp = self.GetVecOfVecDouble("prong_part_E")
        scale = self.GetEMEnergyShift()/electronp[0][3] if (electronp[0][3]>0) else 0
        electron_p = ROOT.TVector3(*tuple(list(electronp[0])[:3]))*(1+scale)
        
        protonX = self.MasterAnaDev_proton_Px_fromdEdx
        protonY = self.MasterAnaDev_proton_Py_fromdEdx
        protonZ = self.MasterAnaDev_proton_Pz_fromdEdx
        proton_p = ROOT.TVector3(*tuple([protonX,protonY,protonZ]))
        dotprod = proton_p.Dot(electron_p)
        mag = proton_p.Mag() * electron_p.Mag()
        if mag > 0:
            dotprod /= mag
        if protonX > -9998 and protonY > -9998 and protonZ > -9998:
            return math.degrees(math.acos(dotprod))
        else:
            return -5


    # --- NEW: vertex-based incident neutrino direction unit vector ---
    def NuDirUnitFromVertex(self, Lmm=900000.0):
        """
        Estimate incident neutrino direction at this event using the reconstructed vertex.

        Returns a unit vector in the SAME coordinate system as ElectronP3D()
        (i.e. after RotationX(BEAM_ANGLE) is applied).

        Lmm is the effective upstream distance to the source point (mm),
        interpreted along the NOMINAL BEAM axis (beam-coordinate -z).
        """
        # Reco vertex in detector coordinates (mm)
        vx = self.GetVecElem("vtx", 0)
        vy = self.GetVecElem("vtx", 1)
        vz = self.GetVecElem("vtx", 2)

        v_det = ROOT.Math.XYZVector(vx, vy, vz)

        # Rotate vertex into the "beam" coordinate convention, matching ElectronP3D()
        r = ROOT.Math.RotationX(SystematicsConfig.BEAM_ANGLE)
        v_beam = r(v_det)

        # Effective source point upstream along -z in BEAM coordinates
        # source_beam = (0, 0, -Lmm)
        # dir_beam    = v_beam - source_beam = (x, y, z + Lmm)
        nu = ROOT.Math.XYZVector(v_beam.X(), v_beam.Y(), v_beam.Z() + Lmm)

        mag = nu.R()
        if mag <= 0:
            return ROOT.Math.XYZVector(0, 0, 1)  # safe fallback (beam +z)
        return nu * (1.0 / mag)


    # --- NEW: build ν-aligned transverse basis (beam coords) ---
    def _NuAlignedBasis(self, Lmm=900000.0):
        """
        Construct an orthonormal basis (ex, ey, ez) in BEAM coordinates such that:
        ez = inferred neutrino direction at this vertex
        ex, ey span the plane transverse to ez

        This is used to define vertex-corrected thetaX/thetaY consistently.
        """
        ez = self.NuDirUnitFromVertex(Lmm=Lmm)  # unit vector in beam coords

        # Pick a reference axis not (nearly) parallel to ez
        yref = ROOT.Math.XYZVector(0, 1, 0)
        if abs(ez.Dot(yref)) > 0.95:
            yref = ROOT.Math.XYZVector(1, 0, 0)

        ex = yref.Cross(ez)
        ex_mag = ex.R()
        if ex_mag <= 0:
            # extreme fallback
            ex = ROOT.Math.XYZVector(1, 0, 0)
            ex_mag = ex.R()
        ex = ex * (1.0 / ex_mag)

        ey = ez.Cross(ex)
        ey_mag = ey.R()
        if ey_mag > 0:
            ey = ey * (1.0 / ey_mag)

        return ex, ey, ez


    # --- NEW: theta using dot product with vertex-based ν direction ---
    def ElectronTheta_vtxcorr(self, prong=0, Lmm=900000.0):
        """
        Vertex-corrected electron theta: opening angle between reconstructed lepton direction
        and estimated incident ν direction at this vertex.
        """
        # Reco lepton momentum (already rotated into beam coords inside ElectronP3D)
        p = self.ElectronP3D()
        pmag = p.R()
        if pmag <= 0:
            return 0.0

        ulep = p * (1.0 / pmag)
        unu  = self.NuDirUnitFromVertex(Lmm=Lmm)

        cosang = _clamp(ulep.Dot(unu))
        return math.acos(cosang)


    def ElectronThetaX_vtxcorr(self, prong=0, Lmm=900000.0):
        """
        Vertex-corrected thetaX: lepton angular component in ν-aligned frame.
        Defined as atan2(u·ex, u·ez).
        """
        p = self.ElectronP3D()
        pmag = p.R()
        if pmag <= 0:
            return 0.0
        ulep = p * (1.0 / pmag)

        ex, ey, ez = self._NuAlignedBasis(Lmm=Lmm)
        ux = ulep.Dot(ex)
        uz = ulep.Dot(ez)

        return math.atan2(ux, uz)


    def ElectronThetaY_vtxcorr(self, prong=0, Lmm=900000.0):
        """
        Vertex-corrected thetaY: lepton angular component in ν-aligned frame.
        Defined as atan2(u·ey, u·ez).
        """
        p = self.ElectronP3D()
        pmag = p.R()
        if pmag <= 0:
            return 0.0
        ulep = p * (1.0 / pmag)

        ex, ey, ez = self._NuAlignedBasis(Lmm=Lmm)
        uy = ulep.Dot(ey)
        uz = ulep.Dot(ez)

        return math.atan2(uy, uz)


    def ElectronTheta2D_vtxcorr(self, prong=0, Lmm=900000.0):
        """
        Vertex-corrected 2D electron angle, consistent with the ν-aligned thetaX/thetaY:
            theta2D = sqrt(thetaX^2 + thetaY^2)
        """
        tx = self.ElectronThetaX_vtxcorr(prong=prong, Lmm=Lmm)
        ty = self.ElectronThetaY_vtxcorr(prong=prong, Lmm=Lmm)
        return math.sqrt(tx*tx + ty*ty)


    def ElectronPhi_vtxcorr(self, prong=0, Lmm=900000.0):
        """
        (Optional) Vertex-corrected azimuth phi of the lepton around the inferred ν axis.
        This is the angle in the ν-transverse plane: atan2(u·ey, u·ex).
        """
        p = self.ElectronP3D()
        pmag = p.R()
        if pmag <= 0:
            return 0.0
        ulep = p * (1.0 / pmag)

        ex, ey, ez = self._NuAlignedBasis(Lmm=Lmm)
        return math.atan2(ulep.Dot(ey), ulep.Dot(ex))

    def PrimaryProtonTheta(self):
        theta = self.MasterAnaDev_proton_theta
        return math.degrees(theta) if theta>=0 else -1

    def ConeInsideE(self): 
        e = self.GetVecElem("prong_TrimmedCaloEnergy",0)/1e3
        psi = self.GetVecElem("Psi",0)
        #print e, 'inside', psi, 'psi'
        return e*psi

    def ConeOutsideE(self):
        localVec = list(self.GetVecDouble("ExtraEnergyVis"))
        e = sum(localVec)
	#print e, 'out'
        return e 


    def NeighborhoodE(self, shell_thickness_mm=50.0, _debug=False):
        def dbg(*args):
            if _debug:
                print("[NeighborhoodE DEBUG]", *args)

        R0_mm = 40.0
        tan_alpha = math.tan(math.radians(10.0))
        ip = 0

        # clusters
        X  = list(self.GetVecDouble("ExtraEnergyClusters_X"))
        Y  = list(self.GetVecDouble("ExtraEnergyClusters_Y"))
        Z  = list(self.GetVecDouble("ExtraEnergyClusters_Z"))
        Ecl = list(self.GetVecDouble("ExtraEnergyClusters_energy"))
        if len(Ecl) == 0:
            dbg("RETURN 0: no clusters")
            return 0.0

        # cone start (mm)
        vtx_all = self.GetVecOfVecDouble("prong_axis_vertex")
        x0, y0, z0 = float(vtx_all[ip][0]), float(vtx_all[ip][1]), float(vtx_all[ip][2])

        # axis: use ElectronP3D (already beam-rotated)
        p3 = self.ElectronP3D()
        px, py, pz = float(p3.X()), float(p3.Y()), float(p3.Z())
        pnorm = math.sqrt(px*px + py*py + pz*pz)
        if pnorm == 0.0:
            dbg("RETURN 0: |p|==0")
            return 0.0
        ux, uy, uz = px/pnorm, py/pnorm, pz/pnorm

        neigh = 0.0
        pass_s = pass_delta = pass_shell = 0

        for j in range(len(Ecl)):
            dx = X[j] - x0
            dy = Y[j] - y0
            dz = Z[j] - z0

            s = dx*ux + dy*uy + dz*uz
            if s <= 0.0:
                continue
            pass_s += 1

            rperp2 = dx*dx + dy*dy + dz*dz - s*s
            if rperp2 < 0.0:
                rperp2 = 0.0
            rperp = math.sqrt(rperp2)

            rcone = R0_mm + s * tan_alpha
            delta = rperp - rcone

            if delta <= 0.0:
                continue
            pass_delta += 1

            if delta < shell_thickness_mm:
                pass_shell += 1
                neigh += Ecl[j]

        dbg("cutflow: s>0", pass_s, "delta>0", pass_delta,
            f"delta<{shell_thickness_mm}", pass_shell, "neigh", neigh)

        return neigh

    def DistToHexEdge(self, apothem_mm=881.25):
        x = self.GetVecElem("vtx", 0)
        y = self.GetVecElem("vtx", 1)
        return _hex_edge_distance_mm(x, y, apothem_mm)

    def MuonP3D_beam(self):
        """
        Reco muon 3-momentum (px,py,pz) in BEAM coordinates, MeV.

        This matches the components used in MuonFunctions.h:GetPmu_nominal()
        via: GetAnaToolName() + "_leptonE".
        """
        lepton_branch = self.GetAnaToolName() + "_leptonE"
        px = self.GetVecElem(lepton_branch, 0)
        py = self.GetVecElem(lepton_branch, 1)
        pz = self.GetVecElem(lepton_branch, 2)

        # print("Muon beam p:", px, py, pz)
        # print("lepton_branch: ", lepton_branch)

        return ROOT.Math.XYZVector(px, py, pz)

    def MuonP3D_det(self):
        """
        Muon 3-momentum rotated into the same DET coordinate convention
        used by ElectronP3D() (RotationX(-BEAM_ANGLE)).
        """
        p = self.MuonP3D_beam()
        r = ROOT.Math.RotationX(-SystematicsConfig.BEAM_ANGLE)
        p_det = r(p)

        # print("Muon det vec:", p_det.X(), p_det.Y(), p_det.Z(), " |p|=", p_det.R())

        return p_det


    @staticmethod
    def Scale4VectorMomentum(vet,scale):
        return vet


    def PrimaryProtonScore(self):
        score = self.MasterAnaDev_proton_score 
        return score if score>=0 else -1

    def PrimaryProtonTheta(self):
        theta = self.MasterAnaDev_proton_theta
        return math.degrees(theta) if theta>=0 else -1


    # =============== collection of all recoil energy stuff.========
    def GetEAvail(self):
        if self.GetAnaToolName() == "MasterAnaDev" or self.GetAnaToolName() == "CCQENu":
            return self.blob_recoil_E_tracker+self.blob_recoil_E_ecal 
        else:
            correction = AvailableECorrection.GetSmoothVisECorrection(self)*1e3 #Eavail correction 
            #return self.recoile_passive_tracker+self.recoile_passive_ecal
            return self.recoile_passive_tracker+self.recoile_passive_ecal+correction

    def GetNuEFuzz(self):
        try:
            fuzz = self.blob_nuefuzz_E_tracker+self.blob_nuefuzz_E_ecal
        except:
            fuzz = 0
        #return 0
        return max(0,fuzz) * SystematicsConfig.AVAILABLE_E_CORRECTION


    def GetLeakageCorrection(self): 
        return SystematicsConfig.LEAKAGE_CORRECTION(self.ElectronEnergyRaw()) + (SystematicsConfig.LEAKAGE_BIAS if self.nsigma is None else 0) 

    def Pi0_Additional_Leakage(self):
        return random.uniform(0,20) if self.mc and (self.mc_intType == 4 or (self.mc_current==2 and self.mc_intType==10)) else 0

    # A function wrapper that maps negative values to zero
    # returning MeV
    @Utilities.decorator_ReLU
    def AvailableEnergy(self):
        #return self.GetEAvail()* SystematicsConfig.AVAILABLE_E_CORRECTION
        if abs(SystematicsConfig.AnaNuPDG) == 12:
            return self.GetEAvail()* SystematicsConfig.AVAILABLE_E_CORRECTION -self.GetCorrection()
        else:
            return self.GetEAvail()

    @Utilities.decorator_ReLU
    def RecoilEnergy(self): 
        result = None
        if self.GetAnaToolName() == "MasterAnaDev" or self.GetAnaToolName() == "CCQENu":            
            return sum(self.__getattr__("blob_recoil_E_{}".format(i)) for i in ["tracker","ecal","hcal","od","nucl"])
        else:
            result = sum(self.__getattr__("recoile_passive_{}".format(i)) for i in ["tracker","ecal","hcal","odclus","nucl"])
        return result - self.GetCorrection() 

    def GetFuzzCorrection(self):
        return 0

    def GetEMEnergyShift(self):
        return self.prong_ECALCalibE[0]*SystematicsConfig.EM_ENERGY_SCALE_SHIFT_ECAL if self.nsigma is None else 0

    def dedxDelta(self):
        dedx = self.GetVecDouble("dEdXs")
        dedx_dz  = self.GetVecDouble("dEdXs_dz")
        dedx_proj  = self.GetVecDouble("dEdXs_proj")
        total = 0
        outer = 0
        if dedx.size() == 0:
            return 0.9999
        for i in range(0,dedx.size()):
            if dedx_proj[i]>400:
                outer+=dedx[i]
            total+=dedx[i]
        return outer/total 

class CVSystematicUniverse(ROOT.PythonMinervaUniverse, CVPythonUniverse):
    def __init__(self,chain,nsigma):
        #self.weight = None
        #self.tuning_weight = None
        super(CVSystematicUniverse,self).__init__(chain,0 if nsigma is None else nsigma)
        super(ROOT.PythonMinervaUniverse,self).InitWithoutSuper(chain, 1 if nsigma is not None else None)

    def SetEntry(self,n_entry):
        #go to another event, discard weight calculated.
        self.weight = None
        self.tuning_weight = None
        CVPythonUniverse.LLR = None
        super(CVSystematicUniverse,self).SetEntry(n_entry)
        super(ROOT.PythonMinervaUniverse,self).ResetWeight()

    def __repr__(self):
        return self.LatexName()+" Universe at sigma: "+ str(self.GetSigma())

class FluxUniverse( ROOT.PlotUtils.FluxUniverse(ROOT.PythonMinervaUniverse),CVPythonUniverse):
    def __init__(self,chain,universe_number):
        super(FluxUniverse,self).__init__(chain,1,universe_number)
        super(ROOT.PlotUtils.FluxUniverse(ROOT.PythonMinervaUniverse),self).InitWithoutSuper(chain,1)

    @staticmethod
    def GetSystematicsUniverses(chain):
        return [FluxUniverse(chain,  i ) for i in range(0,SystematicsConfig.NUM_FLUX_UNIVERSE)]

class GenieUniverse(ROOT.PlotUtils.GenieUniverse(ROOT.PythonMinervaUniverse), CVPythonUniverse):
    def __init__(self,chain,nsigma,universe_name):
        super(GenieUniverse,self).__init__(chain,nsigma,universe_name)
        super(ROOT.PlotUtils.GenieUniverse(ROOT.PythonMinervaUniverse),self).InitWithoutSuper(chain,nsigma)

    @staticmethod
    def GetSystematicsUniverses(chain ):
        return [GenieUniverse(chain,i,j) for j in SystematicsConfig.GENIE_UNIVERSES for i in OneSigmaShift]

class GenieRvx1piUniverse(ROOT.PlotUtils.GenieRvx1piUniverse(ROOT.PythonMinervaUniverse),CVPythonUniverse):
    def __init__(self,chain,nsigma,universe_name):
        super(GenieRvx1piUniverse,self).__init__(chain,nsigma,universe_name)
        super(ROOT.PlotUtils.GenieRvx1piUniverse(ROOT.PythonMinervaUniverse),self).InitWithoutSuper(chain,nsigma)

    @staticmethod
    def GetSystematicsUniverses(chain ):
        return [GenieRvx1piUniverse(chain,i,j) for j in ["Rvn1pi","Rvp1pi"] for i in OneSigmaShift]

class GenieFaCCQEUniverse(ROOT.PlotUtils.GenieFaCCQEUniverse(ROOT.PythonMinervaUniverse),CVPythonUniverse):
    def __init__(self,chain,nsigma,universe_number):
        super(GenieFaCCQEUniverse,self).__init__(chain,nsigma,universe_number)
        super(ROOT.PlotUtils.GenieFaCCQEUniverse(ROOT.PythonMinervaUniverse),self).InitWithoutSuper(chain,nsigma)

    @staticmethod
    def GetSystematicsUniverses(chain ):
        if ROOT.PythonMinervaUniverse.UseZExpansionFaReweight():
            return [GenieFaCCQEUniverse(chain,1.0,i) for i in range(0,SystematicsConfig.NumZExpansionUniverses)]
        else:
            return [GenieFaCCQEUniverse(chain,i,-1) for i in OneSigmaShift]

class GenieNormCCResUniverse(ROOT.PlotUtils.GenieNormCCResUniverse(ROOT.PythonMinervaUniverse),CVPythonUniverse):
    def __init__(self,chain,nsigma):
        super(GenieNormCCResUniverse,self).__init__(chain,nsigma)
        super(ROOT.PlotUtils.GenieNormCCResUniverse(ROOT.PythonMinervaUniverse),self).InitWithoutSuper(chain,nsigma)

    @staticmethod
    def GetSystematicsUniverses(chain ):
        return [GenieNormCCResUniverse(chain,i) for i in OneSigmaShift]

class GenieMaResUniverse(ROOT.PlotUtils.GenieMaResUniverse(ROOT.PythonMinervaUniverse),CVPythonUniverse):
    def __init__(self,chain,nsigma):
        super(GenieMaResUniverse,self).__init__(chain,nsigma)
        super(ROOT.PlotUtils.GenieMaResUniverse(ROOT.PythonMinervaUniverse),self).InitWithoutSuper(chain,nsigma)

    @staticmethod
    def GetSystematicsUniverses(chain ):
        return [GenieMaResUniverse(chain,i) for i in OneSigmaShift]

class GenieMvResUniverse(ROOT.PlotUtils.GenieMvResUniverse(ROOT.PythonMinervaUniverse),CVPythonUniverse):
    def __init__(self,chain,nsigma):
        super(GenieMvResUniverse,self).__init__(chain,nsigma)
        super(ROOT.PlotUtils.GenieMvResUniverse(ROOT.PythonMinervaUniverse),self).InitWithoutSuper(chain,nsigma)

    @staticmethod
    def GetSystematicsUniverses(chain ):
        return [GenieMvResUniverse(chain,i) for i in OneSigmaShift]

# I don't think nu_e analysis need minos shift
# class MinosUniverse(ROOT.PlotUtils.GenieUniverse(ROOT.PythonMinervaUniverse), CVPythonUniverse):

class Universe2p2h(ROOT.PlotUtils.Universe2p2h(ROOT.PythonMinervaUniverse), CVPythonUniverse):
    def __init__(self,chain,universe_number):
        super(Universe2p2h,self).__init__(chain,1,universe_number)
        super(ROOT.PlotUtils.Universe2p2h(ROOT.PythonMinervaUniverse),self).InitWithoutSuper(chain,1)
        self.universe_number = universe_number

    @staticmethod
    def GetSystematicsUniverses(chain ):
        return [Universe2p2h(chain, i) for i in SystematicsConfig.UNIVERSES_2P2H]

class RPAUniverse(ROOT.PlotUtils.RPAUniverse(ROOT.PythonMinervaUniverse), CVPythonUniverse):
    def __init__(self,chain, universe_number,q2_region):
        super(RPAUniverse,self).__init__(chain,1,universe_number,q2_region)
        super(ROOT.PlotUtils.RPAUniverse(ROOT.PythonMinervaUniverse),self).InitWithoutSuper(chain, 1)

    @staticmethod
    def GetSystematicsUniverses(chain ):
        return [RPAUniverse(chain, i,j) for j in SystematicsConfig.RPA_UNIVERSES for i in SystematicsConfig.RPA_UNIVERSES[j] ]

class ResponseUniverse(ROOT.PlotUtils.ResponseUniverse(ROOT.PythonMinervaUniverse), CVPythonUniverse):
    def __init__(self,chain,nsigma,name):
        super(ResponseUniverse,self).__init__(chain,nsigma,name)
        super(ROOT.PlotUtils.ResponseUniverse(ROOT.PythonMinervaUniverse),self).InitWithoutSuper(chain,1)
        self.re = "^blob_.*_E_(tracker|ecal|od|nucl|hcal)$"

    def __getattr__(self,attrName):
        if re.match(self.re,attrName) is None:
            return super(ResponseUniverse,self).__getattr__(attrName)
        else:
            r = super(ResponseUniverse,self).__getattr__(attrName)
            shift = super(ResponseUniverse,self).__getattr__("{}_{}".format(attrName,self.m_name))* self.m_frac_unc
            return r+self.nsigma*shift

    def RecoilClustersCalorimetry(self,splined):
        return super(ResponseUniverse,self).RecoilClustersCalorimetry(splined)+self.GetRecoilShift()

    @staticmethod
    def GetSystematicsUniverses(chain):
        return [ResponseUniverse(chain,i,j) for j in SystematicsConfig.RESPONSE_BRANCHES for i in OneSigmaShift]

class LowQ2PionUniverse(ROOT.PlotUtils.LowQ2PionUniverse(ROOT.PythonMinervaUniverse),CVPythonUniverse):
    def __init__(self,chain,nsigma,channel):
        super(LowQ2PionUniverse,self).__init__(chain,nsigma)
        super(ROOT.PlotUtils.LowQ2PionUniverse(ROOT.PythonMinervaUniverse),self).InitWithoutSuper(chain,nsigma)
        self.channel = channel.upper() if channel is not None else None

    # def GetMyLowQ2PiWeight(self):
    #     if self.channel is None:
    #         return 1.0

    #     return super(LowQ2PionUniverse, self).GetLowQ2PiWeight(self.channel)
    def GetMyLowQ2PiWeight(self, channel=None):

        if channel is None:
            channel = self.channel

        if channel is None:
            return 1.0

        return super(LowQ2PionUniverse, self).GetLowQ2PiWeight(
            channel.upper()
        )


    # def GetMyLowQ2PiWeight(self):
    #     if self.nsigma == 0 or self.channel is None:
    #         return 1
    #     else :
    #         return super(LowQ2PionUniverse,self).GetLowQ2PiWeight(self.channel)

    # @staticmethod
    # def GetSystematicsUniverses(chain):
    #     cvshifts = [LowQ2PionUniverse(chain,i,SystematicsConfig.LowQ2PiWeightChannel) for i in OneSigmaShift] if SystematicsConfig.LowQ2PiWeightChannel is not None else []
    #     #cvshifts.extend([LowQ2PionUniverse(chain,i,j) for j in SystematicsConfig.LowQ2PiWeightSysChannel for i in ([-1,1] if j is not None else [0])])
    #     return cvshifts
    @staticmethod
    def GetSystematicsUniverses(chain):

        if USE_MNVTUNE_V1:
            channel = (
                V1_LOWQ2PI_CHANNEL
                if V1_ADD_LOWQ2PI
                else None
            )
        elif USE_MNVTUNE_V4P4P1:
            channel = (
                V4P4P1_LOWQ2PI_CHANNEL
                if V4P4P1_ADD_LOWQ2PI
                else None
            )
        else:
            channel = SystematicsConfig.LowQ2PiWeightChannel

        if channel is None:
            return []

        return [
            LowQ2PionUniverse(chain, i, channel)
            for i in OneSigmaShift
        ]


class LowQ2PionUniverseAlt(CVPythonUniverse):
    def __init__(self,chain,channel):
        super(LowQ2PionUniverseAlt,self).__init__(chain, 0)
        self.channel = channel

    def ShortName(self):
        return "LowQ2Pi_{}".format(self.channel)

    def LatexName(self):
        return "Low Q2 Pion Suppression Channel {}".format(self.channel)

    def GetMyLowQ2PiWeight(self):
        return super(LowQ2PionUniverseAlt,self).GetMyLowQ2PiWeight(self.channel)

    def IsVerticalOnly(self):
        return True

    @staticmethod
    def GetSystematicsUniverses(chain):
        return [LowQ2PionUniverseAlt(chain,i) for i in SystematicsConfig.LowQ2PiWeightSysChannel]

class MuonAngleXResolutionUniverse(ROOT.PlotUtils.MuonAngleXResolutionUniverse(ROOT.PythonMinervaUniverse),CVPythonUniverse):
    def __init__(self,chain,nsigma):
        super(MuonAngleXResolutionUniverse,self).__init__(chain,nsigma)
        super(ROOT.PlotUtils.MuonAngleXResolutionUniverse(ROOT.PythonMinervaUniverse),self).InitWithoutSuper(chain,nsigma)


    @staticmethod
    def GetSystematicsUniverses(chain):
        return [MuonAngleXResolutionUniverse(chain,i) for i in OneSigmaShift]

class MuonResolutionUniverse(ROOT.PlotUtils.MuonResolutionUniverse(ROOT.PythonMinervaUniverse),CVPythonUniverse):
    def __init__(self,chain,nsigma):
        super(MuonResolutionUniverse,self).__init__(chain,nsigma)
        super(ROOT.PlotUtils.MuonResolutionUniverse(ROOT.PythonMinervaUniverse),self).InitWithoutSuper(chain,nsigma)

    @staticmethod
    def GetSystematicsUniverses(chain):
        return [MuonResolutionUniverse(chain,i) for i in OneSigmaShift]

class MuonAngleYResolutionUniverse(ROOT.PlotUtils.MuonAngleYResolutionUniverse(ROOT.PythonMinervaUniverse),CVPythonUniverse):
    def __init__(self,chain,nsigma):
        super(MuonAngleYResolutionUniverse,self).__init__(chain,nsigma)
        super(ROOT.PlotUtils.MuonAngleYResolutionUniverse(ROOT.PythonMinervaUniverse),self).InitWithoutSuper(chain,nsigma)


    @staticmethod
    def GetSystematicsUniverses(chain):
        return [MuonAngleYResolutionUniverse(chain,i) for i in OneSigmaShift]

class MuonUniverseMinerva(ROOT.PlotUtils.MuonUniverseMinerva(ROOT.PythonMinervaUniverse),CVPythonUniverse):
    def __init__(self,chain,nsigma):
        super(MuonUniverseMinerva,self).__init__(chain,nsigma)
        super(ROOT.PlotUtils.MuonUniverseMinerva(ROOT.PythonMinervaUniverse),self).InitWithoutSuper(chain,nsigma)


    @staticmethod
    def GetSystematicsUniverses(chain):
        return [MuonUniverseMinerva(chain,i) for i in OneSigmaShift]


class MuonUniverseMinos(ROOT.PlotUtils.MuonUniverseMinos(ROOT.PythonMinervaUniverse),CVPythonUniverse):
    def __init__(self,chain,nsigma):
        super(MuonUniverseMinos,self).__init__(chain,nsigma)
        super(ROOT.PlotUtils.MuonUniverseMinos(ROOT.PythonMinervaUniverse),self).InitWithoutSuper(chain,nsigma)


    @staticmethod
    def GetSystematicsUniverses(chain):
        return [MuonUniverseMinos(chain,i) for i in OneSigmaShift]

class MinosEfficiencyUniverse(ROOT.PlotUtils.MinosEfficiencyUniverse(ROOT.PythonMinervaUniverse),CVPythonUniverse):
    def __init__(self,chain,nsigma):
        super(MinosEfficiencyUniverse,self).__init__(chain,nsigma)
        super(ROOT.PlotUtils.MinosEfficiencyUniverse(ROOT.PythonMinervaUniverse),self).InitWithoutSuper(chain,nsigma)

    def GetMyMinosEfficiencyWeight(self):
         if self._is_truth_mode():
             return 1.0
         if self.HasNoBackExitingTracks:
             return 1
         else:
             return super(MinosEfficiencyUniverse,self).GetMinosEfficiencyWeight()


    @staticmethod
    def GetSystematicsUniverses(chain):
        return [MinosEfficiencyUniverse(chain,i) for i in OneSigmaShift]

###########################################################################
class ElectronEnergyShiftUniverse(CVSystematicUniverse):
    def __init__(self,chain, nsigma,region):
        super(ElectronEnergyShiftUniverse,self).__init__(chain, nsigma)
        self.region = region

    def GetEMEnergyShift(self):
        return self.nsigma*SystematicsConfig.EM_ENERGY_SCALE_UNCERTAINTY[self.region]*self.GetVecElem("prong_"+self.region+"CalibE",0)

    def ShortName(self):
        return "elE_"+self.region

    def LatexName(self):
        return "EM Energy Scale in "+ self.region

    @staticmethod
    def GetSystematicsUniverses(chain):
        return [ElectronEnergyShiftUniverse(chain, i,region) for region in SystematicsConfig.EM_ENERGY_SCALE_UNCERTAINTY for i in OneSigmaShift]

###########################################################################
# # Ryan's electron energy scale class
# class ElectronEnergyScaleUniverse(CVSystematicUniverse):
#     def __init__(self,chain, nsigma):
#         super(ElectronEnergyScaleUniverse,self).__init__(chain, nsigma)

#     def ElectronEnergyRaw(self):
#         return self.nsigma*SystematicsConfig.ELECTRON_ENERGY_SCALE* super(ElectronEnergyScaleUniverse,self).ElectronEnergyRaw()

#     def ShortName(self):
#         return "electron_scale"

#     def LatexName(self):
#         return "Electron Energy Scale"

#     @staticmethod
#     def GetSystematicsUniverses(chain):
#         return [ElectronEnergyScaleUniverse(chain, i)  for i in OneSigmaShift]

# # My old electron energy scale class
# class ElectronEnergyScaleUniverse(CVSystematicUniverse):
#     def __init__(self, chain, nsigma):
#         super(ElectronEnergyScaleUniverse, self).__init__(chain, nsigma)

#     def _scale_factor(self):
#         return 1.0 + self.nsigma * SystematicsConfig.ELECTRON_ENERGY_SCALE

#     def ElectronP3D_det(self):
#         """
#         Apply a global ±1σ electron energy/momentum scale shift
#         to the corrected detector-coordinate electron momentum.

#         This makes ElectronEnergy(), ElectronP3D(), and EN4 all see
#         the same shifted reconstructed lepton energy.
#         """
#         p = super(ElectronEnergyScaleUniverse, self).ElectronP3D_det()
#         return p * self._scale_factor()

#     def ElectronP3D_det_noEMshift(self):
#         """
#         Apply the same global scale to the no-EM-shift momentum used by
#         ElectronEnergyRaw() and leakage-related correction paths.
#         """
#         p = super(ElectronEnergyScaleUniverse, self).ElectronP3D_det_noEMshift()
#         return p * self._scale_factor()

#     def ShortName(self):
#         return "electron_scale"

#     def LatexName(self):
#         return "Electron Energy Scale"

#     @staticmethod
#     def GetSystematicsUniverses(chain):
#         return [ElectronEnergyScaleUniverse(chain, i) for i in OneSigmaShift]
class ElectronEnergyScaleUniverse(CVSystematicUniverse):
    def __init__(self, chain, nsigma):
        super(ElectronEnergyScaleUniverse, self).__init__(chain, nsigma)

    def _scale_factor(self):
        return 1.0 + self.nsigma * SystematicsConfig.ELECTRON_ENERGY_SCALE

    def ElectronEnergyRaw(self):
        """
        Apply the global ±1σ electron energy-scale shift
        to the reconstructed electron energy.
        """
        e = super(ElectronEnergyScaleUniverse, self).ElectronEnergyRaw()
        return e * self._scale_factor()

    def ElectronP3D_det(self):
        """
        Apply the same global ±1σ scale shift to the magnitude of the
        corrected detector-coordinate electron momentum.

        The direction is unchanged, so angular observables are unaffected
        by this systematic while momentum-dependent observables shift
        consistently with ElectronEnergy().
        """
        p = super(ElectronEnergyScaleUniverse, self).ElectronP3D_det()
        return p * self._scale_factor()

    def ElectronP3D_det_noEMshift(self):
        """
        Apply the same global scale to the no-EM-shift momentum path.
        """
        p = super(ElectronEnergyScaleUniverse, self).ElectronP3D_det_noEMshift()
        return p * self._scale_factor()

    def ShortName(self):
        return "electron_scale"

    def LatexName(self):
        return "Electron Energy Scale"

    @staticmethod
    def GetSystematicsUniverses(chain):
        return [ElectronEnergyScaleUniverse(chain, i) for i in OneSigmaShift]


###########################################################################
class ElectronAngleShiftUniverse(CVSystematicUniverse):
    def __init__(self,chain, nsigma):
        super(ElectronAngleShiftUniverse,self).__init__(chain, nsigma)
        self.axis_angle = random.random()*math.pi
        self.shift_angle = random.gauss(0,self.nsigma*SystematicsConfig.LEPTON_ANGLE_UNCERTAINTY)

    def ShortName(self):
        return "eltheta"

    def LatexName(self):
        return "Election Candidate Angle"

    def ElectronP3D(self):
        p = super(ElectronAngleShiftUniverse,self).ElectronP3D()
        unit_3_vector = p/p.R()
        normal_vector1 = p.Cross(ROOT.Math.XYZVector(1.0,0,0))
        if normal_vector1.R() != 0:
            normal_vector1 *= 1.0/normal_vector1.R()
        r1 = ROOT.Math.AxisAngle(unit_3_vector,self.axis_angle)
        rotation_axis = r1(normal_vector1)
        r2 = ROOT.Math.AxisAngle(rotation_axis,self.shift_angle)
        return r2(p)

    def ElectronTheta2(self):
        theta = super(ElectronAngleShiftUniverse,self).ElectronTheta()
        return theta+self.shift_angle

    @staticmethod
    def GetSystematicsUniverses(chain ):
        return [ElectronAngleShiftUniverse(chain,  1) for i in range(SystematicsConfig.NUM_UNIVERSE)]

    

###########################################################################
class BirksShiftUniverse(CVSystematicUniverse):
    def __init__(self,chain, nsigma):
        super(BirksShiftUniverse,self).__init__(chain, nsigma)

    @property
    def prong_part_score(self):
        return self.prong_part_score_BirkUP if self.nsigma>0 else self.prong_part_score_BirkDOWN

    @property
    def prong_dEdXMeanFrontTracker(self):
        return self.prong_dEdXMeanFrontTracker_BirkUP if self.nsigma>0 else self.prong_dEdXMeanFrontTracker_BirkDOWN

    @property
    def prong_NonMIPClusFrac(self):
        return self.prong_NonMIPClusFrac_BirkUP if self.nsigma>0 else self.prong_NonMIPClusFrac_BirkDOWN

    def ShortName(self):
        return "birks"

    def LatexName(self):
        return "Birk's Constant"

    @staticmethod
    def GetSystematicsUniverses(chain ):
        return [BirksShiftUniverse(chain, i) for i in OneSigmaShift]

###########################################################################

class BeamAngleShiftUniverse(CVSystematicUniverse):
    def __init__(self,chain, nsigma, x):
        super(BeamAngleShiftUniverse,self).__init__(chain, nsigma)
        self.rotation =  ROOT.Math.RotationX(SystematicsConfig.BEAM_XANGLE_UNCERTAINTY*nsigma) if x else ROOT.Math.RotationY(SystematicsConfig.BEAM_YANGLE_UNCERTAINTY*nsigma)

    def ShortName(self):
        return "beam_angle"

    def LatexName(self):
        return "Beam Angle"

    def ElectronP3D(self):
        p = super(BeamAngleShiftUniverse,self).ElectronP3D()
        return self.rotation(p)

    @staticmethod
    def GetSystematicsUniverses(chain ):
        return [BeamAngleShiftUniverse(chain, i, j) for j in [True,False] for i in OneSigmaShift]


class GeantHadronUniverse(ROOT.PlotUtils.GeantHadronUniverse(ROOT.PythonMinervaUniverse), CVPythonUniverse):
    def __init__(self,chain,nsigma,pdg):
        super(GeantHadronUniverse,self).__init__(chain,nsigma,pdg)
        super(ROOT.PlotUtils.GeantHadronUniverse(ROOT.PythonMinervaUniverse),self).InitWithoutSuper(chain,nsigma)

    @staticmethod
    def GetSystematicsUniverses(chain):
        return [GeantHadronUniverse(chain,i,j) for j in SystematicsConfig.GEANT_PARTICLES for i in OneSigmaShift]

class TargetMassUniverse(ROOT.PlotUtils.TargetMassScintillatorUniverse(ROOT.PythonMinervaUniverse),CVPythonUniverse):
    def __init__(self,chain,nsigma):
        super(TargetMassUniverse,self).__init__(chain,nsigma)
        super(ROOT.PlotUtils.TargetMassScintillatorUniverse(ROOT.PythonMinervaUniverse),self).InitWithoutSuper(chain,nsigma)
    def GetStandardWeight(self):
        weight = super(TargetMassUniverse,self).GetStandardWeight()
        weight*= self.GetWeightRatioToCV()
        return weight

    @staticmethod
    def GetSystematicsUniverses(chain):
        return [TargetMassUniverse(chain, i) for i in OneSigmaShift]

class MKModelUniverse(CVSystematicUniverse):
    def __init__(self,chain,nsigma):
        super(MKModelUniverse,self).__init__(chain,nsigma)
        self.reweighter = ROOT.PlotUtils.MKReweighter(ROOT.PythonMinervaUniverse,ROOT.PlotUtils.detail.empty)()

    def IsVerticalOnly(self):
        return True

    def GetStandardWeight(self):
        weight = super(MKModelUniverse,self).GetStandardWeight()
        weight*= self.reweighter.GetWeight(self,ROOT.PlotUtils.detail.empty())
        return weight

    def ShortName(self):
        return "MK_model"

    def LatexName(self):
        return "MK Model"

    @staticmethod
    def GetSystematicsUniverses(chain):
        return [MKModelUniverse(chain,1)]


class FSIWeightUniverse(CVSystematicUniverse):
    def __init__(self,chain,nsigma,iweight):
        super(FSIWeightUniverse,self).__init__(chain,nsigma)
        self.iweight = iweight
        self.reweighter = ROOT.PlotUtils.FSIReweighter(ROOT.PythonMinervaUniverse,ROOT.PlotUtils.detail.empty)((iweight+1)//2,(iweight+1)%2)

    def IsVerticalOnly(self):
        return True

    def GetStandardWeight(self):
        weight = super(FSIWeightUniverse,self).GetStandardWeight()
        weight*= self.reweighter.GetWeight(self,ROOT.PlotUtils.detail.empty())
        return weight

    def DebugExtraWeight(self, label=""):
        base = super(FSIWeightUniverse, self).GetStandardWeight()
        extra = self.reweighter.GetWeight(self, ROOT.PlotUtils.detail.empty())
        total = self.GetStandardWeight()

        print(
            "[FSI_EXTRA] {} iweight={} "
            "intType={} current={} Ev={:.6g} Q2={:.6g} W={:.6g} "
            "base={:.8g} extra={:.8g} total={:.8g} total/base={:.8g}".format(
                label,
                self.__dict__.get("iweight", "NA"),
                self.mc_intType,
                self.mc_current,
                self.mc_incomingE * 1e-3,
                self.mc_Q2 / 1e6,
                self.mc_w / 1e3,
                base,
                extra,
                total,
                total / base if base else -999,
            )
        )

    def ShortName(self):
        return "fsi_weight"

    def LatexName(self):
        return "FSI weight"

    @staticmethod
    def GetSystematicsUniverses(chain):
        return [FSIWeightUniverse(chain,1,i) for i in range(3)]

class SusaValenciaUniverse(CVSystematicUniverse):
    def __init__(self,chain,nsigma):
        super(SusaValenciaUniverse,self).__init__(chain,nsigma)
        self.reweighter = ROOT.PlotUtils.SuSAFromValencia2p2hReweighter(ROOT.PythonMinervaUniverse,ROOT.PlotUtils.detail.empty)()

    def IsVerticalOnly(self):
        return True

    def GetLowRecoil2p2hWeight(self):
        return 1.0

    def GetRPAWeight(self):
        return 1.0

    def ShortName(self):
        return "SuSA_Valencia_Weight"

    def LatexName(self):
        return "SuSA Valencia Weight"

    def GetStandardWeight(self):
        weight = super(SusaValenciaUniverse,self).GetStandardWeight()
        weight*= self.reweighter.GetWeight(self,ROOT.PlotUtils.detail.empty())
        return weight

    @staticmethod
    def GetSystematicsUniverses(chain):
        return [SusaValenciaUniverse(chain,1)]

class LeakageUniverse(CVSystematicUniverse):
    def __init__(self,chain,nsigma):
        super(LeakageUniverse,self).__init__(chain,nsigma)

    def GetLeakageCorrection(self): 
        return super(LeakageUniverse,self).GetLeakageCorrection() + ( self.nsigma*SystematicsConfig.LEAKAGE_SYSTEMATICS if abs(self.mc_primaryLepton) ==11 else 0)

    def ShortName(self):
        return "Leakage_Uncertainty"

    def LatexName(self):
        return "Leakage Uncertainty"

    @staticmethod
    def GetSystematicsUniverses(chain):
        return [LeakageUniverse(chain,i) for i in OneSigmaShift]


def GetAllSystematicsUniverses(chain,is_data,is_pc =False,exclude=None,playlist=None):
    #use list because I want to control the order, process cv first, then vertical shifts, then lateral shifts
    universes = []
    CVPythonUniverse.is_pc = is_pc
    CVPythonUniverse.LeafGetters= {}
    if is_data:
        #data has only cv universe
        universes.append(CVSystematicUniverse(chain, None))
    else:
        #Set Playlist, only MC cares about playlist
        if playlist is None:
            #attempt guess playlist from the chain
            playlist = Utilities.PlaylistLookup(chain.GetValue("mc_run",0))

        if not is_pc:
            CVSystematicUniverse.SetPlaylist(playlist)
            #Set NuE constraint
            CVSystematicUniverse.SetNuEConstraint(SystematicsConfig.USE_NUE_CONSTRAINT)
        else:
            CVSystematicUniverse.SetPlaylist("dummy")
            CVSystematicUniverse.SetNuEConstraint(False)

        #Set nu_type
        CVSystematicUniverse.SetAnalysisNuPDG(SystematicsConfig.AnaNuPDG)


        # ------------------------------------------------------------
        # GENIE CV tune configuration
        # ------------------------------------------------------------
        active_tunes = [
            USE_MNVTUNE_V1,
            USE_MNVTUNE_V4P3P0,
            USE_MNVTUNE_V4P4P1,
        ]

        if sum(active_tunes) > 1:
            raise RuntimeError(
                "Only one MINERvA tune can be active."
            )

        CVSystematicUniverse.SetNonResPiReweight(True)

        if USE_MNVTUNE_V1:
            tune_name = "v1"
            use_deuterium = V1_ADD_DEUTERIUM
            use_zexp = V1_ADD_ZEXP

        elif USE_MNVTUNE_V4P3P0:
            tune_name = "v4.3.0"
            use_deuterium = True
            use_zexp = False

        elif USE_MNVTUNE_V4P4P1:
            tune_name = "v4.4.1"
            use_deuterium = True
            use_zexp = V4P4P1_ADD_ZEXP

        else:
            tune_name = "legacy_v4"
            use_deuterium = True
            use_zexp = V1_ADD_ZEXP

        CVSystematicUniverse.SetDeuteriumGeniePiTune(
            use_deuterium
        )
        CVSystematicUniverse.SetZExpansionFaReweight(
            use_zexp
        )


        CVSystematicUniverse.SetNFluxUniverses(
            SystematicsConfig.NUM_FLUX_UNIVERSE
        )

        CVSystematicUniverse.RPAMaterials(True)

        # MHRW renorm directory
        if not hasattr(CVPythonUniverse, "_mhrw_filename_set"):
            mhrw_dir = os.path.join(
                os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                "custom_mhrwKineRenorm",
            )

            CVSystematicUniverse.SetMHRWeightFilename(
                "Renorm_Kine_Truth",
                mhrw_dir,
            )

            CVPythonUniverse._mhrw_filename_set = True

            # print("[MHRW CONFIG] directory =", mhrw_dir)

        # if chain.GetTree().GetName() == "Truth":
        #     CVSystematicUniverse.SetTruth(True)
        tree_name = chain.GetTree().GetName()
        is_truth_tree = (tree_name == "Truth")
        CVSystematicUniverse.SetTruth(is_truth_tree)
        # print(
        #     "[UNIVERSE MODE]",
        #     "tree =", tree_name,
        #     "SetTruth =", is_truth_tree,
        # )


        #append cv universe
        universes.append(CVSystematicUniverse(chain,0))

        if exclude is None or "all" not in exclude:
            # Vertical shift first to skip some cut calculation

            # #Electron momentum universe
            if abs(SystematicsConfig.AnaNuPDG)==12:
                universes.extend(ElectronEnergyShiftUniverse.GetSystematicsUniverses(chain ))
            elif abs(SystematicsConfig.AnaNuPDG)==14:
                universes.extend(MuonUniverseMinerva.GetSystematicsUniverses(chain ))
                universes.extend(MuonUniverseMinos.GetSystematicsUniverses(chain ))
                universes.extend(MuonResolutionUniverse.GetSystematicsUniverses(chain ))
                universes.extend(MuonAngleXResolutionUniverse.GetSystematicsUniverses(chain ))
                universes.extend(MuonAngleYResolutionUniverse.GetSystematicsUniverses(chain ))
                universes.extend(MinosEfficiencyUniverse.GetSystematicsUniverses(chain ))
            else:
                raise ValueError ("AnaNuPDG should be \pm 12 or 14, but you set {}".format(SystematicsConfig.AnaNuPDG))

            #Electron angle universe
            universes.extend(ElectronAngleShiftUniverse.GetSystematicsUniverses(chain ))

            #Electron momentum universe
            #universes.extend(ElectronEnergyShiftUniverse.GetSystematicsUniverses(chain ))

            #Electron energy scale universe
            universes.extend(ElectronEnergyScaleUniverse.GetSystematicsUniverses(chain))

            #beam angle shift universe
            universes.extend(BeamAngleShiftUniverse.GetSystematicsUniverses(chain ))

            #particle response shift universe
            universes.extend(ResponseUniverse.GetSystematicsUniverses(chain ))

            #Flux universe
            universes.extend(FluxUniverse.GetSystematicsUniverses(chain ))

            #Genie universe
            universes.extend(GenieUniverse.GetSystematicsUniverses(chain ))
            universes.extend(GenieRvx1piUniverse.GetSystematicsUniverses(chain ))
            universes.extend(GenieFaCCQEUniverse.GetSystematicsUniverses(chain ))
            universes.extend(GenieMaResUniverse.GetSystematicsUniverses(chain ))
            universes.extend(GenieMvResUniverse.GetSystematicsUniverses(chain ))
            universes.extend(GenieNormCCResUniverse.GetSystematicsUniverses(chain ))

            #2p2h universes
            universes.extend(Universe2p2h.GetSystematicsUniverses(chain ))

            #RPA universe:
            universes.extend(RPAUniverse.GetSystematicsUniverses(chain ))

            #Non resonant pion universe
            # #universes.extend(NonResonantPionUniverse.GetSystematicsUniverses(chain ))

            # Low-Q2 pion systematic only when that tune component is active
            if (not USE_MNVTUNE_V1) or V1_ADD_LOWQ2PI:
                universes.extend(
                    LowQ2PionUniverse.GetSystematicsUniverses(chain)
                )
            #universes.extend(LowQ2PionUniverseAlt.GetSystematicsUniverses(chain )) used for warping study variant

            # #birk shift universe
            ##universes.extend(BirksShiftUniverse.GetSystematicsUniverses(chain ))

            #MKModelUniverse
            universes.extend(MKModelUniverse.GetSystematicsUniverses(chain ))

            #FSIWeighUniverse
            universes.extend(FSIWeightUniverse.GetSystematicsUniverses(chain ))

            #SuSAValenciaUniverse
            universes.extend(SusaValenciaUniverse.GetSystematicsUniverses(chain ))

            #hadron reweight shifting universe
            universes.extend(GeantHadronUniverse.GetSystematicsUniverses(chain ))

            #leakage universe
            universes.extend(LeakageUniverse.GetSystematicsUniverses(chain ))

            #target mass universe
            universes.extend(TargetMassUniverse.GetSystematicsUniverses(chain ))


    # Group universes in dict.
    univ_dict = OrderedDict()
    for univ in universes:
        univ_dict.setdefault(univ.ShortName(), []).append(univ)

    if exclude is not None and len(univ_dict)>1:
        for i in exclude:
            if i in univ_dict:
                del univ_dict[i]

    return univ_dict

#code for testing new cv universe.
if __name__ == "__main__":

    #path = "/exp/minerva/data/users/hsu/merged_files/merged-1558038217.root"
    #path = "/pnfs/minerva/persistent/users/jyhan/NuECCQE/v21r1p1_recoilE_final/mcme1A/grid/central_value/minerva/ana/v21r1p1/00/11/00/17/SIM_minerva_00110017_Subruns_0325_NuECCQE_Ana_Tuple_v21r1p1.root"
    path = "/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me1A_DualVertex_p4/MasterAnaDev_mc_AnaTuple_run00110011_Playlist.root"
    
    chainWrapper = PlotUtils.ChainWrapper("MasterAnaDev")
    chainWrapper.Add(path)
    chainWrapper.GetValue("prong_part_score",0)
    univs = GetAllSystematicsUniverses(chainWrapper,False)
    for u in univs.keys():
        print(u)
        for univ in univs[u]:
            univ.SetEntry(5)
            weight = univ.GetWeight(False)
            print(weight)
        print("")