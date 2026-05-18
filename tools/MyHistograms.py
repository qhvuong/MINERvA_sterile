import os
import math
import ROOT
import PlotUtils
from PlotUtils.HistWrapper import HistWrapper
from array import array
from functools import partial
import copy
from tools.PlotLibrary import TranslateSettings,VariantPlotsNamingScheme
from config.SignalDef import TRUTH_CATEGORIES,EXTRA_OTHER
from config.CutConfig import SAMPLE_CUTS,KINEMATICS_CUTS
from tools.CutLibrary import CUTS
from tools.SystematicsUniverse import GetNeutrinoTravelledLengthPDF
from tools.PlotLibrary import is_pion_parent, is_muon_parent, is_kaon_parent

def DebugFluxBandCV(hist, label="", max_bins=20, max_univ=5):
    if hist is None:
        return

    if not hist.HasVertErrorBand("Flux"):
        return

    band = hist.GetVertErrorBand("Flux")

    print("\n[{}] {}".format(label, hist.GetName()))
    print("  title =", hist.GetTitle())
    print("  MAIN CV integral =", hist.Integral())
    print("  FLUX BAND CV integral =", band.Integral())
    print("  bandCV/main integral ratio =", band.Integral() / hist.Integral() if hist.Integral() else "NA")

    print("  bin-by-bin main CV vs Flux band CV:")
    nprinted = 0
    for b in range(0, hist.GetNbinsX() + 2):
        main = hist.GetBinContent(b)
        bcv = band.GetBinContent(b)

        if main != 0 or bcv != 0:
            ratio = bcv / main if main != 0 else "NA"
            print(
                "    bin {:2d} main={:.12g} fluxBandCV={:.12g} ratio={}".format(
                    b, main, bcv, ratio
                )
            )
            nprinted += 1
            if nprinted >= max_bins:
                break

    print("  first few Flux universes:")
    for i in range(min(max_univ, band.GetNHists())):
        hu = band.GetHist(i)
        print(
            "    universe {:3d} integral={:.12g} univ/main={}".format(
                i,
                hu.Integral(),
                hu.Integral() / hist.Integral() if hist.Integral() else "NA"
            )
        )

        nprinted = 0
        for b in range(0, hist.GetNbinsX() + 2):
            uval = hu.GetBinContent(b)
            main = hist.GetBinContent(b)
            if uval != 0 or main != 0:
                print(
                    "      bin {:2d} main={:.12g} univ={:.12g} univ/main={}".format(
                        b,
                        main,
                        uval,
                        uval / main if main != 0 else "NA",
                    )
                )
                nprinted += 1
                if nprinted >= max_bins:
                    break

def PlotBandCVAndUniverses(hist, bandname="Flux", tag="", sync_first=False):
    """
    Plot parent/main CV, band CV, and all universes for one vertical error band.
    If sync_first=True, copy the parent CV into the band CV before plotting.
    """
    if hist is None:
        return

    if not hist.HasVertErrorBand(bandname):
        print(f"[PlotBandCVAndUniverses] {hist.GetName()} has no band {bandname}")
        return

    if sync_first:
        # sync copied band CV to parent CV
        for bname in hist.GetVertErrorBandNames():
            band = hist.GetVertErrorBand(bname)
            if not band:
                continue
            for ibin in range(0, hist.GetNbinsX() + 2):
                band.SetBinContent(ibin, hist.GetBinContent(ibin))
                band.SetBinError(ibin, hist.GetBinError(ibin))

    band = hist.GetVertErrorBand(bandname)

    c = ROOT.TCanvas(f"c_{hist.GetName()}_{bandname}_{tag}", "", 1200, 900)

    # parent/main CV
    h_main = ROOT.TH1D(hist)
    h_main.SetDirectory(0)
    h_main.SetName(f"{hist.GetName()}_{bandname}_mainCV_{tag}")
    h_main.SetStats(0)
    h_main.SetLineColor(ROOT.kBlack)
    h_main.SetLineWidth(4)
    h_main.SetLineStyle(2)   # dashed

    # band CV
    h_bandcv = ROOT.TH1D(band)
    h_bandcv.SetDirectory(0)
    h_bandcv.SetName(f"{hist.GetName()}_{bandname}_bandCV_{tag}")
    h_bandcv.SetStats(0)
    h_bandcv.SetLineColor(ROOT.kBlue)
    h_bandcv.SetLineWidth(2)
    h_bandcv.SetMarkerColor(ROOT.kBlue)
    h_bandcv.SetMarkerStyle(20)
    h_bandcv.SetMarkerSize(1.1)

    # universes
    univs = []
    ymax = max(h_main.GetMaximum(), h_bandcv.GetMaximum())

    for i in range(band.GetNHists()):
        hu = ROOT.TH1D(band.GetHist(i))
        hu.SetDirectory(0)
        hu.SetName(f"{hist.GetName()}_{bandname}_univ_{i}_{tag}")
        hu.SetStats(0)
        hu.SetLineColor(ROOT.kGray + 1)
        hu.SetLineWidth(1)
        univs.append(hu)
        ymax = max(ymax, hu.GetMaximum())

    h_main.SetMaximum(1.25 * ymax if ymax > 0 else 1.0)
    h_main.Draw("HIST")

    for hu in univs:
        hu.Draw("HIST SAME")

    h_bandcv.Draw("HIST SAME")
    h_bandcv.Draw("P SAME")
    h_main.Draw("HIST SAME")  # redraw on top

    leg = ROOT.TLegend(0.55, 0.68, 0.88, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.AddEntry(h_main, "Main CV (dashed)", "l")
    leg.AddEntry(h_bandcv, "Band CV (blue markers)", "lp")
    if univs:
        leg.AddEntry(univs[0], "Universes", "l")
    leg.Draw()

    c.Modified()
    c.Update()

    # keep refs alive in PyROOT
    c._h_main = h_main
    c._h_bandcv = h_bandcv
    c._univs = univs
    c._leg = leg

    outname = f"{hist.GetName()}_{bandname}_{tag}.png"
    print("[PlotBandCVAndUniverses] saving", outname)
    c.SaveAs(outname)

def CheckBadBins1D(hist, label="", check_all_universes=True, max_bad=100):
    if hist is None:
        return False

    found_bad = False
    nbad = 0

    def bad(x):
        try:
            return not math.isfinite(float(x))
        except Exception:
            return True

    def print_bad(where, b, c, e):
        nonlocal nbad

        lo = hist.GetXaxis().GetBinLowEdge(b) if 1 <= b <= hist.GetNbinsX() else "UF/OF"
        hi = hist.GetXaxis().GetBinUpEdge(b) if 1 <= b <= hist.GetNbinsX() else "UF/OF"

        print(
            "[EN4_BAD_BIN] "
            "label={} hist={} where={} bin={} range=[{},{}] content={} error={}".format(
                label, hist.GetName(), where, b, lo, hi, c, e
            ),
            flush=True,
        )

        nbad += 1

    # Main CV
    for b in range(0, hist.GetNbinsX() + 2):
        c = hist.GetBinContent(b)
        e = hist.GetBinError(b)

        if bad(c) or bad(e):
            found_bad = True
            print_bad("MAIN_CV", b, c, e)

            if nbad >= max_bad:
                return found_bad

    # Error bands
    try:
        band_names = list(hist.GetVertErrorBandNames())
    except Exception:
        band_names = []

    for band_name in band_names:
        band = hist.GetVertErrorBand(str(band_name))
        if not band:
            continue

        # Band CV
        for b in range(0, hist.GetNbinsX() + 2):
            c = band.GetBinContent(b)
            e = band.GetBinError(b)

            if bad(c) or bad(e):
                found_bad = True
                print_bad("{}_BAND_CV".format(band_name), b, c, e)

                if nbad >= max_bad:
                    return found_bad

        # Universes
        if check_all_universes:
            for iu in range(band.GetNHists()):
                hu = band.GetHist(iu)

                for b in range(0, hu.GetNbinsX() + 2):
                    c = hu.GetBinContent(b)
                    e = hu.GetBinError(b)

                    if bad(c) or bad(e):
                        found_bad = True
                        print_bad("{}_UNIV_{}".format(band_name, iu), b, c, e)

                        if nbad >= max_bad:
                            return found_bad

    return found_bad

class HistWrapper1D(HistWrapper):
    def __init__(self,title,bins):
        # Note: bins is an array of lower edges of each bin
        super(HistWrapper1D,self).__init__(title,len(bins)-1,array('d',bins),{})

    @property
    def name(self):
        return self.hist.GetName()

    @name.setter
    def name(self,name):
        self.hist.SetName(name)

    # def AddErrorBands(self,cv_universes):
    #     self.AddUniverses(cv_universes)

    # def Fill(self, value, wgt, universe):
    #     self.FillUniverse(universe,value,wgt)
    def Clone(self,name):
        clone = copy.copy(self)
        clone.hist = self.hist.Clone(name)
        clone.eventToUnivMap = {}
        return clone

    # def Write(self):
    #     if not hasattr(self, "hist") or self.hist is None:
    #         print("[WARNING] HistWrapper1D.Write called but self.hist is missing", flush=True)
    #         return

    #     if self.hist.GetName() == "EN4":
    #         CheckBadBins1D(self.hist, "BEFORE_SYNC", check_all_universes=True, max_bad=20)

    #     self.SyncCVHistos()

    #     if self.hist.GetName() == "EN4":
    #         CheckBadBins1D(self.hist, "AFTER_SYNC", check_all_universes=True, max_bad=20)

    #     self.hist.Write()
    #     del self.hist
    def Write(self):
        if not hasattr(self, "hist") or self.hist is None:
            print("[WARNING] HistWrapper1D.Write called but self.hist is missing")
            return
        self.SyncCVHistos()
        self.hist.Write()
        del self.hist   # keep disabled for diagnostic

class HistWrapper2D(HistWrapper):
    def __init__(self,title,Xbins,Ybins):
        # Note: bins is an array of lower edges of each bin
        super(HistWrapper2D,self).__init__(title,len(Xbins)-1,array('d',Xbins),len(Ybins)-1,array('d',Ybins),{})

    @property
    def name(self):
        return self.hist.GetName()

    @name.setter
    def name(self,name):
        self.hist.SetName(name)

    def Clone(self,name):
        clone = copy.copy(self)
        clone.hist = self.hist.Clone(name)
        clone.eventToUnivMap = {}
        return clone

    def Write(self):
        # print("BEFORE SYNC", self.hist.GetName())
        # for i in range(0, self.hist.GetNbinsX()+2):
        #     c = self.hist.GetBinContent(i)
        #     if c != 0:
        #         print("  pre", i, c)

        self.SyncCVHistos()

        # print("AFTER SYNC", self.hist.GetName())
        # for i in range(0, self.hist.GetNbinsX()+2):
        #     c = self.hist.GetBinContent(i)
        #     if c != 0:
        #         print("  post", i, c)

        self.hist.Write()
        del self.hist

class MnvResponseWrapper(object):
    def __init__(self,title, bins):
        self.reco_hist = HistWrapper2D(title,bins[0],bins[1])
        self.truth_hist = HistWrapper2D(title,bins[2],bins[3])

        # +2 for overflow/underflow bins of reco/truth histogram
        migrationNbinsX = (self.reco_hist.hist.GetNbinsX()+2)*(self.reco_hist.hist.GetNbinsY()+2)
        migrationNbinsY = (self.truth_hist.hist.GetNbinsX()+2)*(self.truth_hist.hist.GetNbinsY()+2)
        self.migration_hist = HistWrapper2D(title, list(range(migrationNbinsX+1)),list(range(migrationNbinsY+1)))

    @property
    def name(self):
        return self.migration_hist.name

    @name.setter
    def name(self,name):
        self.reco_hist.name = name+"_reco"
        self.truth_hist.name = name+"_truth"
        self.migration_hist.name = name

    def AddUniverses(self,cv_universes):
        self.reco_hist.AddUniverses(cv_universes)
        self.truth_hist.AddUniverses(cv_universes)
        self.migration_hist.AddUniverses(cv_universes)

    def FillUniverse(self,universe, reco_x,reco_y,truth_x,truth_y, wgt):
        mig_x = self.migration_hist.hist.GetXaxis().GetBinCenter(super(PlotUtils.MnvH2D, self.reco_hist.hist).FindBin(reco_x,reco_y))
        mig_y = self.migration_hist.hist.GetYaxis().GetBinCenter(super(PlotUtils.MnvH2D, self.truth_hist.hist).FindBin(truth_x,truth_y))

        #print(mig_x,mig_y)

        self.reco_hist.FillUniverse(universe,reco_x,reco_y,wgt)
        self.truth_hist.FillUniverse(universe,truth_x,truth_y,wgt)
        self.migration_hist.FillUniverse(universe,mig_x,mig_y,wgt)

    def Write(self):
        self.reco_hist.Write()
        self.truth_hist.Write()
        self.migration_hist.Write()



def _is_bad_number(x):
    try:
        return not math.isfinite(float(x))
    except Exception:
        return True


def _safe_existing_attr(obj, attr, default="NA"):
    """
    Only read Python attributes already stored on the object.
    This avoids triggering CVUniverse dynamic branch lookup.
    """
    try:
        return object.__getattribute__(obj, "__dict__").get(attr, default)
    except BaseException:
        return default


def _debug_universe_context(universe):
    print("  universe type =", type(universe))
    print("  universe str  =", str(universe))

    # This requires adding: universe._debug_entry = counter
    # right after universe.SetEntry(counter) in eventSelection.py.
    print("  debug entry   =", _safe_existing_attr(universe, "_debug_entry"))

    clf = _safe_existing_attr(universe, "classifier", None)
    if clf is not None and clf != "NA":
        print("  side_band     =", _safe_existing_attr(clf, "side_band"))
        print("  truth_class   =", _safe_existing_attr(clf, "truth_class"))
        print("  is_signal     =", _safe_existing_attr(clf, "is_true_signal"))

    kc = _safe_existing_attr(universe, "kin_cal", None)
    if kc is not None and kc != "NA":
        print("  reco_E_lep    =", _safe_existing_attr(kc, "reco_E_lep"))
        print("  reco_visE     =", _safe_existing_attr(kc, "reco_visE"))
        print("  reco_E_nu_cal =", _safe_existing_attr(kc, "reco_E_nu_cal"))

EN4_BAD_FILL_PRINT_LIMIT = 20
EN4_BAD_FILL_PRINT_COUNT = 0

class PlotProcessor():
    def __init__(self,name,title,binning,value_getter = None, cuts = None, weight_function = None,**kwargs):
        
        if value_getter is None:
            print("You didn't tell me what to fill for plot: " + name)
            exit(1)
        
        if len(binning)!=len(value_getter):
            print("inconsistent number of binning and value getter.")
            exit(1)

        for l in binning:
            if len(l)<3 or l[-1]<=l[0]:
                print("wrong binning of plot {}: {}".format(name,l))
                exit(1)

        if len(value_getter) == 1:
            self.histwrapper = HistWrapper1D(title,binning[0])
        elif len (value_getter) == 2:
            self.histwrapper = HistWrapper2D(title,binning[0],binning[1])

        elif len(value_getter) == 4:
            self.histwrapper = MnvResponseWrapper(title, binning)
        else :
            print("Do not support more than 2D histograms yet.")
            exit(1)
        self.value_getter = value_getter
        self.cuts = []
        if cuts is not None:
            for i in cuts:
                self.AddCut(i)

        self.weight_function =  weight_function
        self.histwrapper.name = name

    def Clone(self, name):
        #shalow copy of attrs
        clone = copy.copy(self)
        #deep copy of histwrapper and cuts.
        clone.histwrapper = self.histwrapper.Clone(name)
        clone.cuts = list(self.cuts)
        return clone

    def AddErrorBands(self,universes):
        self.histwrapper.AddUniverses(universes)

    def FillHistWithLPDF(self, universe):
        name = self.histwrapper.name

        if "pionParent" in name and not is_pion_parent(universe):
            return
        if "muonParent" in name and not is_muon_parent(universe):
            return
        if "kaonParent" in name and not is_kaon_parent(universe):
            return

        enu_gev = universe.mc_incomingE * 1e-3
        if enu_gev <= 0:
            return

        w_base = self.weight_function(universe) if self.weight_function else universe.GetWeight()
        hL = GetNeutrinoTravelledLengthPDF(enu_gev, universe.mc_incoming)

        if name.startswith("drawnL_EReco_LE"):
            yval = universe.kin_cal.reco_E_lep + universe.kin_cal.reco_visE
            if yval is None:
                return

            for i in range(1, hL.GetNbinsX() + 1):
                pk = hL.GetBinContent(i)
                if pk <= 0:
                    continue
                Lk = hL.GetXaxis().GetBinCenter(i)
                loe = Lk / enu_gev
                self.histwrapper.FillUniverse(universe, loe, yval, w_base * pk)

        elif name.startswith("drawnL_ElepReco_LE"):
            yval = universe.kin_cal.reco_E_lep
            if yval is None:
                return

            for i in range(1, hL.GetNbinsX() + 1):
                pk = hL.GetBinContent(i)
                if pk <= 0:
                    continue
                Lk = hL.GetXaxis().GetBinCenter(i)
                loe = Lk / enu_gev
                self.histwrapper.FillUniverse(universe, loe, yval, w_base * pk)

        elif name.startswith("drawnL_nu_length"):
            for i in range(1, hL.GetNbinsX() + 1):
                pk = hL.GetBinContent(i)
                if pk <= 0:
                    continue
                Lk = hL.GetXaxis().GetBinCenter(i)
                self.histwrapper.FillUniverse(universe, Lk, w_base * pk)

    def Process(self, universe):

        if all(cut(universe) for cut in self.cuts[::-1]):
            try:
                value = [_(universe) for _ in self.value_getter]
            except Exception as e:
                print(self.histwrapper.name)
                print(("Error", e))
                return None
        else:
            return None

        # Special handling for template-weighted drawn-L plots
        name = self.histwrapper.name
        if (
            name.startswith("drawnL_EReco_LE")
            or name.startswith("drawnL_ElepReco_LE")
            or name.startswith("drawnL_nu_length")
        ):
            self.FillHistWithLPDF(universe)
            return

        wgt = self.weight_function(universe) if self.weight_function else universe.GetWeight()

        if isinstance(value[0], list):
            if not isinstance(wgt, list):
                wgt = len(value[0]) * [wgt]
            for v in map(lambda wgt,*args: (list(args),wgt), wgt, *value):
                self.FillHist(universe, *v)
        else:
            self.FillHist(universe, value, wgt)

    # def Process(self,universe):

    #     if all(cut(universe) for cut in self.cuts[::-1]):
    #         try:
    #             value = [_(universe) for _ in self.value_getter]
    #         except Exception as e:
    #             print(self.histwrapper.name)
    #             print(("Error", e))
    #             return None
    #     else :
    #         return None

    #     wgt = self.weight_function(universe) if self.weight_function else universe.GetWeight()
    #     # base = universe.GetWeight()
    #     # extra = self.weight_function(universe) if self.weight_function else 1.0
    #     # wgt = base * extra
    #     if isinstance(value[0],list):
    #         if not isinstance(wgt,list):
    #             wgt = len(value[0])*[wgt]
    #         #filling the multiple entries per event.
    #         for v in map(lambda wgt,*args: (list(args),wgt), wgt, *value):
    #             self.FillHist(universe,*v)
    #     else:
    #         self.FillHist(universe,value,wgt)


    # def FillHist(self, universe, value, wgt):
    #     global EN4_BAD_FILL_PRINT_COUNT

    #     if value.count(None) != 0:
    #         return

    #     name = self.histwrapper.name
    #     args = value + [wgt]

    #     # Only debug inclusive EN4.
    #     debug_en4 = (name == "EN4")

    #     if not debug_en4:
    #         self.histwrapper.FillUniverse(universe, *args)
    #         return

    #     h = self.histwrapper.hist

    #     try:
    #         x = float(value[0])
    #         b = h.FindBin(x)
    #     except Exception:
    #         b = None

    #     bad_input = any(_is_bad_number(v) for v in value) or _is_bad_number(wgt)

    #     before_c = h.GetBinContent(b) if b is not None else "NA"
    #     before_e = h.GetBinError(b) if b is not None else "NA"

    #     # Fill exactly once.
    #     self.histwrapper.FillUniverse(universe, *args)

    #     after_c = h.GetBinContent(b) if b is not None else "NA"
    #     after_e = h.GetBinError(b) if b is not None else "NA"

    #     made_main_bad = (
    #         b is not None
    #         and (
    #             not math.isfinite(after_c)
    #             or not math.isfinite(after_e)
    #         )
    #     )

    #     bad_band_entries = []

    #     if b is not None:
    #         for bandname in h.GetVertErrorBandNames():
    #             band = h.GetVertErrorBand(str(bandname))
    #             if not band:
    #                 continue

    #             # Band CV at touched bin.
    #             bc = band.GetBinContent(b)
    #             be = band.GetBinError(b)

    #             if not math.isfinite(bc) or not math.isfinite(be):
    #                 bad_band_entries.append(
    #                     ("BAND_CV", str(bandname), None, b, bc, be)
    #                 )

    #             # Universes at touched bin.
    #             for iu in range(band.GetNHists()):
    #                 hu = band.GetHist(iu)
    #                 uc = hu.GetBinContent(b)
    #                 ue = hu.GetBinError(b)

    #                 if not math.isfinite(uc) or not math.isfinite(ue):
    #                     bad_band_entries.append(
    #                         ("UNIVERSE", str(bandname), iu, b, uc, ue)
    #                     )

    #     if not (bad_input or made_main_bad or bad_band_entries):
    #         return

    #     # Print detailed standard-weight breakdown only once, for the first bad input weight.
    #     # This should tell us whether genie/flux/2p2h/rpa/lowq2pi/geant/minos/coh is nan.
    #     if bad_input and EN4_BAD_FILL_PRINT_COUNT == 0:
    #         try:
    #             universe.DebugStandardWeightPieces("EN4_BAD_WEIGHT")
    #         except Exception as e:
    #             print("[EN4_BAD_WEIGHT_PIECES_ERROR] {}".format(e), flush=True)

    #     # Stop flooding after N bad-fill reports.
    #     if EN4_BAD_FILL_PRINT_COUNT >= EN4_BAD_FILL_PRINT_LIMIT:
    #         return

    #     EN4_BAD_FILL_PRINT_COUNT += 1

    #     dbg_entry = _safe_existing_attr(universe, "_debug_entry")
    #     univ_name = str(universe)

    #     base = (
    #         "[EN4_BAD_FILL] "
    #         "n={n}/{limit} entry={entry} bin={bin} value={value} weight={weight} "
    #         "before=({bc},{be}) after=({ac},{ae}) bad_input={bad_input} "
    #         "universe='{univ}'"
    #     ).format(
    #         n=EN4_BAD_FILL_PRINT_COUNT,
    #         limit=EN4_BAD_FILL_PRINT_LIMIT,
    #         entry=dbg_entry,
    #         bin=b,
    #         value=value,
    #         weight=wgt,
    #         bc=before_c,
    #         be=before_e,
    #         ac=after_c,
    #         ae=after_e,
    #         bad_input=bad_input,
    #         univ=univ_name,
    #     )

    #     if made_main_bad:
    #         print(base + " where=MAIN_CV", flush=True)

    #     # Print only first few bad band entries for this fill.
    #     for kind, bandname, iu, ibin, c, e in bad_band_entries[:10]:
    #         if kind == "BAND_CV":
    #             print(
    #                 base
    #                 + " where=BAND_CV band={} content={} error={}".format(
    #                     bandname, c, e
    #                 ),
    #                 flush=True,
    #             )
    #         else:
    #             print(
    #                 base
    #                 + " where=UNIVERSE band={} univ={} content={} error={}".format(
    #                     bandname, iu, c, e
    #                 ),
    #                 flush=True,
    #             )

    #     if len(bad_band_entries) > 10:
    #         print(
    #             "[EN4_BAD_FILL_MORE] entry={} bin={} suppressed_extra_bad_entries={}".format(
    #                 dbg_entry, b, len(bad_band_entries) - 10
    #             ),
    #             flush=True,
    #         )
    def FillHist(self,universe,value,wgt):
        if value.count(None) == 0:
            args = value+[wgt]
            self.histwrapper.FillUniverse(universe, *args)
            # if self.histwrapper.name in ["nu_length", "drawnL_nu_length", "Eel"]:
            #     print("FILL", self.histwrapper.name, value, wgt)

    def Finalize(self):
        # if self.histwrapper.name in ["nu_length", "drawnL_nu_length", "Eel"]:
        #     h = self.histwrapper.hist
        #     print("HIST", h.GetName(), "nbins=", h.GetNbinsX())
        #     for i in range(0, h.GetNbinsX()+2):
        #         c = h.GetBinContent(i)
        #         if c != 0:
        #             print(i, c)
        self.histwrapper.Write()

    def AddCut(self, func):
        self.cuts.append(func)

    

def PlotProcessorProliferater(plot, name_func):
    plots = []
    for k,v  in name_func.items():
        plots.append(plot.Clone(VariantPlotsNamingScheme(plot.histwrapper.name,k)))
        plots[-1].AddCut(v)
    return plots

# complecated function for interpret tags of a plot. Got to revisit after sometime and find out how to organize multiple tags.
def MakePlotProcessors(**kwargs):
    # import inspect
    # print("Calling TranslateSettings from:", inspect.getfile(TranslateSettings))
    # print("TranslateSettings object:", TranslateSettings)
    # print("kwargs['key'] =", repr(kwargs["key"]), type(kwargs["key"]))
    settings = TranslateSettings(kwargs["key"])
    #print settings
    plots = []
    tags = settings["tags"]
    settings.pop("tags")

    
    if "suffix" in kwargs["key"]:
        settings["name"]+=kwargs["key"]["suffix"]

    if "mc_only" in tags and not kwargs["mc"]:
        return plots
    else:
        plots.append( PlotProcessor(**settings))

    if "ignore_selection" not in tags:
        name_func = {}
        if "sideband" in tags:
            for region in kwargs["region"]:
                name_func[VariantPlotsNamingScheme(region)] = MakePlotProcessors.func.setdefault(region,partial(lambda universe,_: universe.classifier.side_band == _, _=region))

        tmp = []
        for plot in plots:
            tmp.extend(PlotProcessorProliferater(plot,name_func))
            plot.AddCut(MakePlotProcessors.func.setdefault("Signal",lambda universe: universe.classifier.side_band =="Signal"))
        plots.extend(tmp)

    if "signal_only" in tags:
        if kwargs["mc"]:
            for plot in plots:
                plot.AddCut(lambda universe: universe.classifier.is_true_signal)
                plot.weight_function = lambda event:event.GetWeight(False)
        else:
            print("cant't make signal only plots for data")
            return plots

    if "truth_class" in tags and kwargs["mc"]:
        name_func = {}
        for cate in list(TRUTH_CATEGORIES.keys())+["Other"]:
            if cate in EXTRA_OTHER:
                continue
            name_func[cate] = MakePlotProcessors.func.setdefault(cate,partial(lambda universe,_: universe.classifier.truth_class == _, _=cate))


        tmp = []
        for plot in plots:
            tmp.extend(PlotProcessorProliferater(plot,name_func))
        plots.extend(tmp)

    return plots

MakePlotProcessors.func = {}

def MakeSelectionPlotsByCut(**kwargs):
    settings = TranslateSettings(kwargs["key"])
    plots = []
    cuts = [ partial(lambda event, cut : event.classifier.reco_cuts_passed[cut] , cut = _ ) for _ in SAMPLE_CUTS["Signal"]]
    cuts.extend( [ partial(lambda event, cut : event.classifier.reco_cuts_passed[cut] , cut = "Reco{}".format(_) ) for _ in KINEMATICS_CUTS] )

    for i in range(len(cuts)):
        local_settings = settings.copy()
        local_settings["name"] = "{}_{}{}".format(settings["name"],i,"cut")
        local_settings["cuts"] = cuts[:i+1]
        plots.append(PlotProcessor(**local_settings))
        local_settings["name"] = "{}_{}{}".format(settings["name"],i,"cut_signal")
        local_settings["cuts"] = [lambda event:event.classifier.is_true_signal]+cuts[:i+1]
        plots.append(PlotProcessor(**local_settings))

    return plots


def MakeCutVariablePlot(cutname,value_getter=None):
    plots = []
    settings={
        "name": cutname,
        "title":";{};NEvents".format(cutname),
    }

    print (binning,value_getter,name)

    if value_getter is None:
        settings["value_getter"]=[CUTS[cutname].Values]
    else :
        settings["value_getter"]=[value_getter]

    settings["binning"]=[CUTS[cutname]._variable_range]

    settings["cuts"] = [lambda event: all(event.classifier.reco_cuts_passed[_] for _ in SAMPLE_CUTS["Signal"] if _ != cutname)]
    settings["cuts"].extend([lambda event: all(event.classifier.reco_cuts_passed["Reco{}".format(_)] for _ in KINEMATICS_CUTS)])

    plots.append(PlotProcessor(**settings))
    settings["cuts"] = settings["cuts"]+[lambda event: event.classifier.is_true_signal]
    settings["name"] = "{}_signal".format(cutname)
    plots.append(PlotProcessor(**settings))

    return plots

