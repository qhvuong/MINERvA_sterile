import os
import copy
from collections import OrderedDict
import json
import ROOT
import PlotUtils
import numpy as np
from array import array

from itertools import compress
from config.SystematicsConfig import CONSOLIDATED_ERROR_GROUPS, ERROR_GROUPS_CONFIG
ccnueroot = os.environ.get('CCNUEROOT')

from tools.Helper import *
from config.AnalysisConfig import AnalysisConfig

minBinCont = 1
errorbandDict = { #keys are the errorbands that need to be renamed, values are what to rename them to
        "GENIE_D2_MaRES":"GENIE_MaRES",
        "GENIE_D2_NormCCRES":"GENIE_NormCCRES",
        "EnergyScale":"ElectronScale",
        "HCALScale":"elE_HCAL",
        "Other_Response":"response_other",
        "Pion_Response":"response_meson",
        "Proton_Response":"response_p",
        "Crosstalk":"response_xtalk",
        "Numu Coh Scale 1":"numucoh1",
        "Numu Coh Scale 2":"numucoh2",
        "Numu Coh Scale 3":"numucoh3",
        "Numu Coh Scale 4":"numucoh4",
        "Numu Coh Scale 5":"numucoh5",
        "Numu Coh Scale 6":"numucoh6",
        "Numu Scale":"Numu CC Scale",
        "RPAHigh":"RPA_HighQ2",
        "RPALow":"RPA_LowQ2",
        "Target_Mass":"Target_Mass_CH",
        }
sample_pdgs = ["numu","nue","anumu","anue"]

class StitchedHistogram:
    def __init__(self,name,dirty=False):
        self.name = name
        self.dirty = dirty

        self.keys = ["fhc_elastic","fhc_imd","fhc_numu_selection","fhc_nue_selection","rhc_elastic","rhc_imd","rhc_numu_selection","rhc_nue_selection"]
        self.titles = {"fhc_elastic":"FHC #nu+e","fhc_imd":"FHC IMD","fhc_numu_selection":"FHC CC #nu_{#mu}","fhc_nue_selection":"FHC CC #nu_{e}",
                "rhc_elastic":"RHC #nu+e","rhc_imd":"RHC IMD","rhc_numu_selection":"RHC CC anti #nu_{#mu}","rhc_nue_selection":"RHC CC anti #nu_{e}"}

        self.stitchKeys = self.keys.copy()

        self.data_hists = OrderedDict()
        self.mc_hists = OrderedDict()

        self.bin_dictionary = {}

        self.nue_hists = {}
        self.numu_hists = {}
        self.swap_hists = {}
        
        self.beam_ids = {}
        self.elastic_ids = {}
        self.ratio_ids = {}

        self.nueel_flavors = {}

        self.nue_templates = {}
        self.numu_templates = {}
        self.swap_templates = {}

        self.preservation_hist_container = {}
        self.preservation_hists = {}

        self.mc_flux_universes = []
        self.A = None # matrix of universes - MC CV

        self.inv_covariance = None
        self.inv_covariance_sans_flux = None
        self.covariance = None
        self.covariance_sans_flux = None
        self.external_covariances = {}

        self.mc_hist = None
        self.data_hist = None
        self.pseudo_hist = None

        self.nue_template = None
        self.numu_template = None
        self.swap_template = None
        
        self.nue_hist = None # nue flavor component
        self.numu_hist = None # numu flavor component
        self.swap_hist = None # flavor swapped numu -> nue sample

        self.osc_hist = None # histogram with an oscillation hypothesis

        self.beam_id = None # 1 for FHC, 0 for RHC
        self.elastic_id = None # 1 for scattering, 0 otherwise
        self.ratio_id = None # 1 for ratio, 0 otherwise

        self.is_processed = False

        self.n_flux_universes = 1000

        # self.use_1000_flux_universes = True

    def __del__(self):
        self.data_hists = {}
        self.mc_hists = {}
        self.nueel_flavors = {}
        self.nue_hists = {}
        self.numu_hists = {}
        self.swap_hists = {}
        self.numu_templates = {}
        self.nue_templates = {}
        self.swap_templates = {}
        self.mc_hist = None
        self.data_hist = None
        self.pseudo_hist = None

    def SetNFluxUniverses(self,n):
        self.n_flux_universes = n

    def AddTemplates(self,name,nue=None,numu=None,swap=None):
        if name not in self.keys:
            raise ValueError("{} does not correspond to a histogram in this object".format(name))
        l = [nue,numu,swap]
        filt = [temp is not None for temp in l]
        if any(filt):
            holder = list(compress(l,filt))[0].Clone()
        else:
            raise ValueError("No templates passed to AddTemplates")
        holder.Reset()
        self.nue_templates[name] = nue.Clone() if nue is not None else holder.Clone()
        self.numu_templates[name] = numu.Clone() if numu is not None else holder.Clone()
        self.swap_templates[name] = swap.Clone() if swap is not None else holder.Clone()

    def AddPreservationHists(self,sample,nuPdg,hist_dict):
        # Store preservation hists keyed by plotName:
        # hist_dict = {"plot name": hist}
        # preservation_hist_container["plot name"] = { "fhc_numu_selection": {"numu" : hist}}
        for plotName in hist_dict.keys():
            if plotName in self.preservation_hist_container.keys():
                self.preservation_hist_container[plotName][sample] = {nuPdg:hist_dict[plotName].Clone()}
            else:
                self.preservation_hist_container[plotName] = {sample:{nuPdg:hist_dict[plotName].Clone()}}
                
    def AddScatteringPreservationHists(self,sample,hist_dict):
        # Store preservation hists keyed by plotName:
        # hist_dict = {"plot name": {"numu":hist1,"nue":hist2,"anumu":hist3,"anue":hist4}
        # preservation_hist_container["plot name"] = { "fhc_elastic": {"numu" : hist1,...}}
        for plotName in hist_dict.keys():
            if plotName in self.preservation_hist_container.keys():
                self.preservation_hist_container[plotName][sample] = hist_dict[plotName]
            else:
                self.preservation_hist_container[plotName] = {sample:hist_dict[plotName]}

    def AddScatteringFlavors(self,name,hist):
        if name in self.nueel_flavors.keys():
            print("{} has already been added to scattering sample dictionary. Doing nothing.".format(name))
        else:
            self.nueel_flavors[name] = hist.Clone()

    # def AddSwappedSample(self, name, hist):
    #     if name in self.swap_hists.keys():
    #         print("{} has already been added to swapped sample dictionary. Doing nothing.".format(name))
    #         return

    #     beam = name[:3]
    #     numu_key = beam + "_numu_selection"

    #     if numu_key not in self.numu_hists:
    #         raise KeyError(
    #             "Cannot add swapped sample {} before {} exists in self.numu_hists. "
    #             "Call AddHistograms('{}', ...) before AddSwappedSample('{}', ...)."
    #             .format(name, numu_key, numu_key, name)
    #         )

    #     scaleIntegral = self.numu_hists[numu_key].Integral()
    #     swapIntegral = hist.Integral()

    #     if swapIntegral == 0 or not np.isfinite(swapIntegral):
    #         raise ValueError(
    #             "Bad swapped-sample integral for {}: {}".format(name, swapIntegral)
    #         )

    #     if not np.isfinite(scaleIntegral):
    #         raise ValueError(
    #             "Bad numu scale integral for {}: {}".format(numu_key, scaleIntegral)
    #         )

    #     print("\n===== AddSwappedSample diagnostic =====")
    #     print("sample                 =", name)
    #     print("numu reference         =", numu_key)
    #     print("numu integral          =", scaleIntegral)
    #     print("raw swap integral      =", swapIntegral)
    #     print("applied scale          =", scaleIntegral / swapIntegral)

    #     for i in range(1, min(hist.GetNbinsX(), 12) + 1):
    #         print(
    #             "bin {:2d}: numu={:12.6g} rawSwap={:12.6g} scaledSwap={:12.6g}".format(
    #                 i,
    #                 self.numu_hists[numu_key].GetBinContent(i),
    #                 hist.GetBinContent(i),
    #                 hist.GetBinContent(i) * scaleIntegral / swapIntegral,
    #             )
    #         )

    #     h = hist.Clone()
    #     h.Scale(scaleIntegral / swapIntegral)
    #     self.swap_hists[name] = h
    def AddSwappedSample(self, name, hist):
        if name in self.swap_hists:
            print(
                "{} has already been added to swapped sample dictionary. "
                "Doing nothing.".format(name)
            )
            return

        swap_integral = hist.Integral()

        if swap_integral == 0 or not np.isfinite(swap_integral):
            raise ValueError(
                "Bad swapped-sample integral for {}: {}".format(
                    name,
                    swap_integral,
                )
            )

        # The swapped sample has already been normalized to the analysis POT
        # in stitch.py. Preserve that normalization.
        h = hist.Clone()
        h.SetDirectory(0)

        self.swap_hists[name] = h

        print("\n===== AddSwappedSample =====")
        print("sample             =", name)
        print("stored swap integral =", h.Integral())
        print("additional scale     = none")
    # def AddSwappedSample(self,name,hist):
    #     if name in self.swap_hists.keys():
    #         print("{} has already been added to swapped sample dictionary. Doing nothing.".format(name))
    #     else:
    #         beam = name[:3]
    #         scaleIntegral = self.numu_hists[beam+"_numu_selection"].Integral()
    #         swapIntegral = hist.Integral()
    #         hist.Scale(scaleIntegral/swapIntegral)
    #         self.swap_hists[name] = hist.Clone()

    def SetDataHistogram(self,hist):
        self.data_hist = hist.Clone()

    def SetPseudoHistogram(self,hist):
        self.pseudo_hist = hist.Clone()

    def SetMCHistogram(self,hist):
        self.mc_hist = hist.Clone()
        self.SetAMatrix()

    def SetAMatrix(self):
        ##### Store flux universes for marginalization procedure #####
        if "Flux" in self.mc_hist.GetVertErrorBandNames():
            band = self.mc_hist.GetVertErrorBand("Flux")
            nhists = band.GetNHists()
            self.mc_flux_universes = np.array([np.array(band.GetHist(i))[1:-1] for i in range(nhists)])

            cv = np.array(self.mc_hist)[1:-1]
            self.A = self.mc_flux_universes - np.array([cv for i in range(nhists)])

    def LoadExternalCovariances(self, root_file):
        """
        Load any external covariance matrices stored in the stitched ROOT file.

        For now this is used for the published MINERvA nu+e elastic table.
        """
        cov = root_file.Get("fhc_elastic_paper_covariance_matrix")

        if cov:
            cov = cov.Clone("fhc_elastic_paper_covariance_matrix_loaded")
            self.external_covariances["fhc_elastic"] = cov

            print("Loaded external covariance: fhc_elastic")
            print("  Nrows =", cov.GetNrows())
            print("  Ncols =", cov.GetNcols())
        # else:
        #     print("No external covariance matrix found for fhc_elastic.")

    def GetExternalCovarianceBinIndices(self, sample_name):
        """
        Return numpy indices for the bins of a sample inside the stitched histogram.

        Uses HIST_CONFIG.json, which is written by TranslateBins().
        The JSON stores zero-based start/end indices, matching numpy indexing.
        """
        try:
            inds = GetSliceIndices("HIST_CONFIG.json", "", [sample_name])
        except Exception as e:
            print("WARNING: could not read HIST_CONFIG.json for external covariance:", e)
            inds = []

        return inds

    def ApplyExternalCovariances(self, covariance, label="covariance"):
        """
        Replace selected blocks of the stitched covariance matrix with
        externally supplied covariance matrices.

        Important: this REPLACES the existing block, not adds to it,
        to avoid double-counting the paper uncertainties.
        """
        cov_out = np.array(covariance, copy=True)

        for sample_name, root_cov in self.external_covariances.items():
            inds = self.GetExternalCovarianceBinIndices(sample_name)

            if len(inds) == 0:
                print("WARNING: no stitched bins found for external covariance sample:", sample_name)
                continue

            ext_cov = TMatrix_to_Numpy(root_cov)

            if ext_cov.shape[0] != len(inds) or ext_cov.shape[1] != len(inds):
                raise RuntimeError(
                    "External covariance shape {} does not match stitched bins {} for sample {}".format(
                        ext_cov.shape, len(inds), sample_name
                    )
                )

            print("\nApplying external covariance to", label)
            print("  sample      =", sample_name)
            print("  global bins =", [i + 1 for i in inds])
            print("  shape       =", ext_cov.shape)

            for a, ia in enumerate(inds):
                for b, ib in enumerate(inds):
                    cov_out[ia, ib] = ext_cov[a, b]

            print("  external covariance block:")
            print(ext_cov)

        return cov_out

    def ApplyExternalCovarianceErrorsToHist(self, hist, fractional=False):
        """
        For plotting only.

        Overwrite the displayed bin errors for bins covered by external
        covariance matrices using sqrt(diag(cov)).

        If fractional=True, set bin errors to sqrt(cov_ii)/CV_i and
        bin contents to 1, suitable for Data/MC ratio error bands.
        """
        h_out = hist.Clone(hist.GetName() + "_external_cov_errors")
        h_out.SetDirectory(0)

        for sample_name, root_cov in self.external_covariances.items():
            inds = self.GetExternalCovarianceBinIndices(sample_name)
            if len(inds) == 0:
                continue

            cov = TMatrix_to_Numpy(root_cov)

            if cov.shape[0] != len(inds):
                raise RuntimeError(
                    "External covariance shape {} does not match {} bins for {}".format(
                        cov.shape, len(inds), sample_name
                    )
                )

            for a, idx0 in enumerate(inds):
                # numpy idx0 is zero-based, ROOT bin is one-based
                ibin = idx0 + 1
                sigma = np.sqrt(max(cov[a, a], 0.0))

                if fractional:
                    cv = self.mc_hist.GetBinContent(ibin)
                    frac = sigma / cv if cv != 0 else 0.0
                    h_out.SetBinContent(ibin, 1.0)
                    h_out.SetBinError(ibin, frac)
                else:
                    h_out.SetBinError(ibin, sigma)

        return h_out


    def SetCovarianceMatrices(self):
        if type(self.mc_hist) == type(self.data_hist) and type(self.mc_hist) == type(self.pseudo_hist) and type(self.mc_hist) == PlotUtils.MnvH1D:
            h_test = self.data_hist.Clone()
            h_test.Add(self.mc_hist, -1)

            covariance = TMatrix_to_Numpy(h_test.GetTotalErrorMatrix(True, False, False))[1:-1, 1:-1]
            flux_covariance = TMatrix_to_Numpy(h_test.GetSysErrorMatrix("Flux"))[1:-1, 1:-1]
            cov_sans_flux = covariance - flux_covariance

            # Diagnostic: check stat/sys decomposition before external replacement
            stat = TMatrix_to_Numpy(h_test.GetStatErrorMatrix())[1:-1, 1:-1]

            stat = TMatrix_to_Numpy(h_test.GetStatErrorMatrix())[1:-1, 1:-1]

            # print("\n===== covariance stat check before external replacement =====")
            # print("trace stat      =", np.trace(stat))
            # print("trace total     =", np.trace(covariance))
            # print("trace sansFlux  =", np.trace(cov_sans_flux))
            # print("trace flux      =", np.trace(flux_covariance))
            # print("trace total-stat =", np.trace(covariance - stat))
            # print("trace sans-stat  =", np.trace(cov_sans_flux - stat))
            # print("trace full - sans - flux =", np.trace(covariance - cov_sans_flux - flux_covariance))
            # print("max |full-sans-flux| =", np.max(np.abs(covariance - cov_sans_flux - flux_covariance)))

            # Replace external covariance blocks, e.g. published FHC nu+e table covariance.
            covariance = self.ApplyExternalCovariances(covariance, label="full covariance")
            cov_sans_flux = self.ApplyExternalCovariances(cov_sans_flux, label="sans-flux covariance")

            self.covariance = covariance
            self.covariance_sans_flux = cov_sans_flux

            self.inv_covariance = np.linalg.inv(covariance)
            self.inv_covariance_sans_flux = np.linalg.inv(cov_sans_flux)
        else:
            raise ValueError("MC and Data histograms must be defined before setting inv_covariance matrix")
    # def SetCovarianceMatrices(self):
    #     if type(self.mc_hist) == type(self.data_hist) and type(self.mc_hist) == type(self.pseudo_hist) and type(self.mc_hist) == PlotUtils.MnvH1D:
    #         h_test = self.data_hist.Clone()
    #         h_test.Add(self.mc_hist,-1)

    #         covariance = TMatrix_to_Numpy(h_test.GetTotalErrorMatrix(True,False,False))[1:-1,1:-1]
    #         flux_covariance = TMatrix_to_Numpy(h_test.GetSysErrorMatrix("Flux"))[1:-1,1:-1]
    #         cov_sans_flux = covariance - flux_covariance

    #         # Replace external covariance blocks, e.g. published FHC nu+e table covariance.
    #         covariance = self.ApplyExternalCovariances(covariance, label="full covariance")
    #         cov_sans_flux = self.ApplyExternalCovariances(cov_sans_flux, label="sans-flux covariance")

    #         self.covariance = covariance
    #         self.covariance_sans_flux = cov_sans_flux

    #         self.inv_covariance = np.linalg.inv(covariance)
    #         self.inv_covariance_sans_flux = np.linalg.inv(cov_sans_flux)
    #         # self.covariance = covariance
    #         # self.covariance_sans_flux = cov_sans_flux

    #         # self.inv_covariance = np.linalg.inv(covariance)
    #         # self.inv_covariance_sans_flux = np.linalg.inv(cov_sans_flux)
    #     else:
    #         raise ValueError("MC and Data histograms must be defined before setting inv_covariance matrix")

    def GetInverseCovarianceMatrix(self,sansFlux=False):
        if sansFlux:
            return(np.copy(self.inv_covariance_sans_flux))
        else:
            return(np.copy(self.inv_covariance))

    def GetCovarianceMatrix(self,sansFlux=False):
        if sansFlux:
            return(np.copy(self.covariance_sans_flux))
        else:
            return(np.copy(self.covariance))

    def SetOscHistogram(self,hist):
        self.osc_hist = hist.Clone()

    def GetMCHistogram(self):
        return(self.mc_hist.Clone())

    def GetDataHistogram(self):
        return(self.data_hist.Clone())

    def GetPseudoHistogram(self):
        return(self.pseudo_hist.Clone())

    def GetOscillatedHistogram(self):
        if self.osc_hist is not None:
            return(self.osc_hist.Clone())
        else:
            raise ValueError("Oscillation histogram not set.")

    def GetFluxUniverses(self):
        return(np.copy(self.mc_flux_universes))

    def GetAMatrix(self):
        return(np.copy(self.A))

    def EmptyHist(self,h):
        h_ret = h.Clone()
        if type(h_ret) == ROOT.TH1D or type(h_ret) == PlotUtils.MnvH1D:
            for i in range(0,h.GetNbinsX()+1):
                h_ret.SetBinContent(i,0)
                h_ret.SetBinError(i,0)
                for name in h_ret.GetVertErrorBandNames():
                    h_ret.GetVertErrorBand(name).SetBinContent(i,0)
                    for univ in range(h_ret.GetVertErrorBand(name).GetNHists()):
                        h_ret.GetVertErrorBand(name).GetHist(univ).SetBinContent(i,0)
        elif type(h_ret) == ROOT.TH2D or type(h_ret) == PlotUtils.MnvH2D:
            for i in range(0,h.GetNbinsX()+1):
                for j in range(0,h.GetNbinsY()+1):
                    h_ret.SetBinContent(i,j,0)
                    h_ret.SetBinError(i,j,0)
                    for name in h_ret.GetVertErrorBandNames():
                        h_ret.GetVertErrorBand(name).SetBinContent(i,j,0)
                        for univ in range(h_ret.GetVertErrorBand(name).GetNHists()):
                            h_ret.GetVertErrorBand(name).GetHist(univ).SetBinContent(i,j,0)

        return(h_ret)

    def UnityHist(self,h):
        h_ret = h.Clone()
        for i in range(0,h.GetNbinsX()+1):
            h_ret.SetBinContent(i,1)
            h_ret.SetBinError(i,0)
            for name in h_ret.GetVertErrorBandNames():
                h_ret.GetVertErrorBand(name).SetBinContent(i,1)
                for univ in range(h_ret.GetVertErrorBand(name).GetNHists()):
                    h_ret.GetVertErrorBand(name).GetHist(univ).SetBinContent(i,1)
        return(h_ret)

    def AddHistograms(self,name,mc_hist,data_hist):
        if len(self.data_hists.keys()) > 0 and type(data_hist) != type(list(self.data_hists.values())[0]):
            raise ValueError("Cannot add {} to data histogram dictionary of {}".format(type(data_hist),type(list(self.data_hists.values())[0])))
        if len(self.mc_hists.keys()) > 0 and type(mc_hist) != type(list(self.mc_hists.values())[0]):
            raise ValueError("Cannot add {} to mc histogram dictionary of {}".format(type(mc_hist),type(list(self.mc_hists.values())[0])))

        self.data_hists[name] = data_hist.Clone()
        self.mc_hists[name] = mc_hist.Clone()

        if not self.dirty:
            beam = name[:4]
            # ----- Set Electron Neutrino Histograms ----- #
            if "nue_selection" in name:
                self.nue_hists[name] = mc_hist.Clone()
            elif "elastic" in name:
                if 'electron_'+name not in self.nueel_flavors.keys():
                    raise ValueError("{} has not been added to nueel flavor dictionary. Do this before AddingHistogram()".format('electron_'+name))
                self.nue_hists[name] = self.nueel_flavors['electron_'+name].Clone()
            elif "ratio" in name:
                if beam+'nue_selection' in self.nue_hists.keys():
                    self.nue_hists[name] = self.nue_hists[beam+'nue_selection'].Clone()
            else:
                self.nue_hists[name] = self.EmptyHist(mc_hist)

            # ----- Set Muon Neutrino Histograms ----- #
            if "numu_selection" in name or "imd" in name:
                self.numu_hists[name] = mc_hist.Clone()
            elif "elastic" in name:
                if 'muon_'+name not in self.nueel_flavors.keys():
                    raise ValueError("{} has not been added to nueel flavor dictionary. Do this before AddingHistogram()".format('muon_'+name))
                self.numu_hists[name] = self.nueel_flavors['muon_'+name].Clone()
            elif "ratio" in name:
                if beam+'numu_selection' in self.numu_hists.keys():
                    self.numu_hists[name] = self.numu_hists[beam+'numu_selection'].Clone()
            else:
                self.numu_hists[name] = self.EmptyHist(mc_hist)

            # ----- Set Flavor Swapped Neutrino Histograms ----- #
            if "nue_selection" in name:
                if name in self.swap_hists.keys():
                    self.swap_hists[name] = self.swap_hists[name]
            elif "ratio" in name:
                if beam+'nue_selection' in self.swap_hists.keys():
                    self.swap_hists[name] = self.swap_hists[beam+'nue_selection'].Clone()
            elif "elastic" in name:
                self.swap_hists[name] = self.nueel_flavors['muon_'+name].Clone()
            else:
                self.swap_hists[name] = self.EmptyHist(mc_hist)

            # ----- Set ID Histograms ----- #
            if "elastic" in name or "imd" in name:
                self.elastic_ids[name] = self.UnityHist(mc_hist)
            else:
                self.elastic_ids[name] = self.EmptyHist(mc_hist)
            if "ratio" in name:
                self.ratio_ids[name] = self.UnityHist(mc_hist)
            else:
                self.ratio_ids[name] = self.EmptyHist(mc_hist)
            if "fhc" in name:
                self.beam_ids[name] = self.UnityHist(mc_hist)
            else:
                self.beam_ids[name] = self.EmptyHist(mc_hist)

        self.UpdateKeys()
        self.is_processed = False

    def RemoveHistograms(self,name):
        if name not in list(self.data_hists.keys()):
            raise ValueError("{} not in data histogram dictionary".format(name))
        if name not in list(self.mc_hists.keys()):
            raise ValueError("{} not in mc histogram dictionary".format(name))

        self.stitchKeys.remove(name)

    def UpdateKeys(self):
        self.keys = list(self.mc_hists.keys())
        self.stitchKeys = self.keys.copy()
        if self.keys != list(self.data_hists.keys()):
            raise ValueError("MC dictionary incompatable with data dictionary")

        if not self.dirty:
            if self.keys != list(self.nue_hists.keys()):
                raise ValueError("MC dictionary incompatable with nue hists")
            if self.keys != list(self.numu_hists.keys()):
                raise ValueError("MC dictionary incompatable with numu hists")

    def CheckRatioFluxUniverses(self, beam, ratio_hist, numu_hist, nue_hist, bins_to_check=[1, 9, 10], universes_to_check=[0, 1, 2, 10]):
        print("\n===== CHECK RATIO FLUX UNIVERSES: {} =====".format(beam))

        if not ratio_hist.HasVertErrorBand("Flux"):
            print("ratio_hist has no Flux band")
            return
        if not numu_hist.HasVertErrorBand("Flux"):
            print("numu_hist has no Flux band")
            return
        if not nue_hist.HasVertErrorBand("Flux"):
            print("nue_hist has no Flux band")
            return

        ratio_band = ratio_hist.GetVertErrorBand("Flux")
        numu_band = numu_hist.GetVertErrorBand("Flux")
        nue_band = nue_hist.GetVertErrorBand("Flux")

        print("NHists ratio =", ratio_band.GetNHists())
        print("NHists numu  =", numu_band.GetNHists())
        print("NHists nue   =", nue_band.GetNHists())

        for b in bins_to_check:
            cv_ratio_hist = ratio_hist.GetBinContent(b)
            cv_numu = numu_hist.GetBinContent(b)
            cv_nue = nue_hist.GetBinContent(b)
            cv_ratio_manual = cv_numu / cv_nue if cv_nue != 0 else 0.0

            print("\nBin", b)
            print("  CV ratio hist   =", cv_ratio_hist)
            print("  CV numu/nue     =", cv_ratio_manual)
            print("  CV diff         =", cv_ratio_hist - cv_ratio_manual)

            for u in universes_to_check:
                if u >= ratio_band.GetNHists():
                    continue

                ru = ratio_band.GetHist(u).GetBinContent(b)
                mu = numu_band.GetHist(u).GetBinContent(b)
                eu = nue_band.GetHist(u).GetBinContent(b)
                manual = mu / eu if eu != 0 else 0.0

                print(
                    "  univ {:3d}: ratio_hist={:12.6g}  manual={:12.6g}  diff={:12.6g}  A={:12.6g}".format(
                        u,
                        ru,
                        manual,
                        ru - manual,
                        ru - cv_ratio_hist
                    )
                )

    def MakeRatio(self, beam):
        beam = beam.lower()
        numuname = beam + "_numu_selection"
        elecname = beam + "_nue_selection"

        if numuname in self.keys and elecname in self.keys:
            mc_clone = self.mc_hists[numuname].Clone()
            mc_clone.Divide(mc_clone, self.mc_hists[elecname])

            data_clone = self.data_hists[numuname].Clone()
            data_clone.Divide(data_clone, self.data_hists[elecname])

            self.CheckRatioFluxUniverses(
                beam,
                mc_clone,
                self.mc_hists[numuname],
                self.mc_hists[elecname],
                bins_to_check=[1, 9, 10],
                universes_to_check=[0, 1, 2, 10]
            )

            self.AddHistograms("{}_ratio".format(beam), mc_clone, data_clone)

            if beam + "_numu_selection" in self.numu_templates:
                self.AddTemplates(
                    "{}_ratio".format(beam),
                    numu=self.numu_templates[beam + "_numu_selection"].Clone(),
                    nue=self.nue_templates[beam + "_nue_selection"].Clone(),
                    swap=self.swap_templates[beam + "_nue_selection"].Clone()
                )

                self.titles["{}_ratio".format(beam)] = "%s CC #nu_{#mu}/#nu_{e} Ratio" % (
                    beam.upper()
                )

        else:
            print("WARNING: cannot make {} ratio.".format(beam))
            print("  missing numu sample?", numuname not in self.keys)
            print("  missing nue sample? ", elecname not in self.keys)

        self.CleanErrorBands()

    def Copy(self):
        return(copy.deepcopy(self))

    def ApplyExclusion(self,exclude):
        for h in self.keys:
            if "fhc" in exclude:
                if "fhc" in h and "selection" in h:
                    self.RemoveHistograms(h)
            if "rhc" in exclude:
                if "rhc" in h and "selection" in h:
                    self.RemoveHistograms(h)
            if "numu" in exclude and "numu" in h:
                self.RemoveHistograms(h)
            if "nue" in exclude and "nue" in h:
                self.RemoveHistograms(h)
            if "elastic" in exclude and "elastic" in h:
                self.RemoveHistograms(h)
            if "imd" in exclude:
                self.RemoveHistograms(h)

    def CleanErrorBands(self,names=[]):
        self.LateralToVertical()
        self.SeparateBeamAngle()
        self.SyncErrorBands()

        for errname in names:
            for h in self.keys:
                if errname in self.data_hists[h].GetVertErrorBandNames():
                    self.data_hists[h].PopVertErrorBand(errname)
                if errname in self.mc_hists[h].GetVertErrorBandNames():
                    self.mc_hists[h].PopVertErrorBand(errname)
                if errname in self.nue_hists[h].GetVertErrorBandNames():
                    self.nue_hists[h].PopVertErrorBand(errname)
                if errname in self.numu_hists[h].GetVertErrorBandNames():
                    self.numu_hists[h].PopVertErrorBand(errname)
                if errname in self.swap_hists[h].GetVertErrorBandNames():
                    self.swap_hists[h].PopVertErrorBand(errname)

        if "fhc_elastic" in self.keys and "rhc_elastic" in self.keys:
            for i in range(self.data_hists["fhc_elastic"].GetNbinsX()+1):
                h_temp = self.data_hists['rhc_elastic']
                h_cont = h_temp.GetBinContent(i)
                h_cont = h_cont if h_cont != 0 else 1
                ratio1 = h_temp.GetVertErrorBand("ElectronScale").GetHist(0).GetBinContent(i)/h_cont
                ratio2 = h_temp.GetVertErrorBand("ElectronScale").GetHist(1).GetBinContent(i)/h_cont
                h_fix = self.data_hists['fhc_elastic']
                h_fix.GetVertErrorBand("ElectronScale").GetHist(0).SetBinContent(i,h_fix.GetBinContent(i) * ratio1)
                h_fix.GetVertErrorBand("ElectronScale").GetHist(1).SetBinContent(i,h_fix.GetBinContent(i) * ratio2)

    def SeparateBeamAngle(self):
        def separate(hist):
            if "beam_angle" in hist.GetVertErrorBandNames():
                univ_hist = hist.GetVertErrorBand("beam_angle").GetHists()
                hist.PopVertErrorBand("beam_angle")

                if "BeamAngleX" not in hist.GetVertErrorBandNames():
                    xuniv_hist = univ_hist[:2]
                    hist.AddVertErrorBand("BeamAngleX",xuniv_hist)
                if "BeamAngleY" not in hist.GetVertErrorBandNames():
                    yuniv_hist = univ_hist[2:]
                    hist.AddVertErrorBand("BeamAngleY",yuniv_hist)

        for h in self.keys:
            separate(self.data_hists[h])
            separate(self.mc_hists[h])
            if not self.dirty:
                separate(self.nue_hists[h])
                separate(self.numu_hists[h])
                separate(self.swap_hists[h])

    def LateralToVertical(self):
        def lateral(hist):
            for name in hist.GetLatErrorBandNames():
                universe_hists = hist.GetLatErrorBand(name).GetHists()
                hist.PopLatErrorBand(name)
                if name not in hist.GetVertErrorBandNames(): 
                    hist.AddVertErrorBand(name,universe_hists)

        for h in self.keys:
            lateral(self.data_hists[h])
            lateral(self.mc_hists[h])
            if not self.dirty:
                lateral(self.nue_hists[h])
                lateral(self.numu_hists[h])
                lateral(self.swap_hists[h])

    def TranslateBins(self):
        histogram_config = {}
        n_bins_new = 0

        for h in self.stitchKeys:
            histogram_config[h] = {"start" : n_bins_new}
            for i in range(1,self.mc_hists[h].GetNbinsX()+1): # skip under and overflow bins
                if self.mc_hists[h].GetBinContent(i) > minBinCont:
                    #print(h,i,self.mc_hists[h].GetBinLowEdge(i),self.mc_hists[h].GetBinLowEdge(i)+self.mc_hists[h].GetBinWidth(i))
                    n_bins_new+=1
                    bin_width = self.mc_hists[h].GetBinWidth(i)
                    if "elastic" in h:
                        bin_norm = 2
                    else:
                        bin_norm = self.mc_hists[h].GetBinWidth(i)

                    bin_center = self.mc_hists[h].GetBinCenter(i)
                    self.bin_dictionary[n_bins_new] = {"sample":h,"bin":i}

            histogram_config[h]["end"] = n_bins_new-1

        with open("HIST_CONFIG.json","w") as file:
            json.dump(histogram_config,file,indent=4)

        return(n_bins_new)




    def GetCommonLEBinEdges(self):
        reference = None
        reference_name = None

        template_groups = [
            ("nue", self.nue_templates),
            ("numu", self.numu_templates),
            ("swap", self.swap_templates),
        ]

        for component_name, template_dict in template_groups:
            for sample_name, template in template_dict.items():
                if template is None:
                    continue

                # Unstitched templates use:
                # x = L/E
                # y = reconstructed bin.
                xaxis = template.GetXaxis()
                nbins = template.GetNbinsX()

                edges = [
                    float(xaxis.GetBinLowEdge(i))
                    for i in range(1, nbins + 1)
                ]
                edges.append(float(xaxis.GetBinUpEdge(nbins)))

                edges = np.asarray(edges, dtype=float)

                if reference is None:
                    reference = edges
                    reference_name = "{}:{}".format(
                        component_name,
                        sample_name,
                    )
                    continue

                if len(edges) != len(reference) or not np.allclose(
                    edges,
                    reference,
                    rtol=0.0,
                    atol=1e-12,
                ):
                    raise RuntimeError(
                        "Inconsistent L/E template axes: {} differs from {}".format(
                            "{}:{}".format(component_name, sample_name),
                            reference_name,
                        )
                    )

        if reference is None:
            raise RuntimeError("No L/E templates were found")

        return reference




    def Stitch(self):
        # ----- Create empty ROOT histograms to fill with stitched content ----- #
        n_bins_new = self.TranslateBins()

        self.mc_hist   = ROOT.TH1D(self.name+"_mc",self.name+"_mc",n_bins_new,0,n_bins_new)
        self.data_hist = ROOT.TH1D(self.name+"_data",self.name+"_data",n_bins_new,0,n_bins_new)
        self.pseudo_hist = ROOT.TH1D(self.name+"_pseudo",self.name+"_pseudo",n_bins_new,0,n_bins_new)
        
        self.nue_hist = ROOT.TH1D(self.name+"_nue",self.name+"_nue",n_bins_new,0,n_bins_new)
        self.numu_hist = ROOT.TH1D(self.name+"_numu",self.name+"_numu",n_bins_new,0,n_bins_new)
        self.swap_hist = ROOT.TH1D(self.name+"_swap",self.name+"_swap",n_bins_new,0,n_bins_new)

        self.beam_id = ROOT.TH1D(self.name+"_beamID",self.name+"_beamID",n_bins_new,0,n_bins_new)
        self.elastic_id = ROOT.TH1D(self.name+"_nueelID",self.name+"_nueelID",n_bins_new,0,n_bins_new)
        self.ratio_id = ROOT.TH1D(self.name+"_ratioID",self.name+"_ratioID",n_bins_new,0,n_bins_new)

        # self.nue_template  = ROOT.TH2D(self.name+"_nue_template",self.name+"_nue_template",n_bins_new,0,n_bins_new,34,0,0.495)
        # self.numu_template  = ROOT.TH2D(self.name+"_numu_template",self.name+"_numu_template",n_bins_new,0,n_bins_new,34,0,0.495)
        # self.swap_template  = ROOT.TH2D(self.name+"_swap_template",self.name+"_swap_template",n_bins_new,0,n_bins_new,34,0,0.495)
        le_edges = self.GetCommonLEBinEdges()
        le_edges_root = array("d", le_edges.tolist())
        n_le_bins = len(le_edges) - 1

        self.nue_template = ROOT.TH2D(
            self.name + "_nue_template",
            self.name + "_nue_template",
            n_bins_new,
            0,
            n_bins_new,
            n_le_bins,
            le_edges_root,
        )

        self.numu_template = ROOT.TH2D(
            self.name + "_numu_template",
            self.name + "_numu_template",
            n_bins_new,
            0,
            n_bins_new,
            n_le_bins,
            le_edges_root,
        )

        self.swap_template = ROOT.TH2D(
            self.name + "_swap_template",
            self.name + "_swap_template",
            n_bins_new,
            0,
            n_bins_new,
            n_le_bins,
            le_edges_root,
        )


        # ----- Combine samples to one histogram ----- #
        print("Filling CV content in stitched histogram...")
        self.StitchThis()   # combine 1D histograms
        if len(list(self.nue_templates.keys())) > 0:
            self.StitchThis2D() # combine 2D templates

        if len(list(self.preservation_hist_container.keys())) > 0:
            # Construct one empty stitched hist for each plot name in preservation hists
            for plotName in self.preservation_hist_container.keys():
                histDict = self.preservation_hist_container[plotName]
                self.preservation_hists[plotName] = {}
                for sample in histDict.keys():
                    refKey = list(histDict[sample].keys())[0]
                    refHist = histDict[sample][refKey]
                    if type(refHist) == PlotUtils.MnvH2D:
                        nBinsY = refHist.GetNbinsY()
                        yLow = refHist.GetYaxis().GetBinLowEdge(1)
                        yUp = refHist.GetYaxis().GetBinLowEdge(nBinsY + 1)
                        stitched_hist = ROOT.TH2D(plotName,plotName,n_bins_new,0,n_bins_new,nBinsY,yLow,yUp)
                        for pdg in ["numu","nue","anumu","anue"]:
                            self.preservation_hists[plotName][pdg] = stitched_hist.Clone()
                        break
                    elif type(refHist) == PlotUtils.MnvH1D:
                        stitched_hist = ROOT.TH1D(plotName,plotName,n_bins_new,0,n_bins_new)
                        for pdg in ["numu","nue","anumu","anue"]:
                            self.preservation_hists[plotName][pdg] = stitched_hist.Clone()
                        break

            self.StitchPreservationHists() # combine histograms used for data preservation

        # ----- Do some errorband cleaning ----- #
        if not self.is_processed:
            print("Processing error systematics...")
            self.SeparateBeamAngle() # make beam angle systematics consistent across samples
            self.LateralToVertical() # convert lateral errorbands (deprecated) from old samples to vertical
            self.SyncErrorBands()  # make sure all samples have the same errorbands, reduce flux universes to 100
            self.is_processed = True

        # ----- Convert ROOT.TH1 to PlotUtils.Mnv1D and fill errorbands ----- #
        self.mc_hist = PlotUtils.MnvH1D(self.mc_hist)
        self.data_hist = PlotUtils.MnvH1D(self.data_hist)
        self.pseudo_hist = PlotUtils.MnvH1D(self.pseudo_hist)

        self.nue_hist = PlotUtils.MnvH1D(self.nue_hist)
        self.numu_hist = PlotUtils.MnvH1D(self.numu_hist)
        self.swap_hist = PlotUtils.MnvH1D(self.swap_hist)

        print("Filling error bands in stitched histogram...")
        self.mc_hist.AddMissingErrorBandsAndFillWithCV(self.mc_hists[self.stitchKeys[0]])
        self.data_hist.AddMissingErrorBandsAndFillWithCV(self.mc_hists[self.stitchKeys[0]])
        self.pseudo_hist.AddMissingErrorBandsAndFillWithCV(self.mc_hists[self.stitchKeys[0]])

        self.nue_hist.AddMissingErrorBandsAndFillWithCV(self.mc_hists[self.stitchKeys[0]])
        self.numu_hist.AddMissingErrorBandsAndFillWithCV(self.mc_hists[self.stitchKeys[0]])
        self.swap_hist.AddMissingErrorBandsAndFillWithCV(self.mc_hists[self.stitchKeys[0]])

        self.FillErrorBandsFromDict()
        self.SetCovarianceMatrices()
        self.SetAMatrix()

        self.SetPlottingStyle()
        if getattr(AnalysisConfig, "debug", False):
            self.DebugPlots()
        # if AnalysisConfig.debug:
        #     self.DebugPlots()

    def Write(self,filename):
        print("Writing histograms to {}...".format(filename))
        f = ROOT.TFile(filename,"RECREATE")
        self.mc_hist.Write()
        self.data_hist.Write()
        self.pseudo_hist.Write()

        self.nue_hist.Write()
        self.numu_hist.Write()
        self.swap_hist.Write()

        self.beam_id.Write()
        self.ratio_id.Write()
        self.elastic_id.Write()

        self.nue_template.Write()
        self.numu_template.Write()
        self.swap_template.Write()

        for h in self.mc_hists:
            self.data_hists[h].SetName("data_"+h)
            self.mc_hists[h].SetName("mc_"+h)

            self.data_hists[h].Write()
            self.mc_hists[h].Write()

        for h in self.keys:
            self.nue_hists[h].SetName("nue_"+h)
            self.numu_hists[h].SetName("numu_"+h)
            self.swap_hists[h].SetName("swap_"+h)
            self.nue_templates[h].SetName("nue_temp_"+h)
            self.numu_templates[h].SetName("numu_temp_"+h)
            self.swap_templates[h].SetName("swap_temp_"+h)

            self.nue_hists[h].Write()
            self.numu_hists[h].Write()
            self.swap_hists[h].Write()
            self.nue_templates[h].Write()
            self.numu_templates[h].Write()
            self.swap_templates[h].Write()

        for name, cov in self.external_covariances.items():
            cov.Write("{}_paper_covariance_matrix".format(name), ROOT.TObject.kOverwrite)

        f.Close()

    # def Load(self,filename):
    #     f = ROOT.TFile.Open(filename)

    #     name = "sample"

    #     hist = f.Get(name+"_mc")
    #     self.SetMCHistogram(hist)

    #     hist = f.Get(name+"_data")
    #     self.SetDataHistogram(hist)

    #     self.pseudo_hist = f.Get(name+"_pseudo")
        
    #     self.nue_hist = f.Get(name+"_nue")
    #     self.numu_hist = f.Get(name+"_numu")
    #     self.swap_hist = f.Get(name+"_swap")

    #     self.beam_id = f.Get(name+"_beamID")
    #     self.elastic_id = f.Get(name+"_nueelID")
    #     self.ratio_id = f.Get(name+"_ratioID")

    #     self.nue_template   = f.Get(name+"_nue_template")
    #     self.numu_template  = f.Get(name+"_numu_template")
    #     self.swap_template  = f.Get(name+"_swap_template")

    #     for h in self.keys:
    #         self.data_hists[h] = f.Get("data_"+h)
    #         self.mc_hists[h] = f.Get("mc_"+h)
    #         if type(self.data_hists[h]) != PlotUtils.MnvH1D:
    #             del self.data_hists[h]
    #             del self.mc_hists[h]

    #         self.nue_hists[h] = f.Get("nue_"+h)
    #         self.numu_hists[h] = f.Get("numu_"+h)
    #         self.swap_hists[h] = f.Get("swap_"+h)
    #         self.nue_templates[h] = f.Get("nue_temp_"+h)
    #         self.numu_templates[h] = f.Get("numu_temp_"+h)
    #         self.swap_templates[h] = f.Get("swap_temp_"+h)

    #     test_hist = f.Get("mc_fhc_ratio")
    #     if type(test_hist) == PlotUtils.MnvH1D:
    #         self.mc_hists["fhc_ratio"] = test_hist
    #         self.mc_hists["rhc_ratio"] = f.Get("mc_rhc_ratio")
    #         self.data_hists["fhc_ratio"] = f.Get("data_fhc_ratio")
    #         self.data_hists["rhc_ratio"] = f.Get("data_rhc_ratio")
    #         self.keys.extend(["fhc_ratio","rhc_ratio"])
    #         self.stitchKeys = self.keys.copy()
    #         self.titles["fhc_ratio"] = "FHC Ratio"
    #         self.titles["rhc_ratio"] = "RHC Ratio"

    #     f.Close()

    #     for h1 in self.nue_hists:
    #         self.nue_hists[h1].AddMissingErrorBandsAndFillWithCV(self.mc_hist)
    #         self.numu_hists[h1].AddMissingErrorBandsAndFillWithCV(self.mc_hist)
    #         self.swap_hists[h1].AddMissingErrorBandsAndFillWithCV(self.mc_hist)
    #         self.data_hists[h1].AddMissingErrorBandsAndFillWithCV(self.mc_hist)
    #         self.nue_hists[h1].AddMissingErrorBandsAndFillWithCV(self.data_hist)
    #         self.numu_hists[h1].AddMissingErrorBandsAndFillWithCV(self.data_hist)
    #         self.swap_hists[h1].AddMissingErrorBandsAndFillWithCV(self.data_hist)
    #         self.data_hists[h1].AddMissingErrorBandsAndFillWithCV(self.data_hist)

    #     self.SetCovarianceMatrices()
    def Load(self, filename):
        f = ROOT.TFile.Open(filename)

        if not f or f.IsZombie():
            raise RuntimeError("Could not open stitched file: {}".format(filename))

        name = "sample"

        hist = f.Get(name + "_mc")
        if not hist or type(hist) != PlotUtils.MnvH1D:
            raise RuntimeError("Could not load {}_mc as MnvH1D".format(name))
        self.SetMCHistogram(hist)

        hist = f.Get(name + "_data")
        if not hist or type(hist) != PlotUtils.MnvH1D:
            raise RuntimeError("Could not load {}_data as MnvH1D".format(name))
        self.SetDataHistogram(hist)

        self.pseudo_hist = f.Get(name + "_pseudo")

        self.nue_hist = f.Get(name + "_nue")
        self.numu_hist = f.Get(name + "_numu")
        self.swap_hist = f.Get(name + "_swap")

        self.beam_id = f.Get(name + "_beamID")
        self.elastic_id = f.Get(name + "_nueelID")
        self.ratio_id = f.Get(name + "_ratioID")

        self.nue_template = f.Get(name + "_nue_template")
        self.numu_template = f.Get(name + "_numu_template")
        self.swap_template = f.Get(name + "_swap_template")

        # Only load sample keys that actually exist in the ROOT file.
        loaded_keys = []

        for h in self.keys:
            data_h = f.Get("data_" + h)
            mc_h = f.Get("mc_" + h)

            if type(data_h) != PlotUtils.MnvH1D or type(mc_h) != PlotUtils.MnvH1D:
                # print("Skipping missing/non-MnvH1D sample:", h)
                continue

            self.data_hists[h] = data_h
            self.mc_hists[h] = mc_h

            nue_h = f.Get("nue_" + h)
            numu_h = f.Get("numu_" + h)
            swap_h = f.Get("swap_" + h)

            if type(nue_h) != PlotUtils.MnvH1D:
                print("  WARNING: missing nue_{}; using empty hist".format(h))
                nue_h = self.EmptyHist(mc_h)

            if type(numu_h) != PlotUtils.MnvH1D:
                print("  WARNING: missing numu_{}; using empty hist".format(h))
                numu_h = self.EmptyHist(mc_h)

            if type(swap_h) != PlotUtils.MnvH1D:
                print("  WARNING: missing swap_{}; using empty hist".format(h))
                swap_h = self.EmptyHist(mc_h)

            self.nue_hists[h] = nue_h
            self.numu_hists[h] = numu_h
            self.swap_hists[h] = swap_h

            self.nue_templates[h] = f.Get("nue_temp_" + h)
            self.numu_templates[h] = f.Get("numu_temp_" + h)
            self.swap_templates[h] = f.Get("swap_temp_" + h)

            loaded_keys.append(h)

        self.keys = loaded_keys
        self.stitchKeys = loaded_keys.copy()

        test_hist = f.Get("mc_fhc_ratio")
        if type(test_hist) == PlotUtils.MnvH1D:
            self.mc_hists["fhc_ratio"] = test_hist
            self.mc_hists["rhc_ratio"] = f.Get("mc_rhc_ratio")
            self.data_hists["fhc_ratio"] = f.Get("data_fhc_ratio")
            self.data_hists["rhc_ratio"] = f.Get("data_rhc_ratio")

            self.keys.extend(["fhc_ratio", "rhc_ratio"])
            self.stitchKeys = self.keys.copy()

            self.titles["fhc_ratio"] = "FHC Ratio"
            self.titles["rhc_ratio"] = "RHC Ratio"

        # print("Loaded stitched samples:")
        # for h in self.keys:
        #     print("  ", h)

        self.LoadExternalCovariances(f)
        # Do not close f before cloning if ROOT owns the objects.
        # Better: clone everything loaded from file before closing.
        # But if your old code worked with f.Close(), leave this for now.
        f.Close()

        for h1 in self.keys:
            if h1 not in self.nue_hists:
                continue

            self.nue_hists[h1].AddMissingErrorBandsAndFillWithCV(self.mc_hist)
            self.numu_hists[h1].AddMissingErrorBandsAndFillWithCV(self.mc_hist)
            self.swap_hists[h1].AddMissingErrorBandsAndFillWithCV(self.mc_hist)
            self.data_hists[h1].AddMissingErrorBandsAndFillWithCV(self.mc_hist)

            self.nue_hists[h1].AddMissingErrorBandsAndFillWithCV(self.data_hist)
            self.numu_hists[h1].AddMissingErrorBandsAndFillWithCV(self.data_hist)
            self.swap_hists[h1].AddMissingErrorBandsAndFillWithCV(self.data_hist)
            self.data_hists[h1].AddMissingErrorBandsAndFillWithCV(self.data_hist)

        self.SetCovarianceMatrices()


    def SetPlottingStyle(self):
        MNVPLOTTER = PlotUtils.MnvPlotter()

        for i,name in enumerate(self.data_hists):
            self.data_hists[name].SetLineColor(ROOT.kBlack)
            self.data_hists[name].SetMarkerStyle(20)
            self.data_hists[name].SetMarkerSize(1)
            self.data_hists[name].SetTitle(self.titles[name])
            self.data_hists[name].SetTitleFont(62)
            self.data_hists[name].SetTitleSize(0.06)
            MNVPLOTTER.ApplyAxisStyle(self.data_hists[name])
            self.data_hists[name].GetXaxis().SetNdivisions(510) #5 minor divisions between 9 major divisions.  I'm trying to match a specific paper here.

            if "elastic" in name:
                self.data_hists[name].GetXaxis().SetTitle("Electron Energy [ GeV ]")
                self.data_hists[name].GetYaxis().SetTitle("NEvents / 2 GeV")
            elif "imd" in name:
                self.data_hists[name].GetXaxis().SetTitle("Muon Energy [ GeV ]")
                self.data_hists[name].GetYaxis().SetTitle("NEvents / GeV")
            else:
                self.data_hists[name].GetXaxis().SetTitle("Neutrino Energy Estimator [ GeV ]")
                self.data_hists[name].GetYaxis().SetTitle("NEvents / GeV")

            self.mc_hists[name].SetLineColor(ROOT.kRed)
            self.mc_hists[name].SetLineWidth(2)
            self.mc_hists[name].SetMarkerStyle(0)
            self.mc_hists[name].SetLineStyle(1)
            self.mc_hists[name].SetTitle(self.titles[name])
            MNVPLOTTER.ApplyAxisStyle(self.mc_hists[name])
            self.mc_hists[name].GetXaxis().SetNdivisions(510) #5 minor divisions between 9 major divisions.  I'm trying to match a specific paper here.

            if "elastic" in name:
                self.mc_hists[name].GetXaxis().SetTitle("Electron Energy [ GeV ]")
                self.mc_hists[name].GetYaxis().SetTitle("NEvents / 2 GeV")
            elif "imd" in name:
                self.mc_hists[name].GetXaxis().SetTitle("Muon Energy [ GeV ]")
                self.mc_hists[name].GetYaxis().SetTitle("NEvents / GeV")
            else:
                self.mc_hists[name].GetXaxis().SetTitle("Neutrino Energy Estimator [ GeV ]")
                self.mc_hists[name].GetYaxis().SetTitle("NEvents / GeV")

        for i,name in enumerate(self.mc_hists):
            self.mc_hists[name].SetLineColor(ROOT.kRed)
            self.mc_hists[name].SetLineWidth(2)
            self.mc_hists[name].SetMarkerStyle(0)
            self.mc_hists[name].SetLineStyle(1)
            self.mc_hists[name].SetTitle(self.titles[name])
            MNVPLOTTER.ApplyAxisStyle(self.mc_hists[name])
            self.mc_hists[name].GetXaxis().SetNdivisions(510) #5 minor divisions between 9 major divisions.  I'm trying to match a specific paper here.
            
            if "elastic" in name:
                self.mc_hists[name].GetXaxis().SetTitle("Electron Energy [ GeV ]")
                self.mc_hists[name].GetYaxis().SetTitle("NEvents / 2 GeV")
            elif "imd" in name:
                self.mc_hists[name].GetXaxis().SetTitle("Muon Energy [ GeV ]")
                self.mc_hists[name].GetYaxis().SetTitle("NEvents / GeV")
            else:
                self.mc_hists[name].GetXaxis().SetTitle("Neutrino Energy Estimator [ GeV ]")
                self.mc_hists[name].GetYaxis().SetTitle("NEvents / GeV")

            self.data_hists[name].SetLineColor(ROOT.kBlack)
            self.data_hists[name].SetMarkerStyle(20)
            self.data_hists[name].SetMarkerSize(1)
            self.data_hists[name].SetTitle(self.titles[name])
            MNVPLOTTER.ApplyAxisStyle(self.data_hists[name])
            self.data_hists[name].GetXaxis().SetNdivisions(510) #5 minor divisions between 9 major divisions.  I'm trying to match a specific paper here.

            if "elastic" in name:
                self.data_hists[name].GetXaxis().SetTitle("Electron Energy [ GeV ]")
                self.data_hists[name].GetYaxis().SetTitle("NEvents / 2 GeV")
            elif "imd" in name:
                self.data_hists[name].GetXaxis().SetTitle("Muon Energy [ GeV ]")
                self.data_hists[name].GetYaxis().SetTitle("NEvents / GeV")
            else:
                self.data_hists[name].GetXaxis().SetTitle("Neutrino Energy Estimator [ GeV ]")
                self.data_hists[name].GetYaxis().SetTitle("NEvents / GeV")

        if type(self.mc_hist) == PlotUtils.MnvH1D:
            self.mc_hist.SetLineColor(ROOT.kRed)
            self.mc_hist.SetLineWidth(2)
            self.mc_hist.SetMarkerStyle(0)
            self.mc_hist.SetLineStyle(1)
            self.mc_hist.GetXaxis().SetTitle("Bin Number")
            self.mc_hist.GetYaxis().SetTitle("Entries")
            MNVPLOTTER.ApplyAxisStyle(self.mc_hist)
            self.mc_hist.GetXaxis().SetNdivisions(510) #5 minor divisions between 9 major divisions.  I'm trying to match a specific paper here.

        if type(self.osc_hist) == PlotUtils.MnvH1D:
            self.osc_hist.SetLineColor(ROOT.kBlue)
            self.osc_hist.SetLineWidth(2)
            self.osc_hist.SetMarkerStyle(0)
            self.osc_hist.SetLineStyle(1)
            self.osc_hist.GetXaxis().SetTitle("Bin number")
            self.osc_hist.GetYaxis().SetTitle("Entries")
            MNVPLOTTER.ApplyAxisStyle(self.osc_hist)
            self.osc_hist.GetXaxis().SetNdivisions(510) #5 minor divisions between 9 major divisions.  I'm trying to match a specific paper here.

        if type(self.data_hist) == PlotUtils.MnvH1D:
            self.data_hist.SetMarkerStyle(20)
            self.data_hist.SetMarkerSize(1)
            self.data_hist.SetMarkerColor(ROOT.kBlack)
            self.data_hist.SetLineColor(ROOT.kBlack)
            self.data_hist.GetXaxis().SetTitle("Bin number")
            self.data_hist.GetYaxis().SetTitle("Entries")
            MNVPLOTTER.ApplyAxisStyle(self.data_hist)
            self.data_hist.GetXaxis().SetNdivisions(510) #5 minor divisions between 9 major divisions.  I'm trying to match a specific paper here.

    # def StitchThis2D(self):
    #     for i in range(1, self.mc_hist.GetNbinsX()+1):
    #         h = self.bin_dictionary[i]['sample']
    #         sample_i = self.bin_dictionary[i]['bin']

    #         h_nue = self.nue_templates[h].Clone()
    #         h_numu = self.numu_templates[h].Clone()
    #         h_swap = self.swap_templates[h].Clone()

    #         # ----- nu_e sample L/E templates ----- #
    #         colInt = 0

    #         # Get integrated column content to normalize
    #         for c in range(1,h_nue.GetNbinsX()+1):
    #             if isinstance(h_nue,ROOT.TH1D):
    #                 colInt+=h_nue.GetBinContent(c)
    #             else:
    #                 colInt+=h_nue.GetBinContent(c,sample_i)

    #         # Set each bin in stiched histogram with fractional column content
    #         for c in range(1,h_nue.GetNbinsX()+1):
    #             if isinstance(h_nue,ROOT.TH1D):
    #                 bin_c = h_nue.GetBinContent(c)
    #             else:
    #                 bin_c = h_nue.GetBinContent(c,sample_i)
    #             self.nue_template.SetBinContent(i, c, bin_c/colInt if colInt > 0 else 0)

    #         # ----- nu_mu sample L/E templates ----- #
    #         colInt = 0

    #         # Get integrated column content to normalize
    #         for c in range(1,h_numu.GetNbinsX()+1):
    #             if isinstance(h_numu,ROOT.TH1D):
    #                 colInt+=h_numu.GetBinContent(c)
    #             else:
    #                 colInt+=h_numu.GetBinContent(c,sample_i)

    #         # Set each bin in stiched histogram with fractional column content
    #         for c in range(1,h_numu.GetNbinsX()+1):
    #             if isinstance(h_numu,ROOT.TH1D):
    #                 bin_c = h_numu.GetBinContent(c)
    #             else:
    #                 bin_c = h_numu.GetBinContent(c,sample_i)
    #             self.numu_template.SetBinContent(i, c, bin_c/colInt if colInt > 0 else 0)

    #         # ----- nu_mu->nu_e sample L/E templates ----- #
    #         colInt = 0

    #         # Get integrated column content to normalize
    #         for c in range(1,h_swap.GetNbinsX()+1):
    #             if isinstance(h_swap,ROOT.TH1D):
    #                 colInt+=h_swap.GetBinContent(c)
    #             else:
    #                 colInt+=h_swap.GetBinContent(c,sample_i)

    #         # Set each bin in stiched histogram with fractional column content
    #         for c in range(1,h_swap.GetNbinsX()+1):
    #             if isinstance(h_swap,ROOT.TH1D):
    #                 bin_c = h_swap.GetBinContent(c)
    #             else:
    #                 bin_c = h_swap.GetBinContent(c,sample_i)
    #             self.swap_template.SetBinContent(i, c, bin_c/colInt if colInt > 0 else 0)
    def StitchThis2D(self):
        def copy_normalized_slice(
            source,
            destination,
            stitched_bin,
            source_reco_bin,
            label,
            sample_name,
        ):
            if source is None:
                raise RuntimeError(
                    "Missing {} template for {}".format(
                        label,
                        sample_name,
                    )
                )

            if source.InheritsFrom("TH2"):
                if source_reco_bin < 1 or source_reco_bin > source.GetNbinsY():
                    raise RuntimeError(
                        "{} reconstructed bin {} is outside template y range "
                        "1..{} for {}".format(
                            label,
                            source_reco_bin,
                            source.GetNbinsY(),
                            sample_name,
                        )
                    )

                n_le_bins = source.GetNbinsX()

                if n_le_bins != destination.GetNbinsY():
                    raise RuntimeError(
                        "{} L/E bin count mismatch for {}: source={}, stitched={}".format(
                            label,
                            sample_name,
                            n_le_bins,
                            destination.GetNbinsY(),
                        )
                    )

                contents = np.asarray(
                    [
                        source.GetBinContent(le_bin, source_reco_bin)
                        for le_bin in range(1, n_le_bins + 1)
                    ],
                    dtype=float,
                )

            elif source.InheritsFrom("TH1"):
                n_le_bins = source.GetNbinsX()

                if n_le_bins != destination.GetNbinsY():
                    raise RuntimeError(
                        "{} 1D L/E bin count mismatch for {}: source={}, stitched={}".format(
                            label,
                            sample_name,
                            n_le_bins,
                            destination.GetNbinsY(),
                        )
                    )

                contents = np.asarray(
                    [
                        source.GetBinContent(le_bin)
                        for le_bin in range(1, n_le_bins + 1)
                    ],
                    dtype=float,
                )

            else:
                raise TypeError(
                    "Unsupported {} template type for {}: {}".format(
                        label,
                        sample_name,
                        type(source),
                    )
                )

            if not np.all(np.isfinite(contents)):
                raise RuntimeError(
                    "Nonfinite {} template content for {}, reco bin {}".format(
                        label,
                        sample_name,
                        source_reco_bin,
                    )
                )

            if np.any(contents < 0):
                raise RuntimeError(
                    "Negative {} template content for {}, reco bin {}".format(
                        label,
                        sample_name,
                        source_reco_bin,
                    )
                )

            normalization = float(np.sum(contents))

            if normalization <= 0:
                # An empty physical component is allowed.
                for le_bin in range(1, destination.GetNbinsY() + 1):
                    destination.SetBinContent(
                        stitched_bin,
                        le_bin,
                        0.0,
                    )
                return

            fractions = contents / normalization

            for le_bin, value in enumerate(fractions, start=1):
                destination.SetBinContent(
                    stitched_bin,
                    le_bin,
                    float(value),
                )

            check = sum(
                destination.GetBinContent(stitched_bin, le_bin)
                for le_bin in range(1, destination.GetNbinsY() + 1)
            )

            if abs(check - 1.0) > 1e-10:
                raise RuntimeError(
                    "{} stitched-template normalization failed for {}, "
                    "stitched bin {}: sum={}".format(
                        label,
                        sample_name,
                        stitched_bin,
                        check,
                    )
                )

        for stitched_bin in range(1, self.mc_hist.GetNbinsX() + 1):
            sample_name = self.bin_dictionary[stitched_bin]["sample"]
            source_reco_bin = self.bin_dictionary[stitched_bin]["bin"]

            copy_normalized_slice(
                self.nue_templates[sample_name],
                self.nue_template,
                stitched_bin,
                source_reco_bin,
                "nue",
                sample_name,
            )

            copy_normalized_slice(
                self.numu_templates[sample_name],
                self.numu_template,
                stitched_bin,
                source_reco_bin,
                "numu",
                sample_name,
            )

            copy_normalized_slice(
                self.swap_templates[sample_name],
                self.swap_template,
                stitched_bin,
                source_reco_bin,
                "swap",
                sample_name,
            )



    def RemoveSystematics(self, exclude):
        for name in self.mc_hist.GetVertErrorBandNames():
            if name in list(ERROR_GROUPS_CONFIG[exclude.upper()].values())[0]:
                print("removing {} from systematics".format(name))
                self.mc_hist.PopVertErrorBand(name)
                self.data_hist.PopVertErrorBand(name)
                self.pseudo_hist.PopVertErrorBand(name)
                for h in self.stitchKeys:
                    self.mc_hists[h].PopVertErrorBand(name)
                    self.data_hists[h].PopVertErrorBand(name)
                    if h in self.nue_hists:
                        self.nue_hists[h].PopVertErrorBand(name)
                        self.numu_hists[h].PopVertErrorBand(name)
                        self.swap_hists[h].PopVertErrorBand(name)
        self.SetCovarianceMatrices()
        self.SetAMatrix()

    def StitchThis(self):
        for i in range(1,self.mc_hist.GetNbinsX()+1):
            h = self.bin_dictionary[i]['sample']
            sample_i = self.bin_dictionary[i]['bin']

            h_mc = self.mc_hists[h].Clone()
            h_data = self.data_hists[h].Clone()


            # ----- do MC stitching ----- #
            bin_mc = h_mc.GetBinContent(sample_i)
            err_mc = h_mc.GetBinError(sample_i)

            self.mc_hist.SetBinContent(i,bin_mc)
            self.pseudo_hist.SetBinContent(i,bin_mc)
            
            self.mc_hist.SetBinError(i,err_mc)
            self.pseudo_hist.SetBinError(i,err_mc)

            # # no statistical error on elastic scattering special production
            # if 'elastic' in h or 'imd' in h:# or "nue" in h:
            #     self.mc_hists[h].SetBinError(sample_i,0)
            #     self.mc_hist.SetBinError(i,0)
            #     self.nue_hist.SetBinError(i,0)
            #     self.numu_hist.SetBinError(i,0)

            # ----- do data stitching ----- #
            bin_data = h_data.GetBinContent(sample_i)
            err_data = h_data.GetBinError(sample_i)

            self.data_hist.SetBinContent(i,bin_data)
            self.data_hist.SetBinError(i,err_data)

            # ----- do other stitching ----- #
            if not self.dirty:
                self.nue_hist.SetBinContent(i,self.nue_hists[h].GetBinContent(sample_i))
                self.numu_hist.SetBinContent(i,self.numu_hists[h].GetBinContent(sample_i))
                self.swap_hist.SetBinContent(i,self.swap_hists[h].GetBinContent(sample_i))

                self.nue_hist.SetBinError(i,self.nue_hists[h].GetBinError(sample_i))
                self.numu_hist.SetBinError(i,self.numu_hists[h].GetBinError(sample_i))
                self.swap_hist.SetBinError(i,self.swap_hists[h].GetBinError(sample_i))

                self.beam_id.SetBinContent(i,self.beam_ids[h].GetBinContent(sample_i))
                self.ratio_id.SetBinContent(i,self.ratio_ids[h].GetBinContent(sample_i))
                self.elastic_id.SetBinContent(i,self.elastic_ids[h].GetBinContent(sample_i))

                self.beam_id.SetBinError(i,self.beam_ids[h].GetBinError(sample_i))
                self.ratio_id.SetBinError(i,self.ratio_ids[h].GetBinError(sample_i))
                self.elastic_id.SetBinError(i,self.elastic_ids[h].GetBinError(sample_i))

    def StitchPreservationHists(self):
        # preservation_hist_container["plot name"] = { "fhc_numu_selection": {"numu" : hist, ...}, ...}
        # preservation_hists = {"plot name": {"numu" : hist1,
        #                                     "nue"  : hist2,
        #                                     "anumu": hist3,
        #                                     "anue" : hist4 }, ...}
        for plotName in self.preservation_hists.keys():
            for pdgKey in self.preservation_hists[plotName].keys():
                stitched_hist = self.preservation_hists[plotName][pdgKey]
                histDict = self.preservation_hist_container[plotName]
                for i in range(1,stitched_hist.GetNbinsX()+1):
                    h = self.bin_dictionary[i]['sample']
                    sample_i = self.bin_dictionary[i]['bin'] 

                    emptyHist = False
                    if h in histDict.keys() and pdgKey in histDict[h].keys():
                        hist = histDict[h][pdgKey]
                    elif 'ratio' in h:
                        beam = h[:3]
                        flavor = pdgKey.replace("a","")
                        newH = beam+"_"+flavor+"_selection"
                        if (newH in histDict.keys() and pdgKey in histDict[newH].keys()):
                            hist = histDict[newH][pdgKey]
                        else:
                            hist = self.EmptyHist(self.mc_hists[h])
                            emptyHist = True
                    else:
                        hist = self.EmptyHist(self.mc_hists[h])
                        emptyHist = True

                    # ----- do stitching ----- #
                    if type(stitched_hist) == ROOT.TH2D:
                        for j in range(1,hist.GetNbinsY()+1):
                            if emptyHist:
                                stitched_hist.SetBinContent(i,j,0)
                            else:
                                bin_mc = hist.GetBinContent(sample_i,j)
                                stitched_hist.SetBinContent(i,j,bin_mc)
                    elif type(stitched_hist) == ROOT.TH1D:
                        bin_mc = hist.GetBinContent(sample_i)
                        stitched_hist.SetBinContent(i,bin_mc)    

            
    def FillErrorBandsFromDict(self):
        for i in range(1,self.mc_hist.GetNbinsX()+1):
            h = self.bin_dictionary[i]['sample']
            sample_i = self.bin_dictionary[i]['bin']
            #print("filling bin {} in stitched histogram with bin {} from {}".format(i,sample_i,h))
            h_mc = self.mc_hists[h]
            h_data = self.data_hists[h]

            if not self.is_processed:
                raise ValueError("Histograms have not been synced")

            errorband_names_vert = h_mc.GetVertErrorBandNames()
            for error_band in errorband_names_vert:
                self.mc_hist.GetVertErrorBand(error_band).SetUseSpreadError(h_mc.GetVertErrorBand(error_band).GetUseSpreadError())
                self.pseudo_hist.GetVertErrorBand(error_band).SetUseSpreadError(h_mc.GetVertErrorBand(error_band).GetUseSpreadError())
                self.data_hist.GetVertErrorBand(error_band).SetUseSpreadError(h_data.GetVertErrorBand(error_band).GetUseSpreadError())

                self.nue_hist.GetVertErrorBand(error_band).SetUseSpreadError(h_mc.GetVertErrorBand(error_band).GetUseSpreadError())
                self.numu_hist.GetVertErrorBand(error_band).SetUseSpreadError(h_mc.GetVertErrorBand(error_band).GetUseSpreadError())
                self.swap_hist.GetVertErrorBand(error_band).SetUseSpreadError(h_mc.GetVertErrorBand(error_band).GetUseSpreadError())

                n_univ = h_mc.GetVertErrorBand(error_band).GetNHists()
                if not self.mc_hist.HasVertErrorBand(error_band) and h_mc.HasVertErrorBand(error_band):
                    raise ValueError("MC histograms were not properly synchronized")
                if not self.data_hist.HasVertErrorBand(error_band) and h_data.HasVertErrorBand(error_band):
                    raise ValueError("Data histograms were not properly synchronized")

                for universe in range(0, n_univ):
                    bin_mc = h_mc.GetBinContent(sample_i)
                    bin_data = h_data.GetBinContent(sample_i)

                    sys_mc = h_mc.GetVertErrorBand(error_band).GetHist(universe).GetBinContent(sample_i)
                    sys_data = h_data.GetVertErrorBand(error_band).GetHist(universe).GetBinContent(sample_i)

                    if not self.dirty:
                        sys_nue = self.nue_hists[h].GetVertErrorBand(error_band).GetHist(universe).GetBinContent(sample_i)
                        sys_numu = self.numu_hists[h].GetVertErrorBand(error_band).GetHist(universe).GetBinContent(sample_i)
                        sys_swap = self.swap_hists[h].GetVertErrorBand(error_band).GetHist(universe).GetBinContent(sample_i)

                    ratio = sys_data/bin_data if bin_data != 0 else 0
                    sys_pseudo = bin_mc*ratio

                    self.mc_hist.GetVertErrorBand(error_band).GetHist(universe).SetBinContent(i,sys_mc)
                    self.data_hist.GetVertErrorBand(error_band).GetHist(universe).SetBinContent(i,sys_data)
                    self.pseudo_hist.GetVertErrorBand(error_band).GetHist(universe).SetBinContent(i,sys_pseudo)

                    if not self.dirty:
                        self.nue_hist.GetVertErrorBand(error_band).GetHist(universe).SetBinContent(i,sys_nue)
                        self.numu_hist.GetVertErrorBand(error_band).GetHist(universe).SetBinContent(i,sys_numu)
                        self.swap_hist.GetVertErrorBand(error_band).GetHist(universe).SetBinContent(i,sys_swap)

    def SyncErrorBands(self,rename_bands=True):
        if len(self.mc_hists) > 0:
            if rename_bands:
                for h in self.keys:
                    self.RenameBands(self.data_hists[h])
                    self.RenameBands(self.mc_hists[h])
                    self.RenameBands(self.nue_hists[h])
                    self.RenameBands(self.numu_hists[h])
                    self.RenameBands(self.swap_hists[h])
            for h1 in self.keys:
                for h2 in self.keys:
                    if h1 != h2:
                        self.mc_hists[h1].AddMissingErrorBandsAndFillWithCV(self.data_hists[h1])
                        self.mc_hists[h1].AddMissingErrorBandsAndFillWithCV(self.data_hists[h2])
                        self.mc_hists[h1].AddMissingErrorBandsAndFillWithCV(self.mc_hists[h2])
                        self.data_hists[h1].AddMissingErrorBandsAndFillWithCV(self.mc_hists[h1])
                        if not self.dirty:
                            self.nue_hists[h1].AddMissingErrorBandsAndFillWithCV(self.mc_hists[h1])
                            self.numu_hists[h1].AddMissingErrorBandsAndFillWithCV(self.mc_hists[h1])
                            self.swap_hists[h1].AddMissingErrorBandsAndFillWithCV(self.mc_hists[h1])

                self.FillOutFlux(self.mc_hists[h1])
                self.FillOutFlux(self.data_hists[h1])
                if not self.dirty:
                    self.FillOutFlux(self.nue_hists[h1])
                    self.FillOutFlux(self.numu_hists[h1])
                    self.FillOutFlux(self.swap_hists[h1])

        if type(self.mc_hist) == PlotUtils.MnvH1D:
            if type(self.mc_hist) == type(self.data_hist) and type(self.mc_hist) == type(self.pseudo_hist):
                self.mc_hist.AddMissingErrorBandsAndFillWithCV(self.data_hist)
                self.data_hist.AddMissingErrorBandsAndFillWithCV(self.mc_hist)
                self.pseudo_hist.AddMissingErrorBandsAndFillWithCV(self.data_hist)
                self.data_hist.AddMissingErrorBandsAndFillWithCV(self.pseudo_hist)
                self.mc_hist.AddMissingErrorBandsAndFillWithCV(self.data_hist)
            else:
                raise ValueError("Histograms are not all set when trying to sync error bands")

    def FillOutFlux(self,h):
        name = "Flux"
        if h.GetVertErrorBand(name).GetNHists() > self.n_flux_universes:
            self.ShortFlux(h)
        elif h.GetVertErrorBand(name).GetNHists() < self.n_flux_universes:
            self.LongFlux(h)

    def ShortFlux(self,h):
        name = "Flux"
        h_hists = h.GetVertErrorBand(name).GetHists()
        h_hists = [h_hists[i] for i in range(self.n_flux_universes)]
        useSpread = h.GetVertErrorBand(name).GetUseSpreadError()
        errband = h.GetVertErrorBand(name)
        h.PopVertErrorBand(name)
        h.AddVertErrorBand(name,h_hists)
        h.GetVertErrorBand(name).SetUseSpreadError(useSpread)
        for i in range(h.GetNbinsX()+1):
            h.GetVertErrorBand(name).SetBinContent(i,errband.GetBinContent(i))


    def LongFlux(self,h):
        name = "Flux"
        cvClone = h.GetCVHistoWithStatError()
        h_hists = h.GetVertErrorBand(name).GetHists()
        h_hists = [h_hists[i] for i in range(h.GetVertErrorBand(name).GetNHists())]
        h_hists.extend([cvClone.Clone() for i in range(self.n_flux_universes-h.GetVertErrorBand(name).GetNHists())])
        useSpread = h.GetVertErrorBand(name).GetUseSpreadError()
        errband = h.GetVertErrorBand(name)
        h.PopVertErrorBand(name)
        h.AddVertErrorBand(name,h_hists)
        h.GetVertErrorBand(name).SetUseSpreadError(useSpread)
        for i in range(h.GetNbinsX()+1):
            h.GetVertErrorBand(name).SetBinContent(i,errband.GetBinContent(i))

    def ReweightCV(self,histogram,fluxSolution,cv=None,mc=None):
        if cv is None:
            cv = np.array(histogram)[1:-1]
        
        band = histogram.GetVertErrorBand("Flux")
        nhists = band.GetNHists()
        universes = np.array([np.array(band.GetHist(l))[1:-1] for l in range(nhists)])
        cv_table = np.array([cv for l in range(len(universes))])
        A = universes - cv_table

        if mc is None:
            mc = cv

        new_cv = mc + fluxSolution @ A

        weights = histogram.GetCVHistoWithStatError()
        for j in range(1,weights.GetNbinsX()+1):
            weight = weights.GetBinContent(j) / new_cv[j-1] if new_cv[j-1] != 0 else weights.GetBinContent(j)
            weights.SetBinContent(j,weight)
            weights.SetBinError(j,0)

        histogram.DivideSingle(histogram,weights)
        return(weights)

    def RenameBands(self,hist):
        for name in hist.GetVertErrorBandNames():
            if str(name) in errorbandDict.keys():
                universes = hist.GetVertErrorBand(name).GetHists()
                useSpread = hist.GetVertErrorBand(name).GetUseSpreadError()
                hist.AddVertErrorBand(errorbandDict[str(name)],universes)
                hist.GetVertErrorBand(errorbandDict[str(name)]).SetUseSpreadError(useSpread)
                for i in range(hist.GetNbinsX()+1):
                    hist.GetVertErrorBand(errorbandDict[str(name)]).SetBinContent(i,hist.GetVertErrorBand(name).GetBinContent(i))
                hist.PopVertErrorBand(name)

    def ReweightFluxToCV(self,name):
        def reweight(name,hist):
            weights = np.loadtxt(name)
            cv = np.array(hist.GetCVHistoWithError())[1:-1]
            summed_dev = np.zeros(cv.shape)
            for i in range(0,100):
                weight = weights[i]
                flux_univ = np.array(hist.GetVertErrorBand("Flux").GetHist(i))[1:-1]
                dev = flux_univ - cv
                dev*= weight
                summed_dev+=dev

            for i in range(hist.GetNbinsX()):
                hist.SetBinContent(i,hist.GetBinContent(i)+summed_dev[i])

            hist.PopVertErrorBand("Flux")

        reweight(name,self.mc_hist)
        reweight(name,self.data_hist)

    def WriteCSVs(self):
        dataprint = np.array(self.data_hist)[1:-1] # store MC bin contents excluding over/underflow bins
        mcprint = np.array(self.mc_hist)[1:-1]
        
        np.savetxt("csvs/mc_cv.csv",mcprint,delimiter=',')
        np.savetxt("csvs/data_cv.csv",dataprint,delimiter=',')
        np.savetxt("csvs/ryan_Cov.csv",self.covariance,delimiter=',')
        np.savetxt("csvs/ryan_CovSansFlux.csv",self.covariance_sans_flux,delimiter=',')
        np.savetxt("csvs/ryan_InvCov.csv",self.inv_covariance,delimiter=',')
        np.savetxt("csvs/ryan_InvCovSansFlux.csv",self.inv_covariance_sans_flux,delimiter=',')
        np.savetxt("csvs/ryan_Amatrix.csv",self.A,delimiter=',')
        print("Finished writing CV files")
        for plotName in self.preservation_hists.keys():
            for pdgKey in self.preservation_hists[plotName].keys():
                hist = self.preservation_hists[plotName][pdgKey]
                x_bins = hist.GetNbinsX()
                y_bins = hist.GetNbinsY()
                bins = np.zeros((x_bins,y_bins))
                for y_bin in range(y_bins):
                    for x_bin in range(x_bins): 
                        bins[x_bin,y_bin] = hist.GetBinContent(x_bin + 1,y_bin + 1)

                fileExt = plotName.replace(" ","_")
                fileExt = fileExt.replace("/","_")
                fileExt+= "_"+pdgKey
                filename = "csvs/"+fileExt+".csv"
                np.savetxt(filename,bins,delimiter=',')

                c1 = ROOT.TCanvas()
                option = "hist"
                if type(hist) == ROOT.TH2D:
                    option = "colz"
                hist.Draw()
                c1.Print("plots/debug_{}.png".format(fileExt))

        print("Finished writing Preservation files")

    def GetSampleGrid(self,fluxSolution,ratio=False):
        if ratio:
            columns = ["elastic","imd","numu_selection","ratio"]
        else:
            columns = ["elastic","imd","numu_selection","nue_selection"]

        rows = ["fhc","rhc"]

        mc_hists = []
        mc_errs = []
        data_hists = []
        ticks = []

        for column in columns:
            row_mc_hists = []
            row_mc_errs = []
            row_data_hists = []
            row_ticks = []

            for row in rows:
                histKey = row+"_"+column
                h_mc = self.mc_hists[histKey].Clone()
                h_data = self.data_hists[histKey].Clone()

                mc = np.array(h_mc)[1:-1]
                band = h_mc.GetVertErrorBand("Flux")
                nhists = band.GetNHists()
                universes = np.array([np.array(band.GetHist(i))[1:-1] for i in range(nhists)])
                cv_table = np.array([mc for i in range(len(universes))])
                A = universes - cv_table

                new_cv = mc + fluxSolution @ A
                weights = h_mc.GetCVHistoWithStatError()
                for i in range(1,weights.GetNbinsX()+1):
                    weight = weights.GetBinContent(i) / new_cv[i-1] if new_cv[i-1] != 0 else weights.GetBinContent(i)
                    weights.SetBinContent(i,weight)
                    weights.SetBinError(i,0)

                h_mc.DivideSingle(h_mc,weights)
                h_mc.PopVertErrorBand("Flux")
                h_mc.AddMissingErrorBandsAndFillWithCV(h_data)

                mc_temp = self.mc_hists[histKey].Clone()
                weights = self.ReweightCV(mc_temp,fluxSolution=fluxSolution)
                mc_temp.PopVertErrorBand("Flux")

                mc_errs_temp = mc_temp.GetCVHistoWithError()
                data_temp = self.data_hists[histKey].Clone()

                if "elastic" in column:
                    mc_temp.Scale(2,"width")
                    mc_errs_temp.Scale(2,"width")
                    data_temp.Scale(2,"width")
                    ticks.append([0,5,10,15,20])
                elif "ratio" not in column:
                    mc_temp.Scale(1,"width")
                    mc_errs_temp.Scale(1,"width")
                    data_temp.Scale(1,"width")
                    if "imd" in column:
                        ticks.append([0,5,50])
                    else:
                        ticks.append([0,5,50,15,20])
                else:
                    ticks.append([0,5,10,15,20])

                row_mc_hists.append(mc_temp)
                row_mc_errs.append(mc_errs_temp)
                row_data_hists.append(data_temp)
                row_ticks.append(ticks)

                if row == "fhc":
                    mc_temp.GetYaxis().SetTitle("#nu-mode    Events/GeV")
                else:
                    mc_temp.GetYaxis().SetTitle("#bar{#nu}-mode    Events/GeV")

                newTitle = mc_temp.GetTitle()
                newTitle = newTitle[4:]
                mc_temp.SetTitle(newTitle)

            mc_hists.append(row_mc_hists)
            mc_errs.append(row_mc_errs)
            data_hists.append(row_data_hists)
            ticks.append(row_ticks)

        return(mc_hists,mc_errs,data_hists,ticks)

    def PlotStitchedHistogram(self,fluxSolution,plotName="sample",chi2=0,penalty=0):
        margin = .12
        bottomFraction = .2
        MNVPLOTTER = PlotUtils.MnvPlotter()
        MNVPLOTTER.draw_normalized_to_bin_width=False
        MNVPLOTTER.error_summary_group_map.clear();
        for k,v in CONSOLIDATED_ERROR_GROUPS.items():
            vec = ROOT.vector("std::string")()
            for vs in v :
                vec.push_back(vs)
            MNVPLOTTER.error_summary_group_map[k]= vec

        h_mc = self.mc_hist.Clone()
        h_data = self.data_hist.Clone()
           
        mc_hists,mc_errs,data_hists,ticks = self.GetSampleGrid(fluxSolution,ratio=False)
        canvas = plot_side_by_side(mc_hists,mc_errs,data_hists,narrow_pads=[1,4],chi2=chi2,penalty=penalty,cName="c1")
        canvas.Print("plots/sample_{}.png".format(plotName))
        err_canvas = plot_errs_side_by_side(mc_hists,[1,4],0.5,MNVPLOTTER,"c2")
        err_canvas.Print("plots/sample_{}_errors.png".format(plotName))

        mc_hists,mc_errs,data_hists,ticks = self.GetSampleGrid(fluxSolution,ratio=True)
        canvas = plot_side_by_side(mc_hists,mc_errs,data_hists,narrow_pads=[1,4],chi2=chi2,penalty=penalty,cName="c3")
        canvas.Print("plots/ratio_{}.png".format(plotName))
        err_canvas = plot_errs_side_by_side(mc_hists,[1,4],0.5,MNVPLOTTER,"c4")
        err_canvas.Print("plots/ratio_{}_errors.png".format(plotName))

    def DebugPlots(self,name=""):
        c1 = ROOT.TCanvas()
        dirname = 'plots/'
        self.nue_hist.Draw("hist")
        c1.Print(dirname+"nue_hist"+name+".png")
        self.numu_hist.Draw("hist")
        c1.Print(dirname+"numu_hist"+name+".png")
        self.swap_hist.Draw("hist")
        c1.Print(dirname+"swap_hist"+name+".png")
        self.beam_id.Draw("hist")
        c1.Print(dirname+"beam_id"+name+".png")
        self.ratio_id.Draw("hist")
        c1.Print(dirname+"ratio_id"+name+".png")
        self.elastic_id.Draw("hist")
        c1.Print(dirname+"elastic_id"+name+".png")
        self.nue_template.Draw("colz")
        c1.Print(dirname+"nue_template"+name+".png")
        self.numu_template.Draw("colz")
        c1.Print(dirname+"numu_template"+name+".png")
        self.swap_template.Draw("colz")
        c1.Print(dirname+"swap_template"+name+".png")

        for plot in self.preservation_hists.keys():
            for pdg in self.preservation_hists[plot].keys():
                hist = self.preservation_hists[plot][pdg]
                if type(hist) == ROOT.TH2D:
                    option = "colz"
                else:
                    option = "hist"

                c1 = ROOT.TCanvas()
                hist.Draw(option)
                c1.Print(dirname+plot.replace("/","")+"_"+pdg+".png")

        MNVPLOTTER = PlotUtils.MnvPlotter()
        MNVPLOTTER.axis_maximum = 1
        MNVPLOTTER.error_summary_group_map.clear();
        for k,v in CONSOLIDATED_ERROR_GROUPS.items():
            vec = ROOT.vector("std::string")()
            for vs in v :
                vec.push_back(vs)
            MNVPLOTTER.error_summary_group_map[k]= vec

        MNVPLOTTER.DrawErrorSummary(self.data_hist,"TR",True,True,0)
        c1.Print(dirname+"data_errsummary.png")
        MNVPLOTTER.DrawErrorSummary(self.mc_hist,"TR",True,True,0)
        c1.Print(dirname+"mc_errsummary.png")

        for key in self.mc_hists:
            MNVPLOTTER.DrawErrorSummary(self.mc_hists[key],"TR",True,True,0)
            c1.Print(dirname+key+"mc_errsummary.png")
            MNVPLOTTER.DrawErrorSummary(self.data_hists[key],"TR",True,True,0)
            c1.Print(dirname+key+"data_errsummary.png")


    def OscillateHistogram(self,m,U_e4,U_mu4,U_tau4,fitPseudodata=False,fitFluxUniverses=False):
        # elastic_id is 1 only for the nueel samples
        # ratio_id is 1 only for the ratio samples
        # if 0 is False, if 1 is True
        hist_nue_energy = self.nue_hist.Clone()
        hist_numu_energy = self.numu_hist.Clone()
        hist_swap_energy = self.swap_hist.Clone()
        hist_nutau_id = self.elastic_id.Clone()
        hist_ratio_id = self.ratio_id.Clone()
        hist_inv_ratio_id = hist_ratio_id.Clone()

        InvertID(hist_inv_ratio_id)

        if fitPseudodata:
            hist = self.GetPseudoHistogram()
        else:
            hist = self.GetMCHistogram()

        nue_weights = hist.GetCVHistoWithStatError().Clone()
        numu_weights = hist.GetCVHistoWithStatError().Clone()
        nuenutau_weights = hist.GetCVHistoWithStatError().Clone()
        numunue_weights = hist.GetCVHistoWithStatError().Clone()
        numunutau_weights = hist.GetCVHistoWithStatError().Clone()

        for weight_hist in [
            nue_weights,
            numu_weights,
            nuenutau_weights,
            numunue_weights,
            numunutau_weights,
        ]:
            weight_hist.Reset()

        # for i in range(0,hist.GetNbinsX() + 1):
        for i in range(1,hist.GetNbinsX() + 1):
            nue_sin = sin_average(i,m,self.nue_template)
            numu_sin = sin_average(i,m,self.numu_template)
            swap_sin = sin_average(i,m,self.swap_template)

            P_ee = float(1 - 4*U_e4*(1-U_e4)*nue_sin)

            P_mue = float(4*(U_e4)*(U_mu4)*swap_sin)

            P_mumu = float(1 - 4*U_mu4*(1-U_mu4)*numu_sin)

            P_mutau = float(4*U_tau4*U_mu4*numu_sin)

            P_etau = float(4*U_e4*U_tau4*nue_sin)

            nue_weights.SetBinContent(i,P_ee)
            numu_weights.SetBinContent(i,P_mumu)
            nuenutau_weights.SetBinContent(i,P_etau)
            numunue_weights.SetBinContent(i,P_mue)
            numunutau_weights.SetBinContent(i,P_mutau)

            nue_weights.SetBinError(i,0)
            numu_weights.SetBinError(i,0)
            nuenutau_weights.SetBinError(i,0)
            numunue_weights.SetBinError(i,0)
            numunutau_weights.SetBinError(i,0)
            hist_ratio_id.SetBinError(i,0)
            hist_nutau_id.SetBinError(i,0)

        numu = hist_numu_energy.Clone()
        numu.MultiplySingle(numu,numu_weights)

        nue = hist_nue_energy.Clone()
        nue.MultiplySingle(nue,nue_weights)

        numunue = hist_swap_energy.Clone()
        numunue.MultiplySingle(numunue,numunue_weights)
        nue.Add(numunue)

        numunutau = hist_numu_energy.Clone()
        numunutau.MultiplySingle(numunutau,numunutau_weights)

        nuenutau  = hist_nue_energy.Clone()
        nuenutau.MultiplySingle(nuenutau,nuenutau_weights)
        nutau = numunutau.Clone()
        nutau.Add(nuenutau)
        nutau.MultiplySingle(nutau,hist_nutau_id)

        ratio = numu.Clone()
        ratio.Divide(ratio,nue)
        ratio.MultiplySingle(ratio,hist_ratio_id)

        nue.MultiplySingle(nue,hist_inv_ratio_id)
        numu.MultiplySingle(numu,hist_inv_ratio_id)

        osc = numu.Clone()
        osc.Add(nue)
        osc.Add(ratio)
        if nutau.Integral() > 0:
            osc.Add(nutau)

        self.SetOscHistogram(osc)


    def OscillateSubHistogram(self,name,m,U_e4,U_mu4,U_tau4):
        # elastic_id is 1 only for the nueel samples
        # ratio_id is 1 only for the ratio samples
        # if 0 is False, if 1 is True
        hist_nueTemp = self.nue_templates[name]
        hist_numuTemp = self.numu_templates[name]
        hist_swapTemp = self.swap_templates[name]

        hist_nue_energy = self.nue_hists[name].Clone()
        hist_numu_energy = self.numu_hists[name].Clone()
        hist_swap_energy = self.swap_hists[name].Clone()

        hist = self.nue_hists[name].Clone()

        nue_weights = hist.GetCVHistoWithStatError().Clone()
        numu_weights = hist.GetCVHistoWithStatError().Clone()
        nuenutau_weights = hist.GetCVHistoWithStatError().Clone()
        numunue_weights = hist.GetCVHistoWithStatError().Clone()
        numunutau_weights = hist.GetCVHistoWithStatError().Clone()

        for weight_hist in [
            nue_weights,
            numu_weights,
            nuenutau_weights,
            numunue_weights,
            numunutau_weights,
        ]:
            weight_hist.Reset()

        # for i in range(0,hist.GetNbinsX() + 1):
        for i in range(1,hist.GetNbinsX() + 1):
            nue_sin = sin_average(i,m,hist_nueTemp,False)
            numu_sin = sin_average(i,m,hist_numuTemp,False)
            swap_sin = sin_average(i,m,hist_swapTemp,False)

            P_ee = float(1 - 4*U_e4*(1-U_e4)*nue_sin)
            
            if ('elastic' in name):
                scale = 1
                P_mue = float(4*(U_e4)*(U_mu4)*swap_sin) * scale
            else:
                P_mue = float(4*(U_e4)*(U_mu4)*swap_sin)

            P_mumu = float(1 - 4*U_mu4*(1-U_mu4)*numu_sin)

            P_mutau = float(4*U_tau4*U_mu4*numu_sin)

            P_etau = float(4*U_e4*U_tau4*nue_sin)

            nue_weights.SetBinContent(i,P_ee)
            numu_weights.SetBinContent(i,P_mumu)
            nuenutau_weights.SetBinContent(i,P_etau)
            numunue_weights.SetBinContent(i,P_mue)
            numunutau_weights.SetBinContent(i,P_mutau)

            nue_weights.SetBinError(i,0)
            numu_weights.SetBinError(i,0)
            nuenutau_weights.SetBinError(i,0)
            numunue_weights.SetBinError(i,0)
            numunutau_weights.SetBinError(i,0)

        numu = hist_numu_energy.Clone()
        numu.MultiplySingle(numu,numu_weights)

        nue = hist_nue_energy.Clone()
        nue.MultiplySingle(nue,nue_weights)

        numunue = hist_swap_energy.Clone()
        numunue.MultiplySingle(numunue,numunue_weights)

        numunutau = hist_numu_energy.Clone()
        numunutau.MultiplySingle(numunutau,numunutau_weights)

        nuenutau  = hist_nue_energy.Clone()
        nuenutau.MultiplySingle(nuenutau,nuenutau_weights)
        nutau = numunutau.Clone()
        nutau.Add(nuenutau)

        total = nue.Clone()
        total.Reset()

        nue.SetFillColor(ROOT.kRed)
        numu.SetFillColor(ROOT.kBlue)
        numunue.SetFillColor(ROOT.kBlue)
        nutau.SetFillColor(ROOT.kGray)

        nue.SetLineColor(ROOT.kRed)
        numu.SetLineColor(ROOT.kBlue)
        numunue.SetLineColor(ROOT.kBlue)
        nutau.SetLineColor(ROOT.kGray)

        nue.SetLineWidth(0)
        numu.SetLineWidth(0)
        numunue.SetLineWidth(0)
        nutau.SetLineWidth(0)
        self.data_hists[name].SetLineWidth(1)

        nue.SetFillStyle(3244)
        numu.SetFillStyle(3744)
        numunue.SetFillStyle(3244)
        nutau.SetFillStyle(3409)

        oscHists = []
        nue.SetTitle("#nu_{e}")
        numu.SetTitle("#nu_{#mu}")
        numunue.SetTitle("#nu_{#mu}#rightarrow #nu_{e}")
        nutau.SetTitle("#nu_{#mu}#rightarrow #nu_{#tau}")
        self.data_hists[name].SetTitle("Oscillated {}".format(name))

        if 'elastic' in name:
            oscHists.append(nue)
            oscHists.append(numunue)
            oscHists.append(numu)
            oscHists.append(nutau)
            total.Add(numu)
            total.Add(numunue)
            total.Add(nutau)
            total.Add(nue)
        elif 'imd' in name or 'numu' in name:
            oscHists.append(numu)
            total.Add(numu)
        # elif 'ratio' in name:
        #     oscHists.append(nue)
        #     oscHists.append(numunue)
        #     oscHists.append(numu)
        #     total.Add(nue)
        #     total.Add(numunue)
        #     total.Add(numu)
        elif 'ratio' in name:
            denominator = nue.Clone(name + "_osc_ratio_denominator")
            denominator.Add(numunue)

            ratio = numu.Clone(name + "_osc_ratio")
            ratio.Reset()

            for i in range(1, ratio.GetNbinsX() + 1):
                numerator_value = numu.GetBinContent(i)
                denominator_value = denominator.GetBinContent(i)

                if denominator_value <= 0:
                    ratio.SetBinContent(i, 0.0)
                    ratio.SetBinError(i, 0.0)
                    continue

                ratio.SetBinContent(
                    i,
                    numerator_value / denominator_value,
                )
                ratio.SetBinError(i, 0.0)

            ratio.SetTitle("CC #nu_{#mu}/#nu_{e}")

            oscHists.append(ratio)
            total.Add(ratio)
        elif 'nue' in name:
            oscHists.append(nue)
            oscHists.append(numunue)
            total.Add(nue)
            total.Add(numunue)
        else:
            print("No self added for {}".format(name))

        return(oscHists,total)
