import os
import copy
from collections import OrderedDict
from itertools import compress
import argparse
import logging, sys
import ROOT
import PlotUtils

import numpy as np
from root_numpy import matrix

import math
from array import array

from config.SignalDef import SWAP_SIGNAL_DEFINITION, SIGNAL_DEFINITION
from config.SystematicsConfig import CONSOLIDATED_ERROR_GROUPS 
from tools import Utilities
from tools.PlotLibrary import HistHolder
ccnueroot = os.environ.get('CCNUEROOT')

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

class StitchedHistogram:
    def __init__(self,name):
        self.name = name

        self.data_hists = OrderedDict()
        self.mc_hists = OrderedDict()

        self.mc_cov = None
        self.data_cov = None
        self.pseudo_cov = None
        self.osc_cov = None

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

        self.keys = []

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

    def __del__(self):
        self.data_hists = {}
        self.mc_hists = {}
        self.nueel_flavors = {}
        self.numu_templates = {}
        self.nue_templates = {}
        self.swap_templates = {}
        self.mc_hist = None
        self.data_hist = None
        self.pseudo_hist = None

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
        self.nue_templates[name] = nue if nue is not None else holder
        self.numu_templates[name] = numu if numu is not None else holder
        self.swap_templates[name] = swap if swap is not None else holder

    def AddScatteringFlavors(self,name,hist):
        if name in self.nueel_flavors.keys():
            print("{} has already been added to scattering sample dictionary. Doing nothing.".format(name))
        else:
            self.nueel_flavors[name] = hist

    def AddSwappedSample(self,name,hist):
        if name in self.swap_hists.keys():
            print("{} has already been added to swapped sample dictionary. Doing nothing.".format(name))
        else:
            self.swap_hists[name] = hist

    def SetDataHistogram(self,hist):
        self.data_hist = hist.Clone()
        self.mc_hist.AddMissingErrorBandsAndFillWithCV(self.data_hist)
        self.data_hist.AddMissingErrorBandsAndFillWithCV(self.mc_hist)
        self.data_cov = np.asarray(matrix(self.data_hist.GetTotalErrorMatrix(True, False, False)))[1:-1,1:-1]

    def SetPseudoHistogram(self,hist):
        self.pseudo_hist = hist.Clone()
        self.mc_hist.AddMissingErrorBandsAndFillWithCV(self.pseudo_hist)
        self.pseudo_hist.AddMissingErrorBandsAndFillWithCV(self.mc_hist)
        self.pseudo_cov = np.asarray(matrix(self.pseudo_hist.GetTotalErrorMatrix(True, False, False)))[1:-1,1:-1]

    def SetMCHistogram(self,hist):
        self.mc_hist = hist.Clone()
        self.mc_hist.AddMissingErrorBandsAndFillWithCV(self.data_hist)
        self.data_hist.AddMissingErrorBandsAndFillWithCV(self.mc_hist)
        self.mc_cov = np.asarray(matrix(self.mc_hist.GetTotalErrorMatrix(True, False, False)))[1:-1,1:-1]

    def SetOscHistogram(self,hist):
        self.osc_hist = hist.Clone()
        self.osc_cov = np.asarray(matrix(self.osc_hist.GetTotalErrorMatrix(True, False, False)))[1:-1,1:-1]

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

    def GetMCCov(self):
        return(self.mc_cov)

    def GetDataCov(self):
        return(self.data_cov)

    def GetOscCov(self):
        return(self.osc_cov)

    def GetPseudoCov(self):
        return(self.pseudo_cov)

    def EmptyHist(self,h):
        h_ret = h.Clone()
        for i in range(0,h.GetNbinsX()+1):
            h_ret.SetBinContent(i,0)
            h_ret.SetBinError(i,0)
            for name in h_ret.GetVertErrorBandNames():
                h_ret.GetVertErrorBand(name).SetBinContent(i,0)
                for univ in range(h_ret.GetVertErrorBand(name).GetNHists()):
                    h_ret.GetVertErrorBand(name).GetHist(univ).SetBinContent(i,0)

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
            raise ValueError("Cannot add {} to data histogram dictionary of {}".format(type(hist),type(list(self.data_hists.values())[0])))
        if len(self.mc_hists.keys()) > 0 and type(mc_hist) != type(list(self.data_hists.values())[0]):
            raise ValueError("Cannot add {} to data histogram dictionary of {}".format(type(hist),type(list(self.mc_hists.values())[0])))

        self.data_hists[name] = data_hist
        self.mc_hists[name] = mc_hist

        # ----- Set Electron Neutrino Histograms ----- #
        if "nue_selection" in name:
            self.nue_hists[name] = mc_hist.Clone()
        elif "elastic" in name:
            if 'electron_'+name not in self.nueel_flavors.keys():
                raise ValueError("{} has not been added to nueel flavor dictionary. Do this before AddingHistogram()".format('electron_'+name))
            self.nue_hists[name] = self.nueel_flavors['electron_'+name]
        elif "ratio" in name:
            self.nue_hists[name] = self.nue_hists[name[:4]+'nue_selection'].Clone()
        else:
            self.nue_hists[name] = self.EmptyHist(mc_hist)

        # ----- Set Muon Neutrino Histograms ----- #
        if "numu_selection" in name or "imd" in name:
            self.numu_hists[name] = mc_hist.Clone()
        elif "elastic" in name:
            if 'muon_'+name not in self.nueel_flavors.keys():
                raise ValueError("{} has not been added to nueel flavor dictionary. Do this before AddingHistogram()".format('muon_'+name))
            self.numu_hists[name] = self.nueel_flavors['muon_'+name]
        elif "ratio" in name:
            self.numu_hists[name] = self.numu_hists[name[:4]+'numu_selection'].Clone()
        else:
            self.numu_hists[name] = self.EmptyHist(mc_hist)

        # ----- Set Flavor Swapped Neutrino Histograms ----- #
        if "nue_selection" in name:
            if name not in self.swap_hists.keys():
                raise ValueError("{} has not been added to swapped flavor dictionary. Do this before AddingHistogram()".format(name))
            self.swap_hists[name] = self.swap_hists[name]
        elif "ratio" in name:
            self.swap_hists[name] = self.swap_hists[name[:4]+'nue_selection'].Clone()
        elif "elastic" in name:
            self.swap_hists[name] = self.nueel_flavors['muon_'+name]
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

        del self.data_hists[name]
        del self.mc_hists[name]
        del self.nue_hists[name]
        del self.numu_hists[name]
        del self.swap_hists[name]
        del self.elastic_ids[name]
        del self.ratio_ids[name]
        del self.beam_ids[name]
        
        if len(list(self.nue_templates.keys())) > 0:
            if name not in list(self.nue_templates.keys()):
                raise ValueError("{} not in nue_templates dictionary".format(name))
            del self.nue_templates[name]

            if name not in list(self.numu_templates.keys()):
                raise ValueError("{} not in numu_templates dictionary".format(name))
            del self.numu_templates[name]

            if name not in list(self.swap_templates.keys()):
                raise ValueError("{} not in swap_templates dictionary".format(name))
            del self.swap_templates[name]

        self.UpdateKeys()

    def UpdateKeys(self):
        self.keys = list(self.mc_hists.keys())
        if self.keys != list(self.data_hists.keys()):
            raise ValueError("MC dictionary incompatable with data dictionary")
        if self.keys != list(self.nue_hists.keys()):
            raise ValueError("MC dictionary incompatable with nue hists")
        if self.keys != list(self.numu_hists.keys()):
            raise ValueError("MC dictionary incompatable with numu hists")

    def MakeRatio(self,beam):
        beam = beam.lower()
        numuname = beam+"_numu_selection"
        elecname = beam+"_nue_selection"

        if numuname in self.keys and elecname in self.keys:
            mc_clone = self.mc_hists[numuname].Clone()
            mc_clone.Divide(mc_clone,self.mc_hists[elecname])
            data_clone = self.data_hists[numuname].Clone()
            data_clone.Divide(data_clone,self.data_hists[elecname])

            self.AddHistograms('{}_ratio'.format(beam),mc_clone,data_clone)
            self.AddTemplates("{}_ratio".format(beam),numu=self.numu_templates[beam+"_numu_selection"],nue=self.nue_templates[beam+"_nue_selection"],swap=self.swap_templates[beam+"_nue_selection"])
            self.RemoveHistograms(elecname)

    def ApplyExclusion(self,exclude):
        for h in self.keys:
            if "fhc" in exclude:
                if "fhc" in h and "selection" in h:
                    self.RemoveHistograms(h)
                    self.RemoveHistograms(h)
            if "rhc" in exclude:
                if "rhc" in h and "selection" in h:
                    self.RemoveHistograms(h)
                    self.RemoveHistograms(h)
            if "numu" in exclude and "numu" in h:
                self.RemoveHistograms(h)
                self.RemoveHistograms(h)
            if "nue" in exclude and "nue" in h:
                self.RemoveHistograms(h)
                self.RemoveHistograms(h)
            if "elastic" in exclude and "elastic" in h:
                self.RemoveHistograms(h)
                self.RemoveHistograms(h)
            if "imd" in exclude:
                self.RemoveHistograms(h)
                self.RemoveHistograms(h)

    def CleanErrorBands(self,names=[]):
        self.LateralToVertical()
        self.SeparateBeamAngle()
        self.SyncErrorBands()

        for errname in names:
            for h in self.keys:
                if errname in self.data_hists[h].GetVertErrorBandNames():
                    self.data_hists[h].PopVertErrorBand(errname)
            for h in self.keys:
                if errname in self.mc_hists[h].GetVertErrorBandNames():
                    self.mc_hists[h].PopVertErrorBand(errname)

        if "fhc_nueel" in self.keys and "rhc_nueel" in self.keys:
            for i in range(self.data_hists["fhc_nueel"].GetNbinsX()+1):
                h_temp = self.data_hists['rhc_nueel']
                h_cont = h_temp.GetBinContent(i)
                h_cont = h_cont if h_cont != 0 else 1
                ratio1 = h_temp.GetVertErrorBand("ElectronScale").GetHist(0).GetBinContent(i)/h_cont
                ratio2 = h_temp.GetVertErrorBand("ElectronScale").GetHist(1).GetBinContent(i)/h_cont
                h_fix = self.data_hists['fhc_nueel']
                h_fix.GetVertErrorBand("ElectronScale").GetHist(0).SetBinContent(i,h_fix.GetBinContent(i) * ratio1)
                h_fix.GetVertErrorBand("ElectronScale").GetHist(1).SetBinContent(i,h_fix.GetBinContent(i) * ratio2)

        self.is_processed = True

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
            lateral(self.nue_hists[h])
            lateral(self.numu_hists[h])
            lateral(self.swap_hists[h])

    def Stitch(self):
        # ----- Create empty ROOT histograms to fill with stitched content ----- #
        n_bins_new = 0
        for h in self.keys:
            for i in range(0,self.mc_hists[h].GetNbinsX()+1):
                if self.mc_hists[h].GetBinContent(i) > minBinCont:
                    n_bins_new+=1

        self.mc_hist   = ROOT.TH1D(self.name+"_mc",self.name+"_mc",n_bins_new,0,n_bins_new)
        self.data_hist = ROOT.TH1D(self.name+"_data",self.name+"_data",n_bins_new,0,n_bins_new)
        self.pseudo_hist = ROOT.TH1D(self.name+"_pseudo",self.name+"_pseudo",n_bins_new,0,n_bins_new)
        
        self.nue_hist = ROOT.TH1D(self.name+"_nue",self.name+"_nue",n_bins_new,0,n_bins_new)
        self.numu_hist = ROOT.TH1D(self.name+"_numu",self.name+"_numu",n_bins_new,0,n_bins_new)
        self.swap_hist = ROOT.TH1D(self.name+"_swap",self.name+"_swap",n_bins_new,0,n_bins_new)

        self.beam_id = ROOT.TH1D(self.name+"_beamID",self.name+"_beamID",n_bins_new,0,n_bins_new)
        self.elastic_id = ROOT.TH1D(self.name+"_nueelID",self.name+"_nueelID",n_bins_new,0,n_bins_new)
        self.ratio_id = ROOT.TH1D(self.name+"_ratioID",self.name+"_ratioID",n_bins_new,0,n_bins_new)

        self.nue_template  = ROOT.TH2D(self.name+"_nue_template",self.name+"_nue_template",n_bins_new,0,n_bins_new,34,0,0.495)
        self.numu_template  = ROOT.TH2D(self.name+"_numu_template",self.name+"_numu_template",n_bins_new,0,n_bins_new,34,0,0.495)
        self.swap_template  = ROOT.TH2D(self.name+"_swap_template",self.name+"_swap_template",n_bins_new,0,n_bins_new,34,0,0.495)

        # ----- Do some errorband cleaning ----- #
        if not self.is_processed:
            self.SeparateBeamAngle() # make beam angle systematics consistent across samples
            self.LateralToVertical() # convert lateral errorbands (deprecated) from old samples to vertical
            self.SyncErrorBands()  # make sure all samples have the same errorbands, reduce flux universes to 100
            self.is_processed = True

        # ----- Combine samples to one histogram ----- #
        self.StitchThis()   # combine 1D histograms
        if len(list(self.nue_templates.keys())) > 0:
            self.StitchThis2D() # combine 2D templates
        
        # ----- Convert ROOT.TH1 to PlotUtils.Mnv1D and fill errorbands ----- #
        self.mc_hist = PlotUtils.MnvH1D(self.mc_hist)
        self.data_hist = PlotUtils.MnvH1D(self.data_hist)
        self.pseudo_hist = PlotUtils.MnvH1D(self.pseudo_hist)

        self.nue_hist = PlotUtils.MnvH1D(self.nue_hist)
        self.numu_hist = PlotUtils.MnvH1D(self.numu_hist)
        self.swap_hist = PlotUtils.MnvH1D(self.swap_hist)

        self.mc_hist.AddMissingErrorBandsAndFillWithCV(self.mc_hists[self.keys[0]])
        self.data_hist.AddMissingErrorBandsAndFillWithCV(self.mc_hists[self.keys[0]])
        self.pseudo_hist.AddMissingErrorBandsAndFillWithCV(self.mc_hists[self.keys[0]])

        self.nue_hist.AddMissingErrorBandsAndFillWithCV(self.mc_hists[self.keys[0]])
        self.numu_hist.AddMissingErrorBandsAndFillWithCV(self.mc_hists[self.keys[0]])
        self.swap_hist.AddMissingErrorBandsAndFillWithCV(self.mc_hists[self.keys[0]])

        self.FillErrorBandsFromDict()

        self.mc_cov = np.asarray(matrix(self.mc_hist.GetTotalErrorMatrix(True, False, False)))[1:-1,1:-1]
        self.data_cov = np.asarray(matrix(self.data_hist.GetTotalErrorMatrix(True, False, False)))[1:-1,1:-1]
        self.pseudo_cov = np.asarray(matrix(self.pseudo_hist.GetTotalErrorMatrix(True, False, False)))[1:-1,1:-1]

    def Write(self,filename):
        f = ROOT.TFile(filename,"UPDATE")
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
        f.Close()

    def Load(self,filename):
        f = ROOT.TFile.Open(filename)

        name = "sample"
        self.mc_hist   = f.Get(name+"_mc")
        self.data_hist = f.Get(name+"_data")
        self.pseudo_hist = f.Get(name+"_pseudo")
        
        self.nue_hist = f.Get(name+"_nue")
        self.numu_hist = f.Get(name+"_numu")
        self.swap_hist = f.Get(name+"_swap")

        self.beam_id = f.Get(name+"_beamID")
        self.elastic_id = f.Get(name+"_nueelID")
        self.ratio_id = f.Get(name+"_ratioID")

        self.nue_template   = f.Get(name+"_nue_template")
        self.numu_template  = f.Get(name+"_numu_template")
        self.swap_template  = f.Get(name+"_swap_template")
        f.Close()

        self.mc_cov = np.asarray(matrix(self.mc_hist.GetTotalErrorMatrix(True, False, False)))[1:-1,1:-1]
        self.data_cov = np.asarray(matrix(self.data_hist.GetTotalErrorMatrix(True, False, False)))[1:-1,1:-1]

    def StitchThis2D(self):
        i_new = 0
        for h in self.keys:
            h_nue = self.nue_templates[h]
            h_numu = self.numu_templates[h]
            h_swap = self.swap_templates[h]
            for x in range(1, self.mc_hists[h].GetNbinsX()+1):
                if self.mc_hists[h].GetBinContent(x) <= minBinCont:
                    continue

                i_new += 1
                colInt = 0
                for c in range(1,h_nue.GetNbinsX()+1):
                    if isinstance(h_nue,ROOT.TH1D):
                        colInt+=h_nue.GetBinContent(c)
                    else:
                        colInt+=h_nue.GetBinContent(c,x)

                for c in range(1,h_nue.GetNbinsX()+1):
                    if isinstance(h_nue,ROOT.TH1D):
                        bin_c = h_nue.GetBinContent(c)
                    else:
                        bin_c = h_nue.GetBinContent(c,x)
                    self.nue_template.SetBinContent(i_new, c, bin_c/colInt if colInt > 0 else 0)
                colInt = 0
                for c in range(1,h_numu.GetNbinsX()+1):
                    if isinstance(h_numu,ROOT.TH1D):
                        colInt+=h_numu.GetBinContent(c)
                    else:
                        colInt+=h_numu.GetBinContent(c,x)

                for c in range(1,h_numu.GetNbinsX()+1):
                    if isinstance(h_numu,ROOT.TH1D):
                        bin_c = h_numu.GetBinContent(c)
                    else:
                        bin_c = h_numu.GetBinContent(c,x)
                    self.numu_template.SetBinContent(i_new, c, bin_c/colInt if colInt > 0 else 0)
                colInt = 0
                for c in range(1,h_swap.GetNbinsX()+1):
                    if isinstance(h_swap,ROOT.TH1D):
                        colInt+=h_swap.GetBinContent(c)
                    else:
                        colInt+=h_swap.GetBinContent(c,x)

                for c in range(1,h_swap.GetNbinsX()+1):
                    if isinstance(h_swap,ROOT.TH1D):
                        bin_c = h_swap.GetBinContent(c)
                    else:
                        bin_c = h_swap.GetBinContent(c,x)
                    self.swap_template.SetBinContent(i_new, c, bin_c/colInt if colInt > 0 else 0)

    def StitchThis(self):
        i_new = 0
        for h in self.keys:
            h_mc = self.mc_hists[h]
            h_data = self.data_hists[h]

            for i in range(1,h_mc.GetNbinsX()+1):
                if h_mc.GetBinContent(i) <= minBinCont:
                    continue # skip empty MC bins
                i_new += 1

                # ----- do MC stitching ----- #
                bin_c = h_mc.GetBinContent(i)
                self.mc_hist.SetBinContent(i_new,bin_c)
                self.pseudo_hist.SetBinContent(i_new,bin_c)

                # no statistical error on elastic scattering special production
                if 'elastic' in h or 'imd':
                    self.mc_hist.SetBinError(i_new,0)
                    self.nue_hist.SetBinError(i_new,0)
                    self.numu_hist.SetBinError(i_new,0)

                # ----- do data stitching ----- #
                bin_c = h_data.GetBinContent(i)
                self.data_hist.SetBinContent(i_new,bin_c)

                # ----- do other stitching ----- #
                self.nue_hist.SetBinContent(i_new,self.nue_hists[h].GetBinContent(i))
                self.numu_hist.SetBinContent(i_new,self.numu_hists[h].GetBinContent(i))
                self.swap_hist.SetBinContent(i_new,self.swap_hists[h].GetBinContent(i))

                self.beam_id.SetBinContent(i_new,self.beam_ids[h].GetBinContent(i))
                self.ratio_id.SetBinContent(i_new,self.ratio_ids[h].GetBinContent(i))
                self.elastic_id.SetBinContent(i_new,self.elastic_ids[h].GetBinContent(i))

    def FillErrorBandsFromDict(self):
        offset = 1
        for h in self.keys:
            h_mc = self.mc_hists[h]
            h_data = self.data_hists[h]
            Nbins = h_mc.GetNbinsX()+1

            if not self.is_processed:
                raise ValueError("Histograms have not been synced")

            errorband_names_vert = h_mc.GetVertErrorBandNames()

            n_univ = 0
            sys_mc = 0.0
            sys_data = 0.0

            for error_band in errorband_names_vert:
                n_univ = h_mc.GetVertErrorBand(error_band).GetNHists()
                if not self.mc_hist.HasVertErrorBand(error_band) and h_mc.HasVertErrorBand(error_band):
                    raise ValueError("MC histograms were not properly synchronized")
                if not self.data_hist.HasVertErrorBand(error_band) and h_data.HasVertErrorBand(error_band):
                    raise ValueError("Data histograms were not properly synchronized")

                for universe in range(0, n_univ):
                    bin_offset = offset
                    for b in range(1, Nbins):
                        bin_mc = h_mc.GetBinContent(b)
                        bin_data = h_data.GetBinContent(b)

                        if bin_mc <= minBinCont:
                            bin_offset += -1
                            continue

                        sys_mc = h_mc.GetVertErrorBand(error_band).GetHist(universe).GetBinContent(b)
                        sys_data = h_data.GetVertErrorBand(error_band).GetHist(universe).GetBinContent(b)

                        sys_nue = self.nue_hists[h].GetVertErrorBand(error_band).GetHist(universe).GetBinContent(b)
                        sys_numu = self.numu_hists[h].GetVertErrorBand(error_band).GetHist(universe).GetBinContent(b)
                        sys_swap = self.swap_hists[h].GetVertErrorBand(error_band).GetHist(universe).GetBinContent(b)

                        ratio = sys_data/bin_data if bin_data != 0 else 0
                        sys_pseudo = bin_mc*ratio

                        bin_new = b + bin_offset - 1

                        self.mc_hist.GetVertErrorBand(error_band).GetHist(universe).SetBinContent(bin_new,sys_mc)
                        self.data_hist.GetVertErrorBand(error_band).GetHist(universe).SetBinContent(bin_new,sys_data)
                        self.pseudo_hist.GetVertErrorBand(error_band).GetHist(universe).SetBinContent(bin_new,sys_pseudo)

                        self.nue_hist.GetVertErrorBand(error_band).GetHist(universe).SetBinContent(bin_new,sys_nue)
                        self.numu_hist.GetVertErrorBand(error_band).GetHist(universe).SetBinContent(bin_new,sys_numu)
                        self.swap_hist.GetVertErrorBand(error_band).GetHist(universe).SetBinContent(bin_new,sys_swap)

            for i in range(1,Nbins):
                if self.mc_hists[h].GetBinContent(i) <= minBinCont:
                    continue
                offset+=1

    def SyncHistograms(self,h_sync):
        if type(h_sync) == StitchedHistogram:
            mc_sync = h_sync.mc_hist
            data_sync = h_sync.data_hist
            if mc_sync.GetVertErrorBandNames() != self.mc_hist.GetVertErrorBandNames():
                self.Sync(mc_sync,self.mc_hist)
            if data_sync.GetVertErrorBandNames() != self.data_hist.GetVertErrorBandNames():
                self.Sync(data_sync,self.data_hist)
        elif type(h_sync) == PlotUtils.MnvH1D:
            if h_sync.GetVertErrorBandNames() != self.mc_hist.GetVertErrorBandNames():
                self.Sync(h_sync,self.mc_hist)
            if h_sync.GetVertErrorBandNames() != self.data_hist.GetVertErrorBandNames():
                self.Sync(h_sync,self.mc_hist)
        else:
            raise ValueError("Cannot sync {} to MnvH1D or StitchedHistogram".format(type(h_sync)))

    def SyncErrorBands(self,rename_bands=True):
        if rename_bands:
            for h in self.keys:
                self.RenameBands(self.mc_hists[h])
                self.RenameBands(self.data_hists[h])
        for h1 in self.keys:
            for h2 in self.keys:
                if h1 != h2:
                    self.mc_hists[h1].AddMissingErrorBandsAndFillWithCV(self.data_hists[h1])
                    self.mc_hists[h1].AddMissingErrorBandsAndFillWithCV(self.data_hists[h2])
                    self.mc_hists[h1].AddMissingErrorBandsAndFillWithCV(self.mc_hists[h2])
                    self.data_hists[h1].AddMissingErrorBandsAndFillWithCV(self.mc_hists[h1])
                    self.nue_hists[h1].AddMissingErrorBandsAndFillWithCV(self.mc_hists[h1])
                    self.numu_hists[h1].AddMissingErrorBandsAndFillWithCV(self.mc_hists[h1])
                    self.swap_hists[h1].AddMissingErrorBandsAndFillWithCV(self.mc_hists[h1])

            self.ShortFlux(self.mc_hists[h1])
            self.ShortFlux(self.data_hists[h1])
            self.ShortFlux(self.nue_hists[h1])
            self.ShortFlux(self.numu_hists[h1])
            self.ShortFlux(self.swap_hists[h1])

    def ShortFlux(self,h):
        name = "Flux"
        if h.GetVertErrorBand(name).GetNHists() > 100:
            h_hists = h.GetVertErrorBand(name).GetHists()
            h_hists = [h_hists[i] for i in range(100)]
            useSpread = h.GetVertErrorBand(name).GetUseSpreadError()
            errband = h.GetVertErrorBand(name)
            h.PopVertErrorBand(name)
            h.AddVertErrorBand(name,h_hists)
            h.GetVertErrorBand(name).SetUseSpreadError(useSpread)
            for i in range(h.GetNbinsX()+1):
                h.GetVertErrorBand(name).SetBinContent(i,errband.GetBinContent(i))

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
        self.mc_cov = np.asarray(matrix(self.mc_hist.GetTotalErrorMatrix(True, False, False)))[1:-1,1:-1]
        self.data_cov = np.asarray(matrix(self.data_hist.GetTotalErrorMatrix(True, False, False)))[1:-1,1:-1]

    def DebugPlots(self):
        c1 = ROOT.TCanvas()

        self.nue_hist.Draw("hist")
        c1.Print("nue_hist.png")
        self.numu_hist.Draw("hist")
        c1.Print("numu_hist.png")
        self.swap_hist.Draw("hist")
        c1.Print("swap_hist.png")
        self.beam_id.Draw("hist")
        c1.Print("beam_id.png")
        self.ratio_id.Draw("hist")
        c1.Print("ratio_id.png")
        self.elastic_id.Draw("hist")
        c1.Print("elastic_id.png")
        self.nue_template.Draw("colz")
        c1.Print("nue_template.png")
        self.numu_template.Draw("colz")
        c1.Print("numu_template.png")
        self.swap_template.Draw("colz")
        c1.Print("swap_template.png")

        MNVPLOTTER = PlotUtils.MnvPlotter()
        MNVPLOTTER.error_summary_group_map.clear();
        for k,v in CONSOLIDATED_ERROR_GROUPS.items():
            vec = ROOT.vector("std::string")()
            for vs in v :
                vec.push_back(vs)
            MNVPLOTTER.error_summary_group_map[k]= vec

        MNVPLOTTER.DrawErrorSummary(self.data_hist,"TR",True,True,0)
        c1.Print("data_errsummary.png")
        MNVPLOTTER.DrawErrorSummary(self.mc_hist,"TR",True,True,0)
        c1.Print("mc_errsummary.png")

