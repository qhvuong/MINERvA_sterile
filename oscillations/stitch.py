import os
import copy
from collections import OrderedDict
import logging, sys
import ROOT
import PlotUtils
from Tools.FitTools import *
from Tools.Histogram import *
from Tools.Helper import *
import numpy as np

import math
from array import array

from config.SignalDef import SWAP_SIGNAL_DEFINITION, SIGNAL_DEFINITION
from config.SystematicsConfig import CONSOLIDATED_ERROR_GROUPS 
from config.AnalysisConfig import AnalysisConfig
from tools import Utilities
from tools.PlotLibrary import HistHolder
ccnueroot = os.environ.get('CCNUEROOT')
MNVPLOTTER = PlotUtils.MnvPlotter()
MNVPLOTTER.error_summary_group_map.clear();
for k,v in CONSOLIDATED_ERROR_GROUPS.items():
    vec = ROOT.vector("std::string")()
    for vs in v :
        vec.push_back(vs)
    MNVPLOTTER.error_summary_group_map[k]= vec
errsToRemove = ["LowQ2Pi"]

# Get This from Rob. Thanks Rob.
# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
# Specifically, w/o this, this script seg faults in the case where I try to instantiate FluxReweighterWithWiggleFit w/ nuE constraint set to False for more than one playlist
ROOT.TH1.AddDirectory(False)

if __name__ == "__main__":
    binwidthScale = AnalysisConfig.binwidth
    if AnalysisConfig.ratio:
        ftag = "ratio"
    else:
        ftag = "stitched"

    fhc_scale_CV = ROOT.TFile.Open("/exp/minerva/data/users/rhowell/nu_e/kin_dist_mcFHC_Electron_Scale_electron_scale_MAD.root").Get("EN4")
    fhc_scale_p1sig = ROOT.TFile.Open("/exp/minerva/data/users/rhowell/nu_e/kin_dist_mcFHC_Electron_Scale_electron_scale_MAD.root").Get("EN4_p1sig")
    fhc_scale_m1sig = ROOT.TFile.Open("/exp/minerva/data/users/rhowell/nu_e/kin_dist_mcFHC_Electron_Scale_electron_scale_MAD.root").Get("EN4_m1sig")
    fhc_scale_p1sig = fhc_scale_p1sig/fhc_scale_CV
    fhc_scale_m1sig = fhc_scale_m1sig/fhc_scale_CV
    rhc_scale_CV = ROOT.TFile.Open("/exp/minerva/data/users/rhowell/antinu_e/kin_dist_mcRHC_Electron_Scale_electron_scale_MAD.root").Get("EN4")
    rhc_scale_p1sig = ROOT.TFile.Open("/exp/minerva/data/users/rhowell/antinu_e/kin_dist_mcRHC_Electron_Scale_electron_scale_MAD.root").Get("EN4_p1sig")
    rhc_scale_m1sig = ROOT.TFile.Open("/exp/minerva/data/users/rhowell/antinu_e/kin_dist_mcRHC_Electron_Scale_electron_scale_MAD.root").Get("EN4_m1sig")
    rhc_scale_p1sig = rhc_scale_p1sig/rhc_scale_CV
    rhc_scale_m1sig = rhc_scale_m1sig/rhc_scale_CV

    type_path_map = {'data':'/exp/minerva/data/users/rhowell/nu_e/kin_dist_dataFHC_Selection_100Univ_thesis_MAD.root','mc':'/exp/minerva/data/users/rhowell/nu_e/kin_dist_mcFHC_Selection_100Univ_thesis_MAD.root'}
    data_file,mc_file,pot_scale,data_pot,mc_pot = Utilities.getFilesAndPOTScale("FHC_Selection_100Univ",type_path_map,"MAD",True)
    standPOT = data_pot if data_pot is not None else mc_pot 
    fhc_nue_selection_mc = ROOT.TFile.Open("/exp/minerva/data/users/rhowell/nu_e/bkgfit_FHC_Selection_100Univ_N4_tune_thesis_MAD.root")
    fhc_nue_selection_mcHold = HistHolder("Predicted MC",fhc_nue_selection_mc,"Signal",True,mc_pot,standPOT)
    h_fhc_nue_selection_data = ROOT.TFile.Open("/exp/minerva/data/users/rhowell/nu_e/bkgfit_FHC_Selection_100Univ_N4_tune_thesis_MAD.root").Get("EN4_data_bkgSubbed")
    fhc_sel = ROOT.TFile.Open("/exp/minerva/data/users/rhowell/nu_e/kin_dist_mcFHC_Selection_100Univ_Genie_thesis_MAD.root")
    fhc_sel_template = HistHolder("Reco Energy vs L/E",fhc_sel,"Signal",True,mc_pot,standPOT)

    type_path_map = {'mc':'/exp/minerva/data/users/rhowell/nu_e_swap/kin_dist_mcFHC_Selection_100Univ_thesis_swap_MAD.root'}
    data_file,mc_file,pot_scale,data_pot,mc_pot = Utilities.getFilesAndPOTScale("FHC_Selection_100Univ",type_path_map,"MAD",True)
    fhc_sel_swap = ROOT.TFile.Open("/exp/minerva/data/users/rhowell/nu_e_swap/kin_dist_mcFHC_Selection_100Univ_thesis_swap_MAD.root")
    fhc_swap_sel_template = HistHolder("Reco Energy vs L/E",fhc_sel_swap,"Signal",True,mc_pot,standPOT)
    fhc_swap_selection_mc = HistHolder("Biased Neutrino Energy",fhc_sel_swap,"Signal",True,mc_pot,standPOT)
    
    type_path_map = {'data':'/exp/minerva/data/users/rhowell/antinu_e/kin_dist_dataRHC_Selection_1000Univ_thesis_MAD.root','mc':'/exp/minerva/data/users/rhowell/antinu_e/kin_dist_mcRHC_Selection_1000Univ_thesis_MAD.root'}
    data_file,mc_file,pot_scale,data_pot,mc_pot = Utilities.getFilesAndPOTScale("RHC_Selection_1000Univ",type_path_map,"MAD",True)
    standPOT = data_pot if data_pot is not None else mc_pot 
    rhc_nue_selection_mc = ROOT.TFile.Open("/exp/minerva/data/users/rhowell/antinu_e/bkgfit_RHC_Selection_1000Univ_N4_tune_thesis_MAD.root")
    rhc_nue_selection_mcHold = HistHolder("Predicted MC",rhc_nue_selection_mc,"Signal",True,mc_pot,standPOT)
    h_rhc_nue_selection_data = ROOT.TFile.Open("/exp/minerva/data/users/rhowell/antinu_e/bkgfit_RHC_Selection_1000Univ_N4_tune_thesis_MAD.root").Get("EN4_data_bkgSubbed")
    rhc_sel = ROOT.TFile.Open("/exp/minerva/data/users/rhowell/antinu_e/kin_dist_mcRHC_Selection_Genie_thesis_MAD.root")
    rhc_sel_template = HistHolder("Reco Energy vs L/E",rhc_sel,"Signal",True,mc_pot,standPOT)

    type_path_map = {'mc':'/exp/minerva/data/users/rhowell/antinu_e_swap/kin_dist_mcRHC_Selection_1000Univ_thesis_swap_MAD.root'}
    data_file,mc_file,pot_scale,data_pot,mc_pot = Utilities.getFilesAndPOTScale("RHC_Selection_1000Univ",type_path_map,"MAD",True)
    rhc_sel_swap = ROOT.TFile.Open("/exp/minerva/data/users/rhowell/antinu_e_swap/kin_dist_mcRHC_Selection_1000Univ_thesis_swap_MAD.root")
    rhc_swap_sel_template = HistHolder("Reco Energy vs L/E",rhc_sel_swap,"Signal",True,mc_pot,standPOT)
    rhc_swap_selection_mc = HistHolder("Biased Neutrino Energy",rhc_sel_swap,"Signal",True,mc_pot,standPOT)

    ### NuMu Selection ###
    cates = ["CCQE","CCDelta","CCDIS","CC2p2h","CCOther","CCWrongSign"]
    
    type_path_map = {'data':'/exp/minerva/data/users/rhowell/nu_mu/kin_dist_dataFHC_Selection_100Univ_thesis_muon_MAD.root','mc':'/exp/minerva/data/users/rhowell/nu_mu/kin_dist_mcFHC_Selection_100Univ_thesis_muon_MAD.root'}
    data_file,mc_file,pot_scale,data_pot,mc_pot = Utilities.getFilesAndPOTScale("FHC_Selection_100Univ",type_path_map,"MAD",True)
    standPOT = data_pot if data_pot is not None else mc_pot 
    fhc_numu_selection_mcHold = HistHolder("Biased Neutrino Energy",mc_file,"Signal",True,mc_pot,standPOT)
    fhc_numu_selection_dataHold = HistHolder("Biased Neutrino Energy",data_file,"Signal",False,data_pot,standPOT)
    fhc_musel_template = HistHolder("Reco Energy vs L/E",mc_file,"Signal",True,mc_pot,standPOT)

    type_path_map = {'data':'/exp/minerva/data/users/rhowell/antinu_mu/kin_dist_dataRHC_Selection_1000Univ_thesis_muon_MAD.root','mc':'/exp/minerva/data/users/rhowell/antinu_mu/kin_dist_mcRHC_Selection_1000Univ_thesis_muon_MAD.root'}
    data_file,mc_file,pot_scale,data_pot,mc_pot = Utilities.getFilesAndPOTScale("RHC_Selection_1000Univ",type_path_map,"MAD",True)
    standPOT = data_pot if data_pot is not None else mc_pot 
    rhc_numu_selection_mcHold = HistHolder("Biased Neutrino Energy",mc_file,"Signal",True,mc_pot,standPOT)
    rhc_numu_selection_dataHold = HistHolder("Biased Neutrino Energy",data_file,"Signal",False,data_pot,standPOT)
    rhc_musel_template = HistHolder("Reco Energy vs L/E",mc_file,"Signal",True,mc_pot,standPOT)

    h2_fhc_elastic_template_nue = ROOT.TFile('/exp/minerva/data/users/rhowell/nu_e/kin_dist_mcFHC_scattering_scatter_MAD.root').Get('nue_scattering_template_NoCuts_nue')
    h2_fhc_elastic_template_numu = ROOT.TFile('/exp/minerva/data/users/rhowell/nu_e/kin_dist_mcFHC_scattering_scatter_MAD.root').Get('nue_scattering_template_NoCuts_numu')
    h2_fhc_elastic_template_anue = ROOT.TFile('/exp/minerva/data/users/rhowell/nu_e/kin_dist_mcFHC_scattering_scatter_MAD.root').Get('nue_scattering_template_NoCuts_antinue')
    h2_fhc_elastic_template_anumu = ROOT.TFile('/exp/minerva/data/users/rhowell/nu_e/kin_dist_mcFHC_scattering_scatter_MAD.root').Get('nue_scattering_template_NoCuts_antinumu')

    h_fhc_imd_template = ROOT.TFile('/exp/minerva/data/users/rhowell/nu_mu/FHC_imd_scattering_template.root').Get('imd_scattering_template')
    h_rhc_imd_template = ROOT.TFile('/exp/minerva/data/users/rhowell/antinu_mu/RHC_imd_scattering_template.root').Get('imd_scattering_template')

    h2_rhc_elastic_template_nue = ROOT.TFile('/exp/minerva/data/users/rhowell/antinu_e/kin_dist_mcRHC_scattering_scatter_MAD.root').Get('nue_scattering_template_NoCuts_nue')
    h2_rhc_elastic_template_numu = ROOT.TFile('/exp/minerva/data/users/rhowell/antinu_e/kin_dist_mcRHC_scattering_scatter_MAD.root').Get('nue_scattering_template_NoCuts_numu')
    h2_rhc_elastic_template_anue = ROOT.TFile('/exp/minerva/data/users/rhowell/antinu_e/kin_dist_mcRHC_scattering_scatter_MAD.root').Get('nue_scattering_template_NoCuts_antinue')
    h2_rhc_elastic_template_anumu = ROOT.TFile('/exp/minerva/data/users/rhowell/antinu_e/kin_dist_mcRHC_scattering_scatter_MAD.root').Get('nue_scattering_template_NoCuts_antinumu')

    h_fhc_swap_selection = fhc_swap_selection_mc.hists["Total"].Clone("sigTotal") 
    h_rhc_swap_selection = rhc_swap_selection_mc.hists["Total"].Clone("sigTotal") 

    h_fhc_numu_selection_mc = fhc_numu_selection_mcHold.GetHist()
    h_rhc_numu_selection_mc = rhc_numu_selection_mcHold.GetHist()
    h_fhc_numu_selection_data = fhc_numu_selection_dataHold.GetHist()
    h_rhc_numu_selection_data = rhc_numu_selection_dataHold.GetHist()
    h2_fhc_musel_template = fhc_musel_template.GetHist()
    h2_rhc_musel_template = rhc_musel_template.GetHist()

    h2_fhc_template = fhc_sel_template.hists["Total"].Clone("sigTotal")
    h2_rhc_template = rhc_sel_template.hists["Total"].Clone("sigTotal")
    h2_fhc_swap_template = fhc_swap_sel_template.hists["Total"].Clone("sigTotal")
    h2_rhc_swap_template = rhc_swap_sel_template.hists["Total"].Clone("sigTotal")

    fhc_numu_selection_mcHold.POTScale(binwidthScale)
    rhc_numu_selection_mcHold.POTScale(binwidthScale)
    fhc_numu_selection_dataHold.POTScale(binwidthScale)
    rhc_numu_selection_dataHold.POTScale(binwidthScale)
    fhc_swap_selection_mc.POTScale(binwidthScale)
    rhc_swap_selection_mc.POTScale(binwidthScale)

    h_fhc_swap_selection.Reset()
    h_rhc_swap_selection.Reset()

    h2_fhc_template.Reset()
    h2_rhc_template.Reset()
    h2_fhc_swap_template.Reset()
    h2_rhc_swap_template.Reset()

    h_fhc_numu_selection_mc.Reset()
    h_rhc_numu_selection_mc.Reset()
    h2_fhc_musel_template.Reset()
    h2_rhc_musel_template.Reset()

    for group in fhc_sel_template.hists:
        if group in SIGNAL_DEFINITION:
            if fhc_sel_template.hists[group]:
                h2_fhc_template.Add(fhc_sel_template.hists[group])

    for group in rhc_sel_template.hists:
        if group in SIGNAL_DEFINITION:
            if rhc_sel_template.hists[group]:
                h2_rhc_template.Add(rhc_sel_template.hists[group])

    for group in fhc_swap_sel_template.hists:
        if group in SWAP_SIGNAL_DEFINITION:
            if fhc_swap_sel_template.hists[group]:
                h_fhc_swap_selection.Add(fhc_swap_selection_mc.hists[group])
                h2_fhc_swap_template.Add(fhc_swap_sel_template.hists[group])

    for group in rhc_swap_sel_template.hists:
        if group in SWAP_SIGNAL_DEFINITION:
            if rhc_swap_sel_template.hists[group]:
                h_rhc_swap_selection.Add(rhc_swap_selection_mc.hists[group])
                h2_rhc_swap_template.Add(rhc_swap_sel_template.hists[group])

    for group in fhc_numu_selection_mcHold.hists:
        if group in cates and type(fhc_numu_selection_mcHold.hists[group]) == PlotUtils.MnvH1D:
            h_fhc_numu_selection_mc.Add(fhc_numu_selection_mcHold.hists[group])
            h2_fhc_musel_template.Add(fhc_musel_template.hists[group])

    for group in rhc_numu_selection_mcHold.hists:
        if group in cates and type(rhc_numu_selection_mcHold.hists[group]) == PlotUtils.MnvH1D:
            h_rhc_numu_selection_mc.Add(rhc_numu_selection_mcHold.hists[group])
            h2_rhc_musel_template.Add(rhc_musel_template.hists[group])

    h_fhc_nue_selection_mc = fhc_nue_selection_mcHold.GetHist()
    h_rhc_nue_selection_mc = rhc_nue_selection_mcHold.GetHist()
    
    h_fhc_nue_selection_mc.AddVertErrorBandAndFillWithCV("ElectronScale", 2)
    h_fhc_swap_selection.AddVertErrorBandAndFillWithCV("ElectronScale", 2)
    sys_p = h_fhc_nue_selection_mc.GetCVHistoWithError() * fhc_scale_p1sig
    sys_m = h_fhc_nue_selection_mc.GetCVHistoWithError() * fhc_scale_m1sig
    swap_p = h_fhc_swap_selection.GetCVHistoWithError() * fhc_scale_p1sig
    swap_m = h_fhc_swap_selection.GetCVHistoWithError() * fhc_scale_m1sig
    for i in range(h_fhc_nue_selection_mc.GetNbinsX()+1):
        h_fhc_nue_selection_mc.GetVertErrorBand("ElectronScale").GetHist(0).SetBinContent(i,sys_p.GetBinContent(i))
        h_fhc_nue_selection_mc.GetVertErrorBand("ElectronScale").GetHist(1).SetBinContent(i,sys_m.GetBinContent(i))

        h_fhc_swap_selection.GetVertErrorBand("ElectronScale").GetHist(0).SetBinContent(i,swap_p.GetBinContent(i))
        h_fhc_swap_selection.GetVertErrorBand("ElectronScale").GetHist(1).SetBinContent(i,swap_m.GetBinContent(i))

    h_rhc_nue_selection_mc.AddVertErrorBandAndFillWithCV("ElectronScale", 2)
    h_rhc_swap_selection.AddVertErrorBandAndFillWithCV("ElectronScale", 2)
    sys_p = h_rhc_nue_selection_mc.GetCVHistoWithError() * rhc_scale_p1sig
    sys_m = h_rhc_nue_selection_mc.GetCVHistoWithError() * rhc_scale_m1sig
    swap_p = h_rhc_swap_selection.GetCVHistoWithError() * rhc_scale_p1sig
    swap_m = h_rhc_swap_selection.GetCVHistoWithError() * rhc_scale_m1sig
    for i in range(h_rhc_nue_selection_mc.GetNbinsX()+1):
        h_rhc_nue_selection_mc.GetVertErrorBand("ElectronScale").GetHist(0).SetBinContent(i,sys_p.GetBinContent(i))
        h_rhc_nue_selection_mc.GetVertErrorBand("ElectronScale").GetHist(1).SetBinContent(i,sys_m.GetBinContent(i))

        h_rhc_swap_selection.GetVertErrorBand("ElectronScale").GetHist(0).SetBinContent(i,swap_p.GetBinContent(i))
        h_rhc_swap_selection.GetVertErrorBand("ElectronScale").GetHist(1).SetBinContent(i,swap_m.GetBinContent(i))

    #get FHC electron energy hist
    h_fhc_elastic_mc = ROOT.TFile('/exp/minerva/data/users/rhowell/nueel/DeepikaMc_LauraData_Mnv.root').Get('mnv_eff_cor_mc')
    h_fhc_elastic_data = ROOT.TFile('/exp/minerva/data/users/rhowell/nueel/DeepikaMc_LauraData_Mnv.root').Get('mnv_eff_cor_data')

    #get RHC electron energy hist
    h_rhc_elastic_mc = ROOT.TFile("/exp/minerva/data/users/rhowell/nueel/electronE_bkgsub_FullSetFinalv2_nobinnorm.root").Get('mc_fluxonly')   #.Get('mc_bkgsub_effcor')
    h_rhc_elastic_data =  ROOT.TFile("/exp/minerva/data/users/rhowell/nueel/electronE_bkgsub_FullSetFinalv2_nobinnorm.root").Get('data_bkgsub_effcor')

    #get IMD histograms
    h_fhc_imd_mc = ROOT.TFile.Open("rootfiles/IMD/PubResult_v18_Combined.root").Get("FHC_MC")
    h_fhc_imd_data = ROOT.TFile.Open("rootfiles/IMD/PubResult_v18_Combined.root").Get("FHC_Data")

    h_rhc_imd_mc = ROOT.TFile.Open("rootfiles/IMD/PubResult_v18_Combined.root").Get("RHC_MC")
    h_rhc_imd_data = ROOT.TFile.Open("rootfiles/IMD/PubResult_v18_Combined.root").Get("RHC_Data")

    fhcnueelnue = ROOT.TFile.Open("rootfiles/NuEnumberEvents/ElectronEnergySpectrum_FHC_withoutconstraint.root").Get('electron_energy_nue')
    fhcnueelnumu = ROOT.TFile.Open("rootfiles/NuEnumberEvents/ElectronEnergySpectrum_FHC_withoutconstraint.root").Get('electron_energy_numu')
    fhcnueelanue = ROOT.TFile.Open("rootfiles/NuEnumberEvents/ElectronEnergySpectrum_FHC_withoutconstraint.root").Get('electron_energy_anue')
    fhcnueelanumu = ROOT.TFile.Open("rootfiles/NuEnumberEvents/ElectronEnergySpectrum_FHC_withoutconstraint.root").Get('electron_energy_anumu')
    rhcnueelnue = ROOT.TFile.Open("rootfiles/NuEnumberEvents/ElectronEnergySpectrum_RHC_withoutconstraint.root").Get('electron_energy_nue')
    rhcnueelnumu = ROOT.TFile.Open("rootfiles/NuEnumberEvents/ElectronEnergySpectrum_RHC_withoutconstraint.root").Get('electron_energy_numu')
    rhcnueelanue = ROOT.TFile.Open("rootfiles/NuEnumberEvents/ElectronEnergySpectrum_RHC_withoutconstraint.root").Get('electron_energy_anue')
    rhcnueelanumu = ROOT.TFile.Open("rootfiles/NuEnumberEvents/ElectronEnergySpectrum_RHC_withoutconstraint.root").Get('electron_energy_anumu')

    # ---------------------- Fix Flavor Weights -----------------------------
    f1 = fhcnueelnue.Clone()
    f1.Add(fhcnueelanue)
    f1.Add(fhcnueelnumu)
    f1.Add(fhcnueelanumu)
    r1 = rhcnueelnue.Clone()
    r1.Add(rhcnueelanue)
    r1.Add(rhcnueelnumu)
    r1.Add(rhcnueelanumu)
    fhcweight = f1.GetCVHistoWithError()
    rhcweight = r1.GetCVHistoWithError()
    for i in range(0,f1.GetNbinsX()+1):
        f_ratio = f1.GetBinContent(i)/h_fhc_elastic_mc.GetBinContent(i) if h_fhc_elastic_mc.GetBinContent(i) != 0 else 1
        r_ratio = r1.GetBinContent(i)/h_rhc_elastic_mc.GetBinContent(i) if h_rhc_elastic_mc.GetBinContent(i) != 0 else 1
        fhcweight.SetBinContent(i,f_ratio)
        fhcweight.SetBinError(i,0)
        rhcweight.SetBinContent(i,r_ratio)
        rhcweight.SetBinError(i,0)

    fhcnueelnumu.DivideSingle(fhcnueelnumu,fhcweight)
    fhcnueelanumu.DivideSingle(fhcnueelanumu,fhcweight)
    fhcnueelanue.DivideSingle(fhcnueelanue,fhcweight)
    fhcnueelnue.DivideSingle(fhcnueelnue,fhcweight)

    rhcnueelnumu.DivideSingle(rhcnueelnumu,rhcweight)
    rhcnueelanumu.DivideSingle(rhcnueelanumu,rhcweight)
    rhcnueelanue.DivideSingle(rhcnueelanue,rhcweight)
    rhcnueelnue.DivideSingle(rhcnueelnue,rhcweight)

    fhcnueelnue.Add(fhcnueelanue)
    rhcnueelnue.Add(rhcnueelanue)
    fhcnueelnumu.Add(fhcnueelanumu)
    rhcnueelnumu.Add(rhcnueelanumu)

    # ---------------------- Elastic scattering flavor universes ----------------------
    for i in range(0,fhcnueelnumu.GetNbinsX()+1):
        for u in range(fhcnueelnumu.GetVertErrorBand("Flux").GetNHists()):
            newBin = h_fhc_elastic_mc.GetVertErrorBand("Flux").GetHist(u).GetBinContent(i)
            ratio = fhcnueelnumu.GetBinContent(i)/h_fhc_elastic_mc.GetBinContent(i) if newBin != 0 else 0
            fhcnueelnumu.GetVertErrorBand("Flux").GetHist(u).SetBinContent(i,newBin*ratio)
            fhcnueelnue.GetVertErrorBand("Flux").GetHist(u).SetBinContent(i,newBin*(1-ratio))

            newBin = h_rhc_elastic_mc.GetVertErrorBand("Flux").GetHist(u).GetBinContent(i)
            ratio = rhcnueelnumu.GetBinContent(i)/h_rhc_elastic_mc.GetBinContent(i) if newBin != 0 else 0
            rhcnueelnumu.GetVertErrorBand("Flux").GetHist(u).SetBinContent(i,newBin*ratio)
            rhcnueelnue.GetVertErrorBand("Flux").GetHist(u).SetBinContent(i,newBin*(1-ratio))

    h2_fhc_elastic_template_nue.Add(h2_fhc_elastic_template_anue)
    h2_fhc_elastic_template_numu.Add(h2_fhc_elastic_template_anumu)
    h2_rhc_elastic_template_nue.Add(h2_rhc_elastic_template_anue)
    h2_rhc_elastic_template_numu.Add(h2_rhc_elastic_template_anumu)
    
    # ---------------------- Create Stitched CV Histograms -----------------------------
    sample_histogram = StitchedHistogram("sample")

    sample_histogram.AddScatteringFlavors("electron_fhc_elastic",fhcnueelnue)
    sample_histogram.AddScatteringFlavors("electron_rhc_elastic",rhcnueelnue)
    sample_histogram.AddScatteringFlavors("muon_fhc_elastic",fhcnueelnumu)
    sample_histogram.AddScatteringFlavors("muon_rhc_elastic",rhcnueelnumu)

    sample_histogram.AddSwappedSample('fhc_nue_selection',h_fhc_swap_selection)
    sample_histogram.AddSwappedSample('rhc_nue_selection',h_rhc_swap_selection)

    # ----- Initialize histogram objects with all samples ----- #
    sample_histogram.AddHistograms('fhc_elastic',h_fhc_elastic_mc,h_fhc_elastic_data)
    sample_histogram.AddHistograms('fhc_imd',h_fhc_imd_mc,h_fhc_imd_data)
    sample_histogram.AddHistograms('fhc_numu_selection',h_fhc_numu_selection_mc,h_fhc_numu_selection_data)
    sample_histogram.AddHistograms('fhc_nue_selection',h_fhc_nue_selection_mc,h_fhc_nue_selection_data)
    sample_histogram.AddHistograms('rhc_elastic',h_rhc_elastic_mc,h_rhc_elastic_data)
    sample_histogram.AddHistograms('rhc_imd',h_rhc_imd_mc,h_rhc_imd_data)
    sample_histogram.AddHistograms('rhc_numu_selection',h_rhc_numu_selection_mc,h_rhc_numu_selection_data)
    sample_histogram.AddHistograms('rhc_nue_selection',h_rhc_nue_selection_mc,h_rhc_nue_selection_data)

    sample_histogram.AddTemplates("fhc_elastic",nue=h2_fhc_elastic_template_nue,numu=h2_fhc_elastic_template_numu,swap=h2_fhc_elastic_template_numu)
    sample_histogram.AddTemplates("fhc_imd",numu=h_fhc_imd_template)
    sample_histogram.AddTemplates("fhc_numu_selection",numu=h2_fhc_musel_template)
    sample_histogram.AddTemplates("fhc_nue_selection",nue=h2_fhc_template,swap=h2_fhc_swap_template)
    sample_histogram.AddTemplates("rhc_elastic",nue=h2_rhc_elastic_template_nue,numu=h2_rhc_elastic_template_numu,swap=h2_rhc_elastic_template_numu)
    sample_histogram.AddTemplates("rhc_imd",numu=h_rhc_imd_template)
    sample_histogram.AddTemplates("rhc_numu_selection",numu=h2_rhc_musel_template)
    sample_histogram.AddTemplates("rhc_nue_selection",nue=h2_rhc_template,swap=h2_rhc_swap_template)

    # ----- Process Systematics and Synchronize across histograms ----- #
    sample_histogram.CleanErrorBands(errsToRemove)
    old_histogram = copy.deepcopy(sample_histogram)

    if AnalysisConfig.ratio: # do we want to replace selection samples with flavor ratios
        if "fhc" not in AnalysisConfig.exclude: # do we care about the fhc component
            sample_histogram.MakeRatio('fhc')
        if "rhc" not in AnalysisConfig.exclude: # do we care about the rhc component
            sample_histogram.MakeRatio('rhc')

    # ----- Remove samples that we want to exclude from analysis ----- #
    sample_histogram.ApplyExclusion(AnalysisConfig.exclude)

    # ----- Stitch histograms together ----- #
    sample_histogram.Stitch()

    mnv_data = sample_histogram.data_hist.Clone()
    mnv_mc   = sample_histogram.mc_hist.Clone()

    dataprint = np.array(mnv_data)[1:-1] # store MC bin contents excluding over/underflow bins
    mcprint = np.array(mnv_mc)[1:-1]
    np.savetxt("mc_cv.csv",mcprint,delimiter=',')
    np.savetxt("data_cv.csv",dataprint,delimiter=',')
    
    filename = "{}/oscillations/rootfiles/NuE_stitched_hists.root".format(ccnueroot)

    sample_histogram.Write(filename)
    #sample_histogram.SetPlottingStyle()
    sample_histogram.DebugPlots()
    
    invCov=sample_histogram.GetInverseCovarianceMatrix(sansFlux=True)
    nullSolution,nullPen = FluxSolution(sample_histogram,invCov=invCov)

    chi2, penalty = Chi2DataMC(sample_histogram,fluxSolution=nullSolution,invCov=invCov,exclude='ratio',lam=1,marginalize=True)
    chi2-=penalty
    #chi2, penalty = Chi2DataMC(sample_histogram,marginalize=False)
    sample_histogram.PlotStitchedHistogram(nullSolution,"bin_width_normalized_ratio",True,chi2,penalty)

    if True:
        old_histogram.Stitch()
        invCov=old_histogram.GetInverseCovarianceMatrix(sansFlux=True)
        chi2, penalty = Chi2DataMC(old_histogram,fluxSolution=nullSolution,invCov=invCov,exclude='ratio',lam=1,marginalize=True)
        chi2-=penalty
        #chi2, penalty = Chi2DataMC(sample_histogram,marginalize=False)
        old_histogram.PlotStitchedHistogram(nullSolution,"bin_width_normalized_noratio",True,chi2,penalty)

    #sample_histogram.PlotSamples(nullSolution)
    #DataMCCVPlot(mnv_data,mnv_mc,"mc_stitched_v2.png")
