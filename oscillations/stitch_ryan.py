import os
import copy
from collections import OrderedDict
import logging, sys
import ROOT
import PlotUtils
from tools.Fitters import *
from tools.StitchedHistogram import *
from tools.Helper import *
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
errsToRemove = ["LowQ2Pi","elETracker"]

# Get This from Rob. Thanks Rob.
# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
# Specifically, w/o this, this script seg faults in the case where I try to instantiate FluxReweighterWithWiggleFit w/ nuE constraint set to False for more than one playlist
ROOT.TH1.AddDirectory(False)

binwidthScale = False

def addSignalHists(hist,cates):
    h_tot = hist.hists["Total"]
    if (h_tot):
        h_tot.Reset()
        for group in hist.hists:
            if group in cates:
                if hist.hists[group]:
                    h_tot.Add(hist.hists[group])

    return(h_tot)

def loadSwapFiles(sample,numuSample,sampleName):
    swapDir = "/exp/minerva/data/users/{}/{}".format(os.environ["USER"],sample["directory_tag"]+"_swap")
    playlist = sample["playlist"]
    selectionTag = sample["swap_tag"]
    cates = sample["signal_categories"]

    AnalysisConfig.input_dir = swapDir
    AnalysisConfig.selection_tag = selectionTag
    AnalysisConfig.playlist = playlist

    type_path_map = { "swap":AnalysisConfig.SelectionHistoPath(AnalysisConfig.playlist,False,False)}

    inputDir = "/exp/minerva/data/users/{}/{}".format(os.environ["USER"],numuSample["directory_tag"])

    AnalysisConfig.input_dir = inputDir
    AnalysisConfig.selection_tag = numuSample["selection_tag"]
    AnalysisConfig.playlist = numuSample["playlist"]

    type_path_map["data"] = AnalysisConfig.SelectionHistoPath(AnalysisConfig.playlist,True,False)
    print("Loading swap files for {}".format(sampleName)) 
    swap_file,mc_file,pot_scale,swap_pot,data_pot = Utilities.getSwapFilesAndPOTScale(type_path_map)
    swap_hist = HistHolder(sample["selection_variable"],swap_file,"Signal",True,swap_pot,data_pot)
    swap_template = HistHolder(sample["selection_template"],swap_file,"Signal",True,swap_pot,data_pot)

    preservation_hists  = {}
    for plotName in sample["preservation_templates"]:
        temp = HistHolder(plotName,swap_file,"Signal",True,data_pot,swap_pot)
        temp.POTScale(binwidthScale)
        temp = addSignalHists(temp, cates)
        preservation_hists[plotName] = temp

    swap_hist.POTScale(binwidthScale)
    swap_hist = addSignalHists(swap_hist,cates)
    print(sampleName, "integral: {}".format(swap_hist.Integral()))
    swap_template = addSignalHists(swap_template,cates)
    return swap_hist, swap_template, preservation_hists 

def loadFiles(sample,sampleName=""):
    inputDir = "/exp/minerva/data/users/{}/{}".format(os.environ["USER"],sample["directory_tag"])
    playlist = sample["playlist"]
    selectionTag = sample["selection_tag"]
    cates = sample["signal_categories"]

    AnalysisConfig.input_dir = inputDir
    AnalysisConfig.selection_tag = selectionTag
    AnalysisConfig.playlist = playlist

    print("Loading files for {}".format(sampleName)) 
    type_path_map = { t:AnalysisConfig.SelectionHistoPath(AnalysisConfig.playlist,t =="data",False) for t in AnalysisConfig.data_types}
    data_file,mc_file,pot_scale,data_pot,mc_pot = Utilities.getFilesAndPOTScale(AnalysisConfig.playlist,type_path_map,"MAD",True)
    standPOT = data_pot if data_pot is not None else mc_pot 
    mc_hist = HistHolder(sample["selection_variable"],mc_file,"Signal",True,mc_pot,standPOT)
    data_hist = HistHolder(sample["selection_variable"],data_file,"Signal",False,data_pot,standPOT)
    template_hist = HistHolder(sample["selection_template"],mc_file,"Signal",True,mc_pot,standPOT)

    preservation_hists  = {}
    for plotName in sample["preservation_templates"]:
        temp = HistHolder(plotName,mc_file,"Signal",True,mc_pot,standPOT)
        temp.POTScale(binwidthScale)
        temp = addSignalHists(temp,cates)
        if (temp):
            preservation_hists[plotName] = temp

    data_hist.POTScale(binwidthScale)
    mc_hist.POTScale(binwidthScale)
    template_hist.POTScale(binwidthScale)

    # Don't POTScale a bkg subtracted data histogram, they are alrady scaled
    if "background_tag" in sample:
        AnalysisConfig.bkgTune_tag = sample["background_tag"]
        filename = AnalysisConfig.BackgroundFitPath(AnalysisConfig.playlist, AnalysisConfig.bkgTune_tag, False)
        data_file =ROOT.TFile.Open(filename,"READ")
        data_hist = HistHolder("Background Subbed Data",data_file,"Signal",False,data_pot,standPOT)

    data_hist = data_hist.GetHist()
    mc_hist = addSignalHists(mc_hist,cates)
    print(sampleName, "integral: {}".format(mc_hist.Integral()))
    template_hist = addSignalHists(template_hist,cates)
    #print("data:",data_hist)
    #print("mc:",mc_hist)
    #print("template:",template_hist)
    #print("preservation:",preservation_hists)
    return data_hist, mc_hist, template_hist, preservation_hists

def ChangeEnergyScaleSystematic(cv,universe,oldScale,newScale):
    for i in range(cv.GetNbinsX()+1):
        difference = universe.GetBinContent(i) - cv.GetBinContent(i)
        newOffset = difference * newScale/oldScale
        newContent = cv.GetBinContent(i) + newOffset
        universe.SetBinContent(i,newContent)

def CreateFromSamples(sampleJson):
    fhc_scale_CV = ROOT.TFile.Open("/exp/minerva/data/users/rhowell/nu_e/kin_dist_mcFHC_Electron_Scale_electron_scale_MAD.root").Get("EN4")
    fhc_scale_p1sig = ROOT.TFile.Open("/exp/minerva/data/users/rhowell/nu_e/kin_dist_mcFHC_Electron_Scale_electron_scale_MAD.root").Get("EN4_p1sig")
    fhc_scale_m1sig = ROOT.TFile.Open("/exp/minerva/data/users/rhowell/nu_e/kin_dist_mcFHC_Electron_Scale_electron_scale_MAD.root").Get("EN4_m1sig")

    rhc_scale_CV = ROOT.TFile.Open("/exp/minerva/data/users/rhowell/antinu_e/kin_dist_mcRHC_Electron_Scale_electron_scale_MAD.root").Get("EN4")
    rhc_scale_p1sig = ROOT.TFile.Open("/exp/minerva/data/users/rhowell/antinu_e/kin_dist_mcRHC_Electron_Scale_electron_scale_MAD.root").Get("EN4_p1sig")
    rhc_scale_m1sig = ROOT.TFile.Open("/exp/minerva/data/users/rhowell/antinu_e/kin_dist_mcRHC_Electron_Scale_electron_scale_MAD.root").Get("EN4_m1sig")

    # Change the 2% electron energy scale to 1%; not perfect for edge case of 0
    ChangeEnergyScaleSystematic(fhc_scale_CV,fhc_scale_p1sig,0.2,0.1)
    ChangeEnergyScaleSystematic(fhc_scale_CV,fhc_scale_m1sig,0.2,0.1)
    ChangeEnergyScaleSystematic(rhc_scale_CV,rhc_scale_p1sig,0.2,0.1)
    ChangeEnergyScaleSystematic(rhc_scale_CV,rhc_scale_m1sig,0.2,0.1)

    fhc_scale_p1sig = fhc_scale_p1sig/fhc_scale_CV
    fhc_scale_m1sig = fhc_scale_m1sig/fhc_scale_CV

    rhc_scale_p1sig = rhc_scale_p1sig/rhc_scale_CV
    rhc_scale_m1sig = rhc_scale_m1sig/rhc_scale_CV

    selectionSamples = {}
    with open(sampleJson, "r") as file:
        selectionSamples = json.load(file)
            
    fhc_numu_selection_data, fhc_numu_selection_mc, fhc_numu_selection_template, fhc_numu_preservation_dict = loadFiles(selectionSamples["fhc_ccnumu"],"fhc_ccnumu")
    rhc_numu_selection_data, rhc_numu_selection_mc, rhc_numu_selection_template, rhc_numu_preservation_dict = loadFiles(selectionSamples["rhc_ccnumu"],"rhc_ccnumu")

    fhc_nue_selection_data, fhc_nue_selection_mc, fhc_nue_selection_template, fhc_nue_preservation_dict = loadFiles(selectionSamples["fhc_ccnue"],"fhc_ccnue")
    rhc_nue_selection_data, rhc_nue_selection_mc, rhc_nue_selection_template, rhc_nue_preservation_dict = loadFiles(selectionSamples["rhc_ccnue"],"rhc_ccnue")

    fhc_nue_selection_swap, fhc_nue_selection_swap_template, fhc_nue_swap_preservation_dict = loadSwapFiles(selectionSamples["fhc_ccnue"],selectionSamples["fhc_ccnumu"],"fhc_ccnue_swap")
    rhc_nue_selection_swap, rhc_nue_selection_swap_template, rhc_nue_swap_preservation_dict = loadSwapFiles(selectionSamples["rhc_ccnue"],selectionSamples["rhc_ccnumu"],"rhc_ccnue_swap")
 
    fhc_elastic_template_nue = ROOT.TFile(selectionSamples["fhc_elastic"]["mc"]["template_file"]).Get(selectionSamples["fhc_elastic"]["mc"]["template_hist_prefix"]+"nue")
    fhc_elastic_template_numu = ROOT.TFile(selectionSamples["fhc_elastic"]["mc"]["template_file"]).Get(selectionSamples["fhc_elastic"]["mc"]["template_hist_prefix"]+"numu")
    fhc_elastic_template_anue = ROOT.TFile(selectionSamples["fhc_elastic"]["mc"]["template_file"]).Get(selectionSamples["fhc_elastic"]["mc"]["template_hist_prefix"]+"antinue")
    fhc_elastic_template_anumu = ROOT.TFile(selectionSamples["fhc_elastic"]["mc"]["template_file"]).Get(selectionSamples["fhc_elastic"]["mc"]["template_hist_prefix"]+"antinumu")

    rhc_elastic_template_nue = ROOT.TFile(selectionSamples["rhc_elastic"]["mc"]["template_file"]).Get(selectionSamples["rhc_elastic"]["mc"]["template_hist_prefix"]+"nue")
    rhc_elastic_template_numu = ROOT.TFile(selectionSamples["rhc_elastic"]["mc"]["template_file"]).Get(selectionSamples["rhc_elastic"]["mc"]["template_hist_prefix"]+"numu")
    rhc_elastic_template_anue = ROOT.TFile(selectionSamples["rhc_elastic"]["mc"]["template_file"]).Get(selectionSamples["rhc_elastic"]["mc"]["template_hist_prefix"]+"antinue")
    rhc_elastic_template_anumu = ROOT.TFile(selectionSamples["rhc_elastic"]["mc"]["template_file"]).Get(selectionSamples["rhc_elastic"]["mc"]["template_hist_prefix"]+"antinumu")
    
    fhc_imd_template = ROOT.TFile(selectionSamples["fhc_imd"]["mc"]["template_file"]).Get(selectionSamples["fhc_imd"]["mc"]["template_hist"])
    rhc_imd_template = ROOT.TFile(selectionSamples["rhc_imd"]["mc"]["template_file"]).Get(selectionSamples["rhc_imd"]["mc"]["template_hist"])

    fhc_nue_selection_mc.AddVertErrorBandAndFillWithCV("ElectronScale", 2)
    fhc_nue_selection_swap.AddVertErrorBandAndFillWithCV("ElectronScale", 2)
    sys_p = fhc_nue_selection_mc.GetCVHistoWithError() * fhc_scale_p1sig
    sys_m = fhc_nue_selection_mc.GetCVHistoWithError() * fhc_scale_m1sig
    swap_p = fhc_nue_selection_swap.GetCVHistoWithError() * fhc_scale_p1sig
    swap_m = fhc_nue_selection_swap.GetCVHistoWithError() * fhc_scale_m1sig
    for i in range(fhc_nue_selection_mc.GetNbinsX()+1):
        fhc_nue_selection_mc.GetVertErrorBand("ElectronScale").GetHist(0).SetBinContent(i,sys_p.GetBinContent(i))
        fhc_nue_selection_mc.GetVertErrorBand("ElectronScale").GetHist(1).SetBinContent(i,sys_m.GetBinContent(i))

        fhc_nue_selection_swap.GetVertErrorBand("ElectronScale").GetHist(0).SetBinContent(i,swap_p.GetBinContent(i))
        fhc_nue_selection_swap.GetVertErrorBand("ElectronScale").GetHist(1).SetBinContent(i,swap_m.GetBinContent(i))

    rhc_nue_selection_mc.AddVertErrorBandAndFillWithCV("ElectronScale", 2)
    rhc_nue_selection_swap.AddVertErrorBandAndFillWithCV("ElectronScale", 2)
    sys_p = rhc_nue_selection_mc.GetCVHistoWithError() * rhc_scale_p1sig
    sys_m = rhc_nue_selection_mc.GetCVHistoWithError() * rhc_scale_m1sig
    swap_p = rhc_nue_selection_swap.GetCVHistoWithError() * rhc_scale_p1sig
    swap_m = rhc_nue_selection_swap.GetCVHistoWithError() * rhc_scale_m1sig
    for i in range(rhc_nue_selection_mc.GetNbinsX()+1):
        rhc_nue_selection_mc.GetVertErrorBand("ElectronScale").GetHist(0).SetBinContent(i,sys_p.GetBinContent(i))
        rhc_nue_selection_mc.GetVertErrorBand("ElectronScale").GetHist(1).SetBinContent(i,sys_m.GetBinContent(i))

        rhc_nue_selection_swap.GetVertErrorBand("ElectronScale").GetHist(0).SetBinContent(i,swap_p.GetBinContent(i))
        rhc_nue_selection_swap.GetVertErrorBand("ElectronScale").GetHist(1).SetBinContent(i,swap_m.GetBinContent(i))

    #get FHC electron energy hist
    fhc_elastic_mc = ROOT.TFile.Open(selectionSamples["fhc_elastic"]["mc"]["file"]).Get(selectionSamples["fhc_elastic"]["mc"]["hist"])
    fhc_elastic_data = ROOT.TFile.Open(selectionSamples["fhc_elastic"]["data"]["file"]).Get(selectionSamples["fhc_elastic"]["data"]["hist"])

    #get RHC electron energy hist
    rhc_elastic_mc = ROOT.TFile.Open(selectionSamples["rhc_elastic"]["mc"]["file"]).Get(selectionSamples["rhc_elastic"]["mc"]["hist"])
    rhc_elastic_data = ROOT.TFile.Open(selectionSamples["rhc_elastic"]["data"]["file"]).Get(selectionSamples["rhc_elastic"]["data"]["hist"])

    #get IMD histograms
    fhc_imd_mc = ROOT.TFile.Open(selectionSamples["fhc_imd"]["mc"]["file"]).Get(selectionSamples["fhc_imd"]["mc"]["hist"])
    fhc_imd_data = ROOT.TFile.Open(selectionSamples["fhc_imd"]["data"]["file"]).Get(selectionSamples["fhc_imd"]["data"]["hist"])

    rhc_imd_mc = ROOT.TFile.Open(selectionSamples["rhc_imd"]["mc"]["file"]).Get(selectionSamples["rhc_imd"]["mc"]["hist"])
    rhc_imd_data = ROOT.TFile.Open(selectionSamples["rhc_imd"]["data"]["file"]).Get(selectionSamples["rhc_imd"]["data"]["hist"])

    fhcnueelnue = ROOT.TFile.Open(selectionSamples["fhc_elastic"]["mc"]["pdg_hist_file"]).Get(selectionSamples["fhc_elastic"]["mc"]["pdg_hist_prefix"]+'nue')
    fhcnueelnumu = ROOT.TFile.Open(selectionSamples["fhc_elastic"]["mc"]["pdg_hist_file"]).Get(selectionSamples["fhc_elastic"]["mc"]["pdg_hist_prefix"]+'numu')
    fhcnueelanue = ROOT.TFile.Open(selectionSamples["fhc_elastic"]["mc"]["pdg_hist_file"]).Get(selectionSamples["fhc_elastic"]["mc"]["pdg_hist_prefix"]+'anue')
    fhcnueelanumu = ROOT.TFile.Open(selectionSamples["fhc_elastic"]["mc"]["pdg_hist_file"]).Get(selectionSamples["fhc_elastic"]["mc"]["pdg_hist_prefix"]+'anumu')
    
    rhcnueelnue = ROOT.TFile.Open(selectionSamples["rhc_elastic"]["mc"]["pdg_hist_file"]).Get(selectionSamples["rhc_elastic"]["mc"]["pdg_hist_prefix"]+'nue')
    rhcnueelnumu = ROOT.TFile.Open(selectionSamples["rhc_elastic"]["mc"]["pdg_hist_file"]).Get(selectionSamples["rhc_elastic"]["mc"]["pdg_hist_prefix"]+'numu')
    rhcnueelanue = ROOT.TFile.Open(selectionSamples["rhc_elastic"]["mc"]["pdg_hist_file"]).Get(selectionSamples["rhc_elastic"]["mc"]["pdg_hist_prefix"]+'anue')
    rhcnueelanumu = ROOT.TFile.Open(selectionSamples["rhc_elastic"]["mc"]["pdg_hist_file"]).Get(selectionSamples["rhc_elastic"]["mc"]["pdg_hist_prefix"]+'anumu')

    fhcnueel_preservation_nue = ROOT.TFile.Open(selectionSamples["fhc_elastic"]["mc"]["preservation_hist_file"]).Get(selectionSamples["fhc_elastic"]["mc"]["preservation_hist_prefix"]+'nue')
    fhcnueel_preservation_numu = ROOT.TFile.Open(selectionSamples["fhc_elastic"]["mc"]["preservation_hist_file"]).Get(selectionSamples["fhc_elastic"]["mc"]["preservation_hist_prefix"]+'numu')
    fhcnueel_preservation_anue = ROOT.TFile.Open(selectionSamples["fhc_elastic"]["mc"]["preservation_hist_file"]).Get(selectionSamples["fhc_elastic"]["mc"]["preservation_hist_prefix"]+'anue')
    fhcnueel_preservation_anumu = ROOT.TFile.Open(selectionSamples["fhc_elastic"]["mc"]["preservation_hist_file"]).Get(selectionSamples["fhc_elastic"]["mc"]["preservation_hist_prefix"]+'anumu')
    
    rhcnueel_preservation_nue = ROOT.TFile.Open(selectionSamples["rhc_elastic"]["mc"]["preservation_hist_file"]).Get(selectionSamples["rhc_elastic"]["mc"]["preservation_hist_prefix"]+'nue')
    rhcnueel_preservation_numu = ROOT.TFile.Open(selectionSamples["rhc_elastic"]["mc"]["preservation_hist_file"]).Get(selectionSamples["rhc_elastic"]["mc"]["preservation_hist_prefix"]+'numu')
    rhcnueel_preservation_anue = ROOT.TFile.Open(selectionSamples["rhc_elastic"]["mc"]["preservation_hist_file"]).Get(selectionSamples["rhc_elastic"]["mc"]["preservation_hist_prefix"]+'anue')
    rhcnueel_preservation_anumu = ROOT.TFile.Open(selectionSamples["rhc_elastic"]["mc"]["preservation_hist_file"]).Get(selectionSamples["rhc_elastic"]["mc"]["preservation_hist_prefix"]+'anumu')

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
        f_ratio = f1.GetBinContent(i)/fhc_elastic_mc.GetBinContent(i) if fhc_elastic_mc.GetBinContent(i) != 0 else 1
        r_ratio = r1.GetBinContent(i)/rhc_elastic_mc.GetBinContent(i) if rhc_elastic_mc.GetBinContent(i) != 0 else 1
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
            newBin = fhc_elastic_mc.GetVertErrorBand("Flux").GetHist(u).GetBinContent(i)
            ratio = fhcnueelnumu.GetBinContent(i)/fhc_elastic_mc.GetBinContent(i) if newBin != 0 else 0
            fhcnueelnumu.GetVertErrorBand("Flux").GetHist(u).SetBinContent(i,newBin*ratio)
            fhcnueelnue.GetVertErrorBand("Flux").GetHist(u).SetBinContent(i,newBin*(1-ratio))

            newBin = rhc_elastic_mc.GetVertErrorBand("Flux").GetHist(u).GetBinContent(i)
            ratio = rhcnueelnumu.GetBinContent(i)/rhc_elastic_mc.GetBinContent(i) if newBin != 0 else 0
            rhcnueelnumu.GetVertErrorBand("Flux").GetHist(u).SetBinContent(i,newBin*ratio)
            rhcnueelnue.GetVertErrorBand("Flux").GetHist(u).SetBinContent(i,newBin*(1-ratio))

    fhc_elastic_template_nue.Add(fhc_elastic_template_anue)
    fhc_elastic_template_numu.Add(fhc_elastic_template_anumu)
    rhc_elastic_template_nue.Add(rhc_elastic_template_anue)
    rhc_elastic_template_numu.Add(rhc_elastic_template_anumu)
    
    # ---------------------- Create Stitched CV Histograms -----------------------------
    sample_histogram = StitchedHistogram("sample")
    sample_histogram.SetNFluxUniverses(110)

    sample_histogram.AddScatteringFlavors("electron_fhc_elastic",fhcnueelnue)
    sample_histogram.AddScatteringFlavors("electron_rhc_elastic",rhcnueelnue)
    sample_histogram.AddScatteringFlavors("muon_fhc_elastic",fhcnueelnumu)
    sample_histogram.AddScatteringFlavors("muon_rhc_elastic",rhcnueelnumu)

    # ----- Initialize histogram objects with all samples ----- #
    sample_histogram.AddHistograms('fhc_elastic',fhc_elastic_mc,fhc_elastic_data)
    sample_histogram.AddHistograms('fhc_imd',fhc_imd_mc,fhc_imd_data)
    sample_histogram.AddHistograms('fhc_numu_selection',fhc_numu_selection_mc,fhc_numu_selection_data)
    sample_histogram.AddHistograms('fhc_nue_selection',fhc_nue_selection_mc,fhc_nue_selection_data)

    sample_histogram.AddHistograms('rhc_elastic',rhc_elastic_mc,rhc_elastic_data)
    sample_histogram.AddHistograms('rhc_imd',rhc_imd_mc,rhc_imd_data)
    sample_histogram.AddHistograms('rhc_numu_selection',rhc_numu_selection_mc,rhc_numu_selection_data)
    sample_histogram.AddHistograms('rhc_nue_selection',rhc_nue_selection_mc,rhc_nue_selection_data)

    sample_histogram.AddTemplates("fhc_elastic",nue=fhc_elastic_template_nue,numu=fhc_elastic_template_numu,swap=fhc_elastic_template_numu)
    sample_histogram.AddTemplates("fhc_imd",numu=fhc_imd_template)
    sample_histogram.AddTemplates("fhc_numu_selection",numu=fhc_numu_selection_template)
    sample_histogram.AddTemplates("fhc_nue_selection",nue=fhc_nue_selection_template,swap=fhc_nue_selection_swap_template)

    sample_histogram.AddTemplates("rhc_elastic",nue=rhc_elastic_template_nue,numu=rhc_elastic_template_numu,swap=rhc_elastic_template_numu)
    sample_histogram.AddTemplates("rhc_imd",numu=rhc_imd_template)
    sample_histogram.AddTemplates("rhc_numu_selection",numu=rhc_numu_selection_template)
    sample_histogram.AddTemplates("rhc_nue_selection",nue=rhc_nue_selection_template,swap=rhc_nue_selection_swap_template)

    sample_histogram.AddSwappedSample('fhc_nue_selection',fhc_nue_selection_swap)
    sample_histogram.AddSwappedSample('rhc_nue_selection',rhc_nue_selection_swap)

    sample_histogram.AddPreservationHists("fhc_numu_selection","numu",fhc_numu_preservation_dict)
    sample_histogram.AddPreservationHists("fhc_nue_selection","nue",fhc_nue_preservation_dict)
    sample_histogram.AddPreservationHists("rhc_numu_selection","anumu",rhc_numu_preservation_dict)
    sample_histogram.AddPreservationHists("rhc_nue_selection","anue",rhc_nue_preservation_dict)

    sample_histogram.AddScatteringPreservationHists("fhc_elastic",{"True Energy vs Biased Neutrino Energy":{
        "numu":fhcnueel_preservation_numu,
        "nue":fhcnueel_preservation_nue,
        "anumu":fhcnueel_preservation_anumu,
        "anue":fhcnueel_preservation_anue},"True Energy vs L/E":{
            "numu":fhc_elastic_template_numu,
            "nue":fhc_elastic_template_nue,
            "anumu":fhc_elastic_template_anumu,
            "anue":fhc_elastic_template_anue}})
    sample_histogram.AddScatteringPreservationHists("rhc_elastic",{"True Energy vs Biased Neutrino Energy":{
        "numu":rhcnueel_preservation_numu,
        "nue":rhcnueel_preservation_nue,
        "anumu":rhcnueel_preservation_anumu,
        "anue":rhcnueel_preservation_anue},"True Energy vs L/E":{
            "numu":rhc_elastic_template_numu,
            "nue":rhc_elastic_template_nue,
            "anumu":rhc_elastic_template_anumu,
            "anue":rhc_elastic_template_anue}})

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

    return(sample_histogram)

if __name__ == "__main__":
    sample_histogram = CreateFromSamples("SAMPLE_CONFIG.json")
    sample_histogram.WriteCSVs()
    exit()

    filename = "{}/oscillations/rootfiles/NuE_stitched_hists.root".format(ccnueroot)
    sample_histogram.Write(filename)
    
    statistic = Statistics(sample_histogram)
    chi2,penalty = statistic.Chi2DataMC()

    fluxSolution = statistic.GetFluxFitter().GetFluxSolution()
    
    sample_histogram.PlotStitchedHistogram(fluxSolution,"bin_width_normalized_ratio",chi2,penalty)
