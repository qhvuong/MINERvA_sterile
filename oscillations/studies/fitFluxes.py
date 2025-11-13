import os
from collections import OrderedDict
import argparse
import logging, sys
import ROOT
import PlotUtils
from Tools.FitTools import *
from Tools.PlotTools import *
from Tools.StitchTools import *
from Tools.Histogram import *

import numpy as np
#np.random.seed(0)
#np.set_printoptions(precision=3)
#np.set_printoptions(linewidth=1520)
#np.set_printoptions(threshold=sys.maxsize)

import math
from array import array

#insert path for modules of this package.
#from config import PlotConfig
#from config.AnalysisConfig import AnalysisConfig
#from config.DrawingConfig import PLOTS_TO_MAKE,Default_Plot_Type,Default_Scale,DefaultPlotters,DefaultSlicer
from config.SignalDef import SWAP_SIGNAL_DEFINITION, SIGNAL_DEFINITION
from config.SystematicsConfig import CONSOLIDATED_ERROR_GROUPS 
from tools import Utilities
from tools.PlotLibrary import HistHolder

MNVPLOTTER = PlotUtils.MnvPlotter()
#config MNVPLOTTER:
MNVPLOTTER.draw_normalized_to_bin_width=False
#MNVPLOTTER.extra_top_margin = -.033# go slightly closer to top of pad
MNVPLOTTER.mc_bkgd_color = 46 
MNVPLOTTER.mc_bkgd_line_color = 46
#MNVPLOTTER.mc_line_width = 0

MNVPLOTTER.data_bkgd_color = 12 #gray
MNVPLOTTER.data_bkgd_style = 24 #circle  

#legend entries are closer
MNVPLOTTER.height_nspaces_per_hist = 1.1
MNVPLOTTER.width_xspace_per_letter = .3
MNVPLOTTER.legend_text_size        = .07
MNVPLOTTER.legend_n_columns        = 1
MNVPLOTTER.axis_draw_grid_x = True
MNVPLOTTER.axis_draw_grid_y = True

MNVPLOTTER.error_summary_group_map.clear();
for k,v in CONSOLIDATED_ERROR_GROUPS.items():
    vec = ROOT.vector("std::string")()
    for vs in v :
        vec.push_back(vs)
    MNVPLOTTER.error_summary_group_map[k]= vec

# Get This from Rob. Thanks Rob.
# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
# Specifically, w/o this, this script seg faults in the case where I try to instantiate FluxReweighterWithWiggleFit w/ nuE constraint set to False for more than one playlist
ROOT.TH1.AddDirectory(False)

muonList =  ["Muon_Energy_Resolution","Muon_Energy_MINOS","Muon_Energy_MINERvA"]
errsToRemove = ["LowQ2Pi"]

def SwapUniverseCV(h_in,err,univ):
    if err in h_in.GetVertErrorBandNames():
        ratio = h_in.GetVertErrorBand(err).GetHist(univ).Clone()
    else:
        print("Histogram doesn't have {} universe in errorband {}, skipping...".format(univ,err))
        return(h_in)

    CV = h_in.GetCVHistoWithStatError()
    ratio.Divide(ratio,CV)
    ratio = PlotUtils.MnvH1D(ratio)
    ratio.AddMissingErrorBandsAndFillWithCV(h_in)
    for i in range(0,ratio.GetNbinsX()+1):
        ratio.SetBinError(i,0)

    h_in.Multiply(h_in,ratio)
    h_in.RenameHistosAndErrorBands(ratio.GetName()+"_{}_{}".format(err,univ))
    del ratio

if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--exclude",
                        dest ="exclude",
                        default = "None",
    )
    parser.add_argument("--ratio",
                        dest ="ratio",
                        default = False,
                        action = "store_true",
    )
    parser.add_argument("--fit_muons",
                        dest ="fit_muons",
                        default = False,
                        action = "store_true",
    )

    args = parser.parse_args()
    exclude = str(args.exclude).lower()
    doratio = args.ratio
    fit_muons = args.fit_muons
    if doratio:
        ftag = "ratio"
    else:
        ftag = "stitched"

    # ---------------------- NuEEl Selection Block -----------------------------
    h_fhc_elastic_mc = ROOT.TFile('/exp/minerva/data/users/rhowell/nueel/DeepikaMc_LauraData_Mnv.root').Get('mnv_eff_cor_mc')
    h_fhc_elastic_data = ROOT.TFile('/exp/minerva/data/users/rhowell/nueel/DeepikaMc_LauraData_Mnv.root').Get('mnv_eff_cor_data')

    h_rhc_elastic_mc = ROOT.TFile("/exp/minerva/data/users/rhowell/nueel/electronE_bkgsub_FullSetFinalv2_nobinnorm.root").Get('mc_fluxonly')   #.Get('mc_bkgsub_effcor')
    h_rhc_elastic_data =  ROOT.TFile("/exp/minerva/data/users/rhowell/nueel/electronE_bkgsub_FullSetFinalv2_nobinnorm.root").Get('data_bkgsub_effcor')

    # ---------------------- IMD Selection Block -----------------------------
    h_fhc_imd_mc = ROOT.TFile.Open("IMD/PubResult_v18_Combined.root").Get("FHC_MC")
    h_fhc_imd_data = ROOT.TFile.Open("IMD/PubResult_v18_Combined.root").Get("FHC_Data")

    h_rhc_imd_mc = ROOT.TFile.Open("IMD/PubResult_v18_Combined.root").Get("RHC_MC")
    h_rhc_imd_data = ROOT.TFile.Open("IMD/PubResult_v18_Combined.root").Get("RHC_Data")


    fhc_sel = ROOT.TFile.Open("/exp/minerva/data/users/rhowell/nu_e/kin_dist_mcFHC_Selection_100Univ_Genie_thesis_MAD.root")
    fhc_sel_template = HistHolder("Reco Energy vs L/E",fhc_sel,"Signal",True,1,1)

    rhc_sel = ROOT.TFile.Open("/exp/minerva/data/users/rhowell/antinu_e/kin_dist_mcRHC_Selection_Genie_thesis_MAD.root")
    rhc_sel_template = HistHolder("Reco Energy vs L/E",rhc_sel,"Signal",True,1,1)
    
    # ---------------------- NuMu Selection Block -----------------------------
    type_path_map = {'data':'/exp/minerva/data/users/rhowell/nu_mu/kin_dist_dataFHC_Selection_100Univ_thesis_muon_MAD.root','mc':'/exp/minerva/data/users/rhowell/nu_mu/kin_dist_mcFHC_Selection_100Univ_thesis_muon_MAD.root'}
    data_file,mc_file,pot_scale,data_pot,mc_pot = Utilities.getFilesAndPOTScale("FHC_Selection_100Univ",type_path_map,"MAD",True)
    standPOT = data_pot if data_pot is not None else mc_pot 
    fhc_numu_selection_mcHold = HistHolder("Biased Neutrino Energy",mc_file,"Signal",True,mc_pot,standPOT)
    fhc_numu_selection_dataHold = HistHolder("Biased Neutrino Energy",data_file,"Signal",False,data_pot,standPOT)
    fhc_musel_template = HistHolder("Reco Energy vs L/E",mc_file,"Signal",True,mc_pot,standPOT)

    type_path_map = {'data':'/exp/minerva/data/users/rhowell/antinu_mu/kin_dist_dataRHC_Selection_1000Univ_thesis_muon_MAD.root','mc':'/exp/minerva/data/users/rhowell/antinu_mu/kin_dist_mcRHC_Selection_1000Univ_thesis_muon_MAD.root'}
    data_file,mc_file,pot_scale,data_pot,mc_pot = Utilities.getFilesAndPOTScale("RHC_Selection_100Univ",type_path_map,"MAD",True)
    standPOT = data_pot if data_pot is not None else mc_pot 
    rhc_numu_selection_mcHold = HistHolder("Biased Neutrino Energy",mc_file,"Signal",True,mc_pot,standPOT)
    rhc_numu_selection_dataHold = HistHolder("Biased Neutrino Energy",data_file,"Signal",False,data_pot,standPOT)
    rhc_musel_template = HistHolder("Reco Energy vs L/E",mc_file,"Signal",True,mc_pot,standPOT)

    type_path_map = {'mc':'/exp/minerva/data/users/rhowell/nu_e_swap/kin_dist_mcFHC_Selection_100Univ_thesis_swap_MAD.root'}
    data_file,mc_file,pot_scale,data_pot,mc_pot = Utilities.getFilesAndPOTScale("FHC_Selection_100Univ",type_path_map,"MAD",True)
    fhc_sel_swap = ROOT.TFile.Open("/exp/minerva/data/users/rhowell/nu_e_swap/kin_dist_mcFHC_Selection_100Univ_thesis_swap_MAD.root")
    fhc_swap_sel_template = HistHolder("Reco Energy vs L/E",fhc_sel_swap,"Signal",True,mc_pot,standPOT)
    fhc_swap_selection_mc = HistHolder("Biased Neutrino Energy",fhc_sel_swap,"Signal",True,mc_pot,standPOT)
    
    type_path_map = {'mc':'/exp/minerva/data/users/rhowell/antinu_e_swap/kin_dist_mcRHC_Selection_1000Univ_thesis_swap_MAD.root'}
    data_file,mc_file,pot_scale,data_pot,mc_pot = Utilities.getFilesAndPOTScale("RHC_Selection_1000Univ",type_path_map,"MAD",True)
    rhc_sel_swap = ROOT.TFile.Open("/exp/minerva/data/users/rhowell/antinu_e_swap/kin_dist_mcRHC_Selection_1000Univ_thesis_swap_MAD.root")
    rhc_swap_sel_template = HistHolder("Reco Energy vs L/E",rhc_sel_swap,"Signal",True,mc_pot,standPOT)
    rhc_swap_selection_mc = HistHolder("Biased Neutrino Energy",rhc_sel_swap,"Signal",True,mc_pot,standPOT)
    
    h_fhc_numu_selection_mc = fhc_numu_selection_mcHold.GetHist()
    h_rhc_numu_selection_mc = rhc_numu_selection_mcHold.GetHist()
    h_fhc_numu_selection_data = fhc_numu_selection_dataHold.GetHist()
    h_rhc_numu_selection_data = rhc_numu_selection_dataHold.GetHist()

    fhc_numu_selection_mcHold.POTScale(False)
    rhc_numu_selection_mcHold.POTScale(False)
    fhc_numu_selection_dataHold.POTScale(False)
    rhc_numu_selection_dataHold.POTScale(False)

    h_fhc_numu_selection_mc.Reset()
    h_rhc_numu_selection_mc.Reset()

    cates = ["CCQE","CCDelta","CCDIS","CC2p2h","CCOther","CCWrongSign"]
    for group in fhc_numu_selection_mcHold.hists:
        if group == "Total":
            continue
        elif group in cates:
            h_fhc_numu_selection_mc.Add(fhc_numu_selection_mcHold.hists[group])

    for group in rhc_numu_selection_mcHold.hists:
        if group == "Total":
            continue
        elif group in cates:
            h_rhc_numu_selection_mc.Add(rhc_numu_selection_mcHold.hists[group])

    h2_fhc_musel_template = fhc_musel_template.GetHist()
    h2_rhc_musel_template = rhc_musel_template.GetHist()

    h2_fhc_template = fhc_sel_template.hists["Total"].Clone("sigTotal")
    h2_rhc_template = rhc_sel_template.hists["Total"].Clone("sigTotal")
    h2_fhc_swap_template = fhc_swap_sel_template.hists["Total"].Clone("sigTotal")
    h2_rhc_swap_template = rhc_swap_sel_template.hists["Total"].Clone("sigTotal")

    # ---------------------- NuE CV Selection Block -----------------------------
    fhc_cv = "/exp/minerva/data/users/rhowell/nu_e/bkgfit_FHC_Selection_100Univ_N4_tune_thesis_MAD.root"
    h_fhc_nue_selection_mc   = ROOT.TFile.Open(fhc_cv).Get("EN4_predicted_Signal")
    h_fhc_nue_selection_data = ROOT.TFile.Open(fhc_cv).Get("EN4_data_bkgSubbed")

    rhc_cv = "/exp/minerva/data/users/rhowell/antinu_e/bkgfit_RHC_Selection_100Univ_N4_tune_thesis_MAD.root"
    h_rhc_nue_selection_mc   = ROOT.TFile.Open(rhc_cv).Get("EN4_predicted_Signal")
    h_rhc_nue_selection_data = ROOT.TFile.Open(rhc_cv).Get("EN4_data_bkgSubbed")

    # ---------------------- Add Electron Energy Scale Systematic to Selection -----------------------------
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

    h_fhc_nue_selection_mc.AddVertErrorBandAndFillWithCV("ElectronScale", 2)
    for i in range(h_fhc_nue_selection_mc.GetNbinsX()+1):
        h_fhc_nue_selection_mc.GetVertErrorBand("ElectronScale").GetHist(0).SetBinContent(i,h_fhc_nue_selection_mc.GetBinContent(i) * fhc_scale_m1sig.GetBinContent(i))
        h_fhc_nue_selection_mc.GetVertErrorBand("ElectronScale").GetHist(1).SetBinContent(i,h_fhc_nue_selection_mc.GetBinContent(i) * fhc_scale_p1sig.GetBinContent(i))

    h_rhc_nue_selection_mc.AddVertErrorBandAndFillWithCV("ElectronScale", 2)
    for i in range(h_rhc_nue_selection_mc.GetNbinsX()+1):
        h_rhc_nue_selection_mc.GetVertErrorBand("ElectronScale").GetHist(0).SetBinContent(i,h_rhc_nue_selection_mc.GetBinContent(i) * rhc_scale_m1sig.GetBinContent(i))
        h_rhc_nue_selection_mc.GetVertErrorBand("ElectronScale").GetHist(1).SetBinContent(i,h_rhc_nue_selection_mc.GetBinContent(i) * rhc_scale_p1sig.GetBinContent(i))

    fhcnueelnue = ROOT.TFile.Open("NuEnumberEvents/ElectronEnergySpectrum_FHC_withoutconstraint.root").Get('electron_energy_nue')
    fhcnueelnumu = ROOT.TFile.Open("NuEnumberEvents/ElectronEnergySpectrum_FHC_withoutconstraint.root").Get('electron_energy_numu')
    fhcnueelanue = ROOT.TFile.Open("NuEnumberEvents/ElectronEnergySpectrum_FHC_withoutconstraint.root").Get('electron_energy_anue')
    fhcnueelanumu = ROOT.TFile.Open("NuEnumberEvents/ElectronEnergySpectrum_FHC_withoutconstraint.root").Get('electron_energy_anumu')
    rhcnueelnue = ROOT.TFile.Open("NuEnumberEvents/ElectronEnergySpectrum_RHC_withoutconstraint.root").Get('electron_energy_nue')
    rhcnueelnumu = ROOT.TFile.Open("NuEnumberEvents/ElectronEnergySpectrum_RHC_withoutconstraint.root").Get('electron_energy_numu')
    rhcnueelanue = ROOT.TFile.Open("NuEnumberEvents/ElectronEnergySpectrum_RHC_withoutconstraint.root").Get('electron_energy_anue')
    rhcnueelanumu = ROOT.TFile.Open("NuEnumberEvents/ElectronEnergySpectrum_RHC_withoutconstraint.root").Get('electron_energy_anumu')

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

    h2_fhc_elastic_template_nue.Add(h2_fhc_elastic_template_anue)
    h2_fhc_elastic_template_numu.Add(h2_fhc_elastic_template_anumu)
    h2_rhc_elastic_template_nue.Add(h2_rhc_elastic_template_anue)
    h2_rhc_elastic_template_numu.Add(h2_rhc_elastic_template_anumu)
    
    fhc_universe = "/exp/minerva/data/users/rhowell/nu_e/bkgfit_FHC_Selection_100Univ_Sigmas_N4_tune_thesis_MAD.root"
    rhc_universe = "/exp/minerva/data/users/rhowell/antinu_e/bkgfit_RHC_Selection_100Univ_Sigmas_N4_tune_thesis_MAD.root"

    chi2filename = "chi2_values_{}_{}.txt".format(ftag,exclude.replace(" ","_"))
    with open(chi2filename,'w') as file:
        file.write("Universe,chi2,osc_chi2,dm^2,U_e4^2,U_mu4^2,U_tau4^2".format(i))
        file.write("\n")
        
        # ---------------------- Create Stitched CV Histograms -----------------------------
        cv_histogram = StitchedHistogram("cv_histogram")

        cv_histogram.AddScatteringFlavors("electron_fhc_elastic",fhcnueelnue)
        cv_histogram.AddScatteringFlavors("electron_rhc_elastic",rhcnueelnue)
        cv_histogram.AddScatteringFlavors("muon_fhc_elastic",fhcnueelnumu)
        cv_histogram.AddScatteringFlavors("muon_rhc_elastic",rhcnueelnumu)

        cv_histogram.AddSwappedSample('fhc_nue_selection',h_fhc_swap_selection)
        cv_histogram.AddSwappedSample('rhc_nue_selection',h_rhc_swap_selection)

        # ----- Initialize histogram objects with all samples ----- #
        cv_histogram.AddHistograms('fhc_elastic',h_fhc_elastic_mc,h_fhc_elastic_data)
        cv_histogram.AddHistograms('fhc_imd',h_fhc_imd_mc,h_fhc_imd_data)
        cv_histogram.AddHistograms('fhc_numu_selection',h_fhc_numu_selection_mc,h_fhc_numu_selection_data)
        cv_histogram.AddHistograms('fhc_nue_selection',h_fhc_nue_selection_mc,h_fhc_nue_selection_data)
        cv_histogram.AddHistograms('rhc_elastic',h_rhc_elastic_mc,h_rhc_elastic_data)
        cv_histogram.AddHistograms('rhc_imd',h_rhc_imd_mc,h_rhc_imd_data)
        cv_histogram.AddHistograms('rhc_numu_selection',h_rhc_numu_selection_mc,h_rhc_numu_selection_data)
        cv_histogram.AddHistograms('rhc_nue_selection',h_rhc_nue_selection_mc,h_rhc_nue_selection_data)

        cv_histogram.AddTemplates("fhc_elastic",nue=h2_fhc_elastic_template_nue,numu=h2_fhc_elastic_template_numu,swap=h2_fhc_elastic_template_numu)
        cv_histogram.AddTemplates("fhc_imd",numu=h_fhc_imd_template)
        cv_histogram.AddTemplates("fhc_numu_selection",numu=h2_fhc_musel_template)
        cv_histogram.AddTemplates("fhc_nue_selection",nue=h2_fhc_template,swap=h2_fhc_swap_template)
        cv_histogram.AddTemplates("rhc_elastic",nue=h2_rhc_elastic_template_nue,numu=h2_rhc_elastic_template_numu,swap=h2_rhc_elastic_template_numu)
        cv_histogram.AddTemplates("rhc_imd",numu=h_rhc_imd_template)
        cv_histogram.AddTemplates("rhc_numu_selection",numu=h2_rhc_musel_template)
        cv_histogram.AddTemplates("rhc_nue_selection",nue=h2_rhc_template,swap=h2_rhc_swap_template)
        # ----- Remove samples that we want to exclude from analysis ----- #
        cv_histogram.ApplyExclusion(exclude)

        # ----- Process Systematics and Synchronize across histograms ----- #
        cv_histogram.CleanErrorBands(errsToRemove)

        if False:
            c1 = ROOT.TCanvas()
            h_fhc_elastic_data.Scale(2,"width")
            band = h_fhc_elastic_data.GetVertErrorBand("ElectronScale")
            univ1 = band.GetHist(0)
            univ2 = band.GetHist(1)
            univ1.SetLineColor(ROOT.kRed)
            univ2.SetLineColor(ROOT.kBlue)
            band.SetLineColor(ROOT.kBlack)
            band.SetLineStyle(7)
            h_fhc_elastic_data.SetLineColor(ROOT.kBlack)
            h_fhc_elastic_data.SetTitle("FHC nu+e Scattering Eff. Corr. Data")
            h_fhc_elastic_data.Draw("p")
            band.Draw("hist same")
            univ1.Draw("hist same")
            univ2.Draw("hist same")
            leg = ROOT.TLegend(.6,.6,.85,.85)
            leg.AddEntry(h_fhc_elastic_data,"CV","l")
            leg.AddEntry(band,"Errorband Hist","l")
            leg.AddEntry(univ1,"Universe 0","l")
            leg.AddEntry(univ2,"Universe 1","l")
            leg.Draw()
            c1.Print("fhc_data.png")
            h_rhc_elastic_data.Scale(2,"width")
            band = h_rhc_elastic_data.GetVertErrorBand("ElectronScale")
            univ1 = band.GetHist(0)
            univ2 = band.GetHist(1)
            univ1.SetLineColor(ROOT.kRed)
            univ2.SetLineColor(ROOT.kBlue)
            band.SetLineColor(ROOT.kBlack)
            band.SetLineStyle(7)
            h_rhc_elastic_data.SetLineColor(ROOT.kBlack)
            h_rhc_elastic_data.SetTitle("RHC nu+e Scattering Eff. Corr. Data")
            h_rhc_elastic_data.Draw("p")
            band.Draw("hist same")
            univ1.Draw("hist same")
            univ2.Draw("hist same")
            leg = ROOT.TLegend(.6,.6,.85,.85)
            leg.AddEntry(h_rhc_elastic_data,"CV","l")
            leg.AddEntry(band,"Errorband Hist","l")
            leg.AddEntry(univ1,"Universe 0","l")
            leg.AddEntry(univ2,"Universe 1","l")
            leg.Draw()
            c1.Print("rhc_data.png")
            exit()

        if doratio: # do we want to replace selection samples with flavor ratios
            if "fhc" not in exclude: # do we care about the fhc component
                cv_histogram.MakeRatio('fhc')
            if "rhc" not in exclude: # do we care about the rhc component
                cv_histogram.MakeRatio('rhc')
            if not fit_muons: # do we want to keep the muon selections in addition to flavor ratios
                cv_histogram.RemoveHistogram('fhc_numu_selection')
                cv_histogram.RemoveHistogram('rhc_numu_selection')

        # ----- Stitch histograms together ----- #
        cv_histogram.Stitch()
        print(Chi2DataMC(cv_histogram.data_hist,cv_histogram.mc_hist))

        # ---------------------- Start Swapping Universes -----------------------------
        err_names = cv_histogram.mc_hist.GetVertErrorBandNames()
        for err in err_names:
            err = str(err)
            if err != "Flux":
                continue
            n_universes = cv_histogram.mc_hist.GetVertErrorBand(err).GetNHists()
            print("loading {} errorband that has {} n_universes".format(err,n_universes))
            for univ in range(n_universes):
                if err != "ElectronScale" and err in rhc_scale_CV.GetVertErrorBandNames():
                    # ---------------------- NuE Universe Selection Block -----------------------------
                    h_fhc_nue_selection_mc_universe   = ROOT.TFile.Open(fhc_universe).Get("EN4_predicted_Signal_{}_{}".format(err,univ))
                    h_fhc_nue_selection_data_universe = ROOT.TFile.Open(fhc_universe).Get("EN4_data_bkgSubbed_{}_{}".format(err,univ))

                    h_rhc_nue_selection_mc_universe   = ROOT.TFile.Open(rhc_universe).Get("EN4_predicted_Signal_{}_{}".format(err,univ))
                    h_rhc_nue_selection_data_universe = ROOT.TFile.Open(rhc_universe).Get("EN4_data_bkgSubbed_{}_{}".format(err,univ))

                    h_fhc_nue_selection_mc_universe.AddVertErrorBandAndFillWithCV("ElectronScale", 2)
                    for i in range(h_fhc_nue_selection_mc_universe.GetNbinsX()+1):
                        h_fhc_nue_selection_mc_universe.GetVertErrorBand("ElectronScale").GetHist(0).SetBinContent(i,h_fhc_nue_selection_mc_universe.GetBinContent(i) * fhc_scale_m1sig.GetBinContent(i))
                        h_fhc_nue_selection_mc_universe.GetVertErrorBand("ElectronScale").GetHist(1).SetBinContent(i,h_fhc_nue_selection_mc_universe.GetBinContent(i) * fhc_scale_p1sig.GetBinContent(i))

                    h_rhc_nue_selection_mc_universe.AddVertErrorBandAndFillWithCV("ElectronScale", 2)
                    for i in range(h_rhc_nue_selection_mc_universe.GetNbinsX()+1):
                        h_rhc_nue_selection_mc_universe.GetVertErrorBand("ElectronScale").GetHist(0).SetBinContent(i,h_rhc_nue_selection_mc_universe.GetBinContent(i) * rhc_scale_m1sig.GetBinContent(i))
                        h_rhc_nue_selection_mc_universe.GetVertErrorBand("ElectronScale").GetHist(1).SetBinContent(i,h_rhc_nue_selection_mc_universe.GetBinContent(i) * rhc_scale_p1sig.GetBinContent(i))

                elif err in muonList:
                    h_fhc_nue_selection_mc_universe = h_fhc_nue_selection_mc.Clone()
                    h_rhc_nue_selection_mc_universe = h_rhc_nue_selection_mc.Clone()
                    h_fhc_nue_selection_data_universe = h_fhc_nue_selection_data.Clone()
                    h_rhc_nue_selection_data_universe = h_rhc_nue_selection_data.Clone()
                else:
                    h_fhc_nue_selection_mc_universe = h_fhc_nue_selection_mc.Clone()
                    h_rhc_nue_selection_mc_universe = h_rhc_nue_selection_mc.Clone()
                    h_fhc_nue_selection_data_universe = h_fhc_nue_selection_data.Clone()
                    h_rhc_nue_selection_data_universe = h_rhc_nue_selection_data.Clone()
                    SwapUniverseCV(h_fhc_nue_selection_mc_universe,err,univ)
                    SwapUniverseCV(h_rhc_nue_selection_mc_universe,err,univ)
                    SwapUniverseCV(h_fhc_nue_selection_data_universe,err,univ)
                    SwapUniverseCV(h_rhc_nue_selection_data_universe,err,univ)


                # ---------------------- NuMu Universe Selection Block -----------------------------
                h_fhc_numu_selection_mc_universe = h_fhc_numu_selection_mc.Clone()
                h_rhc_numu_selection_mc_universe = h_rhc_numu_selection_mc.Clone()
                h_fhc_numu_selection_data_universe = h_fhc_numu_selection_data.Clone()
                h_rhc_numu_selection_data_universe = h_rhc_numu_selection_data.Clone()

                # ---------------------- Scattering Universe Selection Block -----------------------------
                h_fhc_elastic_mc_universe = h_fhc_elastic_mc.Clone()
                h_rhc_elastic_mc_universe = h_rhc_elastic_mc.Clone()
                h_fhc_imd_mc_universe = h_fhc_imd_mc.Clone()
                h_rhc_imd_mc_universe = h_rhc_imd_mc.Clone()

                data = np.array(h_fhc_elastic_data)[1:-1] # store MC bin contents excluding over/underflow bins
                h_fhc_elastic_data_universe = h_fhc_elastic_data.Clone()
                h_rhc_elastic_data_universe = h_rhc_elastic_data.Clone()
                h_fhc_imd_data_universe = h_fhc_imd_data.Clone()
                h_rhc_imd_data_universe = h_rhc_imd_data.Clone()
               
                # ---------------------- Swap Each Systematic Universe with the CV -----------------------------
                # Do monte carlo histograms
                SwapUniverseCV(h_fhc_numu_selection_mc_universe,err,univ)
                SwapUniverseCV(h_rhc_numu_selection_mc_universe,err,univ)
                SwapUniverseCV(h_fhc_elastic_mc_universe,err,univ)
                SwapUniverseCV(h_rhc_elastic_mc_universe,err,univ)
                SwapUniverseCV(h_fhc_imd_mc_universe,err,univ)
                SwapUniverseCV(h_rhc_imd_mc_universe,err,univ)

                # Do data histograms to mimic background subtraction effect, muon selections have no data systematics
                SwapUniverseCV(h_fhc_elastic_data_universe,err,univ)
                SwapUniverseCV(h_rhc_elastic_data_universe,err,univ)
                SwapUniverseCV(h_fhc_imd_data_universe,err,univ)
                SwapUniverseCV(h_rhc_imd_data_universe,err,univ)

                # ---------------------- Create Stitched Universe Histograms -----------------------------
                universe_histogram = StitchedHistogram(err+"_{}".format(univ))

                universe_histogram.AddScatteringFlavors("electron_fhc_elastic",fhcnueelnue)
                universe_histogram.AddScatteringFlavors("electron_rhc_elastic",rhcnueelnue)
                universe_histogram.AddScatteringFlavors("muon_fhc_elastic",fhcnueelnumu)
                universe_histogram.AddScatteringFlavors("muon_rhc_elastic",rhcnueelnumu)

                universe_histogram.AddSwappedSample('fhc_nue_selection',h_fhc_swap_selection)
                universe_histogram.AddSwappedSample('rhc_nue_selection',h_rhc_swap_selection)

                # ----- Initialize histogram objects with all samples ----- #
                universe_histogram.AddHistograms('fhc_elastic',h_fhc_elastic_mc_universe,h_fhc_elastic_data_universe)
                universe_histogram.AddHistograms('fhc_imd',h_fhc_imd_mc_universe,h_fhc_imd_data_universe)
                universe_histogram.AddHistograms('fhc_numu_selection',h_fhc_numu_selection_mc_universe,h_fhc_numu_selection_data_universe)
                universe_histogram.AddHistograms('fhc_nue_selection',h_fhc_nue_selection_mc_universe,h_fhc_nue_selection_data_universe)
                universe_histogram.AddHistograms('rhc_elastic',h_rhc_elastic_mc_universe,h_rhc_elastic_data_universe)
                universe_histogram.AddHistograms('rhc_imd',h_rhc_imd_mc_universe,h_rhc_imd_data_universe)
                universe_histogram.AddHistograms('rhc_numu_selection',h_rhc_numu_selection_mc_universe,h_rhc_numu_selection_data_universe)
                universe_histogram.AddHistograms('rhc_nue_selection',h_rhc_nue_selection_mc_universe,h_rhc_nue_selection_data_universe)

                universe_histogram.AddTemplates("fhc_elastic",nue=h2_fhc_elastic_template_nue,numu=h2_fhc_elastic_template_numu,swap=h2_fhc_elastic_template_numu)
                universe_histogram.AddTemplates("fhc_imd",numu=h_fhc_imd_template)
                universe_histogram.AddTemplates("fhc_numu_selection",numu=h2_fhc_musel_template)
                universe_histogram.AddTemplates("fhc_nue_selection",nue=h2_fhc_template,swap=h2_fhc_swap_template)
                universe_histogram.AddTemplates("rhc_elastic",nue=h2_rhc_elastic_template_nue,numu=h2_rhc_elastic_template_numu)
                universe_histogram.AddTemplates("rhc_elastic",nue=h2_rhc_elastic_template_nue,numu=h2_rhc_elastic_template_numu,swap=h2_rhc_elastic_template_numu)
                universe_histogram.AddTemplates("rhc_imd",numu=h_rhc_imd_template)
                universe_histogram.AddTemplates("rhc_numu_selection",numu=h2_rhc_musel_template)
                universe_histogram.AddTemplates("rhc_nue_selection",nue=h2_rhc_template,swap=h2_rhc_swap_template)
                
                # ----- Remove samples that we want to exclude from analysis ----- #
                universe_histogram.ApplyExclusion(exclude)

                # ----- Process Systematics and Synchronize across histograms ----- #
                universe_histogram.CleanErrorBands(errsToRemove)

                if doratio: # do we want to replace selection samples with flavor ratios
                    if "fhc" not in exclude: # do we care about the fhc component
                        universe_histogram.MakeRatio('fhc')
                    if "rhc" not in exclude: # do we care about the rhc component
                        universe_histogram.MakeRatio('rhc')

                    if not fit_muons: # do we want to keep the muon selections in addition to flavor ratios
                        universe_histogram.RemoveHistogram('fhc_numu_selection')
                        universe_histogram.RemoveHistogram('rhc_numu_selection')

                # ----- Stitch histograms together ----- #
                universe_histogram.Stitch()
                universe_histogram.SyncHistograms(cv_histogram)
                if err == "Flux":
                    chi2_fit, res, templates = universe_histogram.FitOscHypothesis()
                    MakePlot(universe_histogram.mc_hist,universe_histogram.data_hist,templates,res,"{}_{}.png".format(err,univ))

                chi2 = Chi2DataMC(universe_histogram.data_hist,universe_histogram.mc_hist)
                print(chi2)

                file.write("{},{},{},{},{},{},{}".format(univ,chi2,chi2_fit,res['m'],res['ue4'],res['umu4'],res['utau4']))
                file.write("\n")

