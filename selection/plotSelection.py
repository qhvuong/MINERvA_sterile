import os
import sys
import ROOT
import PlotUtils

#insert path for modules of this package.
from config import PlotConfig
from config.AnalysisConfig import AnalysisConfig
from config.DrawingConfig import PLOTS_TO_MAKE,Default_Plot_Type,Default_Scale,DefaultPlotters,DefaultSlicer
from config.SignalDef import SIGNAL_DEFINATION
from tools import Utilities,PlotTools
from tools.PlotLibrary import HistHolder

# Get This from Rob. Thanks Rob.
# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
# Specifically, w/o this, this script seg faults in the case where I try to instantiate FluxReweighterWithWiggleFit w/ nuE constraint set to False for more than one playlist
ROOT.TH1.AddDirectory(False)
SELECTED_SIDEBANDS = AnalysisConfig.sidebands
mnvplotter = PlotUtils.MnvPlotter()
# mnvplotter.draw_normalized_to_bin_width = False
mnvplotter.axis_maximum = 1

def MakePlot(data_hists,mc_hists,config,true_mc=False):
    #scale
    if "scale" in config:
        config["scale"](data_hists)
        config["scale"](mc_hists)
    else: 
        Default_Scale(data_hists)
        Default_Scale(mc_hists)
    CanvasConfig = config.setdefault("canvasconfig",lambda x:True)
    ### want a MC-like pseudodata signal region to avoid preliminary unblinding
    if AnalysisConfig.pseudodata:
        print("using pseudodata")
        data_signal = data_hists.GetHist()
        if data_signal: 
            for q in range(0,data_signal.GetNbinsX()+1):
                data_signal.SetBinContent(q,mc_hists.hists["Total"].GetBinContent(q))

    PlotType = config.setdefault("plot_type",Default_Plot_Type)
    # print(f"\n[MakePlot] Plot '{config['name']}' using plot_type = {PlotType}")

    # if PlotType == "stacked":
    #     # or whatever grouping you use, e.g. CONSOLIDATED_ERROR_GROUPS
    #     mc_ints, colors, titles = mc_hists.GetCateList(CONSOLIDATED_ERROR_GROUPS)

    #     print(f"[MakePlot]   N MC categories for stack: {len(mc_ints)}")
    #     print(f"[MakePlot]   Titles: {titles}")
    
    typeBool = PlotType!="migration" and PlotType!="category_hist" and PlotType!="hist2d"
    slicer = config.setdefault("slicer", DefaultSlicer(data_hists)) if typeBool else PlotTools.IdentitySlicer
    draw_seperate_legend = config.setdefault("draw_seperate_legend",data_hists.dimension!=1 and (PlotType != "migration" or PlotType != "category_hist" or PlotType != "hist2d"))
    try:
        custom_tag = config["tag"]+PlotType if "tag" in config else PlotType
        if PlotType == "custom":
            plotfunction,hists=config["getplotters"](data_hists,mc_hists)
        elif PlotType == "category_hist":
            if "args" in config:
                args = config["args"]
            elif "args" in DefaultPlotters[PlotType]:
                args = DefaultPlotters[PlotType]["args"]
            else:
                args = None
            categories = args[0]
            for category in categories:
                plotfunction,hists = DefaultPlotters[PlotType]["func"](data_hists,mc_hists,categories[category])
                PlotTools.MakeGridPlot(PlotTools.IdentitySlicer,plotfunction,hists,draw_seperate_legend=False,title=category)
                PlotTools.Print(AnalysisConfig.PlotPath(data_hists.plot_name,sideband,category))
                print("plot {} made for category {}.".format(data_hists.plot_name,category))
        else:
            if "args" in config:
                args = config["args"]
            elif "args" in DefaultPlotters[PlotType]:
                args = DefaultPlotters[PlotType]["args"]
            else:
                args = None
            if args is None:
                plotfunction,hists = DefaultPlotters[PlotType]["func"](data_hists,mc_hists)
            else:
                plotfunction,hists = DefaultPlotters[PlotType]["func"](data_hists,mc_hists,*args)

            if PlotType == "2Dstacked":
                PlotTools.SumGridPlots(slicer,plotfunction,hists,draw_seperate_legend=False)
            else:
                PlotTools.MakeGridPlot(slicer,plotfunction,hists,draw_seperate_legend=False)
            PlotTools.Print(AnalysisConfig.PlotPath(data_hists.plot_name,sideband,custom_tag))
            print("plot {} made.".format(data_hists.plot_name))
    except KeyError as e:
        print(e,"plot {} not made.".format(data_hists.plot_name))
        return False
    return True

def MakeRatio(signalHist,sidebandHist,config):
    #scale
    if "scale" in config:
        config["scale"](signalHist)
        config["scale"](sidebandHist)
    else: 
        Default_Scale(signalHist)
        Default_Scale(sidebandHist)

    sig_bkg = signalHist.hists["Total"].Clone()
    sig_bkg.Reset()
    sid_bkg = sidebandHist.hists["Total"].Clone()
    sid_bkg.Reset()

    for group in signalHist.hists:
        if group == "Total":
                continue
        elif group not in SIGNAL_DEFINATION:
            sig_bkg.Add(signalHist.hists[group])
            sid_bkg.Add(sidebandHist.hists[group])

    
    c1 = ROOT.TCanvas()
    sig_bkg.GetVertErrorBand("Flux").DrawAll("hist",True)
    c1.Print("pre_tuneSigBkgFlux.png")
    c1 = ROOT.TCanvas()
    sid_bkg.GetVertErrorBand("Flux").DrawAll("hist",True)
    c1.Print("pre_tuneSidBkgFlux.png")

    sig_errorband = sig_bkg.GetVertErrorBand("Flux")
    sid_errorband = sid_bkg.GetVertErrorBand("Flux")
    for ith in range(sig_errorband.GetNHists()):
        k = 12
        #print(sig_errorband.GetHist(ith).GetBinContent(k),sid_errorband.GetHist(ith).GetBinContent(k))
    c1 = ROOT.TCanvas()
    fluxerr = sig_bkg.GetVertErrorBand("Flux").GetErrorBand(True,False).Clone()
    fluxerr.Draw()
    c1.Print("sig_bkgPreTune_Fluxerrband.png")
    fluxerr = sid_bkg.GetVertErrorBand("Flux").GetErrorBand(True,False).Clone()
    fluxerr.Draw()
    c1.Print("sid_bkgPreTune_Fluxerrband.png")

    sig_bkg.Scale(1/sig_bkg.Integral())
    sid_bkg.Scale(1/sid_bkg.Integral())
    sig_bkg.Divide(sig_bkg,sid_bkg)
    sig_bkg.GetXaxis().SetTitle("E_{available} + E_{lepton}")
    sig_bkg.GetYaxis().SetTitle("Ratio")
    sig_bkg.SetTitle("Signal/Sideband Background Ratio")
    #total.GetYaxis().SetTitle("Background Fraction per Bin")
    mnvplotter.DrawMCWithErrorBand(sig_bkg)

    PlotTools.Print(AnalysisConfig.PlotPath("EN4_ratio","Combined","bkgratioN4_tun"))
    c1 = ROOT.TCanvas()
    mnvplotter.DrawErrorSummary(sig_bkg,"TR",True,True,0)
    c1.Print("{}_pre_tuneSigBkgErrSummary.png".format(signalHist.plot_name))

def MakeEfficiencyPlot(data_hists,mc_hists,true_hists,config,true_mc=False):
    #scale
    if "scale" in config:
        config["scale"](data_hists)
        config["scale"](mc_hists)
        config["scale"](true_hists)
    else: 
        Default_Scale(data_hists)
        Default_Scale(mc_hists)
        Default_Scale(true_hists)

    CanvasConfig = config.setdefault("canvasconfig",lambda x:True)
    PlotType = config.setdefault("plot_type",Default_Plot_Type)
    typeBool = PlotType!="migration"# or PlotType!="2d"
    slicer = config.setdefault("slicer", DefaultSlicer(data_hists)) if typeBool else PlotTools.IdentitySlicer
    draw_seperate_legend = config.setdefault("draw_seperate_legend",data_hists.dimension!=1 and PlotType != "migration")
    try:
        custom_tag = config["tag"]+PlotType if "tag" in config else PlotType
        if PlotType == "custom":
            plotfunction,hists=config["getplotters"](data_hists,mc_hists)
        else:
            if "args" in config:
                args = config["args"]
            elif "args" in DefaultPlotters[PlotType]:
                args = DefaultPlotters[PlotType]["args"]
            else:
                args = None
            if args is None:
                plotfunction,hists = DefaultPlotters[PlotType]["func"](data_hists,mc_hists)
            else:
                plotfunction,hists = DefaultPlotters[PlotType]["func"](data_hists,mc_hists,*args)
            
            PlotTools.MakeEfficiencyGridPlot(slicer,plotfunction,hists,draw_seperate_legend=False)
            PlotTools.Print(AnalysisConfig.PlotPath(data_hists.plot_name,sideband,custom_tag))
            print("plot {} made.".format(data_hists.plot_name))
    except KeyError as e:
        print("plot {} not made.".format(data_hists.plot_name))
        #print(e)
        return False
    return True

if __name__ == "__main__":
    #input knobs
    playlist=AnalysisConfig.playlist
    type_path_map = { t:AnalysisConfig.SelectionHistoPath(playlist,t =="data",False) for t in AnalysisConfig.data_types}
    print("playlist",playlist)
    print("type_path_map",type_path_map)
    print("ntuple_tag",AnalysisConfig.ntuple_tag)
    data_file,mc_file,pot_scale,data_pot,mc_pot = Utilities.getFilesAndPOTScale(playlist,type_path_map,AnalysisConfig.ntuple_tag,True)

    standPOT = data_pot if data_pot is not None else mc_pot 
    sideband_map = {}

    for config in PLOTS_TO_MAKE:
        # if "Front dEdX" in config['name']:
        #     sideband = 'dEdX'
        #     sidebandHist = HistHolder(config["name"] if "name" in config else config,mc_file,"dEdX",True,mc_pot,standPOT)
        #     signalHist = HistHolder(config["name"] if "name" in config else config,mc_file,"Signal",True,mc_pot,standPOT)
        #     signalHist.Add(sidebandHist)
        #     dataSignal = HistHolder(config["name"] if "name" in config else config,data_file,"Signal",False,data_pot,standPOT)
        #     dataSideband = HistHolder(config["name"] if "name" in config else config,data_file,"dEdX",False,data_pot,standPOT)
        #     dataSignal.Add(dataSideband)
        #     MakePlot(dataSignal,signalHist,config)
        #     #MakeRatio(signalHist,sidebandHist,config)
        #     continue

        sideband_group =  config.setdefault("sideband_group",["Signal"]+SELECTED_SIDEBANDS)
        if isinstance(sideband_group,list):
            for sideband in sideband_group:
                data_hists = HistHolder(config["name"] if "name" in config else config,data_file,sideband,False,data_pot,standPOT)
                mc_hists = HistHolder(config["name"] if "name" in config else config,mc_file,sideband,True,mc_pot,standPOT)
                MakePlot(data_hists,mc_hists,config)
        else:
            #assuing sideband_group is a tuple of name, and list of sidebands
            sideband = sideband_group[0]
            sidebands = sideband_group[1]
            data_hists = HistHolder(config["name"] if "name" in config else config,data_file,sidebands[0],False,data_pot,standPOT)
            mc_hists = HistHolder(config["name"] if "name" in config else config,mc_file,sidebands[0],True,mc_pot,standPOT)
            for _ in range(1,len(sidebands)):
                data_hists.Add(HistHolder(config["name"] if "name" in config else config,data_file,sidebands[_],False,data_pot,standPOT))
                mc_hists.Add(HistHolder(config["name"] if "name" in config else config,mc_file,sidebands[_],True,mc_pot,standPOT))
            MakePlot(data_hists,mc_hists,config)
