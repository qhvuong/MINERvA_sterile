import math
import ROOT
import PlotUtils
import heapq
import ctypes
from array import array
from collections import Iterable #for checking whether iterable.
from functools import partial


from config.AnalysisConfig import AnalysisConfig
from config.SignalDef import SIGNAL_DEFINITION
from config.SystematicsConfig import CONSOLIDATED_ERROR_GROUPS 

MNVPLOTTER = PlotUtils.MnvPlotter()
#config MNVPLOTTER:
MNVPLOTTER.draw_normalized_to_bin_width=False
MNVPLOTTER.legend_text_size = 0.02
MNVPLOTTER.extra_top_margin = -0.035# go slightly closer to top of pad
MNVPLOTTER.mc_bkgd_color = 46 
MNVPLOTTER.mc_bkgd_line_color = 46

MNVPLOTTER.data_bkgd_color = 12 #gray
MNVPLOTTER.data_bkgd_style = 24 #circle  

#legend entries are closer
MNVPLOTTER.height_nspaces_per_hist = 1.2
MNVPLOTTER.width_xspace_per_letter = .4
MNVPLOTTER.legend_text_size        = .03
MNVPLOTTER.legend_n_columns = 1

CANVAS = ROOT.TCanvas("c2","c2",1200,1000)
from config.SystematicsConfig import CONSOLIDATED_ERROR_GROUPS 
MNVPLOTTER.error_summary_group_map.clear();
for k,v in CONSOLIDATED_ERROR_GROUPS.items():
    vec = ROOT.vector("std::string")()
    for vs in v :
        vec.push_back(vs)
    MNVPLOTTER.error_summary_group_map[k]= vec

def Logx(canvas):
    canvas.SetLogx(1)
    return True

def Logy(canvas):
    canvas.SetLogy(1)
    return True

def Logz(canvas):
    canvas.SetLogz(1)
    return True

def PrepareSlicer(hist_holder):
    if hist_holder.dimension==1:
        return lambda x:[x.Clone()]
    elif hist_holder.dimension == 2:
        return Make2DSlice
    else:
        return None
        #raise KeyError("Only 1D,2D histograms are supported")

def IdentitySlicer(hist):
    return [hist.Clone()]

def ProjectionX(hist):
    return [hist.ProjectionX()]

def ProjectionY(hist):
    return [hist.ProjectionY()]

def PrepareSignalDecompose(data_hists,mc_hists,cates,sys_on_mc = True, sys_on_data = False,ratio=False,only_cated = False):
    def ReSumHists(h_list):
        if len(h_list) == 0:
            return None
        hnew = h_list[0].Clone()
        for i in range(1,len(h_list)):
            hnew.Add(h_list[i])
        return hnew

    if (data_hists.valid and mc_hists.valid ):
        hists = [data_hists.GetHist().GetCVHistoWithError() if sys_on_data else data_hists.GetHist().GetCVHistoWithStatError()]
        cate_hists,colors,titles  = mc_hists.GetCateList(cates)
        if only_cated :
            totalHist = ReSumHists(cate_hists)
        else:
            totalHist = mc_hists.GetHist()
        hists.append(totalHist.GetCVHistoWithError() if sys_on_mc else totalHist.GetCVHistoWithStatError())
        hists.extend(cate_hists)
        plotfunction = lambda mnvplotter,data_hist, mc_hist, *mc_ints : partial(MakeSignalDecomposePlot,color=colors,title=titles)(data_hist,mc_hist,mc_ints)
    elif mc_hists.valid:
        cate_hists,colors,titles  = mc_hists.GetCateList(cates)
        if only_cated :
            totalHist = ReSumHists(cate_hists)
        else:
            totalHist = mc_hists.GetHist()
        hists = [totalHist.GetCVHistoWithError() if sys_on_mc else totalHist.GetCVHistoWithStatError()]
        hists.extend(cate_hists)
        plotfunction = lambda mnvplotter, mc_hist, *mc_ints : partial(MakeSignalDecomposePlot,color=colors,title=titles)(None,mc_hist,mc_ints)
    else:
        raise KeyError("Non sense making signal decomposition plots without mc")
    return plotfunction,hists

def PrepareSignalDecomposeRatio(data_hists,mc_hists,Grouping,only_cated = False):
    def ReSumHists(h_list):
        if len(h_list) == 0:
            return None
        hnew = h_list[0].Clone()
        for i in range(1,len(h_list)):
            hnew.Add(h_list[i])
        return hnew

    if not mc_hists.valid:
        raise KeyError("Doesn't make sense to plot stacked histogram without MC")

    mc_list,color,title = mc_hists.GetCateList(Grouping)
    if data_hists.valid :
        plotfunction =  lambda mnvplotter, data_hist, *mc_ints: partial(MakeDataMCStackedPlot, color=color,title=title, legend="TR")(data_hist,mc_ints)
        hists = [data_hists.GetHist()]
    else:
        plotfunction =  lambda mnvplotter, mc_hist, *mc_ints: partial(MakeDataMCStackedPlot, color=color,title=title, legend = "TR")(mc_hist,mc_ints)
        tmp = mc_hists.GetHist().Clone()
        tmp.Reset()
        hists =[tmp]
    hists.extend(mc_list)
    if only_cated :
        totalHist = ReSumHists(mc_list)
    else:
        totalHist = mc_hists.GetHist()
    for i in hists:
        i.Divide(i,totalHist)
    return plotfunction,hists

def PrepareComp(data_hists,mc_hists,sys_on_mc = True, sys_on_data = False,as_frac=False):
    if not (data_hists.valid or mc_hists.valid):
        raise KeyError("both data and mc is None")
    if (data_hists.valid and mc_hists.valid):
        plotfunction = lambda mnvplotter,data_hist, mc_hist: mnvplotter.DrawDataMCWithErrorBand(data_hist,mc_hist,1.0,"TR")
        hists = [data_hists.GetHist().GetCVHistoWithError(True,as_frac) if sys_on_data else data_hists.GetHist().GetCVHistoWithStatError(),
                 mc_hists.GetHist().GetCVHistoWithError(True,as_frac) if sys_on_mc else mc_hist.GetHist().GetCVHistoWithStatError()]
    elif data_hists.valid:
        plotfunction = lambda mnvplotter,data_hist: data_hist.Draw("E1X0")
        hists = [data_hists.GetHist().GetCVHistoWithError(True,as_frac) if sys_on_data else data_hists.GetHist().GetCVHistoWithStatError()]
    else:
        plotfunction = lambda mnvplotter,mc_hist: mnvplotter.DrawMCWithErrorBand(mc_hist)
        hists = [mc_hists.GetHist().GetCVHistoWithError(True,as_frac) if sys_on_mc else mc_hist.GetHist().GetCVHistoWithStatError()]

    return plotfunction,hists

def PrepareRatio(data_hists,mc_hists):
    if not (data_hists.valid and mc_hists.valid):
        raise KeyError("both data and mc is Required for ratio")
    plotfunction = lambda mnvplotter,data_hist, mc_hist: mnvplotter.DrawDataMCRatio(data_hist, mc_hist, 1.0 ,True,True,0,2)
    h_data = data_hists.GetHist()
    h_mc = mc_hists.GetHist()
    hists = [h_data,
             h_mc]
    return plotfunction,hists

def PrepareBkgRatio(data_hists,mc_hists):
    if not (mc_hists.valid):
        raise KeyError("mc is Required for BKG ratio")

    out_bkg = mc_hists.hists["Total"].Clone("bkgTotal")
    out_bkg.Reset()
    out_sig = mc_hists.hists["Total"].Clone("sigTotal")
    out_sig.Reset()

    for group in mc_hists.hists:
        if group == "Total":
                continue
        elif group not in SIGNAL_DEFINITION and mc_hists.hists[group]:
            out_bkg.Add(mc_hists.hists[group])
        elif group in SIGNAL_DEFINITION and mc_hists.hists[group]:
            out_sig.Add(mc_hists.hists[group])
            #print("SIGNAL ",group, out_sig.Integral())

    plotfunction = lambda mnvplotter, out_sig, out_bkg: MakeBkgRatioPlot(out_sig, out_bkg, mnvplotter)
    hists = [out_sig,
             out_bkg]

    return plotfunction,hists

def PrepareErr(data_hists,mc_hists,bkgs=None,sys_on_mc=True,sys_on_data=False,grouping=None):
    def ReSumHists(h_list):
        if len(h_list) == 0:
            return None
        hnew = h_list[0].Clone()
        for i in range(1,len(h_list)):
            hnew.Add(h_list[i])
        return hnew
    MNVPLOTTER.axis_maximum = 0.4
    MNVPLOTTER.legend_n_columns = 1
    plotfunction = lambda mnvplotter,data_hist: mnvplotter.DrawErrorSummary(data_hist,"TR",True,True,0.07)
    if data_hists.valid and sys_on_data :
        hist = data_hists.GetHist()
        hist.SetTitle("Background Subtracted Data")
        hist.GetXaxis().SetTitle("Energy_{estimator}")
        hists =[hist]
    elif mc_hists.valid and sys_on_mc:
        if bkgs:
            cate_hists,colors,titles  = mc_hists.GetCateList(bkgs)
            hist = ReSumHists(cate_hists)
            hist.SetTitle("Predicted Signal")
            hist.GetXaxis().SetTitle("Energy_{estimator}")
        else:
            hist = mc_hists.GetHist()

        #hist.PopVertErrorBand("SuSA_Valencia_Weight")
        hist.PopVertErrorBand("LowQ2Pi_None")
        hist.PopVertErrorBand("LowQ2Pi")
        #hist.PopVertErrorBand("MK_model")
        #hist.PopVertErrorBand("fsi_weight")
        hists =[hist]
        #hists =[hist]
    else:
        raise KeyError("Can't make error plot for systematics config mismatch")
    #updatePlotterErrorGroup(grouping)
    #AdaptivePlotterErrorGroup(hists[0],7)
    return plotfunction,hists

def PrepareErrorBand(data_hists,mc_hists, name, sys_on_mc=True,sys_on_data=False):
    plotfunction = lambda mnvplotter,hist: MakeErrorBandPlot(hist,name,mnvplotter)
    if data_hists.valid and sys_on_data:
        hists =[data_hists.GetHist()]
    elif mc_hists.valid and sys_on_mc:
        hists =[mc_hists.GetHist()]
    else:
        raise KeyError("Can't make error plot for systematics config mismatch")
    return plotfunction,hists

def Prepare2DStack(data_hists,mc_hists,Grouping = None):
    if not mc_hists.valid:
        raise KeyError("Doesn't make sense to plot stacked histogram without MC")
    mc_list,color,title = mc_hists.GetCateList(Grouping)

    if data_hists.valid:
        plotfunction =  lambda mnvplotter, data_hist, *mc_ints: partial(MakeDataMCStackedPlot, color=color,title=title,legend="TR")(data_hist,mc_ints)
        hists = [data_hists.GetHist()]
    else:
        plotfunction =  lambda mnvplotter, mc_hist, *mc_ints: partial(MakeDataMCStackedPlot, color=color,title=title, legend = "TR")(mc_hist,mc_ints)
        tmp = mc_hists.GetHist().Clone()
        tmp.Reset()
        hists =[tmp]
    hists.extend(mc_list)
    return plotfunction,hists

def PrepareStack(data_hists,mc_hists,Grouping = None):
    if not mc_hists.valid:
        raise KeyError("Doesn't make sense to plot stacked histogram without MC")
    mc_list,color,title,raw_counts = mc_hists.GetCateList(Grouping, with_raw=True)
    if data_hists.valid:
        plotfunction =  lambda mnvplotter, data_hist, *mc_ints: partial(MakeDataMCStackedPlot, color=color,title=title, raw_counts=raw_counts, legend="TR")(data_hist,mc_ints)
        if "frontdedx" in data_hists.GetHist().GetName():
            hist = data_hists.GetHist()
            #hist.GetXaxis().SetTitle("dE/dx (MeV/cm)")
            hist.GetYaxis().SetTitle("dNEvents/d(dE/dx)")
        hists = [data_hists.GetHist()]
    else:
        plotfunction =  lambda mnvplotter, mc_hist, *mc_ints: partial(MakeDataMCStackedPlot, color=color,title=title, raw_counts=raw_counts, legend = "TR")(mc_hist,mc_ints)
        tmp = mc_hists.GetHist().Clone()
        tmp.Reset()
        hists =[tmp]
    for hist in mc_list:
        if "frontdedx" in hist.GetName():
            #hist.GetXaxis().SetTitle("dE/dx (MeV/cm)")
            hist.GetYaxis().SetTitle("dNEvents/d(dE/dx)")
        if hist == None:
            del hist
    hists.extend(mc_list)
    return plotfunction,hists

def PrepareStackNew(data_hists,mc_hists,Grouping = None):
    if not mc_hists.valid:
        raise KeyError("Doesn't make sense to plot stacked histogram without MC")
    mc_list,color,title = mc_hists.GetCateList(Grouping)
    if data_hists.valid :
        plotfunction =  lambda mnvplotter, data_hist, *mc_ints: partial(MakeDataMCStackedPlot, color=color,title=title, legend="TR")(data_hist,mc_ints)
        hists = [data_hists.GetHist()]
    else:
        plotfunction =  lambda mnvplotter, mc_hist, *mc_ints: partial(MakeDataMCStackedPlot, color=color,title=title, legend = "TR")(mc_hist,mc_ints)
        tmp = mc_hists.GetHist().Clone()
        tmp.Reset()
        hists =[tmp]
    hists.extend(mc_list)
    return plotfunction,hists

def RebinPrepareStack(data_hists,mc_hists,Grouping = None):
    if not mc_hists.valid:
        raise KeyError("Doesn't make sense to plot stacked histogram without MC")
    mc_list,color,title = mc_hists.GetCateList(Grouping)
    if data_hists.valid :
        plotfunction =  lambda mnvplotter, data_hist, *mc_ints: partial(MakeDataMCStackedPlot, color=color,title=title, legend="TR")(data_hist,mc_ints)
        data_hists.GetHist().Rebin(5)
        hists = [data_hists.GetHist()]
    else:
        plotfunction =  lambda mnvplotter, mc_hist, *mc_ints: partial(MakeDataMCStackedPlot, color=color,title=title, legend = "TR")(mc_hist,mc_ints)
        tmp = mc_hists.GetHist().Clone()
        tmp.Reset()
        hists =[tmp]
    for hist in mc_list:
        hist.Rebin(5)
    hists.extend(mc_list)
    return plotfunction,hists

def CategoryHist(data_hists,mc_hists,category):
    if not(mc_hists.valid):
        raise KeyError("No MC histogram to plot migration")
    if mc_hists.valid and category:
        hist = mc_hists.GetHist().Clone()
        hist.Reset()
        for cate in category["cate"]:
            if cate in mc_hists.hists:
                tmp = mc_hists.hists[cate]
                hist.Add(tmp)
        #hist.SetMaximum(100)
        hists = [hist]
        #plotfunction = lambda mnvplotter,mc_hist: mnvplotter.DrawNormalizedMigrationHistogram(mc_hist,False,False,True,True)
        plotfunction = lambda mnvplotter,mc_hist: MakeCategoryHist(mc_hist)
        return plotfunction,hists
    else:
        hist = mc_hists.GetHist()
        #hist.SetMaximum(100)
        hists = [hist]
        #plotfunction = lambda mnvplotter,mc_hist: mnvplotter.DrawNormalizedMigrationHistogram(mc_hist,False,False,True,True)
        plotfunction = lambda mnvplotter,mc_hist: MakeCategoryHist(mc_hist)
        #return plotfunction,hists

def MakeCategoryHist(mc_hist):
    #tcanvas = ROOT.TCanvas()
    #SetMargin(tcanvas)
    mc_hist.DrawCopy("colz")

def PrepareDiff(data_hists,mc_hists):
    if not (data_hists.valid and mc_hists.valid):
        raise KeyError("both data and mc is Required for differece")
    hists = [data_hists.GetHist(),mc_hists.GetHist()]
    def plotDifference(mnvplotter,data_hist,mc_hist):
        #print data_hist
        #ndf = ROOT.Long()
        ndf = ctypes.c_double()
        chi2mnv=mnvplotter.Chi2DataMC(data_hist,mc_hist,ndf)
        #print(chi2mnv,ndf,data_hist.Integral("width"),mc_hist.Integral("width"))
        chi2,ndf=CalChi2(data_hist,mc_hist)
        tmp = data_hist.Clone()
        tmp.Add(mc_hist,-1)
        #print(data_hist.GetName(),tmp.Integral(0,-1,"width"))
        tmp.Draw("E1")
        size = 0.035
        #align = ROOT.Long()
        #xLabel = ROOT.Double()
        #yLabel = ROOT.Double()
        align = ctypes.c_int(1)
        xLabel = ctypes.c_double(1.)
        yLabel = ctypes.c_double(1.)
        mnvplotter.DecodePosition("TR", size, align, xLabel, yLabel )
        mnvplotter.AddPlotLabel("chi2/ndf: {:.4f}/{:d}".format(chi2,ndf),xLabel, yLabel, size, 4, 112, 22)#align)
        print(("chi2/ndf: {:.4f}/{:d}".format(chi2,ndf)))

    return plotDifference,hists

def PrepareMigration(data_hists,mc_hists,Grouping):
    def ReSumHists(h_list):
        if len(h_list) == 0:
            return None
        hnew = h_list[0].Clone()
        for i in range(1,len(h_list)):
            hnew.Add(h_list[i])
        return hnew
    if not(mc_hists.valid):
        raise KeyError("No MC histogram to plot migration")
    mc_list,color,title = mc_hists.GetCateList(Grouping)
    totalHist = ReSumHists(mc_list)
    totalHist.SetTitle("Combined Signal")
    totalHist.SetMaximum(100)
    hists = [totalHist]
    plotfunction = lambda mnvplotter,mc_hist: mnvplotter.DrawNormalizedMigrationHistogram(mc_hist,False,False,True,True)
    return plotfunction,hists

def PrepareHist2D(data_hists,mc_hists,Grouping):
    def ReSumHists(h_list):
        if len(h_list) == 0:
            return None
        hnew = h_list[0].Clone()
        for i in range(1,len(h_list)):
            hnew.Add(h_list[i])
        return hnew
    if not(mc_hists.valid):
        raise KeyError("No MC histogram to plot migration")
    mc_list,color,title = mc_hists.GetCateList(Grouping)
    totalHist = ReSumHists(mc_list)
    totalHist.SetTitle("Combined Signal")
    hists = [totalHist]
    plotfunction = lambda mnvplotter,mc_hist: mc_hist.DrawCopy("colz")
    return plotfunction,hists

def updatePlotterErrorGroup(group,mnvplotter=MNVPLOTTER):
    mnvplotter.error_summary_group_map.clear();
    for k,v in group.items():
        vec = ROOT.vector("std::string")()
        for vs in v :
            vec.push_back(vs)
        mnvplotter.error_summary_group_map[k]= vec
        
def TopNErrorBand(hist,topN):
    #find N largest errorband given histogram hist
    heap = []
    for errorband_name in hist.GetVertErrorBandNames():
        sum_error = hist.GetVertErrorBand(errorband_name).GetErrorBand(False,False).Integral()
        #print(errorband_name,hist.GetVertErrorBand(errorband_name).GetErrorBand(False,False).Integral())
        if len(heap)<topN:
            heapq.heappush(heap, (sum_error,errorband_name))
        elif sum_error>heap[0][0]:
            heapq.heappushpop(heap,(sum_error,errorband_name))
    #print(heap)
    result = []
    while heap:
        result.append(heapq.heappop(heap)[1])
    return result

#def AdaptivePlotterErrorGroup(hist,topN,mnvplotter=MNVPLOTTER):
#    result = TopNErrorBand(hist,topN)
#    rest = [i for i in hist.GetVertErrorBandNames() if i not in result]
#    updatePlotterErrorGroup({"Rest":rest},mnvplotter)

def AdaptivePlotterErrorGroup(hist,result,mnvplotter=MNVPLOTTER):
    rest = [i for i in hist.GetVertErrorBandNames() if i not in result]
    updatePlotterErrorGroup({"Rest":rest},mnvplotter)

def CalMXN(N_plots,max_horizontal_plot = 4):
    height = math.ceil(1.0*N_plots/max_horizontal_plot)  #hard code max 3 plots per line
    width = math.ceil(N_plots/height) #number of plots per line
    return int(width),int(height)

def Make2DSlice(hist2D_o, X_slice=True, bin_start = 1 , bin_end = 0,interval = 1):
    # x slice true produce 1d histograms of Y variable for each x variable bins.
    # bin_end = 0 : all except overflow, = -1: all including overflow, other: bins before(i.e. *not* including ) bin_end

    slicing_hists = []
    hist2D = hist2D_o.Clone()
    axis = hist2D.GetXaxis() if X_slice else hist2D.GetYaxis()
    Nbins = axis.GetNbins()
    start = max(0,bin_start)
    end = Nbins - bin_end + 1 if bin_end <=0 else bin_end
    for i in range(start,end,interval):
        slicing_hists.append(hist2D.ProjectionX(hist2D.GetName()+str(i),i,i+interval-1,"o") if not X_slice else hist2D.ProjectionY(hist2D.GetName()+str(i),i,i+interval-1,"o"))
        slicing_hists[-1].SetTitle("%.2f<%s<%.2f"%(axis.GetBinLowEdge(i),axis.GetTitle(),axis.GetBinUpEdge(i+interval-1)))
        slicing_hists[-1].GetYaxis().SetTitle(hist2D.GetZaxis().GetTitle())
        if "frontdedx" in hist2D.GetName():
            slicing_hists[-1].Rebin(2)

    del hist2D
    return slicing_hists

def MakeDataMCPlot(data_hist, mc_hist, pot_scale=1, sys_on_mc = True, sys_on_data = False, mnvplotter=MNVPLOTTER,canvas=CANVAS):
    local_mc = (mc_hist.GetCVHistoWithError() if sys_on_mc else mc_hist.GetCVHistoWithStatError()) if mc_hist else None
    local_data = (data_hist.GetCVHistoWithError() if sys_on_data else data_hist.GetCVHistoWithStatError()) if data_hist else None
    if not (local_mc or local_data):
        raise KeyError("both data and mc is None")
    if local_mc and local_data:
        mnvplotter.DrawDataMCWithErrorBand(local_data,local_mc,pot_scale,"TR")
    elif local_mc:
        mnvplotter.DrawMCWithErrorBand(local_mc,pot_scale)
    else:
        local_data.Draw("E1X0")

def MakeModelVariantPlot(data_hist, mc_hists, color=None, title=None,legend ="TR",pot_scale=1.0,mnvplotter=MNVPLOTTER,canvas=CANVAS):
    if not mc_hists:
        raise KeyError("Doesn't make sense to plot model variation without MC")
    TArray = ROOT.TObjArray()
    for i in range(len(mc_hists)):
        if color:
            mc_hists[i].SetLineColor(color[i])
        if title:
            mc_hists[i].SetTitle(title[i])
        TArray.Add(mc_hists[i])
    #mnvplotter.DrawDataMCVariations(data_hist,TArray,pot_scale,legend,True,True,False,False,False)
    mnvplotter.DrawDataMCVariations(data_hist,TArray,pot_scale,legend,True,True,False,False)


#============================================================================
## THIS IS DUPLICATE FUNCTION
# def MakeDataMCStackedPlot(data_hist, mc_hists, legend = "TR", pot_scale=1, mnvplotter=MNVPLOTTER,canvas=CANVAS):
#     if not mc_hists:
#         raise KeyError("Doesn't make sense to plot stacked histogram without MC")
#     TArray = ROOT.TObjArray()
#     for i in range(len(mc_hists)):
#         if color:
#             mc_hists[i].SetFillColor(color[i])
#         if title:
#             mc_hists[i].SetTitle(title[i])

#         if mc_hists[i]:
#             TArray.Add(mc_hists[i])
#     if data_hist is not None:
#         mnvplotter.DrawDataStackedMC(data_hist,TArray,pot_scale,legend,"Data",0,0,1001)
#     #else:
#     elif TArray is not None:
#         mnvplotter.DrawStackedMC(TArray,pot_scale,legend,0,0,1001)
#============================================================================



def MakeSignalDecomposePlot(data_hist, mc_hist, mc_hists, title, color, pot_scale = 1.0, mnvplotter=MNVPLOTTER,canvas=CANVAS):
    #if data_hist is not None:
    #    mnvplotter.DrawDataMCWithErrorBand(data_hist,mc_hist,pot_scale,"TR")
    #else:
    #    mnvplotter.DrawMCWithErrorBand(mc_hist,pot_scale)
    mnvplotter.axis_minimum = 0
    tcanvas = canvas.GetPad(1)
    #tcanvas.SetLogy()
    mc_hists[0].SetLineWidth(4)
    mc_hists[0].SetLineColor(color[0])
    #mc_hists[0].Scale(pot_scale)
    mc_hists[0].SetMaximum(1)
    mc_hists[0].DrawCopy("HIST")
    for i in range(1,len(mc_hists)):
        mc_hists[i].SetLineWidth(4)
        mc_hists[i].SetLineColor(color[i])
        #mc_hists[i].Scale(pot_scale)
        mc_hists[i].DrawCopy("HIST SAME")

def MakeRatioPlot(data_hist,mc_hist,mnvplotter=MNVPLOTTER,canvas=CANVAS):
    if not (data_hist and mc_hist):
        raise KeyError("both data and mc is Required for ratio")
    mnvplotter.DrawDataMCRatio(data_hist, mc_hist, 1.0 ,True,True,0,2)

def MakeErrPlot(hist,mnvplotter=MNVPLOTTER,canvas=CANVAS):
    mnvplotter.DrawErrorSummary(hist)

def MakeErrorBandPlot(hist,name,mnvplotter=MNVPLOTTER,canvas=CANVAS):
    #hist.PopVertErrorBand("SuSA_Valencia_Weight")
    #hist.PopVertErrorBand("LowQ2Pi_None")
    errorband = hist.GetVertErrorBand(name)
    errorband.DrawAll("",True)

def MakeBkgRatioPlot(out_sig, out_bkg, mnvplotter=MNVPLOTTER,canvas=CANVAS):
    signal = out_sig.GetCVHistoWithError().Clone()
    background = out_bkg.GetCVHistoWithError().Clone()
    tune = -1
    best_cut = 0
    s = 0
    b = 0
    S = 0
    B = 0
    x = array('d',[])
    y1 = array('d',[])
    width = array('d',[])
    for i in range(signal.GetNbinsX()):
        s+=signal.GetBinContent(i)
        b+=background.GetBinContent(i)
        ratio = (s/math.sqrt(s+b)) if (s + b > 0) else 0
        x.append(signal.GetBinLowEdge(i)+signal.GetBinWidth(i))
        y1.append(ratio)
        width.append(signal.GetBinWidth(i))
        if ratio > tune:
            tune = ratio
            best_cut = i
            S=s
            B=b

    x[0]=0
    x.append(x[-1]+signal.GetBinWidth(signal.GetNbinsX()))
    h1 = ROOT.TH1F("h1", "Signal Significance", signal.GetNbinsX(), x)
    for i in range(signal.GetNbinsX()):
        h1.Fill(x[i]+width[i]/2,y1[i])
    h1.GetXaxis().SetTitle("E_{available}")
    h1.GetYaxis().SetTitle("s/sqrt(s+b)")

    #mnvplotter.mc_error_color = 0
    #total = background/(signal + background)
    total = out_bkg.Clone()
    #total.PopVertErrorBand("SuSA_Valencia_Weight")
    #total.PopVertErrorBand("LowQ2Pi_None")
    #total.PopVertErrorBand("MK_model")
    total.Scale(1/total.Integral())
    total.GetXaxis().SetTitle("E_{available} + E_{lepton}")
    #total.GetYaxis().SetTitle("Background Fraction per Bin")
    mnvplotter.DrawMCWithErrorBand(total)
    #mnvplotter.DrawMCWithErrorBand(h1)
    #mnvplotter.DrawDataMCRatio(out_sig, out_bkg, 1.0 ,True,True,0,4,"Signal/Background")
    #mnvplotter.AddPlotLabel("Max S/sqrt(S+B) at E_avail = {:.2f}".format(x[best_cut]+width[best_cut]),0.65,0.8,0.15)
    #mnvplotter.AddPlotLabel("{:.3f}/sqrt({:.3f}+{:.3f}) = {:.4f}".format(S,S,B,tune),0.65,0.625,0.15)
    #mnvplotter.AddPlotLabel("Max S/sqrt(S+B) at E_lepton = {:.3f}".format(signal.GetBinLowEdge(best_cut)+signal.GetBinWidth(best_cut)),0.65,0.8,0.08)
    #mnvplotter.AddPlotLabel("{:.3f}/sqrt({:.3f}+{:.3f} = {:.4f}".format(S,S,B,tune),0.65,0.625,0.08)

#def MakeGridPlot(MakingSlice,MakingEachPlot,input_hists,CanvasConfig=lambda canvas:True, draw_seperate_legend = False, mnvplotter=MNVPLOTTER,canvas=CANVAS):
#    if canvas is CANVAS:
#        canvas.Clear()
#    slices = list(map(lambda *args: args, *list(map(MakingSlice,input_hists))))
#    N_plots = len(slices)
#    canvas.Divide(*CalMXN(N_plots+int(draw_seperate_legend)))
#    for i in range(N_plots):
#        canvas.cd(i+1)
#        if not CanvasConfig(canvas.GetPad(i+1)):
#            print("Warning: failed setting canvas.")
#        if N_plots>1:
#            SetMargin(canvas.GetPad(i+1))
#        MakingEachPlot(mnvplotter,*slices[i])
#        mnvplotter.AddHistoTitle(slices[i][0].GetTitle())
#    if draw_seperate_legend:
#        for i in range(N_plots):
#            Tleg = GetTLegend(canvas.GetPad(i+1))
#        if Tleg:
#            canvas.cd(N_plots+1)
#            Tleg.SetX1(0)
#            Tleg.SetX2(1)
#            Tleg.SetY1(0)
#            Tleg.SetY2(1)
#            Tleg.SetTextSize(2*MNVPLOTTER.legend_text_size);
#            Tleg.Draw()

def SumGridPlots(MakingSlice,MakingEachPlot,input_hists,CanvasConfig=lambda canvas:True, draw_seperate_legend = False, mnvplotter=MNVPLOTTER,canvas=CANVAS,outname=None,title=None):
    slices = list(map(lambda *args: args, *list(map(MakingSlice,input_hists)))) # MakingSlice is Make2DSlice, takes argument to slice along x or y axis 
    new_slices = slices[0]
    for i in range(1,len(slices)):
        for j in range(len(slices[0])):
            new_slices[j].Add(slices[i][j])
    new_slices = [new_slices]
    del slices
    N_plots = 1
    canvas=CANVAS
    canvas.Divide(*CalMXN(N_plots+int(draw_seperate_legend)))
    for i in range(N_plots):
        canvas.cd(i+1)
        tcanvas = canvas.GetPad(i+1)
        if not CanvasConfig(tcanvas):
            print("Warning: failed setting canvas.")
        if N_plots <= 1:
            SetMargin(tcanvas)
            #tcanvas.SetLogz()
        else:
            SetMarginMulti(tcanvas)
            #i+=34
        MakingEachPlot(mnvplotter,*new_slices[i])
        #if title:
        #    mnvplotter.AddHistoTitle(title)
        #else:
        #    mnvplotter.AddHistoTitle(new_slices[i][0].GetTitle())

        #code.interact(local=locals())
    if draw_seperate_legend:
        Tleg = None
        for i in range(N_plots):
            Tleg = GetTLegend(canvas.GetPad(i+1)) or Tleg
        if Tleg:
            canvas.cd(N_plots+1)
            Tleg.SetX1(0.16)
            Tleg.SetX2(0.98)
            Tleg.SetY1(0.08)
            Tleg.SetY2(0.92)
            Tleg.SetTextSize(2*MNVPLOTTER.legend_text_size);
            Tleg.Draw()
    if outname:
        mnvplotter.MultiPrint(canvas,outname,"png")


def MakeGridPlot(MakingSlice,MakingEachPlot,input_hists,CanvasConfig=lambda canvas:True, draw_seperate_legend = False, mnvplotter=MNVPLOTTER,canvas=CANVAS,outname=None,title=None):
    if canvas is CANVAS:
        canvas.Clear()
        #canvas.SetLeftMargin(0)
    slices = list(map(lambda *args: args, *list(map(MakingSlice,input_hists))))
    N_plots = len(slices)
    #if N_plots > 1:
    #    N_plots = 9
    canvas.Divide(*CalMXN(N_plots+int(draw_seperate_legend)))
    for i in range(N_plots):
        canvas.cd(i+1)
        tcanvas = canvas.GetPad(i+1)
        if not CanvasConfig(tcanvas):
            print("Warning: failed setting canvas.")
        if N_plots <= 1:
            SetMargin(tcanvas)
            #tcanvas.SetLogz()
        else:
            SetMarginMulti(tcanvas)
            #i+=34
        MakingEachPlot(mnvplotter,*slices[i])
        if title:
            mnvplotter.AddHistoTitle(title)
        else:
            mnvplotter.AddHistoTitle(slices[i][0].GetTitle())

        #code.interact(local=locals())
    if draw_seperate_legend:
        Tleg = None
        for i in range(N_plots):
            Tleg = GetTLegend(canvas.GetPad(i+1)) or Tleg
        if Tleg:
            canvas.cd(N_plots+1)
            Tleg.SetX1(0.16)
            Tleg.SetX2(0.98)
            Tleg.SetY1(0.08)
            Tleg.SetY2(0.92)
            Tleg.SetTextSize(2*MNVPLOTTER.legend_text_size);
            Tleg.Draw()
    if outname:
        mnvplotter.MultiPrint(canvas,outname,"png")


def Print(outname,mnvplotter=MNVPLOTTER,canvas=CANVAS):
    mnvplotter.MultiPrint(canvas,outname)

def SetMargin(pad):
    pad.SetRightMargin(0.15)
    pad.SetLeftMargin(.15)
    pad.SetTopMargin(0.08)
    pad.SetBottomMargin(0.2)
    ROOT.TGaxis.SetExponentOffset(0.04,-0.1,"y")

def SetMarginMulti(pad):
    pad.SetRightMargin(0.1)
    pad.SetLeftMargin(0.1)
    pad.SetTopMargin(0.08)
    pad.SetBottomMargin(0.05)
    ROOT.TGaxis.SetExponentOffset(0.04,-0.1,"y")

def GetTLegend(pad):
    tlist = pad.GetListOfPrimitives()
    for i in tlist:
        if (isinstance(i,ROOT.TLegend)) :
            pad.RecursiveRemove(i)
            pad.Update()
            return i
    return None

# def MakeDataMCStackedPlot(data_hist, mc_hists, color=None, title=None, legend = "TR", pot_scale=1, mnvplotter=MNVPLOTTER,canvas=CANVAS):
#     TArray = ROOT.TObjArray()
#     for i in range(len(mc_hists)):
#         if not mc_hists[i] or mc_hists[i].Integral() <= 0:
#             continue
#         if color is not None:
#             mc_hists[i].SetFillColor(color[i])
#         if title is not None:
#             mc_hists[i].SetTitle(title[i])
#         TArray.Add(mc_hists[i])
#     if data_hist.Integral() > 0:
#         try:
#             mnvplotter.DrawDataStackedMC(data_hist,TArray,pot_scale,legend,"Data",0,0,1001)
#         except:
#             pass
#     else:
#         mnvplotter.DrawStackedMC(TArray,pot_scale,legend,0,0,1001)


def MakeDataMCStackedPlot(data_hist, mc_hists, color=None, title=None, legend = "TR", pot_scale=1, raw_counts=None, mnvplotter=MNVPLOTTER,canvas=CANVAS):
    if not mc_hists:
        raise KeyError("Doesn't make sense to plot stacked histogram without MC")

    TArray = ROOT.TObjArray()

    for i, h in enumerate(mc_hists):
        if not h:
            continue

        # ---- NEW: skip tiny raw-count components ----
        if raw_counts is not None:
            if i >= len(raw_counts) or raw_counts[i] is None or raw_counts[i] < 0.1:
                continue
        # --------------------------------------------

        if color is not None and i < len(color):
            h.SetFillColor(color[i])

        if title is not None and i < len(title):
            if raw_counts is not None and i < len(raw_counts):
                h.SetTitle(f"{title[i]} ({raw_counts[i]:.1f})")
            else:
                h.SetTitle(f"{title[i]}")

        TArray.Add(h)

    if data_hist is not None and data_hist.Integral() > 0:
        try:
            mnvplotter.DrawDataStackedMC(
                data_hist, TArray, pot_scale, legend, "Data", 0, 0, 1001
            )
        except Exception as e:
            print("DrawDataStackedMC failed:", e)
    else:
        mnvplotter.DrawStackedMC(TArray, pot_scale, legend, 0, 0, 1001)


def MakeMigrationPlots(hist, output, no_text = False, fix_width = True, mnvplotter= MNVPLOTTER, canvas = CANVAS):
    if canvas is CANVAS:
        canvas.Clear()
    #make sure migration matrix has fix bin width
    if fix_width:
        hist.GetXaxis().Set(hist.GetNbinsX(),0,hist.GetNbinsX())
        hist.GetXaxis().SetTitle("reco bin number")
        hist.GetYaxis().Set(hist.GetNbinsY(),0,hist.GetNbinsY())
        hist.GetYaxis().SetTitle("truth bin number")
    mnvplotter.DrawNormalizedMigrationHistogram(hist,False,False,True,no_text)
    mnvplotter.MultiPrint(canvas,output)

def Make2DPlot(hist,output,mnvplotter= MNVPLOTTER, canvas = CANVAS):
    if canvas is CANVAS:
        canvas.Clear()

    plotfunction = lambda mnvplotter,mc_hist: mc_hist.DrawCopy("colz")
    return plotfunction,hists

def CalChi2(data_hist,mc_hist,pot_scale=1.0):
    chi2 = 0
    ndf = 0
    for i in range (0,data_hist.GetSize()):
        data = data_hist.GetBinContent(i)
        mc = mc_hist.GetBinContent(i)*pot_scale
        sig = data_hist.GetBinError(i)
        #print (i,data,mc,chi2)
        chi2 += (data-mc)**2/sig/sig if data>0 else 0
        ndf += 1 if data>0 else 0
    print("data/mc integral: {}/{}".format(data_hist.Integral(),mc_hist.Integral()))
    print("chi2/ndf of histogram {} is {}/{}.".format(data_hist.GetName(),chi2,ndf))
    return chi2,ndf

