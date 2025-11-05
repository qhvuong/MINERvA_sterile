import ROOT
import PlotUtils
import numpy as np
import ctypes
import json
from tools.PlotLibrary import HistHolder
ROOT.TH1.AddDirectory(False)
ROOT.SetMemoryPolicy(ROOT.kMemoryStrict)

legend_text_size = .025

def plot_osc_side_by_side(mc_hists,null_hists,osc_hists,data_hists,titles,plot_text,MNVPLOTTER,narrow_pads=None,narrow_factor=0.5):
    if len(osc_hists) != 8:
        raise ValueError("Need exactly 8 histograms for 2x4 layout")

    ROOT.gStyle.SetOptTitle(1)

    if narrow_pads is None:
        narrow_pads = []

    # Global min/max for log scale
    global_max = max(h.GetMaximum() for h in data_hists)
    global_min = min(h.GetBinContent(b)
                     for h in osc_hists
                     for b in range(1, h.GetNbinsX()+1)
                     if h.GetBinContent(b) > 0)

    # Margins
    left_margin = 0.2
    right_margin = 0.05
    bottom_margin = 0.
    top_margin = 0.

    # Column widths (relative)
    ncols, nrows = 4, 2
    col_widths = []
    for c in range(ncols):
        w = narrow_factor if c in narrow_pads else 1.0
        col_widths.append(w)
    total_width = sum(col_widths)
    norm_col_widths = [w / total_width for w in col_widths]
    row_height = 1.0 / nrows

    # Create canvas
    c = ROOT.TCanvas("cnorm", "2x4 histograms", 1200, 800)

    # Make pads
    pads = []
    for i, h in enumerate(osc_hists, start=1):
        row = (i - 1) // ncols  # 0 = top, 1 = bottom
        col = (i - 1) % ncols   # 0..3

        # Compute pad coords
        x1 = sum(norm_col_widths[:col])
        x2 = sum(norm_col_widths[:col+1])
        y2 = 1.0 - row * row_height
        y1 = y2 - row_height

        pad = ROOT.TPad(f"pad{i}", "", x1, y1, x2, y2)

        # Margins
        pad.SetLeftMargin(left_margin if col == 0 else 0.0)       # only leftmost has y labels
        pad.SetRightMargin(right_margin if col == ncols-1 else 0.0)
        pad.SetTopMargin(0.0)        # top row has titles
        pad.SetBottomMargin(0.0)  # bottom row has x labels

        pad.Draw()
        pads.append(pad)

    # Font sizes (uniform across pads)
    label_size = 14
    title_size = 18

    # Draw hists
    for i, (pad, h) in enumerate(zip(pads, osc_hists), start=1):
        row = (i - 1) // ncols
        col = (i - 1) % ncols
        null_hist = null_hists[i-1].Clone()

        pad.cd()
        fraction = .6
        top_y1 = 0.4 if row == 1 else 0.28
        top = ROOT.TPad(f"DATAMC{i}", "", 0, top_y1, 1, 1)
        bottom = ROOT.TPad(f"Ratio{i}", "", 0, 0, 1, top_y1)

        top.Draw()
        bottom.Draw()

        top.cd()
        top.SetLogy(True)

        h.SetMaximum(global_max * 1.2)
        h.SetMinimum(1)

        ROOT.gStyle.SetTitleFont(43,"");
        ROOT.gStyle.SetTitleFontSize(28);

        h.GetXaxis().SetTitleFont(43)
        h.GetYaxis().SetTitleFont(43)
        h.GetXaxis().SetLabelFont(43)
        h.GetYaxis().SetLabelFont(43)
        h.GetXaxis().SetLabelSize(label_size)
        h.GetXaxis().SetTitleSize(title_size)
        h.GetYaxis().SetLabelSize(label_size)
        h.GetYaxis().SetTitleSize(title_size)
        null_hist.GetXaxis().SetTitleFont(43)
        null_hist.GetYaxis().SetTitleFont(43)
        null_hist.GetXaxis().SetLabelFont(43)
        null_hist.GetYaxis().SetLabelFont(43)
        null_hist.GetXaxis().SetLabelSize(label_size)
        null_hist.GetXaxis().SetTitleSize(title_size)
        null_hist.GetYaxis().SetLabelSize(label_size)
        null_hist.GetYaxis().SetTitleSize(title_size)

        top.SetBottomMargin(0.0)
        top.SetTopMargin(0.18 if row == 0 else 0.02)   # enough for title
        bottom.SetTopMargin(0.0)
        bottom.SetBottomMargin(0.3 if row == 1 else 0) # enough for axis labels/title

        top.SetLeftMargin(left_margin if col == 0 else 0.0)       # only leftmost has y labels
        top.SetRightMargin(right_margin if col == ncols-1 else 0.0)

        bottom.SetLeftMargin(left_margin if col == 0 else 0.0)       # only leftmost has y labels
        bottom.SetRightMargin(right_margin if col == ncols-1 else 0.0)

        # Hide y-axis except first col
        if col != 0:
            h.GetYaxis().SetLabelSize(0)
            h.GetYaxis().SetTitleSize(0)
            null_hist.GetYaxis().SetLabelSize(0)
            null_hist.GetYaxis().SetTitleSize(0)

        # Hide x-axis labels except bottom row
        if row != 1:
            h.GetXaxis().SetLabelSize(0)
            h.GetXaxis().SetTitleSize(0)
            null_hist.GetXaxis().SetLabelSize(0)
            null_hist.GetXaxis().SetTitleSize(0)

        # Titles: only show on top row
        if row != 0:
            h.SetTitle("")
            null_hist.SetTitle("")

        # Force ticks at multiples of 5
        if i in [2, 6]:
            h.GetXaxis().SetNdivisions(303)
            null_hist.GetXaxis().SetNdivisions(303)
        else:
            h.GetXaxis().SetNdivisions(505)
            null_hist.GetXaxis().SetNdivisions(505)

        TArray = ROOT.THStack("{i}","")
        for j in mc_hists[i-1]:
            j.SetFillStyle(1001)
            TArray.Add(j)

        h.SetLineColor(ROOT.kBlue)
        h.SetLineWidth(0)
        h.Draw("HIST")
        data_hists[i-1].Draw("Same")
        TArray.DrawClone("same hist")

        osc_ratio = osc_hists[i-1].Clone(f"null_{i}")
        osc_ratio.SetLineColor(ROOT.kBlue)
        osc_ratio.SetLineWidth(2)

        if row == 1 and col == 0:
            leg = ROOT.TLegend(0.5,0.4,0.95,0.95)
            leg.SetBorderSize(0)
            leg.SetTextSize(legend_text_size*3)
            leg.AddEntry(null_hist,"Null Hypothesis","l")
            leg.AddEntry(osc_ratio,"Osc. Hypothesis","l")
            for j in mc_hists[0]:
                if "rightarrow" not in j.GetTitle() or "tau" in j.GetTitle():
                    leg.AddEntry(j,j.GetTitle(),"f")
            leg.AddEntry(data_hists[i-1],"Data","p")
            leg.DrawClone()

        if row == 0 and col == 0:
            latex = ROOT.TLatex()
            latex.SetNDC()                # Use normalized device coordinates (0–1)
            latex.SetTextSize(legend_text_size*2.5)       # Adjust text size
            latex.SetTextAlign(13)        # Align left, vertically centered

            x, y = 0.25, 0.8             # Starting coordinates (NDC)
            line_spacing = 0.12           # Vertical spacing between lines
            
            for t,text in enumerate(plot_text[:2]):
                latex.DrawLatex(x, y-t*line_spacing,text)

        if row == 0 and col == 1:
            latex = ROOT.TLatex()
            latex.SetNDC()                # Use normalized device coordinates (0–1)
            latex.SetTextSize(legend_text_size*4)       # Adjust text size
            latex.SetTextAlign(13)        # Align left, vertically centered

            x, y = 0.15, 0.8             # Starting coordinates (NDC)
            line_spacing = 0.12           # Vertical spacing between lines
            
            for t,text in enumerate(plot_text[2:]):
                latex.DrawLatex(x, y-t*line_spacing,text)

        for obj in top.GetListOfPrimitives():
            if type(obj) == ROOT.TH1D:
                obj.SetTitle(titles[i-1])

        h.GetXaxis().SetTickSize()
        h.GetXaxis().SetTickLength()
        h.GetYaxis().SetTickSize()
        h.GetYaxis().SetTickLength()

        bottom.cd()
        
        nullErrors = null_hist.GetTotalError(False, True, False) #The second "true" makes this fractional error, the third "true" makes this cov area normalized
        for whichBin in range(0, nullErrors.GetXaxis().GetNbins()+1): 
            nullErrors.SetBinError(whichBin, max(nullErrors.GetBinContent(whichBin), 1e-9))
            nullErrors.SetBinContent(whichBin, 1)

        #Error envelope for the MC
        nullErrors.SetLineWidth(0)
        nullErrors.SetMarkerStyle(0)
        nullErrors.SetFillColorAlpha(ROOT.kPink + 1, 0.3)
        nullErrors.SetMinimum(.5)
        nullErrors.SetMaximum(1.5)

        nullErrors.SetTitle("")
        nullErrors.GetYaxis().SetTitle("Null Hyp. Ratio")
        nullErrors.GetYaxis().SetNdivisions(505)

        osc_ratio.Divide(osc_ratio,null_hist)
        ratio = data_hists[i-1].Clone(f"ratio_{i}")
        ratio.Divide(ratio,null_hist)

        straightLine = nullErrors.Clone()
        straightLine.SetLineColor(ROOT.kRed)
        straightLine.SetLineWidth(2)
        straightLine.SetFillColor(0)
        straightLine.DrawClone("HIST SAME")

        nullErrors.DrawClone("SAME E2")
        ratio.DrawClone("SAME E1")
        osc_ratio.SetCanExtend(ROOT.TH1.kAllAxes)
        for obin in range(osc_ratio.GetNbinsX()+1):
            if osc_ratio.GetBinContent(obin) <= 0:
                osc_ratio.SetBinContent(obin,1)
        osc_ratio.DrawClone("SAME HIST")

    c.Update()
    return c

def plot_errs_side_by_side(hists,narrow_pads,narrow_factor,MNVPLOTTER):
    if len(hists) != 8:
        raise ValueError("Need exactly 8 histograms for 2x4 layout")

    ROOT.gStyle.SetOptTitle(1)

    if narrow_pads is None:
        narrow_pads = []

    # Margins
    left_margin = 0.2
    right_margin = 0.05
    bottom_margin = 0.18
    top_margin = 0.18

    # Column widths (relative)
    ncols, nrows = 4, 2
    col_widths = []
    for c in range(ncols):
        w = narrow_factor if c in narrow_pads else 1.0
        col_widths.append(w)
    total_width = sum(col_widths)
    norm_col_widths = [w / total_width for w in col_widths]
    row_height = 1.0 / nrows

    # Create canvas
    c = ROOT.TCanvas("cerr", "2x4 histograms", 1200, 800)

    # Make pads
    pads = []
    for i, h in enumerate(hists, start=1):
        ROOT.gStyle.SetTitleFont(43,"");
        ROOT.gStyle.SetTitleFontSize(28);

        label_size = 14
        title_size = 18

        h.GetXaxis().SetTitleFont(43)
        h.GetYaxis().SetTitleFont(43)
        h.GetXaxis().SetLabelFont(43)
        h.GetYaxis().SetLabelFont(43)
        h.GetXaxis().SetLabelSize(label_size)
        h.GetXaxis().SetTitleSize(title_size)
        h.GetYaxis().SetLabelSize(label_size)
        h.GetYaxis().SetTitleSize(title_size)

        row = (i - 1) // ncols  # 0 = top, 1 = bottom
        col = (i - 1) % ncols   # 0..3

        # Hide y-axis except first col
        if col != 0:
            h.GetYaxis().SetLabelSize(0)
            h.GetYaxis().SetTitleSize(0)

        # Hide x-axis labels except bottom row
        if row != 1:
            h.GetXaxis().SetLabelSize(0)
            h.GetXaxis().SetTitleSize(0)

        # Titles: only show on top row
        if row != 0:
            h.SetTitle("")

        # Force ticks at multiples of 5
        if i in [2, 6]:
            h.GetXaxis().SetNdivisions(303)
        else:
            h.GetXaxis().SetNdivisions(505)

        # Compute pad coords
        x1 = sum(norm_col_widths[:col])
        x2 = sum(norm_col_widths[:col+1])
        y2 = 1.0 - row * row_height
        y1 = y2 - row_height

        pad = ROOT.TPad(f"pad{i}", "", x1, y1, x2, y2)

        # Margins
        pad.SetLeftMargin(left_margin if col == 0 else 0.0)       # only leftmost has y labels
        pad.SetRightMargin(right_margin if col == ncols-1 else 0.0)
        pad.SetTopMargin(top_margin if row == 0 else 0)        # top row has titles
        pad.SetBottomMargin(bottom_margin if row == 1 else 0)  # bottom row has x labels

        pad.Draw()
        pads.append(pad)

    # Font sizes (uniform across pads)
    label_size = 14
    title_size = 18

    # Draw hists
    for i, (pad, h) in enumerate(zip(pads, hists), start=1):
        row = (i - 1) // ncols  # 0 = top, 1 = bottom
        col = (i - 1) % ncols   # 0..3

        pad.cd()
        MNVPLOTTER.axis_maximum = 0.25
        MNVPLOTTER.axis_title_font_x = 43
        MNVPLOTTER.axis_title_font_y = 43
        MNVPLOTTER.axis_label_font = 43
        MNVPLOTTER.axis_title_size_x = 18
        MNVPLOTTER.axis_label_size = 14
        MNVPLOTTER.axis_title_size_y = 18
        MNVPLOTTER.title_font = 43
        MNVPLOTTER.extra_top_margin = 0
        MNVPLOTTER.extra_bottom_margin = 0
        MNVPLOTTER.extra_right_margin = 0
        MNVPLOTTER.extra_left_margin = 0
        MNVPLOTTER.headroom = 0
        MNVPLOTTER.footroom = 0
        MNVPLOTTER.legend_n_columns = 1
        MNVPLOTTER.legend_text_size = legend_text_size*1.8
        MNVPLOTTER.legend_offset_y = -0.065
        MNVPLOTTER.legend_offset_x = 0.05

        leg = "TR" if row == 1 and col == 0 else "N"
        MNVPLOTTER.DrawErrorSummary(h,leg,True,True,0)

    c.Update()
    return c

def plot_side_by_side(hists,errs,data,narrow_pads=None, narrow_factor=0.5,chi2=0,penalty=0):
    if len(hists) != 8:
        raise ValueError("Need exactly 8 histograms for 2x4 layout")

    ROOT.gStyle.SetOptTitle(1)

    if narrow_pads is None:
        narrow_pads = []

    # Global min/max for log scale
    #global_max = max(h.GetMaximum() for h in hists) TODO
    global_max = 2300
    [h.GetXaxis().SetRangeUser(1,6) for h in hists]
    global_min = min(h.GetBinContent(b)
                     for h in hists
                     for b in range(1, h.GetNbinsX()+1)
                     if h.GetBinContent(b) > 0)

    # Margins
    left_margin = 0.2
    right_margin = 0.05
    bottom_margin = 0.
    top_margin = 0.

    # Column widths (relative)
    ncols, nrows = 4, 2
    col_widths = []
    for c in range(ncols):
        w = narrow_factor if c in narrow_pads else 1.0
        col_widths.append(w)
    total_width = sum(col_widths)
    norm_col_widths = [w / total_width for w in col_widths]
    row_height = 1.0 / nrows

    # Create canvas
    c = ROOT.TCanvas("cnorm", "2x4 histograms", 1200, 800)

    # Make pads
    pads = []
    for i, h in enumerate(hists, start=1):
        row = (i - 1) // ncols  # 0 = top, 1 = bottom
        col = (i - 1) % ncols   # 0..3

        # Compute pad coords
        x1 = sum(norm_col_widths[:col])
        x2 = sum(norm_col_widths[:col+1])
        y2 = 1.0 - row * row_height
        y1 = y2 - row_height

        pad = ROOT.TPad(f"pad{i}", "", x1, y1, x2, y2)

        # Margins
        pad.SetLeftMargin(left_margin if col == 0 else 0.0)       # only leftmost has y labels
        pad.SetRightMargin(right_margin if col == ncols-1 else 0.0)
        pad.SetTopMargin(0.0)        # top row has titles
        pad.SetBottomMargin(0.0)  # bottom row has x labels

        pad.Draw()
        pads.append(pad)

    # Font sizes (uniform across pads)
    label_size = 14
    title_size = 18

    # Draw hists
    for i, (pad, h) in enumerate(zip(pads, hists), start=1):
        row = (i - 1) // ncols
        col = (i - 1) % ncols

        pad.cd()
        fraction = .6
        top_y1 = 0.4 if row == 1 else 0.28
        top = ROOT.TPad(f"DATAMC{i}", "", 0, top_y1, 1, 1)
        bottom = ROOT.TPad(f"Ratio{i}", "", 0, 0, 1, top_y1)

        top.Draw()
        bottom.Draw()

        top.cd()
        top.SetLogy(False)

        h.SetMaximum(global_max * 1.2)
        h.SetMinimum(0)

        ROOT.gStyle.SetTitleFont(43,"");
        ROOT.gStyle.SetTitleFontSize(28);

        err = errs[i-1]
        d = data[i-1]

        h.GetXaxis().SetTitleFont(43)
        h.GetYaxis().SetTitleFont(43)
        h.GetXaxis().SetLabelFont(43)
        h.GetYaxis().SetLabelFont(43)
        h.GetXaxis().SetLabelSize(label_size)
        h.GetXaxis().SetTitleSize(title_size)
        h.GetYaxis().SetLabelSize(label_size)
        h.GetYaxis().SetTitleSize(title_size)
        d.GetXaxis().SetTitleFont(43)
        d.GetYaxis().SetTitleFont(43)
        d.GetXaxis().SetLabelFont(43)
        d.GetYaxis().SetLabelFont(43)
        d.GetXaxis().SetLabelSize(label_size)
        d.GetXaxis().SetTitleSize(title_size)
        d.GetYaxis().SetLabelSize(label_size)
        d.GetYaxis().SetTitleSize(title_size)

        top.SetBottomMargin(0.0)
        top.SetTopMargin(0.18 if row == 0 else 0.02)   # enough for title
        bottom.SetTopMargin(0.0)
        bottom.SetBottomMargin(0.3 if row == 1 else 0) # enough for axis labels/title

        top.SetLeftMargin(left_margin if col == 0 else 0.0)       # only leftmost has y labels
        top.SetRightMargin(right_margin if col == ncols-1 else 0.0)

        bottom.SetLeftMargin(left_margin if col == 0 else 0.0)       # only leftmost has y labels
        bottom.SetRightMargin(right_margin if col == ncols-1 else 0.0)

        # Hide y-axis except first col
        if col != 0:
            h.GetYaxis().SetLabelSize(0)
            h.GetYaxis().SetTitleSize(0)
            d.GetYaxis().SetLabelSize(0)
            d.GetYaxis().SetTitleSize(0)

        # Hide x-axis labels except bottom row
        if row != 1:
            h.GetXaxis().SetLabelSize(0)
            h.GetXaxis().SetTitleSize(0)
            d.GetXaxis().SetLabelSize(0)
            d.GetXaxis().SetTitleSize(0)

        # Titles: only show on top row
        if row != 0:
            h.SetTitle("")

        # Force ticks at multiples of 5
        if i in [2, 6]:
            h.GetXaxis().SetNdivisions(303)
            d.GetXaxis().SetNdivisions(303)
        else:
            h.GetXaxis().SetNdivisions(505)
            d.GetXaxis().SetNdivisions(505)

        h.GetXaxis().SetTickSize()
        h.GetXaxis().SetTickLength()
        d.GetXaxis().SetTickSize()
        d.GetXaxis().SetTickLength()
        h.GetYaxis().SetTickSize()
        h.GetYaxis().SetTickLength()
        d.GetYaxis().SetTickSize()
        d.GetYaxis().SetTickLength()
        d.SetLineWidth(2)
        h.Draw("HIST")

        err.SetLineColor(ROOT.kRed)
        err.SetLineWidth(2)
        err.SetMarkerStyle(0)
        err.SetFillColorAlpha(ROOT.kPink + 1, 0.3)
        err.Draw("E2 SAME")
        d.SetLineColor(ROOT.kBlack)
        d.Draw("SAME")

        if row == 0 and col == 0:
            leg = ROOT.TLegend(.6,.4)
            leg.SetBorderSize(0)
            leg.SetTextSize(legend_text_size*3)
            leg.AddEntry(h,"Monte Carlo Prediction","l")
            leg.AddEntry(d,"Data","p")
            leg.DrawClone()

        if row == 1 and col == 0:
            latex = ROOT.TLatex()
            latex.SetNDC()                # Use normalized device coordinates (0–1)
            latex.SetTextSize(legend_text_size*3)       # Adjust text size
            latex.SetTextAlign(13)        # Align left, vertically centered

            x, y = 0.45, 0.65             # Starting coordinates (NDC)
            if penalty == 0:
                latex.DrawLatex(x, y, "#chi^{2} = "+"{:.2f}".format(chi2))
            else:
                latex.DrawLatex(x, y, "#chi^{2} = "+"{:.2f} + {:.2f}".format(chi2,penalty))

        bottom.cd()
        
        nullErrors = h.GetTotalError(False, True, False) #The second "true" makes this fractional error, the third "true" makes this cov area normalized
        for whichBin in range(0, nullErrors.GetXaxis().GetNbins()+1): 
            nullErrors.SetBinError(whichBin, max(nullErrors.GetBinContent(whichBin), 1e-9))
            nullErrors.SetBinContent(whichBin, 1)

        #Error envelope for the MC
        nullErrors.SetLineWidth(0)
        nullErrors.SetMarkerStyle(0)
        nullErrors.SetFillColorAlpha(ROOT.kPink + 1, 0.3)
        #nullErrors.SetMinimum(.5) TODO
        nullErrors.SetMinimum(0.6)
        #nullErrors.SetMaximum(1.5) TODO
        nullErrors.SetMaximum(1.85)

        nullErrors.SetTitle("")
        nullErrors.GetYaxis().SetTitle("Data/MC")
        nullErrors.GetYaxis().SetNdivisions(505)

        ratio = d.Clone(f"ratio_{i}")
        ratio.Divide(ratio,h)

        straightLine = nullErrors.Clone()
        straightLine.SetLineColor(ROOT.kRed)
        straightLine.SetLineWidth(2)
        straightLine.SetFillColor(0)
        straightLine.DrawClone("HIST SAME")

        nullErrors.DrawClone("SAME E2")
        ratio.DrawClone("SAME E1")

    c.Update()
    return c

def GetSliceIndices(fname,exclude,keys):
    bin_config = {}
    with open(fname, "r") as file:
        bin_config = json.load(file)

    sliceInds = []
    for h in keys:
        if h not in bin_config.keys():
            continue
        if not checkRemove(exclude,h):
            sliceInds.extend(list(range(bin_config[h]["start"],bin_config[h]["end"]+1)))
    return(sliceInds)

def slicer(arr,inds,axis=None):
    if arr.ndim == 1: # 1D array
        return arr[inds]
    elif arr.ndim == 2: # 2D array
        if axis == 1:
            return arr[:,inds]
        elif axis == 0:
            return arr[inds,:]
        else:
            ret = arr[:,inds]
            ret = ret[inds,:]
            return ret

def checkRemove(exclude,h):
    if "fhc" in exclude:
        if "fhc" in h and "selection" in h:
            return(True)
    if "rhc" in exclude:
        if "rhc" in h and "selection" in h:
            return(True)
    if "numu" in exclude and "numu" in h:
        return(True)
    if "nue" in exclude and "nue" in h:
        return(True)
    if "elastic" in exclude and "elastic" in h:
        return(True)
    if "imd" in exclude and "imd" in h:
        return(True)
    if "ratio" in exclude and "ratio" in h:
        return(True)

    return(False)

def TMatrix_to_Numpy(matrix):
    rows = matrix.GetNrows()
    cols = matrix.GetNcols()
    arr = np.empty((rows, cols), dtype=np.float64)
    for i in range(rows):
        for j in range(cols):
            arr[i, j] = matrix(i, j)
    return arr

def PlotWithRatio(MNVPLOTTER,plotName,h_cv,hists=[],titles=[],colors=[],styles=[]):
    canvas = ROOT.TCanvas()
    margin = .12
    bottomFraction = .2
    top = ROOT.TPad("DATAMC", "DATAMC", 0, bottomFraction, 1, 1)
    bottom = ROOT.TPad("Ratio", "Ratio", 0, 0, 1, bottomFraction+margin)

    top.Draw()
    bottom.Draw()

    top.cd()
    
    h_cv.Draw("hist")

    null = h_cv.GetCVHistoWithError()
    null.SetLineColor(ROOT.kBlack)
    null.SetLineWidth(2)
    null.SetMarkerStyle(0)
    null.SetFillColorAlpha(ROOT.kGray+1, 0.4)
    null.Draw("E2 SAME")

    leg = ROOT.TLegend(.4,.3)
    leg.SetTextSize(legend_text_size);
    leg.AddEntry(h_cv,"CV Prediction","l")

    for i,hist in enumerate(hists):
        if i < len(titles):
            hist.SetTitle(titles[i])
        if i < len(colors):
                hist.SetLineColor(colors[i])
        if i < len(styles):
            hist.SetLineStyle(styles[i])
        hist.Draw("hist same")
        leg.AddEntry(hist,hist.GetTitle(),"l")

    leg.Draw()

    nullRatio =  h_cv.Clone()
    nullRatio.Divide(nullRatio,h_cv)

    bottom.cd()
    bottom.SetTopMargin(0)
    bottom.SetBottomMargin(0.3)

    nullErrors = h_cv.GetTotalError(False, True, False) #The second "true" makes this fractional error, the third "true" makes this cov area normalized
    for whichBin in range(0, nullErrors.GetXaxis().GetNbins()+1): 
        nullErrors.SetBinError(whichBin, max(nullErrors.GetBinContent(whichBin), 1e-9))
        nullErrors.SetBinContent(whichBin, 1)

    nullRatio.SetLineColor(ROOT.kBlack)
    nullRatio.SetLineWidth(3)

    #Error envelope for the MC
    nullErrors.SetLineWidth(0)
    nullErrors.SetMarkerStyle(0)
    nullErrors.SetFillColorAlpha(ROOT.kGray+1, 0.4)
    nullErrors.GetYaxis().SetTitle("Ratio to CV")
    RatioAxis(nullErrors,MNVPLOTTER)
    nullErrors.SetMinimum(0.7)
    nullErrors.SetMaximum(1.3)
    nullErrors.GetXaxis().SetTitle("Neutrino Energy")
    nullErrors.Draw("E2")

    #Draw the data ratios
    nullRatio.Draw("same hist")

    for hist in hists:
        i_ratio = hist.Clone()
        i_ratio.Divide(i_ratio,h_cv)
        i_ratio.DrawClone("same hist")

    straightLine = nullErrors.Clone()
    straightLine.SetLineColor(ROOT.kBlack)
    straightLine.SetLineWidth(2)
    straightLine.SetFillColor(0)
    straightLine.Draw("HIST SAME")
    canvas.Print(plotName)

def RatioAxis(hist,MNVPLOTTER):
    #MNVPLOTTER.ApplyAxisStyle(hist)

    hist.SetTitle("")
    hist.GetYaxis().SetNdivisions(304) #5 minor divisions between 5 major divisions.  I'm trying to match a specific paper here.
    hist.GetYaxis().SetTitleOffset(1.5)
    hist.GetYaxis().SetTitleFont(43)
    hist.GetYaxis().SetTitleSize(26)
    hist.GetYaxis().SetLabelFont(43)
    hist.GetYaxis().SetLabelSize(34)
    hist.GetYaxis().CenterTitle()
    hist.GetXaxis().SetNdivisions(510) #5 minor divisions between 9 major divisions.  I'm trying to match a specific paper here.
    hist.GetXaxis().SetTitleOffset(2.5)
    hist.GetXaxis().SetTitleFont(43)
    hist.GetXaxis().SetTitleSize(26)
    hist.GetXaxis().SetLabelFont(43)
    hist.GetXaxis().SetLabelSize(34)
    hist.GetXaxis().CenterTitle()

    #hist.GetYaxis().SetLabelSize(.13)
    #hist.GetYaxis().SetTitleSize(0.1)
    #hist.GetYaxis().SetTitleOffset(0.6)
    #hist.GetYaxis().SetNdivisions(505) #5 minor divisions between 5 major divisions.  I'm trying to match a specific paper here.
    #hist.GetXaxis().SetTitleSize(0.16)
    #hist.GetXaxis().SetTitleOffset(0.9)
    #hist.GetXaxis().SetLabelSize(.15)

def CanvasPartition(C, Nx, Ny, lMargin, rMargin, bMargin, tMargin):
    vSpacing = 0.0
    vStep = (1. - bMargin - tMargin - (Ny - 1) * vSpacing) / Ny

    hSpacing = 0.0
    hStep = (1. - lMargin - rMargin - (Nx - 1) * hSpacing) / Nx

    vposd, vposu, vmard, vmaru, vfactor = 0, 0, 0, 0, 0
    hposl, hposr, hmarl, hmarr, hfactor = 0, 0, 0, 0, 0

    for i in range(Nx):
        if (i == 0):
            hposl = 0.0
            hposr = lMargin + hStep
            hfactor = hposr - hposl
            hmarl = lMargin / hfactor
            hmarr = 0.0
        elif (i == Nx - 1):
            hposl = hposr + hSpacing
            hposr = hposl + hStep + rMargin
            hfactor = hposr - hposl
            hmarl = 0.0
            hmarr = rMargin / (hposr - hposl)
        else:
            hposl = hposr + hSpacing
            hposr = hposl + hStep
            hfactor = hposr - hposl
            hmarl = 0.0
            hmarr = 0.0

        for j in range(Ny):
            if (j == 0):
                vposd = 0.0
                vposu = bMargin + vStep
                vfactor = vposu - vposd
                vmard = bMargin / vfactor
                vmaru = 0.0
            elif (j == Ny - 1):
                vposd = vposu + vSpacing
                vposu = vposd + vStep + tMargin
                vfactor = vposu - vposd
                vmard = 0.0
                vmaru = tMargin / (vposu - vposd)
            else:
                vposd = vposu + vSpacing
                vposu = vposd + vStep
                vfactor = vposu - vposd
                vmard = 0.0
                vmaru = 0.0

            C.cd(0)

            name = "pad_{}_{}".format(i, j)
            pad = ROOT.TPad(name, "", hposl, vposd, hposr, vposu)
            pad.SetLeftMargin(hmarl)
            pad.SetRightMargin(hmarr)
            pad.SetBottomMargin(vmard)
            pad.SetTopMargin(vmaru)

            pad.SetFrameBorderMode(0)
            pad.SetBorderMode(0)
            pad.SetBorderSize(0)

            pad.Draw()
 
def XtoPad(x):
   xl = ctypes.c_double(0)
   yl = ctypes.c_double(0)
   xu = ctypes.c_double(0)
   yu = ctypes.c_double(0)
   ROOT.gPad.GetPadPar(xl, yl, xu, yu)
   pw = ctypes.c_double(xu.value - xl.value)
   lm = ROOT.gPad.GetLeftMargin()
   rm = ROOT.gPad.GetRightMargin()
   fw = pw.value - pw.value * lm - pw.value * rm
   return (x * fw + pw.value * lm) / pw.value
 
def YtoPad(y):
   xl = ctypes.c_double(0)
   yl = ctypes.c_double(0)
   xu = ctypes.c_double(0)
   yu = ctypes.c_double(0)
   ROOT.gPad.GetPadPar(xl, yl, xu, yu)
   ph = ctypes.c_double(yu.value - yl.value)
   tm = ROOT.gPad.GetTopMargin()
   bm = ROOT.gPad.GetBottomMargin()
   fh = ph.value - ph.value * bm - ph.value * tm
   return (y * fh + bm * ph.value) / ph.value
