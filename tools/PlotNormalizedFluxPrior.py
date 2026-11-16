"""Plot normalized-prior fits with fitData's stitched sample and energy axes.

Public entry point: PlotNormalizedFluxPrior(...). Predictions and the prior
are supplied by the fitter. This module only evaluates that same prior to
retrieve the fitted shifts; it does not optimize or choose flux modes.
ROOT and the existing PlotHistogram helpers are imported only when drawing.
"""
from pathlib import Path
import numpy as np


def fit_plot_values(prior, data, null_pred, best_pred):
    """Use the fitted normalized prior, never the legacy raw-universe solver."""
    null_result = prior.evaluate(data - null_pred)
    best_result = prior.evaluate(data - best_pred)
    return {
        "null": np.asarray(null_pred).copy(),
        "null_profiled": null_pred + null_result[3],
        "best_profiled": best_pred + best_result[3],
        "null_result": null_result,
        "best_result": best_result,
        # C = L L^T includes data/MC statistics and non-flux systematics.
        # This is a residual-error band, not a fitted-model uncertainty.
        "nonflux_sigma": np.sqrt(np.sum(prior.L**2, axis=1)),
    }


def PlotNormalizedFluxPrior(histogram, prior, data, null_pred, best_pred, params,
                     hist_config, output_dir, plot_tag, ntuple_tag, flux_source,
                     flux_scale):
    """fitData-style stitched spectrum and ratio using this fit's corrections.

    Uses the same layout and physical-energy-axis helpers as fitData.
    The supplied prior must be the objective used by the fit; this function
    never calls the legacy raw-universe profiling solver.
    No input histogram is modified, and no fit is repeated here.
    """
    import ROOT
    from tools.PlotHistogram import (GetStitchedSampleLayout,
                                     DrawStitchedSampleLayout,
                                     DrawStitchedEnergyAxes)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    values = fit_plot_values(prior, data, null_pred, best_pred)
    layout = GetStitchedSampleLayout(str(hist_config))
    template = histogram.GetMCHistogram().GetCVHistoWithStatError()
    data_hist = histogram.GetDataHistogram().GetCVHistoWithStatError()
    data_hist.SetDirectory(0)
    data_hist.SetName("normalized_prior_plot_data")
    data_hist.SetTitle("")
    data_hist.SetStats(False)
    data_hist.SetMarkerStyle(20)
    data_hist.SetMarkerSize(.65)
    data_hist.SetMarkerColor(ROOT.kBlack)
    data_hist.SetLineColor(ROOT.kBlack)

    def make_hist(name, contents, color):
        hist = template.Clone(name)
        hist.SetDirectory(0)
        hist.Reset()
        hist.SetTitle("")
        hist.GetXaxis().SetTitle("")
        hist.SetStats(False)
        for i, content in enumerate(contents, 1):
            hist.SetBinContent(i, float(content))
            hist.SetBinError(i, 0.)
        hist.SetLineColor(color)
        hist.SetLineWidth(3)
        hist.SetFillStyle(0)
        hist.SetMarkerStyle(0)
        return hist

    raw = make_hist("normalized_prior_raw_null", values["null"], ROOT.kRed + 1)
    prof_null = make_hist("normalized_prior_profiled_null", values["null_profiled"],
                          ROOT.kGreen + 2)
    prof_best = make_hist("normalized_prior_profiled_best", values["best_profiled"],
                          ROOT.kBlue + 1)
    canvas = ROOT.TCanvas("normalized_prior_fit", "Normalized flux-prior fit", 1500, 900)
    top = ROOT.TPad("normalized_prior_top", "", 0., .31, 1., 1.)
    bottom = ROOT.TPad("normalized_prior_bottom", "", 0., 0., 1., .31)
    top.SetLeftMargin(.085)
    top.SetRightMargin(.025)
    top.SetTopMargin(.15)
    top.SetBottomMargin(.025)
    bottom.SetLeftMargin(.085)
    bottom.SetRightMargin(.025)
    bottom.SetTopMargin(.025)
    bottom.SetBottomMargin(.34)
    top.Draw()
    bottom.Draw()
    top.cd()

    plotted = np.r_[data, values["null"], values["null_profiled"],
                    values["best_profiled"]]
    if not np.isfinite(plotted).all():
        raise ValueError("Cannot plot non-finite fitted spectra")
    # Keep genuine nonpositive contents visible instead of clipping the fit.
    log_y = bool(np.all(plotted > 0))
    top.SetLogy(log_y)
    errors = np.array([data_hist.GetBinError(i) for i in range(1, len(data)+1)])
    ymax = max(float(np.max(plotted)), float(np.max(data + errors)))
    ymin = min(float(np.min(plotted)), float(np.min(data - errors)))
    # Reserve space below the spectra for the two larger fit-result lines.
    raw.SetMinimum(float(np.min(plotted)) * .08 if log_y else min(0., ymin) - .05*max(ymax-ymin, 1.))
    raw.SetMaximum(ymax * (4. if log_y else 1.35))
    raw.GetYaxis().SetTitle("Stitched-bin content")
    raw.GetYaxis().SetTitleSize(.045)
    raw.GetYaxis().SetTitleOffset(1.05)
    raw.GetYaxis().SetLabelSize(.04)
    raw.GetXaxis().SetLabelSize(0.)
    raw.GetXaxis().SetTickLength(0.)
    raw.Draw("HIST")
    prof_null.Draw("HIST SAME")
    prof_best.Draw("HIST SAME")
    data_hist.Draw("E1 SAME")
    # Hold all drawn ROOT objects until SaveAs completes.
    objects = [template, data_hist, raw, prof_null, prof_best, top, bottom]
    objects += DrawStitchedSampleLayout(top, raw, layout)
    # Keep the shared sample mapping, but give this plot larger header labels.
    xaxis = raw.GetXaxis()
    xmin, xmax = xaxis.GetXmin(), xaxis.GetXmax()
    for sample in layout:
        xpos = top.GetLeftMargin() + ((sample["center"]-xmin)/(xmax-xmin)) * (
            1.-top.GetLeftMargin()-top.GetRightMargin())
        label = ROOT.TLatex()
        label.SetNDC(True)
        label.SetTextAlign(22)
        label.SetTextFont(62)
        label.SetTextSize(.042)
        objects.append(label.DrawLatex(xpos, .91, sample["label"]))
        objects.append(label)

    null_q, null_res, null_pen, _ = values["null_result"]
    best_q, best_res, best_pen, _ = values["best_result"]
    legend = ROOT.TLegend(.12, .51, .64, .84)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(.038)
    legend.AddEntry(data_hist, "Data (statistical errors)", "lep")
    legend.AddEntry(raw, "Null before flux profiling", "l")
    legend.AddEntry(prof_null,
                    "Profiled null: #chi^{2} = %.2f + %.2f penalty" % (null_res, null_pen), "l")
    legend.AddEntry(prof_best,
                    "Profiled best fit: #chi^{2} = %.2f + %.2f penalty" % (best_res, best_pen), "l")
    legend.Draw()
    objects.append(legend)
    text = ROOT.TLatex()
    text.SetNDC(True)
    text.SetTextSize(.031)
    text.DrawLatex(.10, .965, "%s: normalized %s flux prior, K=%d, scale=%g" %
                   (ntuple_tag, flux_source, prior.K, flux_scale))
    text.SetTextSize(.038)
    text.DrawLatex(.12, .105,
                   "#Delta m^{2}=%.3f eV^{2}, |U_{e4}|^{2}=%.4f, |U_{#mu4}|^{2}=%.4f, |U_{#tau4}|^{2}=%.4f" %
                   tuple(params[key] for key in ("m", "ue4", "umu4", "utau4")))
    text.DrawLatex(.12, .045, "#Delta#chi^{2}_{null-best} = %.3f" % (null_q-best_q))
    objects.append(text)

    bottom.cd()
    axis = make_hist("normalized_prior_ratio_axis", np.ones(len(data)), ROOT.kRed + 1)
    axis.SetMinimum(.5)
    axis.SetMaximum(1.5)
    axis.GetYaxis().SetTitle("Ratio to raw null")
    axis.GetYaxis().SetTitleSize(.085)
    axis.GetYaxis().SetTitleOffset(.48)
    axis.GetYaxis().SetLabelSize(.075)
    axis.GetYaxis().SetNdivisions(505)
    axis.GetXaxis().SetLabelSize(0.)
    axis.GetXaxis().SetTickLength(0.)
    axis.Draw("HIST")
    objects.append(axis)
    valid = np.flatnonzero(values["null"] != 0.)
    if len(valid) != len(data):
        print("Ratio plot omits bins with a zero raw-null prediction.")

    def ratio_graph(contents, color, errors=None, half_width=False):
        graph = ROOT.TGraphErrors(len(valid))
        graph.SetLineColor(color)
        graph.SetLineWidth(2)
        graph.SetMarkerColor(color)
        graph.SetMarkerStyle(20)
        graph.SetMarkerSize(.55)
        for j, i in enumerate(valid):
            denom = values["null"][i]
            graph.SetPoint(j, template.GetXaxis().GetBinCenter(int(i)+1), float(contents[i]/denom))
            graph.SetPointError(j, template.GetXaxis().GetBinWidth(int(i)+1)/2 if half_width else 0.,
                                float(errors[i]/abs(denom)) if errors is not None else 0.)
        objects.append(graph)
        return graph

    ratio_graph(values["null_profiled"], ROOT.kGreen + 2).Draw("LP SAME")
    ratio_graph(values["best_profiled"], ROOT.kBlue + 1).Draw("LP SAME")
    ratio_graph(data, ROOT.kBlack, errors).Draw("P SAME")
    objects += DrawStitchedSampleLayout(bottom, axis, layout)
    objects += DrawStitchedEnergyAxes(bottom, axis, layout)
    canvas.cd()
    canvas.Update()
    scale_tag = format(flux_scale, "g").replace(".", "p").replace("-", "m")
    name = (f"{ntuple_tag}_stitched_{plot_tag}_normalizedFluxPrior_{flux_source}"
            f"_K{prior.K}_s{scale_tag}_{params['m']:.3f}_{params['ue4']:.4f}"
            f"_{params['umu4']:.4f}_{params['utau4']:.4f}.png")
    path = output_dir / name
    canvas.SaveAs(str(path))
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"ROOT did not create the plot: {path}")
    canvas.Close()
    print(f"Saved normalized-prior fit plot: {path}")
    return path

