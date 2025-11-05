void sum_flux_by_parents(TFile* f) {
  // TFile* f = new TFile(fname, "UPDATE");  // Use "UPDATE" to write new histos
  if (!f || f->IsZombie()) {
    std::cerr << "Failed to open ROOT file." << std::endl;
    return;
  }

  // List of parent PDGs you used when producing the histograms
  std::vector<int> parent_pdgs = {130, 311, 13, -13, 321, -321, 211, -211, 999999};

  PlotUtils::MnvH1D* total = nullptr;

  for (int pdg : parent_pdgs) {
    TString histname = Form("flux_E_cvweighted_parent%d", pdg);
    auto* h = (PlotUtils::MnvH1D*)f->Get(histname);
    if (!h) {
      std::cerr << "Missing: " << histname << std::endl;
      continue;
    }

    if (!total) {
      total = new PlotUtils::MnvH1D(*h);  // First one is the seed
      total->SetName("flux_E_cvweighted");
      total->SetTitle("Total Flux (sum of all parents);E_{#nu} (GeV);Flux / m^{2} / POT / GeV");
    } else {
      total->Add(h);  // Add subsequent histograms
    }
  }

  if (total) {
    f->cd();
    total->Write();  // Will overwrite if exists
    std::cout << "Wrote total flux histogram: flux_E_cvweighted_total" << std::endl;
  } else {
    std::cerr << "No histograms were added. Check input names." << std::endl;
  }

  // f->Close();
}

TH1D* RebinFluxPreservingDensity(TH1D* h1, const TH1D* h_template) {
  if (!h1 || !h_template) return nullptr;

  // Get bin edges from the template histogram
  const TArrayD* binEdges = h_template->GetXaxis()->GetXbins();
  if (binEdges->GetSize() == 0) {
    std::cerr << "Error: Template histogram does not use variable binning." << std::endl;
    return nullptr;
  }

  int nNewBins = binEdges->GetSize() - 1;

  // Create rebinned histogram
  TH1D* h1_rebinned = new TH1D(
    Form("%s_rebinned", h1->GetName()), 
    h1->GetTitle(),
    nNewBins, 
    binEdges->GetArray()
  );

  // Loop over new bins
  for (int i = 1; i <= nNewBins; ++i) {
    double low = h1_rebinned->GetXaxis()->GetBinLowEdge(i);
    double high = h1_rebinned->GetXaxis()->GetBinUpEdge(i);

    int binLow = h1->FindBin(low);
    int binHigh = h1->FindBin(high - 1e-6);  // avoid fencepost issue

    double integral = 0;
    double err2 = 0;
    for (int j = binLow; j <= binHigh; ++j) {
      double content = h1->GetBinContent(j);
      double width   = h1->GetBinWidth(j);
      double error   = h1->GetBinError(j);
      integral += content * width;
      err2     += std::pow(error * width, 2);
    }

    double width = high - low;
    h1_rebinned->SetBinContent(i, integral / width);
    h1_rebinned->SetBinError(i, std::sqrt(err2) / width);
  }

  return h1_rebinned;
}



void plot_two_fluxes() {
  // Load the files
  TFile* f_numu0 = new TFile("flux-gen2thin-pdg14-minerva1.root", "READ");
  // TFile* f_numu  = new TFile("flux-pdg14-minerva1.root", "UPDATE");
  // TFile* f_numu  = new TFile("le1_mu.root", "UPDATE");
  TFile* f_numu = new TFile("flux_numu_wNUE_CONSTRAINT_noErr.root", "UPDATE");
  // TFile* f_numu_nErr = new TFile("flux_numu_wNUE_CONSTRAINT_noErr.root", "UPDATE");
  // TFile* f_nue   = new TFile("flux-pdg12-minerva1.root", "UPDATE");
  TFile* f_nue  = new TFile("flux_nue_wNUE_CONSTRAINT_noErr.root", "UPDATE");
  // TFile* f_nue_nErr  = new TFile("flux_nue_wNUE_CONSTRAINT_noErr.root", "UPDATE");

  // sum_flux_by_parents(f_numu_wErr);
  // sum_flux_by_parents(f_numu_nErr);
  // sum_flux_by_parents(f_nue_wErr);
  // sum_flux_by_parents(f_nue_nErr);

  // Get MnvH1D
  PlotUtils::MnvH1D* h_numu0 = (PlotUtils::MnvH1D*)f_numu0->Get("flux_E_cvweighted");
  PlotUtils::MnvH1D* h_numu  = (PlotUtils::MnvH1D*)f_numu->Get("flux_E_cvweighted");
  // PlotUtils::MnvH1D* h_numu_nErr  = (PlotUtils::MnvH1D*)f_numu_nErr->Get("flux_E_cvweighted");
  // PlotUtils::MnvH1D* h_nue_nErr   = (PlotUtils::MnvH1D*)f_nue_nErr->Get("flux_E_cvweighted");
  PlotUtils::MnvH1D* h_nue   = (PlotUtils::MnvH1D*)f_nue->Get("flux_E_cvweighted");

  // if (!h_numu0 || !h_numu_wErr || !h_numu_nErr || !h_nue_wErr || !h_nue_nErr) {
  //     std::cerr << "❌ Failed to load one of the MnvH1D histograms from file!" << std::endl;
  //     return;
  // }

  // Convert to TH1D
  TH1D* h0 = (TH1D*)h_numu0->GetCVHistoWithError().Clone("h0");
  TH1D* h1 = (TH1D*)h_numu->GetCVHistoWithError().Clone("h1");
  // TH1D* h1n = (TH1D*)h_numu_nErr->GetCVHistoWithError().Clone("h1n");
  TH1D* h2 = (TH1D*)h_nue->GetCVHistoWithError().Clone("h2");
  // TH1D* h2n = (TH1D*)h_nue_nErr->GetCVHistoWithError().Clone("h2n");

  const TArrayD* binEdges = h0->GetXaxis()->GetXbins();

  if (binEdges->GetSize() > 0) {
    for(int i=0; i<binEdges->GetSize(); i++){
      std::cout << (*binEdges)[i] << "\n";
    }
  }

  // TH1D* h1_rebinned = RebinFluxPreservingDensity(h1, h0);
  // TH1D* h2_rebinned = RebinFluxPreservingDensity(h2, h0);

  std::cout << h0->GetNbinsX() << "\t" << h1->GetNbinsX() << "\t" << h2->GetNbinsX() << "\n";

  // Styling
  h0->SetLineColor(kBlack); h0->SetLineWidth(2);
  h1->SetLineColor(kRed);   h1->SetLineWidth(2);
  // h1_rebinned->SetLineColor(kRed);   h1_rebinned->SetLineWidth(2);
  h2->SetLineColor(kBlue);  h2->SetLineWidth(2);
  // h2_rebinned->SetLineColor(kBlue);  h2_rebinned->SetLineWidth(2);

  h0->SetTitle("Fluxes for #nu_{#mu} and #nu_{e};E_{#nu} [GeV];Flux");
  // h0->SetTitle("Fluxes for #bar{#nu_{#mu}} and #bar{#nu_{e}};E_{#nu} [GeV];Flux");

  // // Draw
  // TCanvas* c = new TCanvas("c", "Flux Comparison", 800, 600);
  // c->SetLogy();
  // h0->Draw("HIST");
  // h1->Draw("HIST SAME");
  // h2->Draw("HIST SAME");

  // TLegend* leg = new TLegend(0.6, 0.7, 0.88, 0.88);
  // leg->AddEntry(h0, "existing #nu_{#mu} flux", "l");
  // leg->AddEntry(h1, "#nu_{#mu} flux", "l");
  // leg->AddEntry(h2, "#nu_{e} flux", "l");
  // leg->Draw();

  // c->SaveAs("flux_comparison_logy.png");

  // Create canvas with 2 pads
  TCanvas* c = new TCanvas("c", "Flux Comparison + Ratio", 800, 800);
  c->Divide(1, 2, 0.01, 0.01);  // vertical split

  // === Upper Pad for Fluxes ===
  TPad* topPad = (TPad*)c->cd(1);
  topPad->SetPad(0, 0.3, 1, 1);
  topPad->SetBottomMargin(0.02);
  topPad->SetLogy();
  h0->Draw("HIST");
  h1->Draw("HIST SAME");
  h2->Draw("HIST SAME");

  TLegend* leg = new TLegend(0.6, 0.7, 0.88, 0.88);
  leg->AddEntry(h0, "existing #bar{#nu_{#mu}} flux", "l");
  leg->AddEntry(h1, "new #bar{#nu_{#mu}} flux", "l");
  leg->AddEntry(h2, "#bar{#nu_{e}} flux", "l");
  leg->Draw();

  // === Bottom Pad for Ratio ===
  TPad* bottomPad = (TPad*)c->cd(2);
  bottomPad->SetPad(0, 0.0, 1, 0.3);
  bottomPad->SetTopMargin(0.05);
  bottomPad->SetBottomMargin(0.3);
  bottomPad->SetGridy();

  TH1D* ratio = (TH1D*)h1->Clone("ratio_h1_over_h0");
  ratio->Divide(h0);  // h1 / h0
  // ratio->SetMinimum(0.8);
  // ratio->SetMaximum(1.4);
  ratio->SetLineColor(kRed + 1);
  ratio->SetLineWidth(2);
  ratio->SetTitle("");
  ratio->GetYaxis()->SetTitle("h_{1} / h_{0}");
  ratio->GetYaxis()->SetTitleSize(0.09);
  ratio->GetYaxis()->SetTitleOffset(0.45);
  ratio->GetYaxis()->SetLabelSize(0.08);
  ratio->GetYaxis()->SetNdivisions(505);
  ratio->GetXaxis()->SetTitle("E_{#nu} [GeV]");
  ratio->GetXaxis()->SetTitleSize(0.09);
  ratio->GetXaxis()->SetLabelSize(0.08);
  ratio->SetStats(0);
  ratio->Draw("HIST");

  //Save both
  c->SaveAs("flux_comparison_with_ratio_FHC_nErr.png");
  c->Close();
  f_numu0->Close();
  f_numu->Close();
  f_nue->Close();
}
