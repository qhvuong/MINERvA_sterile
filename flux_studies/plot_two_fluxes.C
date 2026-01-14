void sum_flux_by_parents(TFile* f) {
  // TFile* f = new TFile(fname, "UPDATE");  // Use "UPDATE" to write new histos
  if (!f || f->IsZombie()) {
    std::cerr << "Failed to open ROOT file." << std::endl;
    return;
  }

  // List of parent PDGs you used when producing the histograms
  std::vector<int> parent_pdgs = {130, 311, 13, -13, 321, -321, 211, -211, 999999};

  PlotUtils::MnvH1D* total = nullptr;

  // for (int pdg : parent_pdgs) {
  //   TString histname = Form("flux_E_cvweighted_parent%d", pdg);
  //   auto* h = (PlotUtils::MnvH1D*)f->Get(histname);
  //   if (!h) {
  //     std::cerr << "Missing: " << histname << std::endl;
  //     continue;
  //   }

  //   if (!total) {
  //     total = new PlotUtils::MnvH1D(*h);  // First one is the seed
  //     total->SetName("flux_E_cvweighted");
  //     total->SetTitle("Total Flux (sum of all parents);E_{#nu} (GeV);Flux / m^{2} / POT / GeV");
  //   } else {
  //     total->Add(h);  // Add subsequent histograms
  //   }
  // }

  for (int pdg : parent_pdgs) {
    TString histname = Form("flux_E_unweighted_parent%d", pdg);
    auto* h = (PlotUtils::MnvH1D*)f->Get(histname);
    if (!h) {
      std::cerr << "Missing: " << histname << std::endl;
      continue;
    }

    if (!total) {
      total = new PlotUtils::MnvH1D(*h);  // First one is the seed
      total->SetName("flux_E_unweighted");
      total->SetTitle("Total Flux (sum of all parents);E_{#nu} (GeV);Flux / m^{2} / POT / GeV");
    } else {
      total->Add(h);  // Add subsequent histograms
    }
  }

  // for (int pdg : parent_pdgs) {
  //   TString histname = Form("eventcount_E_unweighted_parent%d", pdg);
  //   auto* h = (PlotUtils::MnvH1D*)f->Get(histname);
  //   if (!h) {
  //     std::cerr << "Missing: " << histname << std::endl;
  //     continue;
  //   }

  //   if (!total) {
  //     total = new PlotUtils::MnvH1D(*h);  // First one is the seed
  //     total->SetName("eventcount_E_unweighted");
  //     total->SetTitle("Event Count unweighted (sum of all parents);E_{#nu} (GeV);Entries");
  //   } 
  //   else {
  //     total->Add(h);  // Add subsequent histograms
  //   }
  // }



  if (total) {
    f->cd();
    total->Write();  // Will overwrite if exists
    std::cout << "Wrote total flux histogram: flux(rate/count)_E_cv(un)weighted_total" << std::endl;
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

TH1D* MakeRatioWithStatError(const TH1D* h1, const TH1D* h0,
                             const char* name = "ratio_h1_over_h0")
{
  if (!h1 || !h0) {
    std::cerr << "[MakeRatioWithStatError] Null histogram pointer!" << std::endl;
    return nullptr;
  }

  TH1D* ratio = (TH1D*)h1->Clone(name);
  ratio->Reset("ICES");  // clear contents, errors, stats but keep binning

  const int nb = h1->GetNbinsX();
  for (int i = 1; i <= nb; ++i)
  {
    const double A  = h1->GetBinContent(i);
    const double B  = h0->GetBinContent(i);
    const double sA = h1->GetBinError(i);

    double R = 0.0, sR = 0.0;
    if (B > 0) {
      R = A / B;
      if (A > 0) sR = sA / B;  // propagate only numerator error
    }

    ratio->SetBinContent(i, R);
    ratio->SetBinError(i, sR);
  }

  // Optional: Copy axis titles from h1
  ratio->GetXaxis()->SetTitle(h1->GetXaxis()->GetTitle());
  ratio->GetYaxis()->SetTitle("h_{1} / h_{0}");

  return ratio;
}


void plot_two_fluxes() {
  // Load the files
  TFile* f_numu0_cvweighted = new TFile("/exp/minerva/app/users/qvuong/MAT_AL9/opt/lib/data/flux/flux-gen2thin-pdg14-minerva1.root", "READ");
  TFile* f_numu0_unweighted = new TFile("/exp/minerva/app/users/qvuong/MAT_AL9/opt/lib/data/flux/flux-g4numiv5-pdg14-minerva1.root", "READ");
  // TFile* f_numu0_cvweighted = new TFile("/exp/minerva/app/users/qvuong/MAT_AL9/opt/lib/data/flux/flux-gen2thin-pdg-14-minerva5.root", "READ");
  // TFile* f_numu0_unweighted  = new TFile("/exp/minerva/app/users/qvuong/MAT_AL9/opt/lib/data/flux/flux-g4numiv5-pdg-14-minerva5.root", "READ");  // TFile* f_numu  = new TFile("flux-pdg14-minerva1.root", "UPDATE");

  // TFile* f_numu    = new TFile("LE5_numu.root", "UPDATE");
  TFile* f_numu = new TFile("LE1_numubar.root", "UPDATE");

  // TFile* f_nue    = new TFile("LE5_nue.root", "UPDATE");
  TFile* f_nue = new TFile("LE1_nuebar.root", "UPDATE");

  // sum_flux_by_parents(f_numu);
  // sum_flux_by_parents(f_numubar);
  // sum_flux_by_parents(f_nue);
  // sum_flux_by_parents(f_nuebar);

  // Get MnvH1D
  PlotUtils::MnvH1D* h_numu0_cvweighted = (PlotUtils::MnvH1D*)f_numu0_cvweighted->Get("flux_E_cvweighted");
  PlotUtils::MnvH1D* h_numu0_unweighted = (PlotUtils::MnvH1D*)f_numu0_unweighted->Get("flux_E_unweighted");
  PlotUtils::MnvH1D* h_numu_cvweighted  = (PlotUtils::MnvH1D*)f_numu->Get("flux_E_cvweighted");
  PlotUtils::MnvH1D* h_numu_unweighted  = (PlotUtils::MnvH1D*)f_numu->Get("flux_E_unweighted");
  PlotUtils::MnvH1D* h_numu_eventcount_cvweighted  = (PlotUtils::MnvH1D*)f_numu->Get("eventcount_E_cvweighted");
  PlotUtils::MnvH1D* h_numu_eventcount_unweighted  = (PlotUtils::MnvH1D*)f_numu->Get("eventcount_E_unweighted");
  PlotUtils::MnvH1D* h_nue_cvweighted   = (PlotUtils::MnvH1D*)f_nue->Get("flux_E_cvweighted");
  PlotUtils::MnvH1D* h_nue_unweighted   = (PlotUtils::MnvH1D*)f_nue->Get("flux_E_unweighted");

  // if (!h_numu0 || !h_numu_wErr || !h_numu_nErr || !h_nue_wErr || !h_nue_nErr) {
  //     std::cerr << "❌ Failed to load one of the MnvH1D histograms from file!" << std::endl;
  //     return;
  // }

  // Convert to TH1D
  TH1D* h0_cv   = (TH1D*)h_numu0_cvweighted->GetCVHistoWithError().Clone("h0_cv");
  TH1D* h1_cv   = (TH1D*)h_numu_cvweighted->GetCVHistoWithError().Clone("h1_cv");
  TH1D* h2_cv   = (TH1D*)h_nue_cvweighted->GetCVHistoWithError().Clone("h2_cv");
  TH1D* hevt_cv = (TH1D*)h_numu_eventcount_cvweighted->GetCVHistoWithError().Clone("hevt_cv");


  TH1D* h0 = (TH1D*)h_numu0_unweighted->GetCVHistoWithError().Clone("h0");
  TH1D* h1 = (TH1D*)h_numu_unweighted->GetCVHistoWithError().Clone("h1");
  TH1D* h2 = (TH1D*)h_nue_unweighted->GetCVHistoWithError().Clone("h2");
  TH1D* hevt = (TH1D*)h_numu_eventcount_unweighted->GetCVHistoWithError().Clone("hevt");

  // const TArrayD* binEdges = h0->GetXaxis()->GetXbins();
  // if (binEdges->GetSize() > 0) {
  //   for(int i=0; i<binEdges->GetSize(); i++){
  //     std::cout << (*binEdges)[i] << "\n";
  //   }
  // }

  // Helper: copy Poisson stat from a counts hist into a flux hist’s bin errors
  auto SetFluxErrorsFromCounts = [](TH1* hFlux, const TH1* hCounts){
    const int nb = hFlux->GetNbinsX();
    for (int i = 1; i <= nb; ++i) {
      const double N    = hCounts->GetBinContent(i);
      const double flux = hFlux->GetBinContent(i);
      if (N > 0 && flux > 0) {
        const double scale = flux / N;           // counts -> flux conversion for this bin
        const double err   = scale * std::sqrt(N);
        // const double err   = std::sqrt(N);
        // std::cout << i << "\t" << err << "\n";
        hFlux->SetBinError(i, err);
      } else {
        hFlux->SetBinError(i, 0.0);
      }
    }
  };
  SetFluxErrorsFromCounts(h1_cv, hevt_cv);
  SetFluxErrorsFromCounts(h1, hevt);

  std::cout << h0_cv->GetNbinsX() << "\t" << h1_cv->GetNbinsX() << "\t" << h2_cv->GetNbinsX() << "\n";
  std::cout << h0->GetNbinsX() << "\t" << h1->GetNbinsX() << "\t" << h2->GetNbinsX() << "\n";

  printf("h0 first bin edge = %f\n", h0_cv->GetXaxis()->GetXmin());
  printf("h0 last bin edge  = %f\n", h0_cv->GetXaxis()->GetXmax());


  // Styling
  h0_cv->SetLineColor(kBlack); h0_cv->SetLineWidth(2);
  h1_cv->SetLineColor(kRed);   h1_cv->SetLineWidth(2);
  h2_cv->SetLineColor(kBlue);  h2_cv->SetLineWidth(2);

  h0_cv->SetTitle("Cvweighted fluxes for #nu_{#mu} and #nu_{e};E_{#nu} [GeV];Flux");
  // h0->SetTitle("Fluxes for #bar{#nu_{#mu}} and #bar{#nu_{e}};E_{#nu} [GeV];Flux");

  // Create canvas with 2 pads
  TCanvas* c_cv = new TCanvas("c_cv", "Flux Comparison + Ratio", 800, 800);
  c_cv->Divide(1, 2, 0.01, 0.01);  // vertical split

  // === Upper Pad for Fluxes ===
  TPad* topPad_cv = (TPad*)c_cv->cd(1);
  topPad_cv->SetPad(0, 0.3, 1, 1);
  topPad_cv->SetBottomMargin(0.02);
  topPad_cv->SetLogx();
  topPad_cv->SetLogy();
  // h0_cv->GetXaxis()->SetRangeUser(0.3, 10.0);
  // h1_cv->GetXaxis()->SetRangeUser(0.3, 10.0);
  // h2_cv->GetXaxis()->SetRangeUser(0.3, 10.0);
  h0_cv->Draw("HIST");
  h1_cv->Draw("HIST SAME");
  h2_cv->Draw("HIST SAME");
  h1_cv->Draw("E SAME");
  h0_cv->GetXaxis()->SetRangeUser(0.5, 10.0);
  gPad->Modified();
  gPad->Update();

  TLegend* leg_cv = new TLegend(0.12, 0.1, 0.45, 0.3);
  leg_cv->AddEntry(h0_cv, "existing #nu_{#mu} flux", "l");
  leg_cv->AddEntry(h1_cv, "new #nu_{#mu} flux", "l");
  leg_cv->AddEntry(h2_cv, "new #nu_{e} flux", "l");
  // leg->AddEntry(h0, "existing #bar{#nu_{#mu}} flux", "l");
  // leg->AddEntry(h1, "new #bar{#nu_{#mu}} flux", "l");
  // leg->AddEntry(h2, "new #bar{#nu_{e}} flux", "l");
  leg_cv->Draw();

  // === Bottom Pad for Ratio ===
  TPad* bottomPad_cv = (TPad*)c_cv->cd(2);
  bottomPad_cv->SetPad(0, 0.0, 1, 0.3);
  bottomPad_cv->SetTopMargin(0.05);
  bottomPad_cv->SetBottomMargin(0.3);
  bottomPad_cv->SetGridy();
  bottomPad_cv->SetLogx();

  // TH1D* ratio_cv = (TH1D*)h1_cv->Clone("ratio_h1_over_h0");
  // ratio_cv->Divide(h0_cv);  // h1 / h0
  TH1D* ratio_cv = MakeRatioWithStatError(h1_cv, h0_cv);
  // ratio_cv->SetMinimum(0.9);
  // ratio_cv->SetMaximum(1.1);
  ratio_cv->SetLineColor(kRed + 1);
  ratio_cv->SetLineWidth(2);
  ratio_cv->SetTitle("");
  ratio_cv->GetYaxis()->SetTitle("h_{1} / h_{0}");
  ratio_cv->GetYaxis()->SetTitleSize(0.09);
  ratio_cv->GetYaxis()->SetTitleOffset(0.45);
  ratio_cv->GetYaxis()->SetLabelSize(0.08);
  ratio_cv->GetYaxis()->SetNdivisions(505);
  ratio_cv->GetXaxis()->SetTitle("E_{#nu} [GeV]");
  ratio_cv->GetXaxis()->SetTitleSize(0.09);
  ratio_cv->GetXaxis()->SetLabelSize(0.08);
  ratio_cv->SetStats(0);
  ratio_cv->GetXaxis()->SetRangeUser(0.5, 10.0);
  ratio_cv->Draw("E1");
  bottomPad_cv->Modified();
  bottomPad_cv->Update();

  //Save both
  c_cv->SaveAs("LE1_nubar_cvweighted.png");
  c_cv->Close();




  // // Styling
  // h0->SetLineColor(kBlack); h0->SetLineWidth(2);
  // h1->SetLineColor(kRed);   h1->SetLineWidth(2);
  // h2->SetLineColor(kBlue);  h2->SetLineWidth(2);

  // h0->SetTitle("Unweighted fluxes for #nu_{#mu} and #nu_{e};E_{#nu} [GeV];Flux");
  // // h0->SetTitle("Fluxes for #bar{#nu_{#mu}} and #bar{#nu_{e}};E_{#nu} [GeV];Flux");

  // TCanvas* c = new TCanvas("c", "Flux Comparison + Ratio", 800, 800);
  // c->Divide(1, 2, 0.01, 0.01);  // vertical split
  // // c->SetLogx();

  // // === Upper Pad for Fluxes ===
  // TPad* topPad = (TPad*)c->cd(1);
  // topPad->SetPad(0, 0.3, 1, 1);
  // topPad->SetBottomMargin(0.02);
  // topPad->SetLogx();
  // topPad->SetLogy();
  // h0->GetXaxis()->SetRangeUser(0.3, 10.0);
  // h1->GetXaxis()->SetRangeUser(0.3, 10.0);
  // h2->GetXaxis()->SetRangeUser(0.3, 10.0);
  // h0->Draw("HIST");
  // h1->Draw("HIST SAME");
  // h2->Draw("HIST SAME");
  // h1->Draw("E SAME");
  // h0->GetXaxis()->SetRangeUser(0.5, 10.0);
  // gPad->Modified();
  // gPad->Update();

  // TLegend* leg = new TLegend(0.12, 0.1, 0.45, 0.3);
  // leg->AddEntry(h0, "existing #nu_{#mu} flux", "l");
  // leg->AddEntry(h1, "new #nu_{#mu} flux", "l");
  // leg->AddEntry(h2, "new #nu_{e} flux", "l");
  // // leg->AddEntry(h0, "existing #bar{#nu_{#mu}} flux", "l");
  // // leg->AddEntry(h1, "new #bar{#nu_{#mu}} flux", "l");
  // // leg->AddEntry(h2, "new #bar{#nu_{e}} flux", "l");
  // leg->Draw();

  // // === Bottom Pad for Ratio ===
  // TPad* bottomPad = (TPad*)c->cd(2);
  // bottomPad->SetPad(0, 0.0, 1, 0.3);
  // bottomPad->SetTopMargin(0.05);
  // bottomPad->SetBottomMargin(0.3);
  // bottomPad->SetGridy();
  // bottomPad->SetLogx();

  // TH1D* ratio = MakeRatioWithStatError(h1, h0);
  // // ratio->SetMinimum(0.9);
  // // ratio->SetMaximum(1.1);
  // ratio->SetLineColor(kRed + 1);
  // ratio->SetLineWidth(2);
  // ratio->SetTitle("");
  // ratio->GetYaxis()->SetTitle("h_{1} / h_{0}");
  // ratio->GetYaxis()->SetTitleSize(0.09);
  // ratio->GetYaxis()->SetTitleOffset(0.45);
  // ratio->GetYaxis()->SetLabelSize(0.08);
  // ratio->GetYaxis()->SetNdivisions(505);
  // ratio->GetXaxis()->SetTitle("E_{#nu} [GeV]");
  // ratio->GetXaxis()->SetTitleSize(0.09);
  // ratio->GetXaxis()->SetLabelSize(0.08);
  // ratio->SetStats(0);
  // ratio->GetXaxis()->SetRangeUser(0.3, 100.0);
  // ratio->GetXaxis()->SetRangeUser(0.5, 10.0);
  // ratio->Draw("E1");
  // bottomPad->Modified();
  // bottomPad->Update();

  // //Save both
  // c->SaveAs("LE1_unweighted.png");
  // c->Close();


  // f_numu0->Close();
  // f_numu->Close();
  // f_nue->Close();
}
