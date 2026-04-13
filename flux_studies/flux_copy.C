{
  TH1::AddDirectory(kFALSE);

  int playlist = 13;
  int pdg = 12;

  TFile* fin = TFile::Open(Form("flux_le%d_nueFHC.root", playlist), "READ");
  if (!fin || fin->IsZombie()) {
    std::cout << "Could not open input file\n";
    return;
  }

  PlotUtils::MnvH1D* h_unw = nullptr;
  PlotUtils::MnvH1D* h_cvw = nullptr;
  fin->GetObject("flux_E_unweighted", h_unw);
  fin->GetObject("flux_E_cvweighted", h_cvw);

  if (!h_unw || !h_cvw) {
    std::cout << "Missing input histograms.\n";
    fin->Close();
    return;
  }

  gROOT->cd();
  PlotUtils::MnvH1D* h1 = (PlotUtils::MnvH1D*)h_unw->Clone("flux_E_unweighted");
  PlotUtils::MnvH1D* h2 = (PlotUtils::MnvH1D*)h_cvw->Clone("flux_E_cvweighted");

  TFile* f1 = TFile::Open(Form("flux-g4numiv5-pdg%d-minerva%d.root", pdg, playlist), "RECREATE");
  if (!f1 || f1->IsZombie()) {
    std::cout << "Could not create g4numi output file\n";
    delete h1;
    delete h2;
    fin->Close();
    return;
  }
  f1->cd();
  h1->Write();
  f1->Close();

  TFile* f2 = TFile::Open(Form("flux-gen2thin-pdg%d-minerva%d.root", pdg, playlist), "RECREATE");
  if (!f2 || f2->IsZombie()) {
    std::cout << "Could not create gen2thin output file\n";
    delete h1;
    delete h2;
    fin->Close();
    return;
  }
  f2->cd();
  h2->Write();
  f2->Close();

  fin->Close();

  delete h1;
  delete h2;
}