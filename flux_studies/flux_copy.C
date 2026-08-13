{
  TH1::AddDirectory(kFALSE);

  const char* dir = "/exp/minerva/data/users/qvuong/flux_studies/producedFluxes_p8_LE";

  std::vector<int> playlists = {1, 5, 13};
  // std::vector<int> playlists = {5};

  struct Species {
    int pdg;
    const char* name;
  };

  std::vector<Species> species_list = {
    {12,  "nue"},
    {-12, "nuebar"},
    {14,  "numu"},
    {-14, "numubar"}
  };

  for (int playlist : playlists) {
    for (const auto& sp : species_list) {

      TString input = Form("%s/LE%d_%s.root", dir, playlist, sp.name);

      TString output_unweighted = Form(
        "%s/flux-g4numiv5-pdg%d-minerva%d.root",
        dir, sp.pdg, playlist
      );

      TString output_cvweighted = Form(
        "%s/flux-gen2thin-pdg%d-minerva%d.root",
        dir, sp.pdg, playlist
      );

      std::cout << "\n==================================================\n";
      std::cout << "Input:  " << input << "\n";
      std::cout << "PDG:    " << sp.pdg << "\n";
      std::cout << "Output unweighted: " << output_unweighted << "\n";
      std::cout << "Output cvweighted: " << output_cvweighted << "\n";
      std::cout << "==================================================\n";

      TFile* fin = TFile::Open(input, "READ");
      if (!fin || fin->IsZombie()) {
        std::cout << "[ERROR] Could not open input file: " << input << "\n";
        if (fin) fin->Close();
        continue;
      }

      PlotUtils::MnvH1D* h_unw = nullptr;
      PlotUtils::MnvH1D* h_cvw = nullptr;

      fin->GetObject("flux_E_unweighted", h_unw);
      fin->GetObject("flux_E_cvweighted", h_cvw);

      if (!h_unw || !h_cvw) {
        std::cout << "[ERROR] Missing input histograms in " << input << "\n";
        std::cout << "  h_unw = " << h_unw << "\n";
        std::cout << "  h_cvw = " << h_cvw << "\n";
        fin->Close();
        continue;
      }

      gROOT->cd();

      PlotUtils::MnvH1D* h1 = (PlotUtils::MnvH1D*)h_unw->Clone("flux_E_unweighted");
      PlotUtils::MnvH1D* h2 = (PlotUtils::MnvH1D*)h_cvw->Clone("flux_E_cvweighted");

      auto band1 = h1->GetVertErrorBand("Flux");
      auto band2 = h2->GetVertErrorBand("Flux");

      std::cout << "  unweighted Flux universes: "
                << (band1 ? band1->GetNHists() : -1) << "\n";
      std::cout << "  cvweighted Flux universes: "
                << (band2 ? band2->GetNHists() : -1) << "\n";

      TFile* f1 = TFile::Open(output_unweighted, "RECREATE");
      if (!f1 || f1->IsZombie()) {
        std::cout << "[ERROR] Could not create output file: " << output_unweighted << "\n";
        delete h1;
        delete h2;
        fin->Close();
        if (f1) f1->Close();
        continue;
      }

      f1->cd();
      h1->Write();
      f1->Close();

      TFile* f2 = TFile::Open(output_cvweighted, "RECREATE");
      if (!f2 || f2->IsZombie()) {
        std::cout << "[ERROR] Could not create output file: " << output_cvweighted << "\n";
        delete h1;
        delete h2;
        fin->Close();
        if (f2) f2->Close();
        continue;
      }

      f2->cd();
      h2->Write();
      f2->Close();

      std::cout << "[OK] Wrote both output files.\n";

      fin->Close();

      delete h1;
      delete h2;
    }
  }

  std::cout << "\nDone.\n";
}