#include "TFile.h"
#include "TH1.h"
#include "TString.h"
#include "TSystem.h"
#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvVertErrorBand.h"

#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

// Convert existing per-playlist ME fluxes to the selection naming convention.
// No merging, rebinning, normalization, or universe reordering is performed.
void flux_copy_ME(
    const char* input_dir = "/exp/minerva/data/users/qvuong/flux_studies/producedFluxes_p8_ME_offsetPPFX",
    const char* output_dir = "") {
  TH1::AddDirectory(kFALSE);
  const std::string destination = output_dir[0] ? output_dir : input_dir;
  const std::vector<std::string> playlists = {"1D", "1M", "1N", "5A", "6A"};
  struct Species { int pdg; const char* name; };
  const std::vector<Species> species = {
      {12, "nue"}, {-12, "nuebar"}, {14, "numu"}, {-14, "numubar"}};
  const std::vector<std::string> hist_names = {
      "flux_E_unweighted", "flux_E_cvweighted"};

  // Check every input before creating or replacing any output.
  for (const auto& playlist : playlists) {
    for (const auto& sp : species) {
      const TString path = Form("%s/ME%s_%s.root", input_dir, playlist.c_str(), sp.name);
      std::unique_ptr<TFile> input(TFile::Open(path, "READ"));
      if (!input || input->IsZombie() || input->TestBit(TFile::kRecovered))
        throw std::runtime_error(std::string("Invalid or recovered input: ") + path.Data());
      for (const auto& name : hist_names) {
        PlotUtils::MnvH1D* hist = nullptr;
        input->GetObject(name.c_str(), hist);
        if (!hist || !hist->HasVertErrorBand("Flux") ||
            hist->GetVertErrorBand("Flux")->GetNHists() != 1000)
          throw std::runtime_error(std::string(path.Data()) + ": " + name +
                                   " missing or not carrying 1000 Flux universes");
      }
      std::cout << "[INPUT OK] " << path << std::endl;
    }
  }
  if (gSystem->mkdir(destination.c_str(), true) != 0 &&
      gSystem->AccessPathName(destination.c_str()))
    throw std::runtime_error("Cannot create output directory: " + destination);

  int written = 0;
  for (const auto& playlist : playlists) {
    for (const auto& sp : species) {
      const TString path = Form("%s/ME%s_%s.root", input_dir, playlist.c_str(), sp.name);
      std::unique_ptr<TFile> input(TFile::Open(path, "READ"));
      if (!input || input->IsZombie())
        throw std::runtime_error(std::string("Cannot reopen input: ") + path.Data());
      for (unsigned int i = 0; i < hist_names.size(); ++i) {
        PlotUtils::MnvH1D* source = nullptr;
        input->GetObject(hist_names[i].c_str(), source);
        if (!source) throw std::runtime_error("Input histogram disappeared");
        std::unique_ptr<PlotUtils::MnvH1D> hist(
            static_cast<PlotUtils::MnvH1D*>(source->Clone(hist_names[i].c_str())));
        hist->SetDirectory(nullptr);
        const TString output_path = Form("%s/flux-%s-pdg%d-minervame%s.root",
            destination.c_str(), i == 0 ? "g4numiv5" : "gen2thin", sp.pdg, playlist.c_str());
        {
          TFile output(output_path, "RECREATE");
          if (output.IsZombie())
            throw std::runtime_error(std::string("Cannot create: ") + output_path.Data());
          output.cd();
          if (hist->Write(hist_names[i].c_str()) <= 0)
            throw std::runtime_error(std::string("Cannot write: ") + output_path.Data());
          output.Close();
        }
        std::unique_ptr<TFile> check(TFile::Open(output_path, "READ"));
        PlotUtils::MnvH1D* saved = nullptr;
        if (check && !check->IsZombie()) check->GetObject(hist_names[i].c_str(), saved);
        if (!saved || !saved->HasVertErrorBand("Flux") ||
            saved->GetVertErrorBand("Flux")->GetNHists() != 1000)
          throw std::runtime_error(std::string("Output verification failed: ") + output_path.Data());
        std::cout << "[OUTPUT OK] " << output_path << " (1000 Flux universes)" << std::endl;
        ++written;
      }
    }
  }
  std::cout << "Done: " << written << " files written for 5 ME playlists and 4 species."
            << std::endl;
}
