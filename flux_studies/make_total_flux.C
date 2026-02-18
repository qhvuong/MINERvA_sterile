  // const char* infile = "flux_me5Ap7_numubarRHC.root";
  // const char* outfile = "flux_numubarRHC.root";

#include "TFile.h"
#include "TKey.h"
#include "TString.h"
#include "TSystem.h"
#include <iostream>

// IMPORTANT: include PlotUtils header AFTER dictionaries are available.
// In ROOT macros, this is usually fine as long as the library is loaded.
#include "PlotUtils/MnvH1D.h"

void make_total_flux(const char* infile = "flux_me5Ap7_numubarRHC.root",
                     const char* outfile = "flux_total.root")
{
  // 1) Make sure PlotUtils is loaded (needed for PlotUtils::MnvH1D dictionary)
  // Try common names; if your env already loads it, this is harmless.
  gSystem->Load("libMAT.so");
  gSystem->Load("libPlotUtils.so");

  TFile* fin = TFile::Open(infile, "READ");
  if (!fin || fin->IsZombie()) {
    std::cerr << "ERROR: cannot open input file: " << infile << "\n";
    return;
  }

  PlotUtils::MnvH1D* h_unw_total = nullptr;
  PlotUtils::MnvH1D* h_cv_total  = nullptr;

  // Loop keys in file
  TIter nextkey(fin->GetListOfKeys());
  TKey* key = nullptr;

  while ((key = (TKey*)nextkey())) {
    TString name = key->GetName();

    // --- sum UNWEIGHTED parents ---
    if (name.BeginsWith("flux_E_unweighted_parent")) {
      auto h = dynamic_cast<PlotUtils::MnvH1D*>(key->ReadObj());
      if (!h) continue;

      if (!h_unw_total) {
        h_unw_total = (PlotUtils::MnvH1D*)h->Clone("flux_E_unweighted_total");
        h_unw_total->Reset();
      }
      h_unw_total->Add(h);
    }

    // --- sum CVWEIGHTED parents ---
    if (name.BeginsWith("flux_E_cvweighted_parent")) {
      auto h = dynamic_cast<PlotUtils::MnvH1D*>(key->ReadObj());
      if (!h) continue;

      if (!h_cv_total) {
        h_cv_total = (PlotUtils::MnvH1D*)h->Clone("flux_E_cvweighted_total");
        h_cv_total->Reset();
      }
      h_cv_total->Add(h);
    }
  }

  if (!h_unw_total || !h_cv_total) {
    std::cerr << "ERROR: did not find both unweighted and cvweighted parent flux histos.\n";
    fin->Close();
    return;
  }

  // keep only total_POT
  auto totalPOT = (TParameter<double>*) fin->Get("total_POT");
  if (!totalPOT) {
    std::cerr << "WARNING: total_POT not found in input.\n";
  }

  // Write a CLEAN output file (only the two totals + total_POT)
  TFile* fout = TFile::Open(outfile, "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "ERROR: cannot open output file: " << outfile << "\n";
    fin->Close();
    return;
  }

  h_unw_total->Write();
  h_cv_total->Write();
  if (totalPOT) totalPOT->Write();

  fout->Write();
  fout->Close();
  fin->Close();

  std::cout << "Wrote: " << outfile << "\n";
}
