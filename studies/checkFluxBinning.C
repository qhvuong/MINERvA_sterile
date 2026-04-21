void checkFluxBinning() {
  TFile* f1 = TFile::Open("/exp/minerva/app/users/qvuong/MAT_AL9/CC-NuE-XSec/custom_plotutils/data/flux/flux-gen2thin-pdg12-minerva1.root");
  TFile* f2 = TFile::Open("/exp/minerva/app/users/qvuong/MAT_AL9/CC-NuE-XSec/custom_plotutils/data/flux/flux-gen2thin-pdg12-minervame1D.root");

  auto h1 = (TH1*)f1->Get("flux_E_cvweighted");
  auto h2 = (TH1*)f2->Get("flux_E_cvweighted");

  bool same = true;

  if (!h1 || !h2) {
    std::cout << "Could not get flux_E_cvweighted from one of the files." << std::endl;
    return;
  }

  if (h1->GetNbinsX() != h2->GetNbinsX()) {
    same = false;
    std::cout << "Different nbins: " << h1->GetNbinsX()
              << " vs " << h2->GetNbinsX() << std::endl;
  } else {
    for (int i = 1; i <= h1->GetNbinsX() + 1; ++i) {
      double e1 = h1->GetXaxis()->GetBinLowEdge(i);
      double e2 = h2->GetXaxis()->GetBinLowEdge(i);
      if (fabs(e1 - e2) > 1e-9) {
        same = false;
        std::cout << "Mismatch at edge " << i << ": "
                  << e1 << " vs " << e2 << std::endl;
        break;
      }
    }
  }

  std::cout << "Same binning? " << (same ? "YES" : "NO") << std::endl;

  if (same) {
    std::cout << "Bin edges: ";
    for (int i = 1; i <= h1->GetNbinsX() + 1; ++i) {
      std::cout << h1->GetXaxis()->GetBinLowEdge(i);
      if (i < h1->GetNbinsX() + 1) std::cout << ", ";
    }
    std::cout << std::endl;
  }
}