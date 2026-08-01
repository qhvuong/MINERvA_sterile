#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <cmath>

#include "TChain.h"
#include "TBranch.h"
#include "TEventList.h"

#include "PlotUtils/MnvH1D.h"
#include "FluxCalculatorLoop.h"

void FluxCalculatorLoop::EventLoop(
  TChain* chain,
  const TEventList* evtList,
  PlotUtils::MnvH1D* histogram,
  std::string branchName,
  double additionalWeight,
  bool cvWeighted,
  unsigned int beamUniverseOffset
)
{
  double cvVal;
  double cvWeight;
  int nuParentID;

  double* combinedWeights = new double[1000];

  std::map<std::string, int*> universeWeightsInt;
  std::map<std::string, double*> universeWeights;
  std::map<std::string, double> cvWeightContributions;

  std::vector<std::string> vertErrorBandNames = histogram->GetVertErrorBandNames();

  for (std::vector<std::string>::const_iterator it_band = vertErrorBandNames.begin();
                                                it_band != vertErrorBandNames.end();
                                                ++it_band) {

    universeWeights[*it_band] = new double[1000];

    const std::string wgtName = "mc_wgt_" + *it_band;
    const bool isIntScaled = (*it_band == "Flux_BeamFocus" || *it_band == "ppfx1_Total");

    if (isIntScaled) {
      universeWeightsInt[*it_band] = new int[1000];
      chain->SetBranchAddress(wgtName.c_str(), universeWeightsInt[*it_band]);
    }
    else {
      chain->SetBranchAddress(wgtName.c_str(), universeWeights[*it_band]);
    }

    std::string name = *it_band;

    if (*it_band == "Flux_BeamFocus") {
      name = "hornCurrent";
    }
    else if (*it_band == "ppfx1_Total") {
      name = "ppfx1";
    }

    const std::string cvWeightContributor = "mc_" + name + "_cvweight";

    if (chain->GetBranch(cvWeightContributor.c_str()) != nullptr) {
      cvWeightContributions[*it_band] = 1.0;
      chain->SetBranchAddress(cvWeightContributor.c_str(), &(cvWeightContributions[*it_band]));
    }
  }

  const unsigned int nCombinedUniv =
    histogram->GetVertErrorBand(histogram->GetVertErrorBandNames().front())->GetNHists();

  if (nCombinedUniv == 0) {
    std::cerr << "ERROR: BeamFocus error band contains zero universes." << std::endl;

    for (auto& it : universeWeightsInt) { delete[] it.second; }
    for (auto& it : universeWeights) { delete[] it.second; }

    delete[] combinedWeights;
    return;
  }

  const unsigned int normalizedBeamOffset = beamUniverseOffset % nCombinedUniv;

  std::cout << "Using BeamFocus universe offset: " << normalizedBeamOffset
            << " out of " << nCombinedUniv << " universes" << std::endl;

  histogram->AddVertErrorBand("Flux", nCombinedUniv);

  if (universeWeights.find("ppfx1_Total") == universeWeights.end() ||
      universeWeights.find("Flux_BeamFocus") == universeWeights.end()) {

    std::cerr << "ERROR: Cannot construct Flux band without both ppfx1_Total and Flux_BeamFocus."
              << std::endl;

    for (auto& it : universeWeightsInt) { delete[] it.second; }
    for (auto& it : universeWeights) { delete[] it.second; }

    delete[] combinedWeights;
    return;
  }

  int nanCount = 0;

  chain->SetBranchAddress(branchName.c_str(), &cvVal);
  chain->SetBranchAddress("mc_cvweight_total", &cvWeight);
  chain->SetBranchAddress("mc_fr_nuParentID", &nuParentID);

  std::map<int, int> parent_freq;

  const unsigned long long N = evtList->GetN();
  const unsigned long long Ndiv10 = static_cast<unsigned long long>(N / 10.0) + 1;

  std::cout << " (using " << N << " events)" << std::endl;

  if (N < 10) { std::cout << "  0/" << N << std::flush; }
  else        { std::cout << "  0%" << std::flush; }

  for (unsigned long long i = 0; i < N; ++i) {
    if (N < 10) { std::cout << "  " << i + 1 << "/" << N << std::flush; }
    else if ((i + 1) % Ndiv10 == 0) {
      std::cout << "  "
                << static_cast<long>(
                     static_cast<double>(i + 1) /
                     static_cast<double>(N) *
                     100.0
                   )
                << "%"
                << std::flush;
    }

    const unsigned long long entrynum = evtList->GetEntry(i);
    chain->GetEntry(entrynum);

    /*
     * Flux_BeamFocus and ppfx1_Total are stored as scaled integers.
     * Convert them back to doubles.
     */
    for (const auto& band : vertErrorBandNames) {
      if (band == "Flux_BeamFocus" || band == "ppfx1_Total") {
        for (unsigned int u = 0; u < nCombinedUniv; ++u) {
          universeWeights[band][u] =
            static_cast<double>(universeWeightsInt[band][u]) * 1.0e-7;
        }
      }
    }

    parent_freq[nuParentID]++;

    if (std::isnan(cvWeight) || cvWeight < 1.0e-6) {
      nanCount++;
      continue;
    }

    const double cv = cvVal * additionalWeight;

    if (cvWeighted) {
      histogram->Fill(cv, cvWeight);
    }
    else {
      histogram->Fill(cv);
    }

    /*
     * Fill standalone flux error bands.
     *
     * PPFX remains unshifted.
     *
     * BeamFocus_saved[u] =
     *   BeamFocus_input[(u + normalizedBeamOffset) % N]
     *
     * P8 stores absolute weights, so standalone universe weights
     * are divided by the corresponding component CV weight.
     */
    for (std::vector<std::string>::const_iterator it_band = vertErrorBandNames.begin();
                                                  it_band != vertErrorBandNames.end();
                                                  ++it_band) {

      if (*it_band != "Flux_BeamFocus" && *it_band != "ppfx1_Total") {
        continue;
      }

      double cvContr = 1.0;

      const auto cvContrIt = cvWeightContributions.find(*it_band);

      if (cvContrIt != cvWeightContributions.end()) {
        cvContr = cvContrIt->second;
      }

      if (std::isnan(cvContr) || std::abs(cvContr) < 1.0e-12) {
        cvContr = 1.0;
      }

      double* reorderedWeights = new double[nCombinedUniv];

      for (unsigned int univ_idx = 0; univ_idx < nCombinedUniv; ++univ_idx) {
        unsigned int sourceIdx = univ_idx;

        if (*it_band == "Flux_BeamFocus") {
          sourceIdx = (univ_idx + normalizedBeamOffset) % nCombinedUniv;
        }

        reorderedWeights[univ_idx] =
          universeWeights[*it_band][sourceIdx] / cvContr;
      }

      if (cvWeighted) {
        histogram->FillVertErrorBand(
          *it_band,
          cv,
          reorderedWeights,
          cvWeight,
          1.0
        );
      }
      else {
        histogram->FillVertErrorBand(
          *it_band,
          cv,
          reorderedWeights
        );
      }

      delete[] reorderedWeights;
    }

    /*
     * Build combined Flux:
     *
     * Flux[u] =
     *   (PPFX[u] / PPFX_CV)
     *   *
     *   (BeamFocus[(u + offset) % N] / BeamFocus_CV)
     */
    double ppfxCVContr = 1.0;
    double beamCVContr = 1.0;

    const auto ppfxCVIt = cvWeightContributions.find("ppfx1_Total");

    if (ppfxCVIt != cvWeightContributions.end()) {
      ppfxCVContr = ppfxCVIt->second;
    }

    const auto beamCVIt = cvWeightContributions.find("Flux_BeamFocus");

    if (beamCVIt != cvWeightContributions.end()) {
      beamCVContr = beamCVIt->second;
    }

    if (std::isnan(ppfxCVContr) || std::abs(ppfxCVContr) < 1.0e-12) {
      ppfxCVContr = 1.0;
    }

    if (std::isnan(beamCVContr) || std::abs(beamCVContr) < 1.0e-12) {
      beamCVContr = 1.0;
    }

    for (unsigned int univ_idx = 0; univ_idx < nCombinedUniv; ++univ_idx) {
      const unsigned int beamIdx =
        (univ_idx + normalizedBeamOffset) % nCombinedUniv;

      const double ppfxRelative =
        universeWeights["ppfx1_Total"][univ_idx] / ppfxCVContr;

      const double beamRelative =
        universeWeights["Flux_BeamFocus"][beamIdx] / beamCVContr;

      combinedWeights[univ_idx] = ppfxRelative * beamRelative;
    }

    static int debugCombinedPrints = 0;

    if (debugCombinedPrints < 5) {
      const unsigned int outIdx = 0;
      const unsigned int ppfxIdx = outIdx;
      const unsigned int beamIdx =
        (outIdx + normalizedBeamOffset) % nCombinedUniv;

      const double ppfxForDebug =
        universeWeights["ppfx1_Total"][ppfxIdx] / ppfxCVContr;

      const double beamForDebug =
        universeWeights["Flux_BeamFocus"][beamIdx] / beamCVContr;

      std::cout << "\nDEBUG combined Flux"
                << " outIdx=" << outIdx
                << " ppfxIdx=" << ppfxIdx
                << " beamIdx=" << beamIdx
                << " ppfxWeight=" << ppfxForDebug
                << " beamWeight=" << beamForDebug
                << " combined=" << combinedWeights[outIdx]
                << std::endl;

      debugCombinedPrints++;
    }

    if (cvWeighted) {
      histogram->FillVertErrorBand(
        "Flux",
        cv,
        combinedWeights,
        cvWeight,
        1.0
      );
    }
    else {
      histogram->FillVertErrorBand(
        "Flux",
        cv,
        combinedWeights
      );
    }
  }

  std::cout << "\nParent PDG Frequency Summary:\n";

  for (const auto& kv : parent_freq) {
    std::cout << "  PDG " << kv.first << " : " << kv.second << " events\n";
  }

  if (N < 10) {
    std::cout << "  " << N << "/" << N << std::endl;
  }
  else {
    std::cout << "  100%" << std::endl;
  }

  std::cout << "Skipped " << nanCount
            << " entries where the CV weight was NaN"
            << std::endl;

  for (auto& it : universeWeightsInt) { delete[] it.second; }
  for (auto& it : universeWeights) { delete[] it.second; }

  delete[] combinedWeights;
}