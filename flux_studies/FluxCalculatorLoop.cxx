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



/*
 * P8 CHANGES:
 *
 * beamUniverseOffset:
 *   Cyclically shifts only the BeamFocus universe index. This is used
 *   to decorrelate the BeamFocus uncertainties between FHC and RHC
 *   while leaving the PPFX universe ordering unchanged.
 *
 * useFlatBeamFocus:
 *   Selects externally supplied flat BeamFocus weights instead of the
 *   event-dependent tuple BeamFocus weights.
 *
 * flatBeamFocusWeights:
 *   Contains one relative BeamFocus weight per universe for species
 *   whose LE tuples do not contain the appropriate focusing variation.
 */
void FluxCalculatorLoop::EventLoop(
  TChain* chain,
  const TEventList* evtList,
  PlotUtils::MnvH1D* histogram,
  std::string branchName,
  double additionalWeight,
  bool cvWeighted,
  unsigned int beamUniverseOffset,
  bool useFlatBeamFocus,
  const std::vector<double>& flatBeamFocusWeights
)
{
  double cvVal;
  double cvWeight;
  int nuParentID;

  double* combinedWeights = nullptr;


  /*
  * P8 stores mc_wgt_Flux_BeamFocus and mc_wgt_ppfx1_Total as
  * scaled integer arrays rather than double arrays. Separate integer
  * branch buffers are therefore required before decoding the weights.
  */
  std::map<std::string, int*> universeWeightsInt;
  std::map<std::string, double*> universeWeights;
  std::map<std::string, double> cvWeightContributions;

  std::vector<std::string> vertErrorBandNames = histogram->GetVertErrorBandNames();

  for (std::vector<std::string>::const_iterator it_band = vertErrorBandNames.begin();
                                                it_band != vertErrorBandNames.end();
                                                ++it_band) {

    universeWeights[*it_band] = new double[1000];

    const std::string wgtName = "mc_wgt_" + *it_band;

    // P8: These two branches are stored as int[1000] with a scale of 1e-7.
    const bool isIntScaled = (*it_band == "Flux_BeamFocus" || *it_band == "ppfx1_Total");

    if (isIntScaled) {
      universeWeightsInt[*it_band] = new int[1000];
      chain->SetBranchAddress(wgtName.c_str(), universeWeightsInt[*it_band]);
    }
    else {
      chain->SetBranchAddress(wgtName.c_str(), universeWeights[*it_band]);
    }


    /*
    * The P8 Flux_BeamFocus and ppfx1_Total universe values are
    * absolute component weights.
    * Read the matching component CV weights so that the universes can
    * later be converted to relative multiplicative weights.
    *
    * Flux_BeamFocus corresponds to mc_hornCurrent_cvweight.
    * ppfx1_Total corresponds to mc_ppfx1_cvweight.
    */
    std::string name = *it_band;

    if (*it_band == "Flux_BeamFocus") { name = "hornCurrent"; }
    else if (*it_band == "ppfx1_Total") { name = "ppfx1"; }

    const std::string cvWeightContributor = "mc_" + name + "_cvweight";

    if (chain->GetBranch(cvWeightContributor.c_str()) != nullptr) {
      cvWeightContributions[*it_band] = 1.0;
      chain->SetBranchAddress(cvWeightContributor.c_str(), &(cvWeightContributions[*it_band]));
    }
  }

  const unsigned int nCombinedUniv = histogram->GetVertErrorBand(histogram->GetVertErrorBandNames().front())->GetNHists();

  /*
  * The P8 tuple branches used here contain 1000 universe entries.
  * Validate the histogram universe count before indexing these
  * fixed-size tuple branch buffers.
  */
  if (nCombinedUniv == 0) {
    std::cerr << "ERROR: BeamFocus error band contains zero universes." << std::endl;

    for (auto& it : universeWeightsInt) { delete[] it.second; }
    for (auto& it : universeWeights) { delete[] it.second; }

    delete[] combinedWeights;
    return;
  }

  if (nCombinedUniv > 1000) {
    std::cerr
      << "ERROR: Requested "
      << nCombinedUniv
      << " universes, but tuple buffers contain only 1000."
      << std::endl;

    for (auto& it : universeWeightsInt) { delete[] it.second; }
    for (auto& it : universeWeights) { delete[] it.second; }

    delete[] combinedWeights;
    return;
  }

  /*
  * Normalize the requested cyclic BeamFocus offset to the available
  * universe count. PPFX is intentionally not shifted.
  */
  const unsigned int normalizedBeamOffset = beamUniverseOffset % nCombinedUniv;
  combinedWeights = new double[nCombinedUniv];


  /*
  * Flat BeamFocus mode requires exactly one externally supplied
  * relative weight for every histogram universe.
  */
  if (useFlatBeamFocus && flatBeamFocusWeights.size() != nCombinedUniv) {
    std::cerr << "ERROR: Flat BeamFocus vector contains "
              << flatBeamFocusWeights.size()
              << " weights, but the histogram contains "
              << nCombinedUniv << " universes." << std::endl;

    for (auto& it : universeWeightsInt) { delete[] it.second; }

    for (auto& it : universeWeights) { delete[] it.second; }

    delete[] combinedWeights;
    return;
  }

  std::cout << "Using BeamFocus universe offset: " << normalizedBeamOffset
            << " out of " << nCombinedUniv << " universes" << std::endl;

  std::cout << "BeamFocus source: "
            << (
                useFlatBeamFocus
                  ? "flat species-dependent weights"
                  : "event-dependent tuple weights"
              )
            << std::endl;


  /*
  * Add the combined Flux error band. Each universe will later contain
  *
  *   PPFX_relative[u] * BeamFocus_relative[u].
  */
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
    * P8 CHANGE:
    * Flux_BeamFocus and ppfx1_Total are stored as integer arrays.
    *
    * The stored value is the physical absolute component weight
    * multiplied by 1e7. Decode it by multiplying by 1e-7.
    *
    * This step only restores the numeric value; the decoded weights
    * must still be divided by their component CV weights to become
    * relative universe weights.
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

    /*
    * Skip events with invalid or effectively zero total CV weight.
    * These events are not safe to use in the flux histograms.
    */
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
      * Tuple PPFX and BeamFocus weights are stored as absolute
      * component weights and are divided by their matching CV
      * contribution.
      *
      * Flat BeamFocus CSV values are already relative weights and
      * are used directly.
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


    /*
    * Protect the absolute-to-relative conversion against missing,
    * invalid, or effectively zero component CV weights. A fallback of
    * 1 preserves the decoded universe value rather than producing NaN
    * or an unbounded ratio.
    */
    if (std::isnan(ppfxCVContr) || std::abs(ppfxCVContr) < 1.0e-12) {
      ppfxCVContr = 1.0;
    }

    if (std::isnan(beamCVContr) || std::abs(beamCVContr) < 1.0e-12) {
      beamCVContr = 1.0;
    }

    /*
    * Return the relative BeamFocus weight for one output universe.
    *
    * The BeamFocus source universe is shifted cyclically:
    *
    *   source = (output + offset) % N.
    *
    * In tuple mode, used for the right-sign species, the BeamFocus
    * weight comes from the event tuple. The tuple value is an absolute
    * component weight, so divide it by mc_hornCurrent_cvweight.
    *
    * In flat mode, used for the other species, use the corresponding
    * species-dependent CSV weight. The CSV values are already relative
    * and must not be divided again.
    */
    const auto getBeamFocusRelative =
      [&](unsigned int outputUniverse) -> double
      {
        const unsigned int sourceUniverse =
          (outputUniverse + normalizedBeamOffset) % nCombinedUniv;

        if (useFlatBeamFocus) {
          // CSV values are already relative BeamFocus weights.
          return flatBeamFocusWeights[sourceUniverse];
        }

        // P8 tuple BeamFocus values are absolute component weights.
        // Convert to relative by dividing by the BeamFocus CV component.
        return
          universeWeights["Flux_BeamFocus"][sourceUniverse] /
          beamCVContr;
      };

    /*
    * Fill the standalone PPFX and BeamFocus error bands using relative
    * universe weights.
    *
    * ppfx1_Total:
    *   Keep universe u unchanged and divide the absolute tuple weight by
    *   mc_ppfx1_cvweight.
    *
    * Flux_BeamFocus:
    *   Apply the cyclic BeamFocus offset and obtain the relative weight
    *   from getBeamFocusRelative(). This may use either the shifted tuple
    *   weight or the shifted flat species-dependent CSV weight.
    */
    for (std::vector<std::string>::const_iterator it_band = vertErrorBandNames.begin();
                                                  it_band != vertErrorBandNames.end();
                                                  ++it_band) {

      if (*it_band != "Flux_BeamFocus" && *it_band != "ppfx1_Total") {
        continue;
      }

      double* reorderedWeights = new double[nCombinedUniv];

      for (unsigned int univ_idx = 0; univ_idx < nCombinedUniv; ++univ_idx) {

        if (*it_band == "Flux_BeamFocus") {
          reorderedWeights[univ_idx] = getBeamFocusRelative(univ_idx);
        }
        else {
          reorderedWeights[univ_idx] = universeWeights["ppfx1_Total"][univ_idx] / ppfxCVContr;
        }
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
    * Build the combined Flux error band:
    *
    *   Flux[u] = PPFX_relative[u] * BeamFocus_relative[u].
    *
    * PPFX remains at universe u. BeamFocus uses the shifted source
    * universe and may come either from the tuple or the flat CSV.
    */
    for (unsigned int univ_idx = 0; univ_idx < nCombinedUniv; ++univ_idx) {

      const double ppfxRelative = universeWeights["ppfx1_Total"][univ_idx] / ppfxCVContr;

      const double beamRelative = getBeamFocusRelative(univ_idx);

      combinedWeights[univ_idx] = ppfxRelative * beamRelative;
    }

    static int debugCombinedPrints = 0;

    /*
    * Print a few representative combined-universe calculations to
    * verify the PPFX index, shifted BeamFocus index, selected BeamFocus
    * source, and final product during validation runs.
    */
    if (debugCombinedPrints < 5) {
      const unsigned int outIdx = 0;
      const unsigned int ppfxIdx = outIdx;
      const unsigned int beamIdx = (outIdx + normalizedBeamOffset) % nCombinedUniv;

      const double ppfxForDebug = universeWeights["ppfx1_Total"][ppfxIdx] / ppfxCVContr;

      const double beamForDebug = getBeamFocusRelative(outIdx);

      std::cout << "\nDEBUG combined Flux"
                << " source="
                << (useFlatBeamFocus ? "flat" : "tuple")
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
            << " entries with invalid or effectively zero total CV weight"
            << std::endl;

  for (auto& it : universeWeightsInt) { delete[] it.second; }
  for (auto& it : universeWeights) { delete[] it.second; }

  delete[] combinedWeights;
}