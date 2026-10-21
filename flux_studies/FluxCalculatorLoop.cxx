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

namespace
{
  const int kBeamFocusTuple    = 0;
  const int kBeamFocusFlat     = 1;
  const int kBeamFocusResponse = 2;
}

/*
 * P8 CHANGES:
 *
 * beamUniverseOffset:
 *   Cyclically shifts only the BeamFocus universe index while leaving
 *   the PPFX universe ordering unchanged.
 *
 * beamFocusMode:
 *   Selects the BeamFocus source:
 *
 *     0 = native tuple BeamFocus
 *     1 = externally supplied flat BeamFocus weights
 *     2 = externally supplied energy-dependent BeamFocus response
 *
 * flatBeamFocusWeights:
 *   One relative BeamFocus weight per universe for LE wrong-sign
 *   species.
 *
 * beamFocusEnergyEdges / beamFocusResponseWeights:
 *   Energy-dependent relative BeamFocus weights for the LE right-sign
 *   electron species. The flattened response is stored as
 *
 *     response[energyBin * nUniverses + universe].
 */

void FluxCalculatorLoop::EventLoop(
  TChain* chain,
  const TEventList* evtList,
  PlotUtils::MnvH1D* histogram,
  std::string branchName,
  double additionalWeight,
  bool cvWeighted,
  unsigned int beamUniverseOffset,
  int beamFocusMode,
  const std::vector<double>& flatBeamFocusWeights,
  const std::vector<double>& beamFocusEnergyEdges,
  const std::vector<double>& beamFocusResponseWeights
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
   * Validate BeamFocus mode.
   */
  if (
    beamFocusMode != kBeamFocusTuple &&
    beamFocusMode != kBeamFocusFlat &&
    beamFocusMode != kBeamFocusResponse
  ) {
    std::cerr
      << "ERROR: Unknown BeamFocus mode: "
      << beamFocusMode
      << std::endl;

    for (auto& it : universeWeightsInt) { delete[] it.second; }
    for (auto& it : universeWeights) { delete[] it.second; }

    delete[] combinedWeights;
    return;
  }


  /*
   * Flat BeamFocus mode requires exactly one relative weight
   * per universe.
   */
  if (
    beamFocusMode == kBeamFocusFlat &&
    flatBeamFocusWeights.size() != nCombinedUniv
  ) {
    std::cerr
      << "ERROR: Flat BeamFocus vector contains "
      << flatBeamFocusWeights.size()
      << " weights, but the histogram contains "
      << nCombinedUniv
      << " universes."
      << std::endl;

    for (auto& it : universeWeightsInt) { delete[] it.second; }
    for (auto& it : universeWeights) { delete[] it.second; }

    delete[] combinedWeights;
    return;
  }


  /*
   * Energy-dependent BeamFocus response mode requires:
   *
   *   Nenergy bins = Nedges - 1
   *
   * and one response weight for every
   *
   *   energy bin x universe.
   */
  if (beamFocusMode == kBeamFocusResponse) {

    if (beamFocusEnergyEdges.size() < 2) {
      std::cerr
        << "ERROR: BeamFocus response requires at least "
        << "two energy edges."
        << std::endl;

      for (auto& it : universeWeightsInt) { delete[] it.second; }
      for (auto& it : universeWeights) { delete[] it.second; }

      delete[] combinedWeights;
      return;
    }

    const unsigned int nBeamEnergyBins =
      beamFocusEnergyEdges.size() - 1;

    const unsigned long long expectedResponseSize =
      static_cast<unsigned long long>(nBeamEnergyBins) *
      static_cast<unsigned long long>(nCombinedUniv);

    if (
      beamFocusResponseWeights.size() != expectedResponseSize
    ) {
      std::cerr
        << "ERROR: BeamFocus response contains "
        << beamFocusResponseWeights.size()
        << " weights, but expected "
        << expectedResponseSize
        << " = "
        << nBeamEnergyBins
        << " energy bins x "
        << nCombinedUniv
        << " universes."
        << std::endl;

      for (auto& it : universeWeightsInt) { delete[] it.second; }
      for (auto& it : universeWeights) { delete[] it.second; }

      delete[] combinedWeights;
      return;
    }

    /*
     * Require strictly increasing energy edges.
     */
    for (
      unsigned int i = 1;
      i < beamFocusEnergyEdges.size();
      ++i
    ) {
      if (
        beamFocusEnergyEdges[i] <=
        beamFocusEnergyEdges[i - 1]
      ) {
        std::cerr
          << "ERROR: BeamFocus response energy edges "
          << "are not strictly increasing."
          << std::endl;

        for (auto& it : universeWeightsInt) { delete[] it.second; }
        for (auto& it : universeWeights) { delete[] it.second; }

        delete[] combinedWeights;
        return;
      }
    }
  }


  std::cout
    << "Using BeamFocus universe offset: "
    << normalizedBeamOffset
    << " out of "
    << nCombinedUniv
    << " universes"
    << std::endl;


  std::string beamFocusSourceName;

  if (beamFocusMode == kBeamFocusTuple) {
    beamFocusSourceName = "event-dependent tuple weights";
  }
  else if (beamFocusMode == kBeamFocusFlat) {
    beamFocusSourceName = "flat species-dependent weights";
  }
  else {
    beamFocusSourceName = "energy-dependent response weights";
  }

  std::cout
    << "BeamFocus source: "
    << beamFocusSourceName
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

  /*
   * Return the zero-based BeamFocus response bin for a neutrino
   * energy in GeV.
   *
   * Bins follow:
   *
   *   [edge_i, edge_{i+1})
   *
   * with the final upper edge included in the last bin.
   *
   * Return -1 outside the response range.
   */
  const auto getBeamFocusEnergyBin =
    [&](double energyGeV) -> int
    {
      if (beamFocusEnergyEdges.size() < 2) {
        return -1;
      }

      if (
        energyGeV < beamFocusEnergyEdges.front() ||
        energyGeV > beamFocusEnergyEdges.back()
      ) {
        return -1;
      }

      if (energyGeV == beamFocusEnergyEdges.back()) {
        return static_cast<int>(
          beamFocusEnergyEdges.size() - 2
        );
      }

      for (
        unsigned int ibin = 0;
        ibin + 1 < beamFocusEnergyEdges.size();
        ++ibin
      ) {
        if (
          energyGeV >= beamFocusEnergyEdges[ibin] &&
          energyGeV <  beamFocusEnergyEdges[ibin + 1]
        ) {
          return static_cast<int>(ibin);
        }
      }

      return -1;
    };

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
      * Flat BeamFocus CSV and regenerated response values are already
      * relative weights and are used directly.
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
     * BeamFocus universe ordering is shifted cyclically:
     *
     *   source = (output + offset) % N.
     *
     * Three possible BeamFocus sources are supported:
     *
     *   tuple:
     *     Native P8 BeamFocus weight. This is an absolute component
     *     weight, so divide by mc_hornCurrent_cvweight.
     *
     *   flat:
     *     Externally supplied relative weight for LE wrong-sign
     *     species. Use directly.
     *
     *   response:
     *     Externally supplied energy-dependent relative weight for
     *     LE right-sign electron species. Use directly.
     */
    const auto getBeamFocusRelative =
      [&](unsigned int outputUniverse, double energyGeV) -> double
      {
        const unsigned int sourceUniverse =
          (outputUniverse + normalizedBeamOffset) %
          nCombinedUniv;


        /*
         * Flat LE wrong-sign BeamFocus.
         */
        if (beamFocusMode == kBeamFocusFlat) {
          return flatBeamFocusWeights[sourceUniverse];
        }


        /*
         * Energy-dependent LE right-sign electron BeamFocus.
         */
        if (beamFocusMode == kBeamFocusResponse) {

          const int energyBin =
            getBeamFocusEnergyBin(energyGeV);

          /*
           * The regenerated response covers 0--100 GeV.
           * If an event nevertheless falls outside that range,
           * leave its BeamFocus weight unchanged.
           */
          if (energyBin < 0) {
            return 1.0;
          }

          const unsigned long long responseIndex =
            static_cast<unsigned long long>(energyBin) *
            static_cast<unsigned long long>(nCombinedUniv) +
            sourceUniverse;

          return beamFocusResponseWeights[responseIndex];
        }


        /*
         * Native tuple BeamFocus.
         *
         * P8 tuple BeamFocus values are absolute component weights.
         * Convert them to relative weights using the corresponding
         * BeamFocus CV component.
         */
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
    *   from getBeamFocusRelative(). This may use the shifted tuple,
    *   flat species-dependent CSV, or energy-dependent response weight.
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
          reorderedWeights[univ_idx] = getBeamFocusRelative(univ_idx, cv);
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
    * universe and may come from the tuple, flat CSV, or regenerated
    * energy-dependent response.
    */
    for (unsigned int univ_idx = 0; univ_idx < nCombinedUniv; ++univ_idx) {

      const double ppfxRelative = universeWeights["ppfx1_Total"][univ_idx] / ppfxCVContr;

      const double beamRelative = getBeamFocusRelative(univ_idx, cv);

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

      const double beamForDebug = getBeamFocusRelative(outIdx, cv);

      std::cout << "\nDEBUG combined Flux"
                << " source=" << beamFocusSourceName
                << " E=" << cv
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