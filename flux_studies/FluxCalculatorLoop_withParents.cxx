#include <iostream>
#include <vector>
#include <string>
#include "TChain.h"
#include "TBranch.h"
#include "TEventList.h"
#include "PlotUtils/MnvH1D.h"
#include "FluxCalculatorLoop.h"
void FluxCalculatorLoop::EventLoop(TChain * chain, const TEventList * evtList, std::vector<PlotUtils::MnvH1D*>& parentHistos, std::string branchName, double additionalWeight, bool cvWeighted)
{
  double cvVal, cvWeight;
  
  // double  mc_lebrun_cv_wgt;
  // double  mc_horntilt_cv_wgt;
  // double  mc_beamspot_cv_wgt;   

  int nuParentID;   
  
  double * combinedWeights = new double[1000];
  double combinedCVContr;
  
  std::map<std::string, double*> universeWeights;
  std::map<std::string, double> cvWeightContributions;
  std::vector<std::string> vertErrorBandNames = parentHistos.front()->GetVertErrorBandNames();
  for (std::vector<std::string>::const_iterator it_band = vertErrorBandNames.begin();
       it_band != vertErrorBandNames.end(); ++it_band) {
    
  //  std::cout << "band names" << "  " << *it_band;
    universeWeights[*it_band] = new double[1000];
    std::string wgtName = "mc_wgt_" + *it_band;
    chain->SetBranchAddress( wgtName.c_str(), universeWeights[*it_band] );
    
    std::string name = *it_band;
    if( *it_band == "Flux_BeamFocus" ) name = "hornCurrent";                
    else if( *it_band == "ppfx1_Total" ) name = "ppfx1";
    std::string cvWeightContributor = "mc_" + name + "_cvweight";
    if (chain->GetBranch(cvWeightContributor.c_str()) != NULL)
      {
        cvWeightContributions[*it_band] = 1;
        chain->SetBranchAddress( cvWeightContributor.c_str(), &(cvWeightContributions[*it_band]) );
      }
  }
  // Use the first histogram in the map to get the number of universes
  unsigned int nCombinedUniv = parentHistos.front()->GetVertErrorBand(vertErrorBandNames.front())->GetNHists();

  // Add "Flux" error band to all histograms
  for (auto* hist : parentHistos) {
      hist->AddVertErrorBand("Flux", nCombinedUniv);
  }
  
  int nanCount = 0;
  
  chain->SetBranchAddress(branchName.c_str(), &cvVal);
  chain->SetBranchAddress("mc_cvweight_total", &cvWeight);
  
  chain->SetBranchAddress("mc_fr_nuParentID", &nuParentID);
  // chain->SetBranchAddress("mc_lebrun_cv_wgt", &mc_lebrun_cv_wgt);
  // chain->SetBranchAddress("mc_horntilt_cv_wgt", &mc_horntilt_cv_wgt);
  // chain->SetBranchAddress("mc_beamspot_cv_wgt", &mc_beamspot_cv_wgt);        
  
  std::map<int, int> parent_freq;
  const int OTHER_PDG = 999999; 

  unsigned long long N = evtList->GetN();
  unsigned long long Ndiv10 = long(N/10.) + 1;
  
  std::cout << " (using " << N << " events)" << std::endl;
  if (N < 10)
    std::cout << "  0/" << N << std::flush;
  else
                std::cout << "  0%" << std::flush;


  std::vector<int> pdg_order = {130, 311, 13, -13, 321, -321, 211, -211};
  const int OTHER_INDEX = pdg_order.size();  // Last index in vector is for "other"

  for (unsigned long long i = 0; i < N; ++i)
  {
    if (N < 10)
      std::cout << "  " << i + 1 << "/" << N << std::flush;
    else if ((i + 1) % Ndiv10 == 0)
      std::cout << "  " << long(float(i + 1) / N * 100) << "%" << std::flush;

    unsigned long long entrynum = evtList->GetEntry(i);
    chain->GetEntry(entrynum);

    if (cvWeight != cvWeight || cvWeight < 1.0E-6) {
      nanCount++;
      continue;
    }

    int pdg = nuParentID;
    parent_freq[pdg]++;

    double cv = cvVal * additionalWeight;

    // Find the correct index for the histogram
    int index = OTHER_INDEX;
    for (size_t j = 0; j < pdg_order.size(); ++j) {
      if (pdg == pdg_order[j]) {
        index = j;
        break;
      }
    }

    PlotUtils::MnvH1D* targetHist = parentHistos[index];

    if (cvWeighted)
      targetHist->Fill(cv, cvWeight);
    else
      targetHist->Fill(cv);

    // Reset combinedWeights
    for (unsigned int univ_idx = 0; univ_idx < nCombinedUniv; ++univ_idx)
      combinedWeights[univ_idx] = 1.0;

    for (const auto& band : vertErrorBandNames)
    {
      if (band == "Flux") continue;

      if (band == "Flux_BeamFocus") {
        if (cvWeighted)
          targetHist->FillVertErrorBand(band, cv, universeWeights[band], cvWeight, 1.0);
        else
          targetHist->FillVertErrorBand(band, cv, universeWeights[band]);
      }
      else if (band == "ppfx1_Total") {
        if (cvWeighted)
          targetHist->FillVertErrorBand(band, cv, universeWeights[band], 1.0, 1.0);
        else
          targetHist->FillVertErrorBand(band, cv, universeWeights[band]);
      }

      for (unsigned int univ_idx = 0; univ_idx < nCombinedUniv; ++univ_idx)
        combinedWeights[univ_idx] *= universeWeights[band][univ_idx];
    }

    // Fill combined "Flux" band
    if (cvWeighted)
      targetHist->FillVertErrorBand("Flux", cv, combinedWeights, 1.0, 1.0);
    else
      targetHist->FillVertErrorBand("Flux", cv, combinedWeights);
  }
      
  std::cout << "\nParent PDG Frequency Summary:\n";
  for (const auto& kv : parent_freq) {
    std::cout << "  PDG " << kv.first << " : " << kv.second << " events\n";
  }

  if (N < 10)
    std::cout << "  " << N << "/" << N << std::endl;
  else
    std::cout << "  100%" << std::endl;
  
  std::cout << "Skipped " << nanCount << " entries where the CV weight was NaN" << std::endl;
  
  // clean up...
  for (auto& it_univ : universeWeights)
  {
    if (it_univ.second)
      delete[] it_univ.second;
  }
  if (combinedWeights)
    delete[] combinedWeights;

}