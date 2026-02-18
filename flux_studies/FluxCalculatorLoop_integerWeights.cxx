#include <iostream>
#include <vector>
#include <string>
#include "TChain.h"
#include "TBranch.h"
#include "TEventList.h"
#include "PlotUtils/MnvH1D.h"
#include "FluxCalculatorLoop.h"
void FluxCalculatorLoop::EventLoop(TChain * chain, const TEventList * evtList, PlotUtils::MnvH1D * histogram, std::string branchName, double additionalWeight, bool cvWeighted)
{
  double cvVal, cvWeight;
  
  // double  mc_lebrun_cv_wgt;
  // double  mc_horntilt_cv_wgt;
  // double  mc_beamspot_cv_wgt;   

  int nuParentID;   
  
  double * combinedWeights = new double[1000];
  double combinedCVContr;
  
  std::map<std::string, int*> universeWeightsInt;
  std::map<std::string, double*> universeWeights;  // converted, scaled
  std::map<std::string, double> cvWeightContributions;
  std::vector<std::string> vertErrorBandNames = histogram->GetVertErrorBandNames();
  for (std::vector<std::string>::const_iterator it_band = vertErrorBandNames.begin();
       it_band != vertErrorBandNames.end(); ++it_band) {
    
  //  std::cout << "band names" << "  " << *it_band;
  universeWeights[*it_band] = new double[1000];

  std::string wgtName = "mc_wgt_" + *it_band;
  const bool isIntScaled = (*it_band == "Flux_BeamFocus" || *it_band == "ppfx1_Total");

  if (isIntScaled) {
    universeWeightsInt[*it_band] = new int[1000];
    chain->SetBranchAddress(wgtName.c_str(), universeWeightsInt[*it_band]);
  } else {
    chain->SetBranchAddress(wgtName.c_str(), universeWeights[*it_band]); // read doubles directly
  }

    
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
  // also add a 'combined' histogram, assuming they all have the same number of universes...
  unsigned int nCombinedUniv = histogram->GetVertErrorBand(histogram->GetVertErrorBandNames().front())->GetNHists();
  histogram->AddVertErrorBand("Flux", nCombinedUniv);
  
  int nanCount = 0;
  
  chain->SetBranchAddress(branchName.c_str(), &cvVal);
  chain->SetBranchAddress("mc_cvweight_total", &cvWeight);
  
  chain->SetBranchAddress("mc_fr_nuParentID", &nuParentID);
  // chain->SetBranchAddress("mc_lebrun_cv_wgt", &mc_lebrun_cv_wgt);
  // chain->SetBranchAddress("mc_horntilt_cv_wgt", &mc_horntilt_cv_wgt);
  // chain->SetBranchAddress("mc_beamspot_cv_wgt", &mc_beamspot_cv_wgt);        
  
  std::map<int, int> parent_freq;

  unsigned long long N = evtList->GetN();
  unsigned long long Ndiv10 = long(N/10.) + 1;
  
  std::cout << " (using " << N << " events)" << std::endl;
  if (N < 10)
    std::cout << "  0/" << N << std::flush;
  else
                std::cout << "  0%" << std::flush;
  for (unsigned long long i = 0; i < N; ++i)
    {
      if (N < 10)
        std::cout << "  " << i+1 << "/" << N << std::flush;
      else if ((i+1) % Ndiv10 == 0)
        std::cout << "  " << long(float(i+1)/N*100) << "%" << std::flush;
      
      unsigned long long entrynum = evtList->GetEntry(i);
      chain->GetEntry(entrynum);

      for (const auto& band : vertErrorBandNames) {
        if (band == "Flux") continue;

        if (band == "Flux_BeamFocus" || band == "ppfx1_Total") {
          for (unsigned int u = 0; u < nCombinedUniv; ++u) {
            universeWeights[band][u] = static_cast<double>(universeWeightsInt[band][u]) * 1.0e-7;
          }
        }
        // else: already read into universeWeights[band] as doubles, no scaling needed
      }


      int pdg = nuParentID;
      parent_freq[pdg]++;
      // std::cout << "Entry " << i << ": E = " << cvVal << ", PDG = " << pdg << std::endl;
      
      if( cvWeight != cvWeight || cvWeight < 1.0E-6 ) {
        nanCount++;
        continue;
      }
      
      double cv = cvVal * additionalWeight;
 //     std::cout << "central value           : " << cv << std::endl;
 //     std::cout << "central value weight    : " << cvWeight << std::endl;
      if (cvWeighted)
        histogram->Fill(cv, cvWeight);
      else
        histogram->Fill(cv);
      
      combinedCVContr = 1.0;
      for (unsigned int univ_idx = 0; univ_idx < nCombinedUniv; univ_idx++)
        combinedWeights[univ_idx] = 1.0;
            
      for (std::vector<std::string>::const_iterator it_band = vertErrorBandNames.begin();
           it_band != vertErrorBandNames.end(); ++it_band) {
        if ((*it_band) == "Flux")
          continue;
        
        double cvWeightFromMe = 1.0;
        //double abc = 1;
        if( *it_band == "Flux_BeamFocus" ) {
          if (cvWeighted)
            histogram->FillVertErrorBand( *it_band, cv, universeWeights[*it_band], cvWeight, 1.0);
          else
            histogram->FillVertErrorBand( *it_band, cv, universeWeights[*it_band] );
        } else if( *it_band == "ppfx1_Total" ) { 
          //abc = mc_lebrun_cv_wgt*mc_horntilt_cv_wgt*mc_beamspot_cv_wgt;
          if (cvWeighted)
            //histogram->FillVertErrorBand( *it_band, cv, universeWeights[*it_band], abc , 1.0);
            histogram->FillVertErrorBand( *it_band, cv, universeWeights[*it_band], 1.0, 1.0);
          else
            histogram->FillVertErrorBand( *it_band, cv, universeWeights[*it_band] );
          
        }
//        int jj = 1;
//        std::cout << "ErrBand: " << *it_band <<" UniId: " <<jj <<" UniWight: "<<  universeWeights[*it_band][jj] << std::endl;                  
        for (unsigned int univ_idx = 0; univ_idx < nCombinedUniv; univ_idx++)
          combinedWeights[univ_idx] *= universeWeights[*it_band][univ_idx];
      }
      
//      std::cout << "universe weight for flux: " << combinedWeights[0] << std::endl;
//      std::cout << "cvweightFromMe for flux : " << combinedCVContr << std::endl;
//      std::cout << "cvweight for flux       : " << cvWeight << std::endl;
//      std::cout << "-------------------------------------------" << std::endl;
      if (cvWeighted)
        histogram->FillVertErrorBand( "Flux", cv, combinedWeights, 1.0, 1.0);
      else
        histogram->FillVertErrorBand( "Flux", cv, combinedWeights );
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
  for (auto& it : universeWeightsInt) {
    delete [] it.second;
  }
  for (auto& it : universeWeights) {
    delete [] it.second;
  }
  delete [] combinedWeights;
}