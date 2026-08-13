#ifndef PYTHONMINERVAUNIVERSE_H
#define PYTHONMINERVAUNIVERSE_H
#include "PlotUtils/MinervaUniverse.h"
#include "PlotUtils/NSFDefaults.h"
#include "PlotUtils/PlotUtilsPhysicalConstants.h"
#include "TVector3.h"
#include "PlotUtils/FluxSystematics.cxx" // PlotUtils::flux_reweighter
#include "PlotUtils/GenieSystematics.cxx" //PlotUtils::IsNonResPi
#include "PlotUtils/MinosMuonEfficiencyCorrection.h" // PlotUtils::MinosMiuonEfficiencyCorrection
#include "PlotUtils/MnvNormalization.h" // PlotUtils::MnvNormalizer
#include "PlotUtils/weight_fsi.h" // PlotUtils::weight_fsi
#include "PlotUtils/weightCoherentPi.h" // PlotUtils::weight_coherent


class PythonMinervaUniverse : public PlotUtils::MinervaUniverse {
public:
  PythonMinervaUniverse(PlotUtils::TreeWrapper* chw=NULL, double nsigma =0):
    PlotUtils::MinervaUniverse(chw,nsigma) {};
#include "PlotUtils/WeightFunctions.h"
#include "PlotUtils/TruthFunctions.h"
#include "PlotUtils/LowRecoilFunctions.h"
#include "PlotUtils/MuonFunctions.h"
  virtual double GetCalRecoilEnergy() const {
  	return 0;
  }
  std::vector<int> m_non_cal_idx;
  ROOT::Math::XYZTVector GetVertex() const {
    ROOT::Math::XYZTVector result;
    result.SetCoordinates(GetVec<double>("vtx").data());
    return result;
  }
  int IsMinosMatchMuon() const {
    return GetInt("isMinosMatchTrack") == 1;
  }

  double GetMuonQP() const {
    double qp = GetDouble("MasterAnaDev_minos_trk_qp");
    return qp==-1? 0:qp;
  }
};

//PythonMinervaUniverse::Q0Spline PythonMinervaUniverse::current_spline;
#endif /* PYTHONMINERVAUNIVERSE_H */
