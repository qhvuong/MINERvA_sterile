#ifndef UnfoldCovmxAsErrorBand_h
#define UnfoldCovmxAsErrorBand_h

#include <set>
#include <vector>

#include <TMatrixD.h>

#include <PlotUtils/MnvH1D.h>

TMatrixD RemoveZeroRows(const TMatrixD& covmx, std::set<int>& zero_index_vec);

std::vector<TMatrixD> GenerateCorrelatedErrors(const TMatrixD& covmx, int n_universe);

// convert the covariance to bin-to-bin correlated error band named 'Unfold_error', which
// has 500 universes
PlotUtils::MnvH1D* UnfoldCovmxAsErrorBand(const PlotUtils::MnvH1D* mnvh1d_unfolded, const TMatrixD& covmx);

#endif
