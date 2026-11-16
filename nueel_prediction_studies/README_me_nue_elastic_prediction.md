ME neutrino–electron prediction
This is the ME equivalent of the existing LE production in this directory.
make_me_nue_elastic_prediction.py imports the cross-section functions directly
from make_nue_elastic_from_flux_TomalakXS_fixedTotalNorm_mnv.py: the couplings,
radiative-correction CSVs, and existing temporary total-cross-section factors
are the same. It uses the new ME flux files without modifying their universes.

The data and measurement covariance come from Table IV of Zazueta et al.,
https://arxiv.org/abs/2209.05540 (Phys. Rev. D 107, 012001).
Both the FHC and RHC six-bin spectra are included. The last published bin is
9 GeV to infinity; prediction integration ends at the supplied flux maximum.
The numerical histogram edge is 100 GeV.
Inputs and calculation
For each of nue, nuebar, numu, and numubar, the calculation is
events = data_POT * target_electrons * integral(flux * differential_cross_section)
Flux is a density per GeV per square metre per POT. Cross sections are in
square centimetres, so the code includes the factor 1e-4. The input histogram
defaults to flux_E_cvweighted. Only its combined Flux band is propagated;
its separate PPFX and BeamFocus bands are not added a second time.
Corresponding universe indices are summed across species and playlists.

Outputs retain the LE histogram names h_nue_elastic_total and
h_nue_elastic_{nue,nuebar,numu,numubar}, each with 1000 Flux universes.
Bin errors on these prediction objects are zero: input MC statistical errors
are not yet propagated. The flux uncertainty is in the universe band.

Published data objects are paper_data_FHC and paper_data_RHC.
Their error bars are total measurement errors, not statistical errors alone.
The covariance objects are:
Object	Contents
cov_nonflux_ME_FHC	Published FHC 6x6 measurement block
cov_nonflux_ME_RHC	Published RHC 6x6 measurement block
cov_nonflux_ME_FHC_RHC	Joint 12x12 block, retaining FHC/RHC correlations


These are aliases of the paper_measurement_covariance_* objects. They supply
the data/non-flux covariance term; do not add a second Poisson data variance.
They are the published measurement errors, including uncertainties from its
background subtraction and efficiency correction. They are not a prediction
flux covariance. The IMD measurement is omitted by taking the upper-left
12x12 block, without conditioning on IMD.
Run on GPVM
Place this script and me_nue_elastic_paper_inputs.json alongside the existing
LE script. Use the existing ROOT/PlotUtils environment. ROOT is loaded before
NumPy/SciPy.

First run a single-source diagnostic; these commands do not claim to reproduce
the full published playlist mixture:
cd "$CCNUEROOT/nueel_prediction_studies"

OPENBLAS_NUM_THREADS=1 python make_me_nue_elastic_prediction.py \
  --beam FHC --playlist ME1D --pot 1.16e21 \
  --check-quadrature --plots \
  --output me_nue_elastic_ME1D_diagnostic.root

OPENBLAS_NUM_THREADS=1 python make_me_nue_elastic_prediction.py \
  --beam RHC --playlist ME5A --pot 1.12e21 \
  --check-quadrature --plots \
  --output me_nue_elastic_ME5A_diagnostic.root
The FHC exposure and 1.99e30 target electrons are stated in the earlier FHC
publication, https://arxiv.org/abs/1906.00111. The RHC command explicitly uses
the 1.12e21 POT quoted in Luis's conclusion. Section II instead quotes 1.2e21:
confirm the exact exposure before precision use. RHC uses the same default
electron count based on its stated 5.99-ton fiducial mass; confirm that
normalization too. There is no fitted normalization or normalization to the
578 selected RHC events: that count is before the efficiency correction.

For the full mixture, repeat --playlist-pot PLAYLIST=DATA_POT for ME1D,
ME1M, and ME1N in FHC, or ME5A and ME6A in RHC. Each POT must be the data
exposure represented by that distinct MC source. The total is their sum.
Do not use MC POT fractions or give the full exposure to every source.
Per-POT fluxes cannot simply be added. The exact data-to-source exposure
mapping still needs to be supplied.

The default flux directory is
/exp/minerva/data/users/qvuong/flux_studies/producedFluxes_p8_ME_offsetPPFX.
The Tomalak CSV directory is the existing LE default,
/exp/minerva/data/users/qvuong/nueel_prediction_studies/tomalak_inputs.
Both can be overridden. Use fluxes that have not already been constrained
by these same elastic-scattering data.
Electron response
As with the base LE flux-integration script, the default prediction is in
true total electron energy. Luis's measured spectrum is efficiency-corrected
reconstructed energy, not an unfolded truth spectrum. The default plot is
therefore a diagnostic. No LE detector response is silently applied to ME.

An optional --response-json PATH accepts true_edges_GeV,
reco_edges_GeV, and matrix_reco_by_true. The matrix has reconstructed bins
as rows and true bins as columns, with nonnegative probabilities and column
sums at most one. It describes migration, not selection efficiency, because
the published data are already efficiency-corrected. Reconstructed edges must
match the six published bins. Include fine true bins below 0.8 GeV if their
migration into the measured range matters.

With this argument the script additionally writes
h_nue_elastic_reco_total and the four reconstructed flavor components,
including the response-folded Flux universes. This is not yet an oscillation
template builder or a change to stitch.
Verification and outputs
python make_me_nue_elastic_prediction.py --self-test runs without ROOT.
It checks an analytic integration example, flux-bin subdivision, exposure
linearity, universe propagation, published totals, and covariance positivity.
The physical integration kernel was also checked for all four flavors against
an independent analytic electron-energy integral followed by adaptive
neutrino-energy integration, with radiative corrections disabled: maximum
relative difference 1.36e-9 for the checked bins. This checks the integration,
not the external cross-section model or production ROOT files.

--check-quadrature doubles both integration orders and records/prints the
changes using the actual fluxes and radiative corrections. ROOT I/O and the
full prediction need execution on GPVM; ROOT is unavailable in the development
environment. Each run writes ROOT histograms, metadata JSON, a truth-prediction
CSV, the joint published covariance CSV, and optionally a comparison PNG.