Read diagnostics in this order:
1. checks.csv: FAIL is a failed invariant; WARN needs interpretation. Exit code is nonzero for FAIL.
2. CASE/input_samples.csv and template_coverage.csv: confirm actual fitted samples and components.
3. CASE/oscillation_bins.csv: actual vs independent oscillation prediction, probabilities and L/E averages.
   Both reconstruction paths use the same input templates; agreement does NOT prove input units/weights are physical.
4. CASE/controlled_profile_tests.csv: actual profiling vs independent whitened least squares.
   Injection recovery is compared to the RIDGE optimum, not the injected coefficients.
5. CASE/chi2_breakdown.csv: exact additive sample contributions q_i = r_i*(V^-1 r)_i.
   Sample contributions can be negative. Null-minus columns sum with the penalty to total improvement.
6. CASE/ue_conditional_slice.csv: fix dm2, Umu, Utau and vary Ue, re-solving flux each time.
   This is a conditional slice, NOT a minimization over the other oscillation parameters.
7. flux_solution_comparisons.csv: within-run null/BF, cross-run identical points, and own best fits.
   Compare induced shifts as well as coefficients. Redundant universe directions can make coefficients misleading.
   shared_basis_agrees checks common sample columns, not external flux-universe identity/provenance.
8. common_sample_shift_comparisons.csv and coefficient/spectra PNGs.
9. CASE/data_side_flux_response.csv: quantify the data-side Flux response; no correction is applied.

A is fixed and lambda/Nprof are unchanged. The script forces relative-A and hybrid-covariance
switches OFF and records previous values. --exclude ratio removes ratio bins from the flux
SOLVE ONLY, following fitData; ratios remain in the final residual chi2.
The three default runs use the same common diagnostic points. Approximate candidates are
never claimed as measured best fits. Supply exact best_fit values in config, or use --refit.
If all inputs and predictions pass, a Ue=0 optimum can still be physical: muon-only disappearance
changes CCnumu, ratios and nu+e while direct CCnue stays unchanged; the combined profiled
objective can prefer this even when the inclusive CCnumu integral is nearly correct.
Plot error bars show only covariance diagonals. CSV chi2 uses the full correlated covariance.
No production ROOT files or fitter source are modified. No HIST_CONFIG.json is overwritten.
