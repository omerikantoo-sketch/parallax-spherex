# Research log (2026-11-14 / 15)

* **Data access.** IRSA SIA returned 502 and the IBE directories 404 throughout. We switched to the public AWS mirror, using HTTP range reads of only the needed HDUs.
  - FLAGS is a RICE-compressed bintable, decoded per tile row.
  - The ePSF table has 441 zones of 33×33 at ×5 oversampling; it sums to 1 on the oversampled grid, so the native response is 25 × the sampled value.
  - WCS-WAVE is a 9×9 control grid, indexed as `vals[iy, ix, (wave, band)]` in 1-based pixels.
* **Coverage without SIA.** The first prior (ecliptic geometry) left 57 candidate periods. A solar-elongation prior (boresight 91–116°) plus sampled headers and exact TAN-SIP containment cut this down. Barnard's Star: 238 exposures in 1122 s. Off-detector `all_world2pix` raised NoConvergence; this was fixed with a TAN pre-test.
* **Cadence measured from real headers.** Bursts of 2–4 exposures 0.0015 d apart, each shifting the source ~115 px (about one spectral channel); detector pairs (1,4), (2,5), (3,6) see the same sky. The simulator was changed to reproduce this.
* **E1 v1.** The bank gave only +0.17 mag overall, and brown-dwarf m50 was undefined on the coarse grid. Rerun as v2 (finer grid, same nulls). Result: brown dwarfs +0.65 mag, solar-like ≈ 0, ordering consistent with the analytic R.
* **Branch-and-bound.** Exact but no pruning: the bound is loose over a 20-day baseline. Recorded as a negative result.
* **Gravity nesting bug.** H1 beat H0 on null data by Δχ² ≈ 39 because H0 had less spectral freedom. Giving the static/blend term a 2K basis in every hypothesis makes the models nested; null Δχ² is now ~10–15.
* **E3.** Joint chromatic vs grey vs single band:
  - the chromatic model gives *no* detection gain;
  - it improves t_E accuracy (0.15 vs 0.42 median |ln| error, deep cadence);
  - with BIC alone, symmetric variables are often labelled as lenses (30%; see the guard entry below).

  This narrows the Gravity claims to characterisation and blend decomposition.
* **Feasibility.** The deep-field yield is ≲0.02 per year. Large bulge numbers are upper bounds, limited by confusion and cadence. Conclusion: SPHEREx is a spectral follow-up instrument for known events, not a discovery engine for microlensing.
* **E2 Barnard's Star (real).** 29 of 85 exposures gave usable PSF centroids; the rest have saturated or flagged cores. The fit gives μα* = −1700 ± 922 and μδ = 11952 ± 853 mas/yr, against Gaia DR3 −802 and 10362 mas/yr. Motion is detected at 14σ in declination, and the result is consistent with Gaia within 1.0σ and 1.9σ. Fixed a bug on the way: simultaneous detectors share a timestamp, which mis-indexed the spacecraft positions.
* **E2b real search.** The first run was dominated by static-star residuals (z ≈ 1058 at Barnard's position).
  - Root cause: `chi2_to_z` saturated at z ≈ 37, so every bright source had the same z. The static catalogue then picked plateau-edge positions instead of the stars, and subtraction missed them.
  - Fix: Wilson–Hilferty above p < 1e-250. Bright sources are also masked, with radius set from the fitted flux, and a minimum velocity of 1 px per visit is imposed.
  - After the fix, the top statistic is 232 against a time-scrambled null maximum of 233 on the same frames, and no candidate exceeds it.
  - Near a 9-mag star the search is limited by residuals (z ≈ 230), not by noise. E1 thresholds (z ≈ 5) do not apply to such fields.
* **WCS-WAVE grid size bug (caught by a real-data unit test).** D3's spectral control grid is 9×14, but the parser assumed 9×9 for every detector. As a result, 33 of 86 cached D3 cutouts carried bandwidth-like values (~0.055 µm) as their wavelength.
  - The parser now reads TFORM/TDIM.
  - `scripts/repair_wavelengths.py` recomputed the wave/band planes of all cached cutouts; no pixel data changed.
  - Derived field caches were rebuilt, and E2 was rerun after the fix.
  - Synthetic benchmarks E1/E3/E4 do not use archive wavelengths and are unaffected.
* **Gravity on real photometry produced a false "lens" (caught in browser testing).** A bright, saturated field star in the Barnard field gave "Microlensing-compatible" with t_E = 22 470 d, a negative source fraction and χ²/dof ≈ 900. BIC alone rewarded the lens model's extra spectral freedom. Guards were added to `GravityModel.compare`:
  1. model adequacy (best reduced χ² ≤ 3, otherwise "Model inadequate, no inference");
  2. physical validity (0 < f_s ≤ 1.05, t_E ≤ 3× span, t₀ near the data);
  3. a 1% relative calibration floor for real photometry.

  After the fix, the bright star gives reduced χ² 419 and a blank or faint position gives 8; both are refused. The real-photometry noise model (archive variance plus forced photometry) underestimates the scatter, and a smooth 6-term SED cannot follow real M-star molecular bands at S/N ≈ 200. E3 was rerun so that its classification uses the same guarded decision. Deep cadence: symmetric variables labelled lens fell from 30/100 (BIC alone) to 0/100; true lenses labelled lens fell from 28 to 19 of 300 injections (all 19 at u₀=1; 9 more at u₀=1 flagged "best lens fit is unphysical"); nulls stayed 0/300. Detection efficiencies and t_E errors are unchanged, since they do not depend on the verdict.
