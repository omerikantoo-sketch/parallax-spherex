# Known limitations

1. **Synthetic benchmarks**:
   - E1 and E3 use noise-limited simulated fields, with a Gaussian PSF (E1) and white noise.
   - Real fields add confusion, PSF variation across the LVF, artefacts (ghosts, persistence, cosmic rays) and correlated noise.
   - Noise σ for D2–D6 in the simulator are placeholders; D1 was measured.
2. **Real-field search significance**: computed under a Gaussian null. It is not calibrated on real null fields, so candidates are *Exploratory* and need vetting.
3. **Coverage without IRSA SIA**:
   - The footprint index samples every 12th header per survey period, then densely scans nearby periods.
   - A very short period near the position could in principle be missed.
   - Building the index for a new position takes ~5–20 min the first time.
4. **Frame wavelength label**: the label is the wavelength at the cutout centre. Inference uses per-pixel wavelengths.
5. **Display resampling**: the display uses north-up bilinear resampling. Inference always uses native pixels.
6. **Astrometry**:
   - Linear motion plus a fixed parallax.
   - Bright stars can be flagged as saturated or nonlinear; such frames are excluded.
   - No absolute astrometric calibration beyond the archive WCS.
7. **Gravity**:
   - Point-source point-lens only: no finite source, binary lens or orbital parallax.
   - BIC with a fixed threshold.
   - BIC alone mis-labels 30/100 symmetric variables as lenses (E3, deep cadence). The physical-validity guard removes all of these in E3, but it also turns 9 of the 28 lens verdicts on true lenses into "unphysical"; it has only been tested against one variable family (an achromatic Gaussian bump).
8. **Gravity on real photometry**: the adequacy guard currently refuses every tested real position.
   - The reduced χ² is ~8 on faint or blank positions and ~400 on bright stars.
   - Causes: archive variance underestimates the forced-photometry scatter; a smooth SED basis cannot fit real stellar molecular bands at high S/N; bright stars saturate.
   - A binned SED basis and an empirical noise model are needed before real-data lens inference.
9. **Feasibility**: order-of-magnitude only. N (resolvable stars) and τ are literature ranges; the bulge and disk numbers are upper bounds.
10. **References** in research/PRIOR_ART.md were recalled without online verification.
11. **Targets**: positions and proper motions are hard-coded from published catalogues. The ε Indi B and WISE 0855 values are approximate (see `parallax/data/targets.py`).
