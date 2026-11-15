# Next steps (ranked by scientific value per effort)

1. **Real-background injection-recovery.** Inject synthetic movers into real SPHEREx frames (`inject_into_frames`, provenance "real+injected") at faint, uncrowded positions. Calibrate thresholds on many time-scrambled realisations of the same field, then repeat E1 there. This is the step that turns the synthetic improvement into a real-data one.
2. **Scene model for bright stars.** Fit a wavelength-dependent ePSF and halo near saturated stars, or use larger masks driven by the flag planes. Near Barnard's Star the real-field search is currently residual-limited (z ≈ 230).
3. **Parallax on real data.** Fit parallax freely for Barnard's Star using three visits and the spacecraft positions, and compare with Gaia.
4. **Brown-dwarf targets.** Run WISE 0855−0714 and ε Indi B once their coverage indices finish. For those spectra the spectral bank should matter most (E1: +0.65 mag).
5. **Real-photometry noise model for Gravity.** Use an empirical scatter from time-scrambled residuals and a binned (per-channel) SED basis, so the adequacy guard can pass on real data.
6. **Gravity discrimination.** Add light-curve shape tests (wing profile, colour constancy of A(t)) so that discrimination does not rest on the physical-validity guard alone (which rejects 9 of the 28 true lenses that BIC alone labelled lens in E3) and is tested against more variable families than a Gaussian bump. Apply the method to known ground-survey events that fall in SPHEREx visits.
7. **Velocity-space pruning.** A heuristic coarse-to-fine search with measured recall loss, since the exact bound does not prune.
8. **Restore IRSA SIA** as the primary discovery path when the service is reachable. The S3 footprint index is the fallback.
