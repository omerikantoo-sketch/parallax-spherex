# Current status (2026-11-15)

| Component | Status | Evidence |
|---|---|---|
| Archive access without IRSA (AWS mirror, range reads, RICE flags, ePSF, WCS-WAVE with per-detector grid size) | Validated, real | Barnard's Star: 238 exposures found; 18/18 tests pass, including real-cutout unit and wavelength checks |
| Observatory: RA/Dec → all exposures, time-ordered playback, λ / detector / UTC / obs-ID labels, epoch comparison, difference, co-add, flags overlay, spectrum / light curve, permalink | Implemented, checked in browser | Barnard's Star field (160 frames) |
| Proper motion of a known star from SPHEREx | Validated, real (with systematic floor) | Motion detected; per-visit positions within 1.2″ of the Gaia prediction; formal errors underestimate ~0.7″ per-visit systematics |
| Spectral-template matched filter for moving objects | Improvement demonstrated (synthetic) | E1: +0.65 mag m50 for cold brown-dwarf SEDs vs constant-flux stacking at equal FAP; ≈0 for solar-like reflectance, as predicted |
| Exact velocity branch-and-bound | Negative result | Exact (test) but no pruning at useful thresholds |
| Search on real fields | Exploratory | Bright-star residual-limited; no candidate |
| Gravity: PSPL + chromatic blend vs variability | Validated, synthetic | E3: no detection gain over grey fit; better t_E; strong gain over single band; guarded verdict: 0/100 symmetric variables and 0/300 nulls labelled lens (BIC alone: 30/100), at the cost of 19 instead of 28 lens verdicts among 300 true lens injections (deep cadence) |
| Gravity on real SPHEREx photometry | Implemented, guarded | Adequacy guard refuses all tested real positions (reduced χ² 8–419); no real-data lens inference yet |
| Microlensing yield for SPHEREx | Exploratory | E4 order-of-magnitude intervals |
| Candidates / confirmed discoveries | None | — |

## Not done

* Real-background injection-recovery.
* Coverage indices are built for all four test targets (Barnard 238, WISE 0855 255, RW Aur 222, ε Indi B 341 exposures). Fields are pre-built for Barnard's Star and WISE 0855 (160 frames each); RW Aur and ε Indi B fields build on first open (a few minutes).
* No motion search or astrometry has yet been run on the WISE 0855 field.
* Free-parallax fit on real data.
