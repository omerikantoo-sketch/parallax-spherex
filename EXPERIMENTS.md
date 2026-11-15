# Experiments

Each experiment has a script, fixed seeds and a JSON report in `reports/`.

| ID | Question | Script | Report | Status |
|---|---|---|---|---|
| E0 | Can we read real SPHEREx L2 without IRSA (headers, flags, ePSF, λ)? | `parallax/data/archive.py`, `scripts/build_coverage.py` | `cache/index.sqlite`, `logs/coverage.log` | Validated, real: Barnard's Star 238 exposures found by exact WCS containment |
| E1 | Does a spectral-template matched filter detect fainter moving objects than constant-flux stacking at equal FAP? | `parallax/validation/motion_benchmark.py` | `motion_benchmark.json` | Validated, synthetic. Improvement demonstrated for non-flat SEDs (+0.65 mag for cold brown dwarfs); negligible for solar-like reflectance |
| E1b | Does exact velocity branch-and-bound prune enough to be faster? | `parallax/inference/matched.py::branch_and_bound` | `tests/test_core.py` (exactness) | Negative result: exact, but no pruning at threshold 6, 60 exposures |
| E2 | Does the pipeline recover the known proper motion of Barnard's Star from real SPHEREx images? | `scripts/real_validation.py barnard` | `real_validation.json` | see CURRENT_STATUS.md |
| E2b | What does the search return on a real field? | `scripts/real_validation.py search` | `real_validation.json` | Exploratory |
| E3 | Joint chromatic PSPL vs grey multi-band vs single-band: detection, timescale accuracy, lens/variable confusion | `parallax/validation/gravity_benchmark.py` | `gravity_benchmark_{deep,wide}.json` | Validated, synthetic. No detection gain over grey; better t_E; single-band clearly worse; guarded verdict: 0/100 symmetric variables and 0/300 nulls labelled lens (BIC alone: 30/100 deep, 4/100 wide); true-lens verdicts 19/300 deep, 3/300 wide |
| E4 | Order-of-magnitude microlensing yield for SPHEREx | `parallax/validation/gravity_feasibility.py` | `gravity_feasibility.json` | Exploratory (model estimate) |

## Protocol notes

* **Seeds:** nulls 1000+ (E1) / 100000+ (E3); injections 5000+ / 200000+; variables 300000+. Calibration seeds are never reused for evaluation.
* **Thresholds:** the (1 − FAP) quantile of each method's *global* maximum on null data, so the look-elsewhere effect is included.
* **Baselines:** each runs at full strength on identical data. The oracle (true-SED template) is reported as an upper bound, not as a method.
* **E1 v1 → v2:**
  - v1 used a magnitude grid of 19.5–22.5 in 0.5 steps, which was too coarse for brown dwarfs (m50 undefined).
  - v2 reuses the same 40 null fields (`reports/_motion_null_v1.json`) and refines the grid to 17.5–22 in 0.25 steps, with 24 injections per bin.
* **E3 baseline fix:** the single-band baseline first used a flat in-band SED, which gave a null Δχ² ≈ 190 because an SED slope was mistaken for a lens. It was fixed before the reported run by modelling a linear in-band slope for both source and blend.
