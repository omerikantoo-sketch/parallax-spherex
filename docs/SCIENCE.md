# PARALLAX: science and methods

Status vocabulary: **Implemented**, **Tested** (unit tests), **Validated, synthetic**,
**Validated, real**, **Improvement demonstrated** (against a named baseline, at a
matched false-alarm rate), **Candidate**, **Confirmed**, **Exploratory**.

## 1. Observation model (SPHEREx L2)

Each SPHEREx exposure is a 2040×2040 image per detector (D1–D6), with these extensions:
- `IMAGE`: MJy/sr, zodiacal light *not* subtracted;
- `VARIANCE`;
- `FLAGS`: a bitmask whose names are read from the `MSKN` header cards;
- `ZODI`: a model of zodiacal light;
- `PSF`: 441 oversampled (×5) ePSFs on a zone grid;
- `WCS-WAVE`: a 9×9 control grid of central wavelength and bandwidth.

A linear variable filter makes the wavelength change along one detector axis. A source therefore moves through wavelength as the pointing steps. At Barnard's Star we measured bursts of 2–4 exposures, 0.0015 d apart, each moving the source by ~115 px (about one spectral channel), and ~25 bursts per ~20-day visit, twice a year.

For pixel *p* in exposure *e*, the measurement model is

    d_ep = F(λ_ep, t_e) · P_e(x_ep − x_e(θ)) + B_ep + n_ep,   n ~ N(0, σ_ep²)

with these terms:
- F: the source spectrum at the pixel's wavelength (from `WCS-WAVE`), in µJy;
- P_e: the zone ePSF (sum = 1 on the native grid);
- B: background (zodi model plus a sigma-clipped residual), with static sources removed by forced fitting;
- θ: the trajectory.

Flagged pixels are excluded and stay missing (never filled in). Conversion: µJy/px = MJy/sr × Ω_px × 2.35e-11 × 1e12 (≈ 890 for Ω ≈ 37.9 arcsec²).

## 2. Motion: joint spectro-orbital inference

For each exposure we precompute ψ_e = Σ_p w P d and φ_e = Σ_p w P² on a ×2 sub-pixel grid (the KBMOD ψ/φ formalism). For a trajectory θ:

| Statistic | Formula | Meaning |
|---|---|---|
| constant (stack) | Σψ / √Σφ | KBMOD-like, flat spectrum |
| template *s* | Σ sψ / √Σ s²φ | known spectral shape |
| bank | max over 10 physical templates | reflected light ± thermal, M/L/T/Y-like dwarfs |
| subspace K | χ² of a K-dim smooth spectral model | shape-agnostic |
| per-detector | χ² with one free flux per detector | baseline |
| best detector | max single-detector stack | baseline |
| detect-then-link | ≥3 exposures with per-frame SNR ≥ 5 within 1.5 px | baseline |

For a source with true spectrum *s*, the constant-flux SNR is lower than the matched SNR by R = Σsφ / √(Σφ · Σs²φ). The expected gain in limiting magnitude is therefore 2.5 log₁₀(1/R).

**Calibration.** For each method, the detection threshold is the 95th percentile of its maximum over the *whole* search grid on 40 noise-only fields with the real burst cadence. Every method therefore operates at the same per-field false-alarm probability, including the look-elsewhere effect.

### E1 result (Validated, synthetic; Improvement demonstrated)

`reports/motion_benchmark.json`: 24 injections per 0.25 mag bin, 17.5–22 AB.

| Method | m50 all | reflected | reflected+thermal | cold brown dwarf |
|---|---|---|---|---|
| detect-then-link | 18.64 | 18.70 | 18.83 | 18.08 |
| best single detector | 19.95 | 20.25 | 20.06 | — |
| constant-flux stack | 20.32 | 20.80 | 20.44 | 18.55 |
| per-detector χ² | 20.30 | 20.55 | 20.44 | 19.10 |
| **template bank** | **20.38** | 20.80 | 20.58 | **19.20** |
| true-SED oracle | 20.42 | 20.75 | 20.58 | 19.20 |

The spectral bank gains **+0.65 mag** over constant-flux stacking for cold brown-dwarf spectra. The R prediction is 1.03 mag; the bank does not contain the exact SED, and the bins have finite statistics. For solar-like reflectance the gain is ≤0.14 mag, consistent with the prediction of 0.07–0.09 mag. Overall (equal mix of classes) the gain is +0.05 mag.

Shift-and-stack of any kind beats detect-then-link by ~1.7 mag. Null fields produced no false tracklets in the link baseline. Each m50 bin has 8 injections per class, so per-class m50 values are uncertain by roughly ±0.15 mag.

### Branch-and-bound (Negative result)

We implemented a rigorous upper bound over velocity cells (max-pooled ψ, min-pooled φ, plus one q-pixel rounding slack). It returns exactly the exhaustive set: 70/70 trajectories in development, and it is enforced by `tests/test_core.py::test_branch_and_bound_is_exact`. At 60 exposures and threshold 6, however, the bound pruned nothing, because the long time baseline inflates the pooling boxes. No speed-up is claimed. KBMOD's region search is prior art.

## 3. Gravity: point lens with chromatic blending

The model is F(t, λ) = A(t; t₀, u₀, t_E) S(λ) + B(λ), with Paczyński A(u). S is expanded on K smooth basis functions and B on 2K. The amplitudes are profiled by weighted least squares, so the nonlinear search is 3-D (a grid, then Nelder–Mead). The hypotheses are nested:
- H0: constant, 2K parameters;
- H1: lens, 3K+3;
- H2: achromatic Gaussian bump, 3K+2.

Comparison uses BIC with |ΔBIC| > 10, after an adequacy check (best reduced χ² ≤ 3) and, for a lens verdict, a physical-validity check.

### E3 result (Validated, synthetic)

`reports/gravity_benchmark_{deep,wide}.json`: 300 nulls, 60 events per u₀, 100 symmetric variables. Noise is 2%, the source is 6000 K and the blend 3200 K, and FAP is 0.01 per light curve.

* **Detection:** the joint chromatic fit is **not** better than a grey (single blend fraction) multi-band fit: −0.08 to +0.02 in efficiency. Its extra freedom raises the null threshold (Δχ² 24.9 vs 17.1). Both multi-band fits beat a single-band fit strongly: deep cadence +0.65 at u₀ = 2, wide cadence +0.28.
* **Characterisation:** for detected u₀ = 1 events, the median |ln(t_E,fit/t_E,true)| is 0.15 (joint) vs 0.42 (grey) in the deep cadence, and 0.96 vs 1.53 in the wide cadence.
* **False lenses:** with BIC alone, 30/100 symmetric achromatic variables in the deep cadence (4/100 in the wide) received a "lens" verdict, so BIC by itself does not separate a symmetric bump from a lens. With the app's guarded rule (adequacy, then physical validity: 0 < f_s ≤ 1.05, t_E ≤ 3× span), **0/100** variables and 0/300 nulls get a lens verdict in either cadence. The guard has a cost: true lenses labelled lens drop from 28 to 19 of 300 (deep) and from 5 to 3 (wide); the rest become "ambiguous". The guard was tested against one variable family only, so any real candidate would still need light-curve shape tests (wings, colour constancy of A) and external data.

### E4 feasibility (Exploratory, model estimate)

The yield is Γ = (2/π) N τ / ⟨t_E⟩ × ε, with ε taken from E3. See `reports/gravity_feasibility.json`. Bulge and disk numbers are large because N τ is large, but they are upper bounds: 6.2″ pixels are confusion-limited, and ε was measured for an isolated bright star. In the wide cadence t_E is constrained only to a factor of ~2.7.

The deep-field scenarios give ≲0.02 events per year. The realistic use is spectral blend decomposition of ground-survey events that fall in a SPHEREx visit. **No real event was searched for or found.**

Mass: t_E alone gives a prior-dominated interval spanning more than a factor of 10 (`mass_posterior`). No object type (black hole, planet, primordial object) can be inferred from a single light curve.

## 4. Real-data checks (E2, Validated, real)

`reports/real_validation.json`, Barnard's Star field: 160 exposures, 54 with usable ePSF centroids, 3 visits over 262 d.

| | SPHEREx (this work) | Gaia DR3 |
|---|---|---|
| μα* [mas/yr] | −2404 ± 601 (stat) ± 1673 (sys) | −802 |
| μδ [mas/yr] | 12183 ± 597 (stat) ± 1673 (sys) | 10362 |

* **Per-visit positions.** Each visit's median position agrees with the Gaia-predicted position (proper motion plus parallax, seen from SPHEREx) to ≤1.2″ (0.2 px). The star moves 7.4″ over the span. Our parallax model was checked independently against astropy's GCRS transformation from the header spacecraft position: <0.1″ difference in ξ, <0.3″ in η.
* **Statistical errors.** The statistical-only errors under-represent per-visit systematics of ~0.7″, which likely come from saturated-core centroiding and the absolute WCS between visits. D1–D3 (few points, saturated) show offsets up to 3.5″.
* **Systematic floor.** The floor is derived from offsets relative to Gaia, so "consistent within 1σ" is not an independent test. The robust results are that motion is detected and that the per-visit agreement is ≤1.2″.
* **Real-field search (Exploratory).** Bright-star residuals limit the search (z ≈ 240) before noise does. One of eight listed trajectories slightly exceeds the maximum from a single time-scrambled null (242.9 vs 239.2); this is not a candidate.

## 5. Reproducibility

```bash
python -m pytest -q tests
python -m parallax.validation.motion_benchmark --null 40 --per-mag 21 --out reports/_motion_null_v1.json   # null calibration (v1)
python -m parallax.validation.motion_benchmark --per-mag 24 --mags 17.5:22.25:0.25 --reuse-null reports/_motion_null_v1.json
python -m parallax.validation.gravity_benchmark --cadence deep
python -m parallax.validation.gravity_benchmark --cadence wide
python -m parallax.validation.gravity_feasibility
python scripts/build_coverage.py barnard
python scripts/real_validation.py barnard search
```

Seeds are fixed in each script: nulls, injections and variables use disjoint seed ranges.
