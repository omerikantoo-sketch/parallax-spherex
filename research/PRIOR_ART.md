# Prior-art matrix

Bibliographic details below were recalled during a time-limited hackathon and
**could not be checked online from this environment** (IRSA/ADS access was
limited). Before any publication, verify every reference against ADS. Items
marked *(verify)* are ones where the recall is uncertain.

## Faint moving objects

| Approach | Representative reference | What it does | Relation to PARALLAX Motion |
|---|---|---|---|
| Shift-and-stack along trial orbits | Bernstein et al. 2004, AJ 128, 1364 (HST TNO search) | Co-adds images along hypothesised linear/orbital motion to detect sources below single-image limit | Our `stack` baseline is this, with PSF-weighted ψ/φ maps |
| GPU shift-and-stack with ψ/φ likelihood maps | Whidden et al. 2019, AJ 157, 119 (KBMOD) | Per-image matched-filter "ψ, φ" maps, max-likelihood trajectory search, region search | The constant-flux statistic and the ψ/φ formalism are KBMOD's. KBMOD region search is prior art for velocity-space pruning |
| Synthetic tracking | Shao et al. 2014, ApJ 782, 1 | Shift-and-add of short exposures for fast NEOs | Same family as the stack baseline |
| Matched-filter asteroid detection | Gural, Larsen & Gleason 2005, AJ 130, 1951 | Velocity-matched filter over image sequences | Same family |
| Detect-then-link | e.g. THOR, Moeyens et al. 2021, AJ 162, 143; HelioLinC (Heinze et al.) *(verify)* | Single-epoch detections are linked into tracklets/orbits | Our `link` baseline (≥3 exposures at 5σ) is a minimal version |
| Shift-stacking for Planet Nine in survey data | Rice & Laughlin 2020 (TESS shift-stacking) *(verify)* | Shift-stack along slow outer-solar-system motions in wide-field data | Same idea in a different survey; motivates the challenge |
| Multi-band co-addition with assumed colours | Common practice in multi-filter TNO surveys | Weights bands by an assumed colour | The template statistic generalises this to SPHEREx's per-exposure wavelength; **this is the closest prior art to our spectral bank** |
| Forward-modelled photometry | Lang, Hogg & Mykytyn 2016, The Tractor (ASCL) | Joint forward modelling of pixels | Conceptual model for `forced_photometry` / `subtract_static` |
| Difference imaging | Alard & Lupton 1998, ApJ 503, 325 | PSF-matched subtraction | We use static-scene model subtraction instead (PSF varies across the LVF) |
| SPHEREx tools | IRSA SPHEREx Spectrophotometry Tool and Source Discovery Tool (listed in the challenge resources) *(not inspected)* | Archive-side spectrophotometry and discovery | Our cutout and forced-photometry steps overlap. We did not inspect their algorithms |

## Microlensing

| Approach | Reference | Relation to PARALLAX Gravity |
|---|---|---|
| Point-source point-lens magnification | Paczyński 1986, ApJ 304, 1 | `magnification`, `u_of_t` |
| Blending in microlensing | Woźniak & Paczyński 1997, ApJ 487, 55 *(verify)* | Blend flux per band is standard; we parametrise source and blend SEDs smoothly in λ |
| Multi-band fits with free blend per band | Standard in OGLE/MOA/KMTNet analyses | Equivalent to our model with an unconstrained SED. The smooth basis only adds a regulariser, so **no novelty is claimed for the lens model itself** |
| Mass–distance degeneracy, θE formalism | Gould 2000, ApJ 542, 785 | `mass_posterior` makes the degeneracy explicit |
| Bulge optical depth | Sumi et al. 2013, ApJ 778, 150 (MOA-II); Mróz et al. 2019, ApJS 244, 29 (OGLE-IV) | Inputs to the feasibility estimate (order of magnitude only) |
| LMC optical depth | Alcock et al. 2000, ApJ 542, 281 (MACHO); Wyrzykowski et al. 2011, MNRAS 416, 2949 (OGLE-III) | Range used for the ecliptic-pole scenario |
| Space-based microlensing surveys | Johnson et al. 2020, AJ 160, 123 (Roman free-floating planet yield) *(verify)* | Contrast: Roman has a 15-min cadence. SPHEREx's cadence is unsuited to discovery outside the deep fields |

## SPHEREx mission

| Reference | Use |
|---|---|
| Doré et al. 2014, arXiv:1412.4872 | Mission concept, LVF spectral sampling |
| SPHEREx Explanatory Supplement / IRSA documentation | FITS layout (IMAGE, FLAGS, VARIANCE, ZODI, PSF, WCS-WAVE); we inferred details from real headers and recorded them in `parallax/data/archive.py` |

## Novelty assessment (honest)

* **Spectral-template matched filter for moving objects on LVF data.** Mathematically this is a matched filter with a known source spectrum. That is standard detection theory, and close to colour-weighted multi-band stacking. What is specific here:
  1. it is applied to SPHEREx's sampling, where *every exposure is a different wavelength*;
  2. we compare at matched, null-calibrated false-alarm probability against strong baselines;
  3. the measured gains follow a closed-form prediction (the SNR ratio R).

  We claim an *application and quantified evaluation*, not a new algorithm.
* **Rigorous velocity-space branch-and-bound.** Exact, but it gave no speed-up in our tests. KBMOD's region search is prior art. This is reported as a negative result.
* **Chromatic-blend PSPL fit.** Standard. The E3 test shows that it **does not** improve detection over a grey multi-band fit. It does improve the timescale accuracy of detected events in our simulations. Any claim is limited to that.
* **Exotic lenses (black holes, primordial objects, "baby universes").** Out of scope. A light curve alone cannot identify a black hole, and no such claim is made.
