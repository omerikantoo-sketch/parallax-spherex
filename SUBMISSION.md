# PARALLAX: Space Apps 2026 submission notes

**Challenge:** Planet X and SPHEREx.

## Summary (≈100 words)

PARALLAX is an observatory for SPHEREx's time dimension. Type any sky position and it finds every SPHEREx exposure that covers it, directly from the public archive. It plays them in time order, labelled with wavelength, detector and time, and lets you compare epochs, difference them and measure motion.

On real data it recovers the proper motion of Barnard's Star. Its search engine treats SPHEREx's frame-by-frame wavelength sampling as a feature: matching a spectrum while stacking along trial orbits finds cold, brown-dwarf-like movers about 0.65 mag fainter than standard stacking at the same false-alarm rate (simulated).

A Gravity mode tests microlensing against variability, with an honest SPHEREx feasibility estimate.

## NASA data and resources used

* SPHEREx Level-2 spectral images (DR1), from the NASA/IPAC Infrared Science Archive via the AWS Open Data Registry `nasa-irsa-spherex`, with IRSA SIA when reachable.
* SPHEREx data-product documentation (FITS layout: IMAGE, FLAGS, VARIANCE, ZODI, PSF, WCS-WAVE).

**Citation text** (also shown in the app):

> This work uses data from SPHEREx, a joint project of the Jet Propulsion Laboratory and the California Institute of Technology, funded by NASA. This research has made use of the NASA/IPAC Infrared Science Archive, which is funded by NASA and operated by the California Institute of Technology.

Please also cite the SPHEREx DR1 dataset DOI listed by IRSA.

No NASA logos or branding are used in the interface or the generated content.

## Use of AI (disclosure)

* **Code and documentation** were written with an AI coding assistant (Anthropic Claude, via Claude Code) under the team's direction. The repository's README, About page and this file acknowledge this.
* **Data values:** every number reported comes from running the code in this repository on archive data or on the clearly labelled simulator. The AI did not invent results.
* **Images:** no AI-generated images or video are used. If a demo video is made with AI-generated footage, it must carry a visible watermark per the rules.
* The project idea and challenge choice must be the team's own (rule: AI should not be used to come up with the idea). The team should confirm this statement before submitting.

## Demo: 7 slides

1. **The hidden dimension.** SPHEREx samples a different wavelength in every exposure, and revisits the sky. Screenshot: Observatory playing Barnard's Star across three visits.
2. **Observatory.** RA/Dec in, every exposure out, labelled with λ, detector, UTC and observation ID; epoch comparison, difference and co-add. Live from the archive (works when IRSA SIA is down).
3. **Real motion.** Barnard's Star proper motion measured from SPHEREx centroids vs Gaia DR3 (numbers from `reports/real_validation.json`).
4. **Spectrum-aware search.** Matched filter over orbit × spectrum. E1 completeness curves: +0.65 mag for cold brown-dwarf spectra, ~0 for solar-like reflectance, as predicted by the SNR ratio R.
5. **What did not work.** Exact branch-and-bound gives no speed-up. On real fields bright-star residuals dominate, and no candidate exceeds the time-scrambled null.
6. **Gravity.** Point lens with chromatic blending vs variability. No detection gain over a grey fit, but better timescales. Feasibility: SPHEREx is a spectral follow-up tool for known events, not a discovery engine.
7. **Honesty and next steps.** Evidence matrix (Results page), real-background injection, brown-dwarf targets.
