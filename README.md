# PARALLAX

A SPHEREx time-domain viewer and physics-driven inference engine, built for the
NASA Space Apps Challenge 2026 challenge “Planet X and SPHEREx”.

* **Observatory**: enter any RA/Dec (or pick a test target) and PARALLAX finds every SPHEREx exposure covering it. It plays them in time order, every frame labelled with wavelength, detector, UTC and observation ID. It also offers epoch-to-epoch differences, co-adds, a per-pixel spectrum or light curve, and proper-motion measurement for a known star.
* **Motion**: joint spectro-orbital matched filtering for faint moving objects, benchmarked against shift-and-stack, per-detector and detect-then-link baselines at a null-calibrated false-alarm rate.
* **Gravity**: point-lens microlensing with chromatic blending, compared with intrinsic variability. Includes a feasibility estimate for SPHEREx.
* **Results** and **Data & Methods**: every claim with its evidence level, plus data credits, SI disclosure and limitations.

## Run

```bash
python -m venv .venv
.venv/Scripts/pip install numpy scipy astropy fastapi uvicorn httpx pytest   # (Linux/macOS: .venv/bin/pip)
.venv/Scripts/python -m uvicorn api.app:app --port 8600
# open http://localhost:8600
```

Data come live from the public SPHEREx archive: IRSA SIA when reachable, otherwise the AWS Open Data mirror `nasa-irsa-spherex`. Everything fetched is cached under `cache/`, so fields opened once (and the pre-built demo fields) work offline.

## Verify

```bash
.venv/Scripts/python -m pytest -q tests
```

See [docs/SCIENCE.md](docs/SCIENCE.md) for methods and results and [research/PRIOR_ART.md](research/PRIOR_ART.md) for prior art and the novelty assessment. Status and limitations are in [CURRENT_STATUS.md](CURRENT_STATUS.md), [KNOWN_LIMITATIONS.md](KNOWN_LIMITATIONS.md), [EXPERIMENTS.md](EXPERIMENTS.md) and [SUBMISSION.md](SUBMISSION.md).

## Layout

```
api/            FastAPI server (JSON + binary image stacks) and static hosting
parallax/data   archive access: remote FITS range reads, RICE flags, ePSF, WCS-WAVE, coverage index
parallax/science PSF, spectra, SPHEREx-like simulator
parallax/inference ψ/φ maps, statistics, search, branch-and-bound, baselines, scene/astrometry
parallax/gravity point lens + chromatic blend, hypothesis comparison, mass degeneracy
parallax/validation benchmarks E1, E3, E4
scripts/        coverage builder, real-data validation (E2)
web/            Observatory frontend (vanilla ES modules, no build step)
reports/        JSON outputs of every experiment
tests/          pytest suite
```

## Credits

The SPHEREx data are from a joint project of JPL and Caltech, funded by NASA. They are served by the NASA/IPAC Infrared Science Archive. Code and documentation were written with an AI coding assistant (Anthropic Claude, via Claude Code) under the team's direction; no AI-generated images are used. Full disclosure in SUBMISSION.md.
