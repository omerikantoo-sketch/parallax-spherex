"""E4: Order-of-magnitude microlensing yield for SPHEREx.

Gamma = (2/pi) N tau / <tE>  (events with u0 < 1 per year; standard relation
between optical depth and event rate for a fixed mean timescale).
Detectable = Gamma * epsilon, epsilon taken from the E3 injection-recovery
benchmark (joint method, u0 = 1.0 -- a conservative proxy for u0 < 1 events)
for the matching cadence.  Every input is a stated assumption with a low/high
range; the output is an interval, not a prediction.

Literature orders of magnitude used for tau (to be checked against the
original papers before any quantitative citation):
  Galactic bulge  ~1e-6 to 3e-6   (MOA-II, OGLE-IV bulge surveys)
  Galactic disk   ~3e-8 to 3e-7   (OGLE-IV disk survey; strong longitude dependence)
  LMC direction   ~1e-8 to 1e-7   (MACHO vs EROS/OGLE disagreement)
  high latitude   ~1e-9 to 1e-8   (disk self-lensing; very uncertain)
N (stars SPHEREx resolves above ~5% photometric precision per visit) is the
largest uncertainty: 6.2" pixels are confusion-limited toward the bulge.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]

SCENARIOS = [
    # name, cadence, N (lo, mid, hi), tau (lo, mid, hi), <tE> days
    ("Galactic bulge, all-sky survey", "wide", (3e6, 1e7, 3e7), (1e-6, 2e-6, 3e-6), 20.0),
    ("Galactic disk |b|<5 deg, all-sky survey", "wide", (3e7, 1e8, 3e8), (3e-8, 1e-7, 3e-7), 40.0),
    ("South ecliptic pole deep field (LMC direction, if overlapping)", "deep", (3e4, 1e5, 1e6), (1e-8, 3e-8, 1e-7), 40.0),
    ("North ecliptic pole deep field (high latitude)", "deep", (1e4, 3e4, 1e5), (1e-9, 3e-9, 1e-8), 30.0),
]


def rate(n, tau, te):
    return 2.0 / np.pi * n * tau * 365.25 / te


def main():
    eff = {}
    for cad in ("deep", "wide"):
        p = ROOT / "reports" / f"gravity_benchmark_{cad}.json"
        if p.exists():
            b = json.loads(p.read_text())
            eff[cad] = float(b["methods"]["joint"]["efficiency"][b["u0_values"].index(1.0)])
    rows = []
    for name, cad, N, tau, te in SCENARIOS:
        e = eff.get(cad)
        g = [rate(N[0], tau[0], te), rate(N[1], tau[1], te), rate(N[2], tau[2], te)]
        rows.append({"name": name, "cadence": cad, "n_stars": N[1], "n_range": [N[0], N[2]], "tau": tau[1],
                     "tau_range": [tau[0], tau[2]], "tE_days": te, "events_per_yr": g[1],
                     "events_range": [g[0], g[2]], "efficiency": e,
                     "detectable_per_yr": g[1] * e if e is not None else None,
                     "detectable_range": [g[0] * e, g[2] * e] if e is not None else None})
    res = {
        "experiment": "E4 gravity feasibility (order of magnitude)", "provenance": "model estimate",
        "summary": ("Yield = (2/pi) N tau / <tE> x detection efficiency from the E3 injection-recovery benchmark. "
                    "Inputs are ranges, not measurements; the central column uses mid values."),
        "scenarios": rows, "efficiency_source": {k: f"reports/gravity_benchmark_{k}.json, joint method, u0=1.0" for k in eff},
        "caveats": ("Efficiency was measured on a bright (16 AB at 2 um), isolated synthetic source with 2% noise; "
                    "real bulge photometry at 6.2 arcsec pixels is confusion-limited, so bulge/disk numbers are upper "
                    "bounds. The most realistic near-term use is spectral blend decomposition of events already "
                    "found by ground surveys (OGLE, MOA, KMTNet) that happen to fall in a SPHEREx visit. "
                    "No event has been searched for or found in real SPHEREx data in this project. "
                    "With the wide-survey cadence a 'detection' is a significant brightening during one visit; the "
                    "timescale is then poorly constrained (E3: median |ln tE error| about 1, a factor ~2.7), so most "
                    "such events could not be characterised from SPHEREx alone."),
    }
    out = ROOT / "reports" / "gravity_feasibility.json"
    out.write_text(json.dumps(res, indent=1))
    for r in rows:
        print(r["name"], f"{r['events_per_yr']:.3g}/yr", r["efficiency"], r["detectable_per_yr"])


if __name__ == "__main__":
    main()
