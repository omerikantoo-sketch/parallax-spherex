"""Build footprint coverage for challenge targets from the AWS SPHEREx mirror.

usage: python scripts/build_coverage.py [slug ...] [--release dr1]
"""
from __future__ import annotations

import argparse
import logging
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from parallax.data import archive as A  # noqa: E402
from parallax.data.targets import TARGETS  # noqa: E402

logging.getLogger("httpx").setLevel(logging.WARNING)
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s %(message)s")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("slugs", nargs="*", default=list(TARGETS))
    ap.add_argument("--release", default="dr1")
    args = ap.parse_args()
    for slug in args.slugs:
        t = TARGETS[slug]
        ra, dec = (float(v) for v in t.at(61000.0))  # mid-2026 position for discovery
        t0 = time.time()
        hits, how = A.find_coverage(ra, dec, args.release,
                                    progress=lambda m, f: logging.info("%s %s %.0f%%", slug, m, 100 * f))
        logging.info("%s: %d exposures via %s in %.0fs", slug, len(hits), how, time.time() - t0)
        for h in hits[:200]:
            logging.info("  %s det%d mjd=%.4f", h["obsid"], h["detector"], h["mjd"])


if __name__ == "__main__":
    main()
