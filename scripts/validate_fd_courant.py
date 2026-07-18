#!/usr/bin/env python3
"""Validate Courant-matched FD propagation against d'Alembert.

Acceptance case
---------------
Sci. Rep. Fig. 3 Cr/GaAs preset at 1.8 ns.  Beyond the monitor plane the
substrate is homogeneous, source-free (after subtracting local quasistatic
thermal strain), and right-going.  Therefore d'Alembert translation is the
continuum reference and the Courant-one FD field must match it.

Writes ``results/fd_courant_acceptance.json`` and exits nonzero if the
configured numerical tolerances are exceeded.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from strain_wave import get_preset, run_simulation
from strain_wave.models.ttm_dalembert_cr_gaas import MONITOR_OFFSET_CELLS

REPO = Path(__file__).resolve().parents[1]
REPORT = REPO / "results" / "fd_courant_acceptance.json"

# These are deliberately much tighter than any physically meaningful strain
# difference.  At C=1 the residual should be interpolation/roundoff only.
MAX_ABS_TOL = 1e-10
NORMALIZED_RMS_TOL = 1e-7


def main() -> int:
    dalembert = run_simulation(
        get_preset("paper_fig3_gaas", model="ttm_dalembert_cr_gaas")
    )
    fd = run_simulation(
        get_preset("paper_fig3_gaas", model="ttm_fd_courant_cr_gaas")
    )

    monitor = dalembert.n_bin_film - 1 + MONITOR_OFFSET_CELLS
    reference = dalembert.strain[monitor + 1 :]
    candidate = fd.strain[monitor + 1 :]
    delta = candidate - reference

    rms = float(np.sqrt(np.mean(delta**2)))
    reference_rms = float(np.sqrt(np.mean(reference**2)))
    normalized_rms = rms / reference_rms
    max_abs = float(np.max(np.abs(delta)))
    correlation = float(np.corrcoef(reference, candidate)[0, 1])

    reference_peak = int(np.argmax(np.abs(reference)))
    candidate_peak = int(np.argmax(np.abs(candidate)))
    dz_nm = float((dalembert.z[1] - dalembert.z[0]) * 1e9)

    passed = max_abs <= MAX_ABS_TOL and normalized_rms <= NORMALIZED_RMS_TOL
    report = {
        "case": "paper_fig3_gaas",
        "reference_model": "ttm_dalembert_cr_gaas",
        "candidate_model": "ttm_fd_courant_cr_gaas",
        "source_free_region_starts_nm": float(dalembert.z[monitor + 1] * 1e9),
        "max_abs_strain_difference": max_abs,
        "rms_strain_difference": rms,
        "normalized_rms_difference": normalized_rms,
        "correlation": correlation,
        "peak_position_difference_nm": float(
            abs(candidate_peak - reference_peak) * dz_nm
        ),
        "tolerances": {
            "max_abs_strain_difference": MAX_ABS_TOL,
            "normalized_rms_difference": NORMALIZED_RMS_TOL,
        },
        "passed": passed,
    }

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    print(f"\nReport saved to {REPORT}")
    print("Acceptance", "PASSED" if passed else "FAILED")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
