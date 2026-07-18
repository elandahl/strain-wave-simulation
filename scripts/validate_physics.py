#!/usr/bin/env python3
"""Run the physics acceptance suite for the Cr/GaAs strain models.

This packages the model-independent, physics-based sanity checks that should be
green *before* any comparison to experimental data (APS 7ID, PLS) or extension
to new substrates. It complements ``validate_fd_courant.py`` (numerical
equivalence of the C=1 FD model to d'Alembert) with checks that tie the models
to physical expectations:

  * acoustic pulse-train echo spacing = 2*L_film*v_GaAs/v_Cr,
  * echo-to-echo decay set by the Cr/GaAs impedance mismatch |r|,
  * echo spacing scales linearly with film thickness,
  * lumped electron-phonon thermalization is picosecond-scale and ~1/G,
  * leapfrog-vs-d'Alembert far-field error grows with propagation distance
    (confirming the numerical-dispersion diagnosis).

Writes ``docs/physics_acceptance.json`` and, unless ``--no-figure`` is given,
``docs/images/physics_acceptance_pulse_train.png``. Exits nonzero if any check
fails so it can be wired into CI.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from strain_wave import run_simulation
from strain_wave.acceptance import (
    _acoustic_far_field,
    _paper_config,
    expected_echo_spacing_nm,
    measure_echo_spacing_nm,
    run_acceptance_suite,
)

REPO = Path(__file__).resolve().parents[1]
REPORT = REPO / "docs" / "physics_acceptance.json"
FIGURE = REPO / "docs" / "images" / "physics_acceptance_pulse_train.png"


def _jsonable(obj):
    """Recursively convert numpy scalars/arrays to plain Python for json.dump."""
    if isinstance(obj, dict):
        return {k: _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, np.generic):
        return obj.item()
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    return obj


def _write_pulse_train_figure(path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from strain_wave.models.ttm_cr_gaas_solver import V_GAAS

    cfg = _paper_config("ttm_dalembert_cr_gaas", t_max=200e-12, L_film=80e-9, L_sub=1.2e-6)
    res = run_simulation(cfg, verbose=False)
    z_nm, acoustic, monitor = _acoustic_far_field(res, cfg)
    front_nm = cfg.t_max * V_GAAS * 1e9
    sel = (np.arange(len(z_nm)) > monitor) & (z_nm < front_nm - 10.0)

    spacing = measure_echo_spacing_nm(res, cfg)
    expected = expected_echo_spacing_nm(cfg)

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(z_nm[sel], acoustic[sel] * 1e3, lw=1.2, color="#1f77b4")
    # Mark expected echo positions stepping back from the leading front.
    n_echoes = int((front_nm - z_nm[monitor]) / expected)
    for k in range(1, n_echoes + 1):
        ax.axvline(
            front_nm - k * expected,
            color="0.6",
            ls="--",
            lw=0.8,
            label="analytic echo spacing" if k == 1 else None,
        )
    ax.set_xlabel("depth into substrate z (nm)")
    ax.set_ylabel("acoustic strain (×10⁻³)")
    ax.set_title(
        "Cr/GaAs acoustic pulse train (d'Alembert far field)\n"
        f"measured spacing {spacing:.1f} nm vs analytic {expected:.1f} nm "
        "(= 2·L_film·v_GaAs/v_Cr)"
    )
    ax.legend(loc="upper left", fontsize=9)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--model",
        default="ttm_dalembert_cr_gaas",
        help="strain model to run the pulse-train checks against",
    )
    parser.add_argument(
        "--no-figure", action="store_true", help="skip writing the diagnostic figure"
    )
    args = parser.parse_args()

    report = run_acceptance_suite(model=args.model)

    print(f"Physics acceptance suite ({args.model})\n")
    for check in report.checks:
        status = "PASS" if check.passed else "FAIL"
        print(f"  [{status}] {check.name}: {check.detail}")

    payload = _jsonable({"model": args.model, **report.to_dict()})
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(payload, indent=2) + "\n")
    print(
        f"\n{payload['n_passed']}/{payload['n_total']} checks passed. "
        f"Report saved to {REPORT}"
    )

    if not args.no_figure:
        _write_pulse_train_figure(FIGURE)
        print(f"Figure saved to {FIGURE}")

    print("Acceptance", "PASSED" if report.passed else "FAILED")
    return 0 if report.passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
