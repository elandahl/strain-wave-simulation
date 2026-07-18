#!/usr/bin/env python3
"""Validate split repos reproduce thermo-elastic-gaas paper outputs."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np

PROJECTS = Path.home() / "Projects"
STRAIN_REPO = PROJECTS / "strain-wave-simulation"
XRD_REPO = PROJECTS / "xrd-strain-simulation"
PAPER_REPO = PROJECTS / "thermo-elastic-gaas"
REFERENCE = PAPER_REPO / "docs" / "reference_outputs.json"
RESULTS = STRAIN_REPO / "results" / "validation_split"


def run(cmd: list[str], cwd: Path) -> None:
    print("$", " ".join(cmd))
    subprocess.run(cmd, cwd=cwd, check=True)


def python_in(repo: Path) -> str:
    venv_python = repo / ".venv" / "bin" / "python"
    if venv_python.exists():
        return str(venv_python)
    return sys.executable


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    strain_file = RESULTS / "strain_profile.npz"
    strain_fig = RESULTS / "strain_figure.png"
    xrd_fig = RESULTS / "xrd_figure.png"

    # Pin the historical leapfrog model: this script checks bit-for-bit
    # equivalence with the frozen thermo-elastic-gaas reference, which uses
    # that solver. The repo-wide default is the d'Alembert model.
    run(
        [
            python_in(STRAIN_REPO),
            "scripts/run.py",
            "--model",
            "ttm_cr_gaas",
            "--no-show",
            "--output",
            str(strain_file),
            "--figure",
            str(strain_fig),
        ],
        STRAIN_REPO,
    )

    run(
        [
            python_in(XRD_REPO),
            "scripts/run.py",
            "--strain-file",
            str(strain_file),
            "--no-show",
            "--output",
            str(xrd_fig),
        ],
        XRD_REPO,
    )

    with np.load(strain_file) as data:
        displacement = data["displacement"]
        strain = data["strain"]

    xrd_python = python_in(XRD_REPO)
    compare_code = """
import json, sys
import numpy as np
from pathlib import Path
from xrd_strain.config import XrdConfig
from xrd_strain.io import load_strain_profile
from xrd_strain.pipeline import run_xrd

strain_file, reference_path = sys.argv[1], sys.argv[2]
with np.load(strain_file) as data:
    displacement, strain = data["displacement"], data["strain"]
profile = load_strain_profile(strain_file)
xrd_result = run_xrd(profile, config=XrdConfig())
with open(reference_path) as f:
    ref = json.load(f)["package_summary"]
checks = {
    "max_displacement_nm": float(np.max(np.abs(displacement)) * 1e9),
    "max_strain": float(np.max(np.abs(strain))),
    "xrd_peak_angle_deg": float(xrd_result.angle_deg[xrd_result.intensity.argmax()]),
    "xrd_peak_log10_intensity": float(np.max(xrd_result.intensity)),
}
print(json.dumps({"checks": checks, "reference": ref}))
"""
    proc = subprocess.run(
        [xrd_python, "-c", compare_code, str(strain_file), str(REFERENCE)],
        cwd=XRD_REPO,
        check=True,
        capture_output=True,
        text=True,
    )
    payload = json.loads(proc.stdout.strip())
    checks = payload["checks"]
    ref = payload["reference"]

    print("\n=== Split pipeline vs paper reference ===")
    passed = True
    for key, expected in ref.items():
        if key not in checks:
            continue
        got = checks[key]
        diff = abs(got - expected)
        ok = diff < 1e-6
        passed &= ok
        status = "OK" if ok else "FAIL"
        print(f"  {key}: got={got:.6g} expected={expected:.6g} diff={diff:.3e} [{status}]")

    report = {"checks": checks, "reference": ref, "passed": passed}
    report_path = RESULTS / "validation_report.json"
    report_path.write_text(json.dumps(report, indent=2))
    print(f"\nReport saved to {report_path}")
    print("Validation", "PASSED" if passed else "FAILED")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
