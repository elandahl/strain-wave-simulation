"""Validation matrix for the Sci. Rep. Fig. 3 case: 2 strain models x 2 instruments.

Runs (all with the paper_fig3_gaas preset):

  strain models : ttm_cr_gaas (leapfrog, paper-faithful)
                  ttm_dalembert_cr_gaas (dispersion-free far field)
  instruments   : notebook (multi-Gaussian detector model from the notebook)
                  aps_7idc (single Gaussian, 1.8 arcsec FWHM per the paper)

Outputs (results/matrix/):
  strain_<model>.npz          portable strain profiles
  xrd_<model>_<instrument>.npz  rocking-curve arrays (via xrd-strain-simulation)
  matrix_strain_far.png       near/far-field strain comparison figure
  matrix_rocking.png          rocking-curve comparison figure
  matrix_summary.json         numeric diagnostics

The XRD step is executed with the xrd-strain-simulation venv via subprocess,
mirroring scripts/validate_split.py.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
XRD_REPO = REPO.parent / "xrd-strain-simulation"
OUT = REPO / "results" / "matrix"

MODELS = ["ttm_cr_gaas", "ttm_dalembert_cr_gaas"]
INSTRUMENTS = ["notebook", "aps_7idc"]

# +/- ~200 arcsec around the GaAs (004) Bragg angle at 10 keV, sampled at
# ~0.7 arcsec so the 1.8 arcsec instrument Gaussian is resolved.
ANGLE_MIN = 25.975
ANGLE_MAX = 26.091
N_POINTS = 581

ARCSEC_PER_DEG = 3600.0


def python_in(repo: Path) -> str:
    venv = repo / ".venv" / "bin" / "python"
    return str(venv) if venv.exists() else sys.executable


def run_strain_models() -> dict[str, Path]:
    sys.path.insert(0, str(REPO / "src"))
    from strain_wave import get_preset, run_simulation, save_strain_profile

    profiles = {}
    for model in MODELS:
        print(f"\n=== strain model: {model} ===")
        config = get_preset("paper_fig3_gaas", model=model)
        result = run_simulation(config=config)
        path = OUT / f"strain_{model}.npz"
        save_strain_profile(result.to_profile(), path)
        profiles[model] = path
        print(f"saved {path}")
    return profiles


def run_xrd(profiles: dict[str, Path]) -> dict[tuple[str, str], Path]:
    xrd_python = python_in(XRD_REPO)
    outputs = {}
    for model, profile_path in profiles.items():
        for instrument in INSTRUMENTS:
            out = OUT / f"xrd_{model}_{instrument}.npz"
            fig = OUT / f"xrd_{model}_{instrument}.png"
            print(f"\n=== xrd: {model} / {instrument} ===")
            subprocess.run(
                [
                    xrd_python,
                    str(XRD_REPO / "scripts" / "run.py"),
                    "--strain-file", str(profile_path),
                    "--instrument", instrument,
                    "--angle-min", str(ANGLE_MIN),
                    "--angle-max", str(ANGLE_MAX),
                    "--n-points", str(N_POINTS),
                    "--output", str(fig),
                    "--save-arrays", str(out),
                    "--no-show",
                ],
                check=True,
                cwd=XRD_REPO,
            )
            outputs[(model, instrument)] = out
    return outputs


def dominant_far_field_wavelength_nm(z: np.ndarray, strain: np.ndarray) -> float:
    """Dominant oscillation wavelength in the 7.9-8.5 um window via FFT."""
    mask = (z >= 7.9e-6) & (z <= 8.5e-6)
    s = strain[mask] - strain[mask].mean()
    if not np.any(np.abs(s) > 1e-12):
        return float("nan")
    dz = z[1] - z[0]
    spec = np.abs(np.fft.rfft(s))
    freqs = np.fft.rfftfreq(len(s), d=dz)
    k = spec[1:].argmax() + 1
    return float(1.0 / freqs[k] * 1e9)


def make_strain_figure(profiles: dict[str, Path]) -> dict:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    data = {m: np.load(p, allow_pickle=True) for m, p in profiles.items()}
    labels = {
        "ttm_cr_gaas": "leapfrog (paper-faithful)",
        "ttm_dalembert_cr_gaas": "d'Alembert far field",
    }
    colors = {"ttm_cr_gaas": "tab:blue", "ttm_dalembert_cr_gaas": "tab:red"}

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)
    windows = [(0, 800), (7900, 8700)]
    titles = ["near field (film + interface)", "far field (~8 um, wavefront)"]
    for ax, (lo, hi), title in zip(axes, windows, titles):
        for m in MODELS:
            z_nm = data[m]["z"] * 1e9
            mask = (z_nm >= lo) & (z_nm <= hi)
            ax.plot(
                z_nm[mask],
                data[m]["strain"][mask] * 1e3,
                color=colors[m],
                lw=1.0,
                alpha=0.8,
                label=labels[m],
            )
        ax.set_xlabel("depth z (nm)")
        ax.set_title(title)
        ax.grid(alpha=0.3)
    axes[0].set_ylabel(r"strain $\times 10^{3}$")
    axes[0].legend()
    fig.suptitle(
        "Sci. Rep. Fig. 3 preset, t = 1.8 ns: leapfrog vs d'Alembert strain models"
    )
    fig.tight_layout()
    path = OUT / "matrix_strain_far.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"saved {path}")

    diag = {}
    for m in MODELS:
        z = data[m]["z"]
        strain = data[m]["strain"]
        diag[m] = {
            "far_field_dominant_wavelength_nm": dominant_far_field_wavelength_nm(
                z, strain
            ),
            "max_strain": float(strain.max()),
            "min_strain": float(strain.min()),
        }
    return diag


def make_rocking_figure(outputs: dict[tuple[str, str], Path]) -> dict:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = {
        "ttm_cr_gaas": "leapfrog strain",
        "ttm_dalembert_cr_gaas": "d'Alembert strain",
    }
    colors = {"ttm_cr_gaas": "tab:blue", "ttm_dalembert_cr_gaas": "tab:red"}

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)
    diag = {}
    for ax, instrument in zip(axes, INSTRUMENTS):
        curves = {}
        for m in MODELS:
            d = np.load(outputs[(m, instrument)], allow_pickle=True)
            th = d["angle_deg"]
            inten = d["intensity"]
            th_bragg = th[inten.argmax()]
            x = (th - th_bragg) * ARCSEC_PER_DEG
            ax.plot(x, inten, color=colors[m], lw=1.0, alpha=0.8, label=labels[m])
            curves[m] = inten
        ax.set_xlabel(r"$\theta - \theta_{Bragg}$ (arcsec)")
        ax.set_xlim(-200, 200)
        ax.set_title(f"instrument: {instrument}")
        ax.grid(alpha=0.3)
        diag[instrument] = {
            "max_abs_delta_log10I_between_models": float(
                np.max(np.abs(curves[MODELS[0]] - curves[MODELS[1]]))
            )
        }
    axes[0].set_ylabel(r"$\log_{10}$ intensity")
    axes[0].legend()
    fig.suptitle(
        "GaAs (004) rocking curves, Fig. 3 preset: strain model x instrument"
    )
    fig.tight_layout()
    path = OUT / "matrix_rocking.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"saved {path}")
    return diag


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    profiles = run_strain_models()
    outputs = run_xrd(profiles)
    strain_diag = make_strain_figure(profiles)
    xrd_diag = make_rocking_figure(outputs)

    summary = {
        "preset": "paper_fig3_gaas",
        "angle_window_deg": [ANGLE_MIN, ANGLE_MAX],
        "n_points": N_POINTS,
        "strain": strain_diag,
        "xrd": xrd_diag,
    }
    summary_path = OUT / "matrix_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))
    print(f"\nsaved {summary_path}")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
