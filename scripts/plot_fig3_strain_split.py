#!/usr/bin/env python3
"""Plot the d'Alembert strain profile on Sci. Rep. Fig. 3's split z-axis.

This produces a direct overlay against the published Fig. 3 strain panel:
the corrected (dispersion-free) ``ttm_dalembert_cr_gaas`` strain drawn on the
paper's exact broken depth axis (0-800 nm and 7900-8700 nm). The paper-faithful
leapfrog ``ttm_cr_gaas`` strain is shown as a faint reference so the removal of
the numerical-dispersion wake is visible.

Reads the saved profiles from ``results/matrix/`` if present (created by
``scripts/validation_matrix.py``); otherwise runs the two paper-preset
simulations itself (~2-3 min).

Output: ``results/fig3_strain_split_dalembert.png`` and a checked-in reference
copy under ``docs/images/``.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

REPO = Path(__file__).resolve().parents[1]

# Depth windows (nm) matching the published split strain axis.
NEAR_WINDOW_NM = (0.0, 800.0)
FAR_WINDOW_NM = (7900.0, 8700.0)


def _load_or_simulate(model: str) -> tuple[np.ndarray, np.ndarray]:
    npz = REPO / "results" / "matrix" / f"strain_{model}.npz"
    if npz.exists():
        data = np.load(npz, allow_pickle=True)
        return data["z"] * 1e9, data["strain"]

    import sys

    sys.path.insert(0, str(REPO / "src"))
    from strain_wave import get_preset, run_simulation

    config = get_preset("paper_fig3_gaas", model=model)
    result = run_simulation(config=config)
    return result.z * 1e9, result.strain


def main() -> None:
    z_da, strain_da = _load_or_simulate("ttm_dalembert_cr_gaas")
    z_lf, strain_lf = _load_or_simulate("ttm_cr_gaas")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for ax, (lo, hi) in zip(axes, (NEAR_WINDOW_NM, FAR_WINDOW_NM)):
        mask_lf = (z_lf >= lo) & (z_lf <= hi)
        ax.plot(
            z_lf[mask_lf],
            strain_lf[mask_lf],
            color="0.7",
            lw=0.9,
            label="leapfrog (paper-faithful solver, dispersive)",
        )
        mask_da = (z_da >= lo) & (z_da <= hi)
        ax.plot(
            z_da[mask_da],
            strain_da[mask_da],
            "-k",
            lw=1.1,
            label="d'Alembert (dispersion-free far field)",
        )
        ax.set_xlim(lo, hi)
        ax.set_xlabel("z (nm)")
        ax.axhline(0.0, color="0.85", lw=0.8, zorder=0)

    axes[0].set_ylabel("strain")
    axes[0].set_title("near field (film + interface)", fontsize=10)
    axes[1].set_title("far field (acoustic wavefront at 1.8 ns)", fontsize=10)
    axes[1].legend(fontsize=8, loc="upper left")
    fig.suptitle(
        "Cr/GaAs strain vs depth, Sci. Rep. Fig. 3 preset (Δt = 1.8 ns) — "
        "paper split z-axis",
        fontsize=11,
    )
    fig.tight_layout()

    out = REPO / "results" / "fig3_strain_split_dalembert.png"
    ref = REPO / "docs" / "images" / "fig3_strain_split_dalembert.png"
    for path in (out, ref):
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved {out}")
    print(f"Saved reference copy {ref}")


if __name__ == "__main__":
    main()
