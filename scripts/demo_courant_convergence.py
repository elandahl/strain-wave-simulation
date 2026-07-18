#!/usr/bin/env python3
"""Demonstrate how the acoustic step size (Courant number) drives agreement.

This is an *analysis / figure* script, not a registered model. It takes the
real recorded acoustic boundary history from the Cr/GaAs paper-preset run and
propagates it through the *same* second-order leapfrog scheme at a sweep of
acoustic Courant numbers ``C = v * dt_acoustic / dz``:

- at the original solver's value (C ~ 0.003, set by thermal stability) the
  finite-difference wave develops the spurious dispersive wake;
- as ``dt_acoustic`` is increased toward the acoustic "magic step"
  ``dt = dz / v`` (C = 1) the wake vanishes and the field converges to the
  exact d'Alembert translation.

That convergence is exactly why ``ttm_fd_courant_cr_gaas`` (which runs at
C = 1) agrees with ``ttm_dalembert_cr_gaas``.

Outputs ``results/courant_convergence.png`` and a checked-in copy under
``docs/images/``.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from numba import njit  # noqa: E402

from strain_wave import get_preset  # noqa: E402
from strain_wave.models.ttm_cr_gaas_solver import V_GAAS, solver  # noqa: E402
from strain_wave.models.ttm_dalembert_cr_gaas import (  # noqa: E402
    BETA_GAAS,
    MONITOR_OFFSET_CELLS,
    T0,
)

REPO = Path(__file__).resolve().parents[1]

# Courant numbers below and immediately approaching the stability boundary.
# The native full-grid solver runs near the smallest value; C = 1 is the
# magic step used by ttm_fd_courant_cr_gaas.
STABLE_COURANTS = [
    0.003,
    0.01,
    0.03,
    0.1,
    0.3,
    0.5,
    0.7,
    0.85,
    0.9,
    0.95,
    0.98,
    0.99,
    0.995,
    0.999,
    1.0,
]

# Values just above and farther above the CFL boundary. These are formally
# unstable, so they are assessed by time-to-divergence rather than presented
# as physically meaningful final strain profiles.
ABOVE_CFL_COURANTS = [1.00001, 1.0001, 1.001, 1.005, 1.01, 1.05, 1.1, 1.3]

# Stop an unstable run once its field is this many times larger than the
# actual boundary pulse. This avoids floating-point overflow and gives a
# useful, comparable time-to-divergence metric.
DIVERGENCE_FACTOR = 100.0
FAR_WINDOW_NM = (7900.0, 8700.0)


@njit
def _propagate_leapfrog(
    boundary_samples: np.ndarray,
    n_space: int,
    c2: float,
    stop_amplitude: float,
) -> tuple[np.ndarray, int, bool]:
    """Boundary-driven leapfrog, optionally stopping an unstable run."""
    previous = np.zeros(n_space)
    current = np.zeros(n_space)
    current[0] = boundary_samples[1] if len(boundary_samples) > 1 else boundary_samples[0]

    for step in range(1, len(boundary_samples) - 1):
        following = np.zeros(n_space)
        following[0] = boundary_samples[step + 1]
        if n_space > 2:
            following[1:-1] = (
                2.0 * current[1:-1]
                - previous[1:-1]
                + c2 * (current[2:] - 2.0 * current[1:-1] + current[:-2])
            )
        if n_space > 1:
            # First-order outgoing (Mur) boundary at Courant number sqrt(c2).
            following[-1] = current[-1] + np.sqrt(c2) * (current[-2] - current[-1])
        if np.max(np.abs(following)) >= stop_amplitude:
            return following, step + 1, True
        previous = current
        current = following
    return current, len(boundary_samples) - 1, False


def _propagate_at_courant(
    boundary_times: np.ndarray,
    boundary_vals: np.ndarray,
    n_space: int,
    dz: float,
    speed: float,
    t_final: float,
    courant: float,
    stop_amplitude: float,
) -> tuple[np.ndarray, int, int, float, bool]:
    dt_ac = courant * dz / speed
    n_steps = int(np.ceil(t_final / dt_ac))
    t_start = t_final - n_steps * dt_ac
    acoustic_times = t_start + np.arange(n_steps + 1) * dt_ac
    samples = np.interp(
        acoustic_times, boundary_times, boundary_vals, left=0.0, right=boundary_vals[-1]
    )
    field, completed_steps, diverged = _propagate_leapfrog(
        samples, n_space, courant * courant, stop_amplitude
    )
    elapsed = completed_steps * dt_ac
    return field, completed_steps, n_steps, elapsed, diverged


def main() -> None:
    config = get_preset("paper_fig3_gaas")
    z = np.linspace(0, config.L_tot, config.n_total)
    monitor = config.n_bin_film - 1 + MONITOR_OFFSET_CELLS

    print("Running Cr/GaAs paper-preset solver to record boundary history...")
    _, _, _, t_s, n_iter, dt, strain_rec, ts_rec = solver(
        config.G,
        config.R_ps,
        config.k_e_factor,
        config.k_s_factor,
        config.J,
        config.cr_density_factor,
        config.cr_v_factor,
        config.L_film,
        config.L_sub,
        config.dz,
        config.t_max,
        monitor,
    )

    t_end = n_iter * dt
    t_grid = (np.arange(n_iter) + 1) * dt
    acoustic_record = strain_rec - 3.0 * BETA_GAAS * (ts_rec - T0)
    boundary_times = np.concatenate(([0.0], t_grid))
    boundary_vals = np.concatenate(([0.0], acoustic_record))

    dz_grid = z[monitor + 1] - z[monitor]
    n_space = len(z) - monitor
    z_far = z[monitor:] * 1e9  # nm, aligned with propagated field index 0 = monitor
    thermal_far = 3.0 * BETA_GAAS * (t_s[monitor:] - T0)

    # Exact d'Alembert translation (the C -> 1 limit / continuum reference).
    t_retarded = t_end - (z[monitor:] - z[monitor]) / V_GAAS
    dalembert = np.interp(t_retarded, t_grid, acoustic_record, left=0.0) + thermal_far

    print("Sweeping Courant numbers...")
    fields = {}
    errors = {}
    mask = (z_far >= FAR_WINDOW_NM[0]) & (z_far <= FAR_WINDOW_NM[1])
    ref_rms = float(np.sqrt(np.mean((dalembert[mask] - thermal_far[mask]) ** 2)))
    for c in STABLE_COURANTS:
        field, _, _, _, diverged = _propagate_at_courant(
            boundary_times,
            boundary_vals,
            n_space,
            dz_grid,
            V_GAAS,
            t_end,
            c,
            stop_amplitude=1e300,
        )
        if diverged:
            raise RuntimeError(f"Unexpected divergence for stable C={c}")
        total = field + thermal_far
        fields[c] = total
        acoustic_only = field
        dalembert_acoustic = dalembert - thermal_far
        errors[c] = float(
            np.sqrt(np.mean((acoustic_only[mask] - dalembert_acoustic[mask]) ** 2))
            / ref_rms
        )
        print(f"  C = {c:6.3f} : normalized far-field RMS error = {errors[c]:.3e}")

    print("Testing values above the CFL boundary...")
    stop_amplitude = DIVERGENCE_FACTOR * float(np.max(np.abs(boundary_vals)))
    above_cfl = {}
    for c in ABOVE_CFL_COURANTS:
        _, completed, total_steps, elapsed, diverged = _propagate_at_courant(
            boundary_times,
            boundary_vals,
            n_space,
            dz_grid,
            V_GAAS,
            t_end,
            c,
            stop_amplitude=stop_amplitude,
        )
        above_cfl[c] = {
            "completed": completed,
            "total": total_steps,
            "elapsed_ps": elapsed * 1e12,
            "diverged": diverged,
        }
        status = (
            f"reached {DIVERGENCE_FACTOR:g}x amplitude"
            if diverged
            else "did not reach threshold within 1.8 ns"
        )
        print(f"  C = {c:8.5f} : {status} at {elapsed * 1e12:.2f} ps")

    # ---- Figure -----------------------------------------------------------
    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(18, 5))

    # d'Alembert reference drawn thick underneath everything.
    axA.plot(
        z_far[mask],
        dalembert[mask] * 1e3,
        color="0.6",
        lw=4.0,
        alpha=0.9,
        label="d'Alembert (exact)",
        zorder=1,
    )
    profile_cs = [0.003, 0.3, 0.7, 0.9, 0.99, 0.999]
    cmap = plt.cm.winter(np.linspace(0.0, 1.0, len(profile_cs)))
    for color, c in zip(cmap, profile_cs):
        axA.plot(
            z_far[mask],
            fields[c][mask] * 1e3,
            color=color,
            lw=1.0,
            alpha=0.85,
            label=f"C = {c:g}",
            zorder=2,
        )
    # The magic step C=1 lands on the reference; draw it bold on top.
    axA.plot(
        z_far[mask],
        fields[1.0][mask] * 1e3,
        color="crimson",
        lw=1.6,
        ls="--",
        label="C = 1 (magic step, used by model)",
        zorder=3,
    )
    axA.set_xlabel("depth z (nm)")
    axA.set_ylabel(r"strain $\times 10^{3}$")
    axA.set_title(
        "Same recorded boundary wave, different acoustic step size\n"
        "wake contracts only very near C=1; C=1 lands on exact"
    )
    axA.legend(fontsize=8, loc="upper left")
    axA.grid(alpha=0.3)

    cs = np.array(STABLE_COURANTS)
    errs = np.array([errors[c] for c in STABLE_COURANTS])
    axB.semilogy(cs, errs, "o-", color="tab:red", lw=1.4, ms=4, zorder=3)
    axB.axvline(0.003, color="0.6", ls=":", lw=1.2)
    axB.annotate(
        "original solver (C ≈ 0.003):\ntiny step, full wake",
        xy=(0.003, errs[0]),
        xytext=(0.0045, errs[0] * 0.06),
        fontsize=8,
        arrowprops=dict(arrowstyle="->", color="0.4"),
    )
    axB.annotate(
        "ttm_fd_courant (C = 1):\nexact, ~1e-12",
        xy=(1.0, errs[-1]),
        xytext=(0.58, errs[-1] * 300),
        fontsize=8,
        arrowprops=dict(arrowstyle="->", color="0.4"),
    )
    axB.text(
        0.05,
        1.6,
        "small-C sinc limit:\nshrinking the step does NOT help",
        fontsize=8,
        color="0.35",
    )
    axB.set_xlim(-0.02, 1.02)
    axB.set_xlabel("acoustic Courant number  C = v·dt/dz")
    axB.set_ylabel("normalized far-field RMS error vs d'Alembert")
    axB.set_title("Stable side: dense sweep just below C = 1")
    axB.grid(alpha=0.3, which="both")

    above_cs = np.array(ABOVE_CFL_COURANTS)
    elapsed_ps = np.array([above_cfl[c]["elapsed_ps"] for c in above_cs])
    diverged = np.array([above_cfl[c]["diverged"] for c in above_cs])
    axC.plot(above_cs[diverged], elapsed_ps[diverged], "rx-", lw=1.2, label="reached 100× amplitude")
    if np.any(~diverged):
        axC.plot(
            above_cs[~diverged],
            elapsed_ps[~diverged],
            "^",
            mec="tab:blue",
            mfc="none",
            ms=7,
            label="survived 1.8 ns (still formally unstable)",
        )
    axC.axhline(t_end * 1e12, color="0.6", ls=":", lw=1.0)
    axC.axvline(1.0, color="0.3", ls="--", lw=1.0)
    axC.set_xlim(0.995, max(ABOVE_CFL_COURANTS) + 0.015)
    axC.set_xlabel("Courant number C above the CFL limit")
    axC.set_ylabel("time to 100× amplitude (ps)")
    axC.set_title("Above C = 1: numerical instability")
    axC.grid(alpha=0.3)
    axC.legend(fontsize=8)

    fig.suptitle(
        "Cr/GaAs Fig. 3 preset (Δt = 1.8 ns): below, at, and above the acoustic CFL boundary",
        fontsize=12,
    )
    fig.tight_layout()

    out = REPO / "results" / "courant_convergence.png"
    ref = REPO / "docs" / "images" / "courant_convergence.png"
    for path in (out, ref):
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\nSaved {out}")
    print(f"Saved reference copy {ref}")


if __name__ == "__main__":
    main()
