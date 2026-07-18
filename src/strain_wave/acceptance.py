"""Physics acceptance diagnostics for the Cr/GaAs strain models.

These are *model-independent, physics-based* checks that can be run before any
comparison to experimental data. They are meant to catch regressions and to
back up the claims made in docs/ACOUSTIC_MODELS.md with quantitative,
reproducible numbers:

1. Analytic acoustic pulse-train diagnostics
   - echo spacing set by the film round-trip time,
     Δz = 2 * L_film * v_GaAs / v_Cr;
   - echo-to-echo amplitude decay set by the Cr/GaAs acoustic impedance
     mismatch, r = (Z_Cr - Z_GaAs) / (Z_Cr + Z_GaAs);
   - echo spacing scales linearly with film thickness.

2. Lumped two-temperature thermalization sanity
   - electron-phonon equilibration happens on a picosecond timescale and the
     time constant scales as ~C_e/G (i.e. inversely with the coupling G).

3. Numerical-dispersion convergence
   - the difference between the raw leapfrog far field and the dispersion-free
     d'Alembert far field grows monotonically with propagation distance,
     confirming the artifact is propagation-accumulated numerical dispersion
     rather than physics.

The heavy lifting (running the solver) lives in the model registry; this module
only provides light post-processing and analytic reference values, plus a
``run_acceptance_suite`` entry point that packages everything into a structured,
serialisable report. ``scripts/validate_physics.py`` and
``tests/test_physics_acceptance.py`` both build on the functions here.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field

import numpy as np

from strain_wave.config import SimulationConfig
from strain_wave.models.ttm_cr_gaas_solver import V_GAAS, source
from strain_wave.models.ttm_dalembert_cr_gaas import (
    BETA_GAAS,
    MONITOR_OFFSET_CELLS,
    T0,
)
from strain_wave.pipeline import run_simulation
from strain_wave.presets import paper_fig3_gaas

# --- Material constants (mirrors the hard-coded values in the solver) -------
# Bulk Cr density and longitudinal sound speed used as the solver base.
_RHO_CR_BASE = 7190.0  # kg m^-3
_V_CR_BASE = 6477.0  # m s^-1
# Bulk GaAs.
_RHO_GAAS = 5320.0  # kg m^-3


def _paper_config(model: str, **overrides) -> SimulationConfig:
    """Return the paper Fig. 3 config with the given model and field overrides.

    Used to build fast, small-domain variants of the paper case for the
    acceptance checks (shorter ``t_max`` / ``L_sub`` keep runtimes to seconds).
    """
    cfg = paper_fig3_gaas(model=model)
    for key, value in overrides.items():
        setattr(cfg, key, value)
    return cfg


def impedance_reflection(cfg: SimulationConfig) -> float:
    """Cr->GaAs acoustic pressure reflection coefficient r.

    r = (Z_Cr - Z_GaAs) / (Z_Cr + Z_GaAs) with Z = rho * v. The magnitude of r
    sets the fraction of the strain pulse retained at each film round trip and
    therefore the echo-to-echo amplitude decay.
    """
    z_cr = (_RHO_CR_BASE * cfg.cr_density_factor) * (_V_CR_BASE * cfg.cr_v_factor)
    z_gaas = _RHO_GAAS * V_GAAS
    return (z_cr - z_gaas) / (z_cr + z_gaas)


def expected_echo_spacing_nm(cfg: SimulationConfig) -> float:
    """Spatial spacing (nm) between successive acoustic echoes in the substrate.

    A pulse makes a film round trip in time 2*L_film/v_Cr; in that time the
    leading front has advanced 2*L_film*(v_GaAs/v_Cr) into the substrate.
    """
    v_cr = _V_CR_BASE * cfg.cr_v_factor
    return 2.0 * cfg.L_film * (V_GAAS / v_cr) * 1e9


def _acoustic_far_field(result, cfg: SimulationConfig):
    """Return (z_nm, acoustic_strain) beyond the d'Alembert monitor plane.

    Subtracts the quasi-static thermal strain 3*beta*(T_s - T0) so that only the
    propagating acoustic component remains.
    """
    monitor = cfg.n_bin_film - 1 + MONITOR_OFFSET_CELLS
    acoustic = result.strain - 3.0 * BETA_GAAS * (result.T_s - T0)
    z_nm = result.z * 1e9
    return z_nm, acoustic, monitor


def measure_echo_spacing_nm(result, cfg: SimulationConfig) -> float:
    """Measure echo spacing (nm) from the acoustic far field via autocorrelation.

    Autocorrelation is robust to the bipolar shape of each echo and to the
    thermal background, unlike naive peak picking.
    """
    z_nm, acoustic, monitor = _acoustic_far_field(result, cfg)
    front_nm = cfg.t_max * V_GAAS * 1e9
    dz_nm = (result.z[1] - result.z[0]) * 1e9

    idx = np.arange(len(z_nm))
    sel = (idx > monitor) & (z_nm > z_nm[monitor] + 20.0) & (z_nm < front_nm - 20.0)
    signal = acoustic[sel]
    if signal.size < 8:
        return float("nan")
    signal = signal - signal.mean()

    corr = np.correlate(signal, signal, mode="full")[signal.size - 1 :]
    for i in range(2, corr.size - 1):
        if corr[i] > corr[i - 1] and corr[i] > corr[i + 1] and corr[i] > 0.2 * corr[0]:
            return i * dz_nm
    return float("nan")


def measure_first_echo_ratio(result, cfg: SimulationConfig) -> float:
    """Amplitude ratio of the first trailing echo to the leading acoustic pulse.

    Windows the acoustic far field backward from the leading front in steps of
    the measured echo spacing and takes the peak |strain| in each window. The
    ratio of the second window (first echo) to the first window (leading pulse)
    should approach the impedance reflection magnitude |r|. Later windows are
    contaminated by the thermal tail near the interface and are not used.
    """
    z_nm, acoustic, _ = _acoustic_far_field(result, cfg)
    spacing = measure_echo_spacing_nm(result, cfg)
    if not np.isfinite(spacing):
        return float("nan")
    dz_nm = (result.z[1] - result.z[0]) * 1e9
    spacing_cells = int(round(spacing / dz_nm))
    front_nm = cfg.t_max * V_GAAS * 1e9

    i = int((front_nm - 30.0) / dz_nm)
    envelope = []
    for _ in range(2):
        lo = max(0, i - spacing_cells)
        segment = np.abs(acoustic[lo:i])
        if segment.size == 0:
            break
        envelope.append(float(segment.max()))
        i -= spacing_cells
    if len(envelope) < 2 or envelope[0] == 0.0:
        return float("nan")
    return envelope[1] / envelope[0]


@dataclass
class ThermalizationResult:
    peak_time_ps: float
    peak_delta_T: float
    tau_1e_ps: float


def lumped_thermalization(
    G: float, J: float, cr_density_factor: float = 0.85
) -> ThermalizationResult:
    """0-D two-temperature electron-phonon thermalization at the film surface.

    Integrates
        C_e dT_e/dt = -G (T_e - T_p) + S(t),   C_e = gamma_Cr * T_e
        C_p dT_p/dt = +G (T_e - T_p)
    with the same laser source, gamma_Cr and C_p as the full solver, ignoring
    diffusion (lumped/0-D). Returns the time of peak (T_e - T_p), the peak value,
    and the 1/e decay time of (T_e - T_p) after the peak. This isolates the film
    heat channel from the acoustics; tau_1e should be picosecond-scale and scale
    as ~C_e/G (inversely with G).
    """
    R_cr = 0.6
    knu_cr = 4 * np.pi / 800e-9 * 3.135
    tp = 50e-15
    rho_cr = _RHO_CR_BASE * cr_density_factor
    c_p = rho_cr * 460.0
    gamma = rho_cr / 51.9961e-3 * 1.4e-3

    dt = 1e-16
    n = int(15e-12 / dt)
    t_e = t_p = T0
    times = np.empty(n)
    diff = np.empty(n)
    for i in range(n):
        t = i * dt
        c_e = gamma * t_e
        s = source(R_cr, knu_cr, J, tp, 0.0, t)
        d_te = (-G * (t_e - t_p) + s) / c_e
        d_tp = (G * (t_e - t_p)) / c_p
        t_e += d_te * dt
        t_p += d_tp * dt
        times[i] = t
        diff[i] = t_e - t_p

    i_peak = int(diff.argmax())
    peak = float(diff[i_peak])
    post = diff[i_peak:]
    below = np.where(post < peak / np.e)[0]
    tau = float(times[i_peak + below[0]] - times[i_peak]) if below.size else float("nan")
    return ThermalizationResult(
        peak_time_ps=float(times[i_peak] * 1e12),
        peak_delta_T=peak,
        tau_1e_ps=tau * 1e12,
    )


@dataclass
class Check:
    """One acceptance check with a machine-readable pass/fail verdict."""

    name: str
    passed: bool
    detail: str
    values: dict = field(default_factory=dict)


@dataclass
class AcceptanceReport:
    checks: list[Check]

    @property
    def passed(self) -> bool:
        return all(c.passed for c in self.checks)

    def to_dict(self) -> dict:
        return {
            "passed": self.passed,
            "n_passed": sum(c.passed for c in self.checks),
            "n_total": len(self.checks),
            "checks": [asdict(c) for c in self.checks],
        }


def check_pulse_train(model: str = "ttm_dalembert_cr_gaas") -> list[Check]:
    """Echo spacing, impedance-set decay, and thickness scaling checks."""
    checks: list[Check] = []

    cfg80 = _paper_config(model, t_max=200e-12, L_film=80e-9, L_sub=1.2e-6)
    res80 = run_simulation(cfg80, verbose=False)
    spacing = measure_echo_spacing_nm(res80, cfg80)
    expected = expected_echo_spacing_nm(cfg80)
    rel = abs(spacing - expected) / expected
    checks.append(
        Check(
            name="echo_spacing",
            passed=rel < 0.12,
            detail=(
                f"measured {spacing:.1f} nm vs analytic "
                f"2*L_film*v_GaAs/v_Cr = {expected:.1f} nm ({rel * 100:.1f}% off)"
            ),
            values={"measured_nm": spacing, "expected_nm": expected, "rel_error": rel},
        )
    )

    r = abs(impedance_reflection(cfg80))
    ratio = measure_first_echo_ratio(res80, cfg80)
    checks.append(
        Check(
            name="echo_reflectivity",
            passed=0.5 * r < ratio < 2.0 * r,
            detail=(
                f"first echo/leading amplitude ratio {ratio:.3f} vs "
                f"impedance |r| = {r:.3f}"
            ),
            values={"measured_ratio": ratio, "impedance_r": r},
        )
    )

    cfg120 = _paper_config(model, t_max=200e-12, L_film=120e-9, L_sub=1.2e-6)
    res120 = run_simulation(cfg120, verbose=False)
    spacing120 = measure_echo_spacing_nm(res120, cfg120)
    scale = spacing120 / spacing if spacing else float("nan")
    expected_scale = 120.0 / 80.0
    rel_scale = abs(scale - expected_scale) / expected_scale
    checks.append(
        Check(
            name="thickness_scaling",
            passed=rel_scale < 0.12,
            detail=(
                f"spacing(120 nm)/spacing(80 nm) = {scale:.3f} vs "
                f"thickness ratio {expected_scale:.3f} ({rel_scale * 100:.1f}% off)"
            ),
            values={
                "spacing_120_nm": spacing120,
                "spacing_80_nm": spacing,
                "measured_scale": scale,
                "expected_scale": expected_scale,
            },
        )
    )
    return checks


def check_thermalization() -> list[Check]:
    """Picosecond timescale and ~1/G scaling of electron-phonon coupling."""
    cfg = paper_fig3_gaas()
    base = lumped_thermalization(cfg.G, cfg.J, cfg.cr_density_factor)
    checks = [
        Check(
            name="thermalization_timescale",
            passed=0.1 < base.tau_1e_ps < 3.0,
            detail=(
                f"lumped electron-phonon 1/e time {base.tau_1e_ps:.3f} ps "
                f"(peak ΔT {base.peak_delta_T:.0f} K at {base.peak_time_ps:.3f} ps)"
            ),
            values=asdict(base),
        )
    ]

    half = lumped_thermalization(cfg.G * 0.5, cfg.J, cfg.cr_density_factor)
    double = lumped_thermalization(cfg.G * 2.0, cfg.J, cfg.cr_density_factor)
    monotonic = half.tau_1e_ps > base.tau_1e_ps > double.tau_1e_ps
    checks.append(
        Check(
            name="thermalization_scales_with_G",
            passed=monotonic,
            detail=(
                f"tau(G/2)={half.tau_1e_ps:.3f} > tau(G)={base.tau_1e_ps:.3f} "
                f"> tau(2G)={double.tau_1e_ps:.3f} ps"
            ),
            values={
                "tau_half_G": half.tau_1e_ps,
                "tau_G": base.tau_1e_ps,
                "tau_double_G": double.tau_1e_ps,
            },
        )
    )
    return checks


def check_dispersion_convergence(
    t_max_values=(80e-12, 160e-12, 320e-12), l_sub: float = 1.8e-6
) -> list[Check]:
    """Leapfrog-vs-d'Alembert far-field error grows with propagation distance."""
    fronts = []
    errors = []
    for t_max in t_max_values:
        lf_cfg = _paper_config("ttm_cr_gaas", t_max=t_max, L_sub=l_sub)
        da_cfg = _paper_config("ttm_dalembert_cr_gaas", t_max=t_max, L_sub=l_sub)
        lf = run_simulation(lf_cfg, verbose=False)
        da = run_simulation(da_cfg, verbose=False)

        monitor = da_cfg.n_bin_film - 1 + MONITOR_OFFSET_CELLS
        z_nm = lf.z * 1e9
        front_nm = t_max * V_GAAS * 1e9
        idx = np.arange(len(z_nm))
        sel = (idx > monitor) & (z_nm < front_nm - 20.0)
        err = np.sqrt(np.mean((lf.strain[sel] - da.strain[sel]) ** 2))
        scale = np.sqrt(np.mean(da.strain[sel] ** 2))
        fronts.append(front_nm)
        errors.append(err / scale if scale else float("nan"))

    monotonic = all(errors[i] < errors[i + 1] for i in range(len(errors) - 1))
    return [
        Check(
            name="dispersion_grows_with_distance",
            passed=monotonic,
            detail=(
                "relative far-field RMS(leapfrog - d'Alembert) increases with "
                "propagation distance: "
                + ", ".join(
                    f"{f:.0f} nm -> {e:.3f}" for f, e in zip(fronts, errors)
                )
            ),
            values={"fronts_nm": fronts, "rel_rms": errors},
        )
    ]


def run_acceptance_suite(model: str = "ttm_dalembert_cr_gaas") -> AcceptanceReport:
    """Run the full physics acceptance suite and return a structured report."""
    checks: list[Check] = []
    checks.extend(check_pulse_train(model))
    checks.extend(check_thermalization())
    checks.extend(check_dispersion_convergence())
    return AcceptanceReport(checks=checks)
