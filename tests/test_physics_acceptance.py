"""Physics acceptance tests for the Cr/GaAs strain models.

These wrap the diagnostics in ``strain_wave.acceptance`` so the physical claims
in docs/ACOUSTIC_MODELS.md and docs/VALIDATION.md are enforced as regressions.
Each test asserts on a physically meaningful quantity (echo spacing, impedance
reflectivity, thickness scaling, thermalization timescale, dispersion growth),
not on a frozen numerical snapshot.
"""

import numpy as np

from strain_wave.acceptance import (
    check_dispersion_convergence,
    check_pulse_train,
    check_thermalization,
    expected_echo_spacing_nm,
    impedance_reflection,
    lumped_thermalization,
    measure_echo_spacing_nm,
    measure_first_echo_ratio,
    _paper_config,
)
from strain_wave import run_simulation


def test_echo_spacing_matches_film_round_trip():
    cfg = _paper_config(
        "ttm_dalembert_cr_gaas", t_max=200e-12, L_film=80e-9, L_sub=1.2e-6
    )
    res = run_simulation(cfg, verbose=False)
    measured = measure_echo_spacing_nm(res, cfg)
    expected = expected_echo_spacing_nm(cfg)
    assert abs(measured - expected) / expected < 0.12


def test_first_echo_decays_like_impedance_reflectivity():
    cfg = _paper_config(
        "ttm_dalembert_cr_gaas", t_max=200e-12, L_film=80e-9, L_sub=1.2e-6
    )
    res = run_simulation(cfg, verbose=False)
    ratio = measure_first_echo_ratio(res, cfg)
    r = abs(impedance_reflection(cfg))
    assert 0.5 * r < ratio < 2.0 * r


def test_pulse_train_suite_passes():
    checks = check_pulse_train("ttm_dalembert_cr_gaas")
    failed = [c.name for c in checks if not c.passed]
    assert not failed, f"failed pulse-train checks: {failed}"


def test_thermalization_is_picosecond_and_scales_inversely_with_G():
    cfg = _paper_config("ttm_dalembert_cr_gaas")
    base = lumped_thermalization(cfg.G, cfg.J, cfg.cr_density_factor)
    assert 0.1 < base.tau_1e_ps < 3.0

    half = lumped_thermalization(cfg.G * 0.5, cfg.J, cfg.cr_density_factor)
    double = lumped_thermalization(cfg.G * 2.0, cfg.J, cfg.cr_density_factor)
    assert half.tau_1e_ps > base.tau_1e_ps > double.tau_1e_ps


def test_thermalization_suite_passes():
    checks = check_thermalization()
    failed = [c.name for c in checks if not c.passed]
    assert not failed, f"failed thermalization checks: {failed}"


def test_numerical_dispersion_grows_with_propagation_distance():
    # Smaller/shorter than the script default to keep the test fast.
    checks = check_dispersion_convergence(
        t_max_values=(60e-12, 180e-12), l_sub=1.2e-6
    )
    assert all(c.passed for c in checks), checks[0].detail
