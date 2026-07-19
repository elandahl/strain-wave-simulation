"""Smoke tests."""

import numpy as np
import pytest

from strain_wave import SimulationConfig, get_preset, run_simulation
from strain_wave.materials import get_substrate


def test_default_model_is_dalembert():
    config = SimulationConfig(t_max=1e-12)
    result = run_simulation(config, verbose=False)

    assert result.model == "ttm_dalembert_cr_gaas"
    assert result.displacement.shape == result.z.shape
    assert result.strain.shape == result.z.shape
    assert result.n_iter > 0


def test_leapfrog_reference_model_still_runs():
    config = SimulationConfig(model="ttm_cr_gaas", t_max=1e-12)
    result = run_simulation(config, verbose=False)

    assert result.model == "ttm_cr_gaas"
    assert result.strain.shape == result.z.shape


def test_dalembert_model_runs_and_matches_near_field():
    leapfrog = run_simulation(
        SimulationConfig(model="ttm_cr_gaas", t_max=1e-12), verbose=False
    )
    dalembert = run_simulation(
        SimulationConfig(model="ttm_dalembert_cr_gaas", t_max=1e-12), verbose=False
    )

    assert dalembert.model == "ttm_dalembert_cr_gaas"
    assert dalembert.strain.shape == dalembert.z.shape

    # Film + near-interface strain must be identical to the leapfrog model;
    # only the far field beyond the monitor plane is reconstructed.
    monitor = leapfrog.n_bin_film - 1 + 10
    np.testing.assert_allclose(
        dalembert.strain[: monitor + 1], leapfrog.strain[: monitor + 1], atol=1e-18
    )


def test_paper_preset():
    config = get_preset("paper_fig3_gaas")
    assert config.model == "ttm_dalembert_cr_gaas"
    assert config.t_max == 1.8e-9
    assert config.L_film == 80e-9
    assert config.J == 80.0


def test_paper_fig2_si_preset_and_short_run():
    config = get_preset("paper_fig2_si")
    assert config.model == "ttm_dalembert_cr_si"
    assert config.substrate == "Si"
    assert config.t_max == 0.34e-9
    assert abs(config.R_ps - 1.0 / 1.1e8) < 1e-15

    short = SimulationConfig(
        model="ttm_dalembert_cr_si",
        substrate="Si",
        t_max=1e-12,
        L_film=80e-9,
        L_sub=500e-9,
        dz=2.67e-9,
    )
    result = run_simulation(short, verbose=False)
    assert result.substrate_material == "Si"
    assert result.model == "ttm_dalembert_cr_si"


@pytest.mark.parametrize(
    "substrate,model,expected_speed",
    [
        ("Ge", "ttm_dalembert_cr_ge", 4910.0),
        ("InSb", "ttm_dalembert_cr_insb", 3395.0),
    ],
)
def test_ge_insb_materials_and_short_runs(substrate, model, expected_speed):
    material = get_substrate(substrate)
    assert material.v == expected_speed
    assert material.rho > 5000.0
    assert material.cp > 0.0
    assert material.k_bulk > 0.0

    config = SimulationConfig(
        model=model,
        # Alias must enforce the intended substrate even when the caller
        # leaves the SimulationConfig default ("GaAs").
        t_max=1e-12,
        L_film=80e-9,
        L_sub=500e-9,
        dz=2.67e-9,
    )
    result = run_simulation(config, verbose=False)
    assert result.substrate_material == substrate
    assert result.model == model
