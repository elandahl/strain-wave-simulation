"""Smoke tests."""

import numpy as np

from strain_wave import SimulationConfig, get_preset, run_simulation


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
