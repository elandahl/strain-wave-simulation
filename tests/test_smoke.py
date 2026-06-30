"""Smoke tests."""

from strain_wave import SimulationConfig, run_simulation


def test_short_simulation_runs():
    config = SimulationConfig(t_max=1e-12)
    result = run_simulation(config, verbose=False)

    assert result.model == "ttm_cr_gaas"
    assert result.displacement.shape == result.z.shape
    assert result.strain.shape == result.z.shape
    assert result.n_iter > 0
