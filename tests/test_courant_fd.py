"""Acceptance tests for the Courant-matched acoustic field solver."""

import numpy as np

from strain_wave import SimulationConfig, run_simulation
from strain_wave.models.courant_fd import propagate_right_going_fd


def test_courant_one_fd_equals_dalembert_translation():
    """A source-free boundary wave must translate exactly at C=1."""
    speed = 5000.0
    dz = 2.5e-9
    dt = dz / speed
    t_final = 40 * dt

    # A smooth bipolar acoustic pulse, sampled much more finely than the
    # acoustic FD clock to exercise boundary interpolation too.
    boundary_times = np.linspace(0.0, t_final, 1001)
    center = 8.0 * dt
    width = 1.5 * dt
    x = (boundary_times - center) / width
    boundary_strain = x * np.exp(-0.5 * x**2)
    boundary_strain[0] = 0.0

    field, dt_acoustic, courant, n_steps = propagate_right_going_fd(
        boundary_times,
        boundary_strain,
        n_space=80,
        dz=dz,
        t_final=t_final,
        speed=speed,
    )

    retarded_times = t_final - np.arange(80) * dz / speed
    expected = np.interp(
        retarded_times,
        boundary_times,
        boundary_strain,
        left=0.0,
        right=boundary_strain[-1],
    )

    assert courant == 1.0
    assert dt_acoustic == dt
    assert n_steps == 40
    np.testing.assert_allclose(field, expected, atol=2e-14, rtol=2e-14)


def test_cr_gaas_fd_matches_dalembert_in_source_free_substrate():
    """Integrated Cr/GaAs acceptance limit for the new FD backend."""
    common = dict(
        t_max=40e-12,
        L_film=80e-9,
        L_sub=400e-9,
        dz=2e-9,
    )
    dalembert = run_simulation(
        SimulationConfig(model="ttm_dalembert_cr_gaas", **common),
        verbose=False,
    )
    fd = run_simulation(
        SimulationConfig(model="ttm_fd_courant_cr_gaas", **common),
        verbose=False,
    )

    monitor = dalembert.n_bin_film - 1 + 10
    np.testing.assert_array_equal(
        fd.strain[: monitor + 1],
        dalembert.strain[: monitor + 1],
    )
    np.testing.assert_allclose(
        fd.strain[monitor + 1 :],
        dalembert.strain[monitor + 1 :],
        atol=2e-12,
        rtol=2e-8,
    )
