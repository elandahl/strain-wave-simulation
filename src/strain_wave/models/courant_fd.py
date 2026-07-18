"""Courant-matched finite-difference propagation for a homogeneous substrate.

This module is deliberately independent of the Cr/GaAs thermal model.  It is
the acoustic field-solver foundation for future models with distributed
substrate sources, reflections, or spatially varying properties.

For the present source-free, homogeneous, right-going GaAs limit, the solver
is driven by a strain history at the left boundary and uses the acoustic
"magic time step" ``dt = dz / v`` (Courant number C = 1).  The second-order
leapfrog dispersion relation is then exact on the grid, so this finite-
difference solution must reproduce d'Alembert translation.  That identity is
used as an acceptance test.
"""

from __future__ import annotations

import numpy as np
from numba import njit


@njit
def _propagate_courant_one(
    boundary_samples: np.ndarray,
    n_space: int,
) -> np.ndarray:
    """Advance a boundary-driven 1D wave at Courant number exactly one."""
    previous = np.zeros(n_space)
    previous[0] = boundary_samples[0]

    if len(boundary_samples) == 1:
        return previous

    current = np.zeros(n_space)
    current[0] = boundary_samples[1]
    if n_space > 1:
        # Exact right-going initialization.  In normal use the simulation
        # starts just before t=0, so previous[0] is zero.
        current[1] = previous[0]

    for step in range(1, len(boundary_samples) - 1):
        following = np.zeros(n_space)
        following[0] = boundary_samples[step + 1]

        if n_space > 2:
            following[1:-1] = (
                2.0 * current[1:-1]
                - previous[1:-1]
                + current[2:]
                - 2.0 * current[1:-1]
                + current[:-2]
            )

        if n_space > 1:
            # First-order outgoing (Mur) boundary.  At C=1 it is exact:
            # the value at the last cell is the prior value one cell inward.
            following[-1] = current[-2]

        previous = current
        current = following

    return current


def propagate_right_going_fd(
    boundary_times: np.ndarray,
    boundary_strain: np.ndarray,
    *,
    n_space: int,
    dz: float,
    t_final: float,
    speed: float,
) -> tuple[np.ndarray, float, float, int]:
    """Propagate a boundary strain history through a homogeneous 1D domain.

    The acoustic clock starts up to one acoustic step before ``t=0`` so that
    it lands exactly on ``t_final`` while retaining ``dt = dz / speed`` and
    therefore Courant number one.  Boundary values before the first supplied
    time are zero.

    Returns
    -------
    strain_final, dt_acoustic, courant, n_steps
    """
    if n_space < 1:
        raise ValueError("n_space must be at least 1")
    if dz <= 0.0 or speed <= 0.0 or t_final < 0.0:
        raise ValueError("dz and speed must be positive; t_final must be nonnegative")
    if boundary_times.ndim != 1 or boundary_strain.ndim != 1:
        raise ValueError("boundary_times and boundary_strain must be one-dimensional")
    if len(boundary_times) != len(boundary_strain) or len(boundary_times) < 2:
        raise ValueError("boundary history arrays must have equal length >= 2")
    if np.any(np.diff(boundary_times) <= 0.0):
        raise ValueError("boundary_times must be strictly increasing")

    dt_acoustic = dz / speed
    n_steps = int(np.ceil(t_final / dt_acoustic))
    if n_steps == 0:
        value = np.interp(
            t_final,
            boundary_times,
            boundary_strain,
            left=0.0,
            right=boundary_strain[-1],
        )
        field = np.zeros(n_space)
        field[0] = value
        return field, dt_acoustic, 1.0, 0

    # Starting slightly before t=0 lets the final step land exactly at the
    # requested physical time without moving away from Courant number one.
    t_start = t_final - n_steps * dt_acoustic
    acoustic_times = t_start + np.arange(n_steps + 1) * dt_acoustic
    boundary_samples = np.interp(
        acoustic_times,
        boundary_times,
        boundary_strain,
        left=0.0,
        right=boundary_strain[-1],
    )

    field = _propagate_courant_one(boundary_samples, n_space)
    return field, dt_acoustic, 1.0, n_steps
