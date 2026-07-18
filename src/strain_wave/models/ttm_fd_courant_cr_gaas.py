"""TTM Cr/GaAs model with Courant-matched FD substrate propagation.

This model is the finite-difference counterpart to
``ttm_dalembert_cr_gaas``.  Both use the same trusted TTM + historical
leapfrog calculation to generate the film/interface response and record the
outgoing acoustic strain at a monitor plane.  They differ only beyond that
plane:

* d'Alembert translates the recorded waveform analytically;
* this model advances a boundary-driven GaAs finite-difference field with
  ``dt_acoustic = dz / v_GaAs`` (Courant number one).

For a homogeneous, source-free GaAs substrate these must agree.  That is the
built-in acceptance limit for this solver.  Unlike direct translation, the FD
field is also the numerical foundation that can later be extended with
distributed carrier/deformation-potential sources, buried interfaces,
reflections, damping, or nonlinear constitutive terms.

The historical ``ttm_cr_gaas`` model remains unchanged for notebook
reproducibility.  This new model does not claim to fix its Cr-film near field;
it fixes and validates substrate propagation, where the long-distance
dispersion artifact accumulated.
"""

import numpy as np

from strain_wave.config import SimulationConfig
from strain_wave.materials import CHROMIUM, GAAS
from strain_wave.models.base import register_model
from strain_wave.models.courant_fd import propagate_right_going_fd
from strain_wave.models.result import StrainSimulationResult
from strain_wave.models.ttm_cr_gaas_solver import (
    V_GAAS,
    simulation_diagnostics,
    solver,
)
from strain_wave.models.ttm_dalembert_cr_gaas import (
    BETA_GAAS,
    MONITOR_OFFSET_CELLS,
    T0,
)


class TtmFdCourantCrGaAsModel:
    """TTM near field plus Courant-one FD propagation in homogeneous GaAs."""

    name = "ttm_fd_courant_cr_gaas"

    def run(
        self, config: SimulationConfig, verbose: bool = True
    ) -> StrainSimulationResult:
        z = np.linspace(0, config.L_tot, config.n_total)
        monitor = config.n_bin_film - 1 + MONITOR_OFFSET_CELLS

        displacement, t_e, t_p, t_s, n_iter, dt, strain_rec, ts_rec = solver(
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

        strain = np.gradient(displacement, z)
        t_grid = (np.arange(n_iter) + 1) * dt
        t_end = n_iter * dt

        # Separate the outgoing acoustic component from the local quasistatic
        # thermal strain at the monitor plane, exactly as in the d'Alembert
        # reference model.
        acoustic_record = strain_rec - 3.0 * BETA_GAAS * (ts_rec - T0)
        boundary_times = np.concatenate(([0.0], t_grid))
        boundary_strain = np.concatenate(([0.0], acoustic_record))

        # Use the actual output-grid spacing.  np.linspace differs minutely
        # from config.dz whenever L_tot/dz is not an integer.
        dz_acoustic = z[monitor + 1] - z[monitor]
        acoustic_field, dt_acoustic, courant, n_acoustic_steps = (
            propagate_right_going_fd(
                boundary_times,
                boundary_strain,
                n_space=len(z) - monitor,
                dz=dz_acoustic,
                t_final=t_end,
                speed=V_GAAS,
            )
        )

        thermal_far = 3.0 * BETA_GAAS * (t_s[monitor + 1 :] - T0)
        strain = strain.copy()
        strain[monitor + 1 :] = acoustic_field[1:] + thermal_far

        displacement = displacement.copy()
        displacement[monitor + 1 :] = displacement[monitor] + np.cumsum(
            strain[monitor + 1 :]
        ) * config.dz

        elastic_wave_reach_nm = simulation_diagnostics(
            config.dz, config.t_max, n_iter, dt
        )

        if verbose:
            print("Info")
            print("model =", self.name)
            print("t_max = ", config.t_max)
            print("thermal n_iter =", n_iter)
            print("acoustic n_iter =", n_acoustic_steps)
            print("acoustic dt =", dt_acoustic)
            print("acoustic Courant number =", courant)
            print("Elastic wave can reach to.. ", elastic_wave_reach_nm, "(nm)")
            print(
                "FD boundary/monitor plane at z =",
                round(z[monitor] * 1e9, 2),
                "nm (interface +",
                MONITOR_OFFSET_CELLS,
                "cells)",
            )

        return StrainSimulationResult(
            z=z,
            displacement=displacement,
            strain=strain,
            T_e=t_e,
            T_p=t_p,
            T_s=t_s,
            n_iter=n_iter,
            elastic_wave_reach_nm=elastic_wave_reach_nm,
            model=self.name,
            film_material=CHROMIUM.name,
            substrate_material=GAAS.name,
            dz=config.dz,
            n_bin_film=config.n_bin_film,
            L_film=config.L_film,
            L_sub=config.L_sub,
            t_max=config.t_max,
        )


register_model(TtmFdCourantCrGaAsModel())
