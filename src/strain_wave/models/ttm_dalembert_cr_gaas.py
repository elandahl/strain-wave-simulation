"""TTM Cr/GaAs model with d'Alembert far-field reconstruction.

Motivation
----------
The leapfrog elastic solver in ``ttm_cr_gaas`` runs at an acoustic Courant
number C = v*dt/dz ~ 0.003 (the time step is set by thermal-diffusion
stability, not acoustics). At small C the leapfrog scheme is numerically
dispersive: short-wavelength components travel slower than v, so a sharp
strain pulse propagating several micrometres develops an oscillatory wake
(~35 nm wavelength for the Fig. 3 parameters) that is a numerical artifact,
not physics. See docs/ACOUSTIC_MODELS.md for the full analysis.

Method
------
This model exploits the fact that dispersion accumulates with *propagation
distance*, while everything near the film/substrate interface (pulse
generation, film round trips over ~80 nm, interface transmission) is computed
accurately by the same trusted finite-difference solver.

1. Run the unmodified TTM + leapfrog solver, recording the strain and
   substrate temperature time-histories at a monitor plane a few grid cells
   inside the substrate.
2. Split the recorded strain into a quasi-static thermal part,
   eta_th(t) = 3*beta_GaAs*(T_s(z_m, t) - T0)  (stress-free equilibrium,
   since B = rho*v^2 in this model), and a purely right-going acoustic part
   eta_ac(t) = eta(z_m, t) - eta_th(t).
3. The substrate beyond the monitor is homogeneous with an absorbing far
   boundary, so the acoustic field there is exactly a right-going d'Alembert
   solution: eta_ac(z, t_max) = eta_ac(z_m, t_max - (z - z_m)/v_GaAs).
   This translation is exact — no numerical dispersion.
4. Reassemble: solver strain up to the monitor plane, then
   eta_ac (translated) + 3*beta*(T_s(z, t_max) - T0) beyond it. Displacement
   beyond the monitor is rebuilt by integrating the strain.

The film, interface, and near-field substrate strain are bit-identical to
``ttm_cr_gaas``; only the far field (z > monitor plane) is replaced by the
dispersion-free reconstruction.
"""

import numpy as np

from strain_wave.config import SimulationConfig
from strain_wave.materials import CHROMIUM, GAAS
from strain_wave.models.base import register_model
from strain_wave.models.result import StrainSimulationResult
from strain_wave.models.ttm_cr_gaas_solver import (
    V_GAAS,
    simulation_diagnostics,
    solver,
)

BETA_GAAS = 5.73e-6  # linear thermal expansion coefficient used by the solver
T0 = 300.0

# Monitor plane offset from the interface, in grid cells. Close enough that
# numerical dispersion has not accumulated, far enough to sit cleanly inside
# the substrate grid.
MONITOR_OFFSET_CELLS = 10


class TtmDalembertCrGaAsModel:
    """TTM + leapfrog near field, d'Alembert-reconstructed far field."""

    name = "ttm_dalembert_cr_gaas"

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

        # --- d'Alembert far-field reconstruction -------------------------
        # Record times: step i stores the state after update i+1.
        t_grid = (np.arange(n_iter) + 1) * dt
        t_end = n_iter * dt

        # Right-going acoustic component at the monitor plane.
        ac_rec = strain_rec - 3.0 * BETA_GAAS * (ts_rec - T0)

        z_m = z[monitor]
        z_far = z[monitor + 1 :]
        # Arrival-time mapping: eta_ac(z, t_end) = eta_ac(z_m, t_end - (z-z_m)/v)
        t_retarded = t_end - (z_far - z_m) / V_GAAS
        eta_ac_far = np.interp(t_retarded, t_grid, ac_rec, left=0.0)
        eta_th_far = 3.0 * BETA_GAAS * (t_s[monitor + 1 :] - T0)

        strain = strain.copy()
        strain[monitor + 1 :] = eta_ac_far + eta_th_far

        # Rebuild far-field displacement consistently with the new strain.
        displacement = displacement.copy()
        displacement[monitor + 1 :] = displacement[monitor] + np.cumsum(
            strain[monitor + 1 :]
        ) * config.dz
        # -----------------------------------------------------------------

        elastic_wave_reach_nm = simulation_diagnostics(
            config.dz, config.t_max, n_iter, dt
        )

        if verbose:
            print("Info")
            print("model =", self.name)
            print("t_max = ", config.t_max)
            print("n_iter =", n_iter)
            print("Elastic wave can reach to.. ", elastic_wave_reach_nm, "(nm)")
            print(
                "d'Alembert monitor plane at z =",
                round(z_m * 1e9, 2),
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


register_model(TtmDalembertCrGaAsModel())
