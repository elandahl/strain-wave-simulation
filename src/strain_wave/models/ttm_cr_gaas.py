"""Cr/GaAs two-temperature model with elastic wave propagation."""

import numpy as np

from strain_wave.config import SimulationConfig
from strain_wave.materials import CHROMIUM, GAAS
from strain_wave.models.base import register_model
from strain_wave.models.result import StrainSimulationResult
from strain_wave.models.ttm_cr_gaas_solver import simulation_diagnostics, solver


class TtmCrGaAsModel:
    """Paper-reproduction TTM + elastic wave model for Cr on GaAs."""

    name = "ttm_cr_gaas"

    def run(
        self, config: SimulationConfig, verbose: bool = True
    ) -> StrainSimulationResult:
        z = np.linspace(0, config.L_tot, config.n_total)

        displacement, t_e, t_p, t_s, n_iter, dt = solver(
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
        )

        elastic_wave_reach_nm = simulation_diagnostics(
            config.dz, config.t_max, n_iter, dt
        )

        if verbose:
            print("Info")
            print("model =", self.name)
            print("t_max = ", config.t_max)
            print("n_iter =", n_iter)
            print("Elastic wave can reach to.. ", elastic_wave_reach_nm, "(nm)")

        strain = np.gradient(displacement, z)

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


register_model(TtmCrGaAsModel())
