"""TTM Cr/substrate model with d'Alembert far-field reconstruction.

Generalizes ``ttm_dalembert_cr_gaas`` to any substrate registered in
``strain_wave.materials``. GaAs behavior is preserved when
``config.substrate == "GaAs"`` (the default).
"""

from __future__ import annotations

import numpy as np

from strain_wave.config import SimulationConfig
from strain_wave.materials import CHROMIUM, get_substrate
from strain_wave.models.base import register_model
from strain_wave.models.result import StrainSimulationResult
from strain_wave.models.ttm_cr_gaas_solver import simulation_diagnostics, solver

T0 = 300.0
MONITOR_OFFSET_CELLS = 10


class TtmDalembertCrSubstrateModel:
    """TTM + leapfrog near field, d'Alembert far field, substrate from config."""

    name = "ttm_dalembert_cr_substrate"

    def run(
        self, config: SimulationConfig, verbose: bool = True
    ) -> StrainSimulationResult:
        sub = get_substrate(config.substrate)
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
            sub.v,
            sub.rho,
            sub.beta,
            sub.cp,
            sub.k_bulk,
        )

        strain = np.gradient(displacement, z)

        t_grid = (np.arange(n_iter) + 1) * dt
        t_end = n_iter * dt
        ac_rec = strain_rec - 3.0 * sub.beta * (ts_rec - T0)

        z_m = z[monitor]
        z_far = z[monitor + 1 :]
        t_retarded = t_end - (z_far - z_m) / sub.v
        eta_ac_far = np.interp(t_retarded, t_grid, ac_rec, left=0.0)
        eta_th_far = 3.0 * sub.beta * (t_s[monitor + 1 :] - T0)

        strain = strain.copy()
        strain[monitor + 1 :] = eta_ac_far + eta_th_far

        displacement = displacement.copy()
        displacement[monitor + 1 :] = displacement[monitor] + np.cumsum(
            strain[monitor + 1 :]
        ) * config.dz

        elastic_wave_reach_nm = simulation_diagnostics(
            config.dz, config.t_max, n_iter, dt, sub.v
        )

        if verbose:
            print("Info")
            print("model =", self.name)
            print("substrate =", sub.name)
            print("t_max = ", config.t_max)
            print("n_iter =", n_iter)
            print("Elastic wave can reach to.. ", elastic_wave_reach_nm, "(nm)")
            print(
                "d'Alembert monitor plane at z =",
                round(z_m * 1e9, 2),
                "nm",
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
            substrate_material=sub.name,
            dz=config.dz,
            n_bin_film=config.n_bin_film,
            L_film=config.L_film,
            L_sub=config.L_sub,
            t_max=config.t_max,
        )


# Convenience alias registered under a Si-specific name for presets/CLI.
class TtmDalembertCrSiModel(TtmDalembertCrSubstrateModel):
    name = "ttm_dalembert_cr_si"

    def run(
        self, config: SimulationConfig, verbose: bool = True
    ) -> StrainSimulationResult:
        if config.substrate != "Si":
            config = SimulationConfig(**{**config.__dict__, "substrate": "Si"})
        return super().run(config, verbose=verbose)


class TtmDalembertCrGeModel(TtmDalembertCrSubstrateModel):
    """Convenience alias for Cr/Ge campaign runs."""

    name = "ttm_dalembert_cr_ge"

    def run(
        self, config: SimulationConfig, verbose: bool = True
    ) -> StrainSimulationResult:
        if config.substrate != "Ge":
            config = SimulationConfig(**{**config.__dict__, "substrate": "Ge"})
        return super().run(config, verbose=verbose)


class TtmDalembertCrInSbModel(TtmDalembertCrSubstrateModel):
    """Convenience alias for Cr/InSb campaign runs."""

    name = "ttm_dalembert_cr_insb"

    def run(
        self, config: SimulationConfig, verbose: bool = True
    ) -> StrainSimulationResult:
        if config.substrate != "InSb":
            config = SimulationConfig(**{**config.__dict__, "substrate": "InSb"})
        return super().run(config, verbose=verbose)


register_model(TtmDalembertCrSubstrateModel())
register_model(TtmDalembertCrSiModel())
register_model(TtmDalembertCrGeModel())
register_model(TtmDalembertCrInSbModel())
