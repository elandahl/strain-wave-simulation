"""Simulation configuration."""

from dataclasses import dataclass


@dataclass
class SimulationConfig:
    """Parameters for the Cr-on-semiconductor two-temperature strain simulation."""

    # Default model: d'Alembert far field (physically realistic; dispersion-
    # free). "ttm_cr_gaas" is the historical notebook-faithful leapfrog
    # reference. See docs/ACOUSTIC_MODELS.md.
    model: str = "ttm_dalembert_cr_gaas"
    # Substrate key into strain_wave.materials.SUBSTRATES ("GaAs", "Si", ...).
    # Historical GaAs models ignore this and always use GaAs props; the
    # generalized d'Alembert model honors it.
    substrate: str = "GaAs"
    t_max: float = 300e-12
    L_film: float = 180e-9
    L_sub: float = 1800e-9
    dz: float = 1e-9

    G: float = 42e16
    R_ps: float = 5e-8
    k_e_factor: float = 10.0
    k_s_factor: float = 0.4
    J: float = 50.0
    cr_density_factor: float = 0.85
    cr_v_factor: float = 1.0

    @property
    def L_tot(self) -> float:
        return self.L_film + self.L_sub

    @property
    def n_total(self) -> int:
        return int(self.L_tot / self.dz)

    @property
    def n_bin_film(self) -> int:
        return int(self.L_film / self.dz)
