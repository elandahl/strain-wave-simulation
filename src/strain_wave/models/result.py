"""Result containers and interchange format."""

from dataclasses import dataclass, field

import numpy as np

FORMAT_VERSION = "1"


@dataclass
class StrainSimulationResult:
    """Output of a strain/displacement simulation."""

    z: np.ndarray
    displacement: np.ndarray
    strain: np.ndarray
    T_e: np.ndarray
    T_p: np.ndarray
    T_s: np.ndarray
    n_iter: int
    elastic_wave_reach_nm: float
    model: str
    film_material: str
    substrate_material: str
    dz: float
    n_bin_film: int
    L_film: float
    L_sub: float
    t_max: float

    def to_profile(self) -> "StrainProfile":
        return StrainProfile.from_simulation_result(self)


@dataclass
class StrainProfile:
    """Portable strain profile for downstream XRD or other analysis."""

    format_version: str
    model: str
    film_material: str
    substrate_material: str
    z: np.ndarray
    displacement: np.ndarray
    strain: np.ndarray
    substrate_strain: np.ndarray
    dz: float
    n_bin_film: int
    L_film: float
    L_sub: float
    t_max: float
    metadata: dict = field(default_factory=dict)

    @classmethod
    def from_simulation_result(cls, result: StrainSimulationResult) -> "StrainProfile":
        return cls(
            format_version=FORMAT_VERSION,
            model=result.model,
            film_material=result.film_material,
            substrate_material=result.substrate_material,
            z=result.z,
            displacement=result.displacement,
            strain=result.strain,
            substrate_strain=result.strain[result.n_bin_film :].copy(),
            dz=result.dz,
            n_bin_film=result.n_bin_film,
            L_film=result.L_film,
            L_sub=result.L_sub,
            t_max=result.t_max,
            metadata={
                "n_iter": result.n_iter,
                "elastic_wave_reach_nm": result.elastic_wave_reach_nm,
            },
        )
