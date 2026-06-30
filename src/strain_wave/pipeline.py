"""End-to-end strain simulation pipeline."""

from strain_wave.config import SimulationConfig
from strain_wave.io import save_strain_profile
from strain_wave.models.base import get_model
from strain_wave.models.result import StrainProfile, StrainSimulationResult


def run_simulation(
    config: SimulationConfig | None = None,
    verbose: bool = True,
) -> StrainSimulationResult:
    config = config or SimulationConfig()
    model = get_model(config.model)
    return model.run(config, verbose=verbose)


def run_and_save_profile(
    output_path: str,
    config: SimulationConfig | None = None,
    verbose: bool = True,
) -> StrainProfile:
    result = run_simulation(config=config, verbose=verbose)
    profile = result.to_profile()
    save_strain_profile(profile, output_path)
    return profile
