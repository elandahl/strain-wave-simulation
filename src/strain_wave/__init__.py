"""Strain and elastic wave simulation."""

from strain_wave import models as _models  # noqa: F401 — registers default models
from strain_wave.config import SimulationConfig
from strain_wave.io import load_strain_profile, save_strain_profile
from strain_wave.models import StrainProfile, StrainSimulationResult
from strain_wave.pipeline import run_simulation
from strain_wave.presets import get_preset, list_presets

__all__ = [
    "SimulationConfig",
    "StrainProfile",
    "StrainSimulationResult",
    "get_preset",
    "list_presets",
    "load_strain_profile",
    "run_simulation",
    "save_strain_profile",
]
