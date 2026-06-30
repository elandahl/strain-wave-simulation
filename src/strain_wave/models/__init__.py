"""Strain simulation result types and model registry."""

from strain_wave.models import ttm_cr_gaas as _ttm_cr_gaas  # noqa: F401
from strain_wave.models.base import StrainModel, get_model, list_models, register_model
from strain_wave.models.result import StrainProfile, StrainSimulationResult

__all__ = [
    "StrainModel",
    "StrainProfile",
    "StrainSimulationResult",
    "get_model",
    "list_models",
    "register_model",
]
