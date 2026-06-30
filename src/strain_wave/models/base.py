"""Base interface for strain simulation models."""

from __future__ import annotations

from typing import Protocol

from strain_wave.config import SimulationConfig
from strain_wave.models.result import StrainSimulationResult


class StrainModel(Protocol):
    """Protocol for strain/displacement simulation backends."""

    name: str

    def run(self, config: SimulationConfig, verbose: bool = True) -> StrainSimulationResult:
        ...


_MODELS: dict[str, StrainModel] = {}


def register_model(model: StrainModel) -> StrainModel:
    _MODELS[model.name] = model
    return model


def get_model(name: str) -> StrainModel:
    if name not in _MODELS:
        available = ", ".join(sorted(_MODELS)) or "(none registered)"
        raise KeyError(f"Unknown strain model {name!r}. Available: {available}")
    return _MODELS[name]


def list_models() -> list[str]:
    return sorted(_MODELS)
