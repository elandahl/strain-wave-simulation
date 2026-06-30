"""Material property definitions.

Additional substrates (Si, Ge, InSb) and film materials can be added here
as new dataclasses without changing the model registry interface.
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class FilmMaterial:
    name: str


@dataclass(frozen=True)
class SubstrateMaterial:
    name: str


CHROMIUM = FilmMaterial("Cr")
GAAS = SubstrateMaterial("GaAs")
