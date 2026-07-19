"""Material property definitions for film/substrate stacks."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class FilmMaterial:
    name: str


@dataclass(frozen=True)
class SubstrateMaterial:
    """Acoustic and thermal properties used by the TTM / wave solver."""

    name: str
    v: float  # longitudinal sound speed (m/s)
    rho: float  # mass density (kg/m^3)
    beta: float  # linear thermal expansion (1/K)
    cp: float  # specific heat (J kg^-1 K^-1)
    k_bulk: float  # bulk thermal conductivity (W m^-1 K^-1)


CHROMIUM = FilmMaterial("Cr")

GAAS = SubstrateMaterial(
    name="GaAs",
    v=4726.5,
    rho=5.32e3,
    beta=5.73e-6,
    cp=350.0,
    k_bulk=55.0,
)

# Si [001] longitudinal speed ~8430 m/s (standard literature); density,
# expansion, and heat capacity are room-temperature values used for the
# Cr/Si Sci. Rep. Fig. 2 forward model. Bulk k = 148 W/m/K matches the
# paper's Table 1 reference (effective k_s = 34 → factor 34/148).
SI = SubstrateMaterial(
    name="Si",
    v=8430.0,
    rho=2329.0,
    beta=2.6e-6,
    cp=700.0,
    k_bulk=148.0,
)

# Ge [001]: v_L = sqrt(C11/rho) with C11 ~= 128.5 GPa. Remaining values
# are room-temperature bulk properties; see docs/MATERIALS_PROVENANCE.md.
GE = SubstrateMaterial(
    name="Ge",
    v=4910.0,
    rho=5323.0,
    beta=5.9e-6,
    cp=320.0,
    k_bulk=60.2,
)

# InSb [001]: v_L = sqrt(C11/rho) with C11 ~= 66.6 GPa. InSb is much softer
# and less thermally conductive than the other campaign substrates.
INSB = SubstrateMaterial(
    name="InSb",
    v=3395.0,
    rho=5775.0,
    beta=5.37e-6,
    cp=199.0,
    k_bulk=18.0,
)

SUBSTRATES = {
    GAAS.name: GAAS,
    GE.name: GE,
    INSB.name: INSB,
    SI.name: SI,
}


def get_substrate(name: str) -> SubstrateMaterial:
    if name not in SUBSTRATES:
        available = ", ".join(sorted(SUBSTRATES))
        raise KeyError(f"Unknown substrate {name!r}. Available: {available}")
    return SUBSTRATES[name]
