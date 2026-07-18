"""Named simulation presets for published and validation cases.

Mirrors ``thermo_elastic.presets`` in the archival thermo-elastic-gaas repo
(tag paper-v1.0) so paper-reproduction runs can be launched from this repo
with any registered strain model.
"""

from __future__ import annotations

from strain_wave.config import SimulationConfig

# Bulk Cr thermal conductivity used in the solver (W m^-1 K^-1).
_K_CR_BULK = 93.7
# Bulk GaAs thermal conductivity used in the solver (W m^-1 K^-1).
_K_GAAS_BULK = 55.0
# Bulk Si thermal conductivity (W m^-1 K^-1); paper Table 1 reference.
_K_SI_BULK = 148.0
# Bulk Cr longitudinal sound speed used as the solver base (m/s).
_V_CR_BASE = 6477.0


def paper_fig3_gaas(model: str = "ttm_dalembert_cr_gaas") -> SimulationConfig:
    """Parameters aimed at reproducing Sci. Rep. Fig. 3 (Cr/GaAs, Δt = 1.8 ns).

    Source: Jo et al., Sci. Rep. 12, 16606 (2022), DOI 10.1038/s41598-022-20715-5.

    Mapping notes
    -------------
    - Fluence 8 mJ/cm² → 80 J/m².
    - Paper σ_TBC (GaAs) = 0.5×10^8 W m^-2 K^-1 → R_ps = 1/σ_TBC = 2×10^-8.
    - Paper k_s (GaAs) = 5.1 → k_s_factor = 5.1 / 55.
    - Paper k_e,eff = 18.7 → k_e_factor ≈ 93.7 / 18.7 (room-temperature limit).
    - Paper v_Cr = 6608 m/s → cr_v_factor = 6608 / 6477.
    - Substrate depth is extended so the acoustic wave does not hit the
      absorbing boundary before 1.8 ns (v_GaAs ≈ 4726.5 m/s → ~8.5 μm).
    - Spatial step matches the paper FD grid (~2.67 nm).
    """
    sigma_tbc = 0.5e8
    k_s_paper = 5.1
    k_e_eff = 18.7
    v_cr_paper = 6608.0

    return SimulationConfig(
        model=model,
        substrate="GaAs",
        t_max=1.8e-9,
        L_film=80e-9,
        L_sub=10e-6,
        dz=2.67e-9,
        G=84e16,
        R_ps=1.0 / sigma_tbc,
        k_e_factor=_K_CR_BULK / k_e_eff,
        k_s_factor=k_s_paper / _K_GAAS_BULK,
        J=80.0,
        cr_density_factor=0.85,
        cr_v_factor=v_cr_paper / _V_CR_BASE,
    )


def paper_fig2_si(model: str = "ttm_dalembert_cr_si") -> SimulationConfig:
    """Parameters aimed at Sci. Rep. Fig. 2 (Cr/Si, Δt = 0.34 ns).

    Source: Jo et al., Sci. Rep. 12, 16606 (2022). Fig. 2 is the early-time
    acoustic / Brillouin-fringe figure used to extract G and k_e on Si.

    Mapping notes
    -------------
    - Same film drive as Fig. 3: 80 nm Cr, F = 8 mJ/cm², G = 84e16, k_e = 18.7.
    - Table 1 Si: σ_TBC = 1.1×10^8 → R_ps = 1/σ_TBC; k_s = 34 → factor 34/148.
    - Delay 0.34 ns; substrate depth ~4 μm (v_Si·t ≈ 2.9 μm plus margin).
    """
    sigma_tbc = 1.1e8
    k_s_paper = 34.0
    k_e_eff = 18.7
    v_cr_paper = 6608.0

    return SimulationConfig(
        model=model,
        substrate="Si",
        t_max=0.34e-9,
        L_film=80e-9,
        L_sub=4e-6,
        dz=2.67e-9,
        G=84e16,
        R_ps=1.0 / sigma_tbc,
        k_e_factor=_K_CR_BULK / k_e_eff,
        k_s_factor=k_s_paper / _K_SI_BULK,
        J=80.0,
        cr_density_factor=0.85,
        cr_v_factor=v_cr_paper / _V_CR_BASE,
    )


PRESETS = {
    "default": lambda model="ttm_dalembert_cr_gaas": SimulationConfig(model=model),
    "paper_fig3_gaas": paper_fig3_gaas,
    "paper_fig2_si": paper_fig2_si,
}


def get_preset(name: str, model: str | None = None) -> SimulationConfig:
    if name not in PRESETS:
        available = ", ".join(sorted(PRESETS))
        raise KeyError(f"Unknown preset {name!r}. Available: {available}")
    if model is None:
        return PRESETS[name]()
    return PRESETS[name](model=model)


def list_presets() -> list[str]:
    return sorted(PRESETS)
