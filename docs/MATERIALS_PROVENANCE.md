# Substrate material properties

Room-temperature bulk properties used by the generalized
`ttm_dalembert_cr_substrate` model.

| substrate | `v_L[001]` (m/s) | density (kg/m³) | linear expansion (K⁻¹) | Cp (J/kg/K) | bulk k (W/m/K) |
|---|---:|---:|---:|---:|---:|
| GaAs | 4726.5 | 5320 | 5.73e-6 | 350 | 55 |
| Si | 8430 | 2329 | 2.60e-6 | 700 | 148 |
| Ge | 4910 | 5323 | 5.90e-6 | 320 | 60.2 |
| InSb | 3395 | 5775 | 5.37e-6 | 199 | 18.0 |

## Ge

- The [001] longitudinal speed is `sqrt(C11/rho)` using `C11 ≈ 128.5 GPa`,
  giving approximately 4.91 km/s.
- Density, thermal expansion, heat capacity, and conductivity are standard
  near-300 K bulk values.

## InSb

- The [001] longitudinal speed is `sqrt(C11/rho)` using `C11 ≈ 66.6 GPa`,
  giving approximately 3.40 km/s.
- InSb is substantially softer and less thermally conductive than the other
  substrates; this changes acoustic impedance, pulse spacing, and heat-flow
  timescales.

## What is deliberately not in the material registry

These are intrinsic substrate properties only. The following remain explicit
per-sample configuration or fit parameters:

- Cr film thickness and density/velocity corrections
- absorbed optical fluence
- Cr/substrate thermal boundary conductance (`R_ps`)
- effective near-surface substrate conductivity (`k_s_factor`)
- pump/probe spot averaging and instrument response

The published Sci. Rep. paper supplies fitted interface/effective-transport
values for its GaAs and Si examples, but not a universal Cr/Ge or Cr/InSb
value. Therefore this repository does **not** create misleading “paper”
presets for Ge/InSb. Campaign presets should be added only when the sample
metadata and fitting assumptions are available.
