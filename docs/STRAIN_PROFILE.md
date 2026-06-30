# Strain Profile Interchange Format (v1)

Portable `.npz` files produced by `strain-wave-simulation` for consumption by XRD or other tools.

## Required fields

| Key | Type | Description |
|-----|------|-------------|
| `format_version` | str | `"1"` |
| `model` | str | e.g. `"ttm_cr_gaas"` |
| `film_material` | str | e.g. `"Cr"` |
| `substrate_material` | str | e.g. `"GaAs"` |
| `z` | float array | depth grid (m) |
| `displacement` | float array | displacement (m) |
| `strain` | float array | strain (dimensionless) |
| `substrate_strain` | float array | strain in substrate only (`strain[n_bin_film:]`) |
| `dz` | float | grid step (m) |
| `n_bin_film` | int | number of film grid bins |
| `L_film` | float | film thickness (m) |
| `L_sub` | float | substrate thickness (m) |
| `t_max` | float | simulation end time (s) |
| `metadata_json` | str | JSON dict (e.g. `n_iter`, `elastic_wave_reach_nm`) |

## Python API

```python
from strain_wave.io import load_strain_profile, save_strain_profile
```
