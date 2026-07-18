# Strain Wave Simulation

Compute displacement and strain profiles from laser-excited film/substrate systems. This repository is split from the published [thermo-elastic-gaas](https://github.com/elandahl/thermo-elastic-gaas) package and is intended for extension to new strain models and material stacks.

## Current models

- **`ttm_dalembert_cr_gaas` (default)** — TTM near field + exact d'Alembert far-field propagation (dispersion-free; recovers the discrete pulse train of Sci. Rep. Fig. 3). Physically realistic choice for new work.
- **`ttm_cr_gaas`** — Two-temperature model + leapfrog elastic wave propagation (notebook/paper-faithful). Kept as the **historical reference**: its far field carries a known numerical-dispersion wake at this solver's tiny acoustic Courant number. Used by `scripts/validate_split.py` for bit-for-bit equivalence with the frozen thermo-elastic-gaas repo.

See `docs/ACOUSTIC_MODELS.md` for the numerical-dispersion analysis and when to use which.

## Quick start

```bash
git clone https://github.com/elandahl/strain-wave-simulation.git
cd strain-wave-simulation
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
python scripts/run.py --no-show

# Paper Fig. 3 parameters (d'Alembert far field is the default):
python scripts/run.py --preset paper_fig3_gaas --no-show

# Historical notebook-faithful leapfrog solver:
python scripts/run.py --preset paper_fig3_gaas --model ttm_cr_gaas --no-show
```

Outputs:

- `results/strain_profile.npz` — portable strain profile for downstream XRD
- `results/strain_figure.png` — displacement and strain plots

## Downstream use

Feed the strain profile to [xrd-strain-simulation](https://github.com/elandahl/xrd-strain-simulation):

```bash
python ../xrd-strain-simulation/scripts/run.py \
  --strain-file results/strain_profile.npz --no-show
```

## Strain profile format (v1)

The `.npz` file includes `z`, `displacement`, `strain`, `substrate_strain`, `dz`, `n_bin_film`, material names, model id, and metadata. See `docs/STRAIN_PROFILE.md`.

## Extension points

| Future work | Where to add |
|-------------|--------------|
| Phonon MFP strain models | `src/strain_wave/models/` + register with `register_model()` |
| Si, Ge, InSb substrates | `src/strain_wave/materials/` + new model or generalized solver |
| Multilayers / variable thickness | extend `SimulationConfig` and add new model class |
| Different films | new material modules + model implementation |

## Paper reproduction and validation

The combined pipeline for the published paper remains frozen at [thermo-elastic-gaas](https://github.com/elandahl/thermo-elastic-gaas) (tag `paper-v1.0`). This repo reproduces the **strain step** only.

- `--preset paper_fig3_gaas` mirrors the Sci. Rep. Fig. 3 parameter set (`src/strain_wave/presets.py`).
- `scripts/validation_matrix.py` runs both strain models through both XRD instrument models and produces comparison figures (`docs/ACOUSTIC_MODELS.md`, `docs/images/matrix_strain_far.png`).
