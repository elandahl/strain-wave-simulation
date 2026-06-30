# Strain Wave Simulation

Compute displacement and strain profiles from laser-excited film/substrate systems. This repository is split from the published [thermo-elastic-gaas](https://github.com/elandahl/thermo-elastic-gaas) package and is intended for extension to new strain models and material stacks.

## Current model

- **`ttm_cr_gaas`** — Two-temperature model + elastic wave propagation for Cr on GaAs (paper reproduction)

## Quick start

```bash
git clone https://github.com/elandahl/strain-wave-simulation.git
cd strain-wave-simulation
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
python scripts/run.py --no-show
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

## Paper reproduction

The combined pipeline for the published paper remains frozen at [thermo-elastic-gaas](https://github.com/elandahl/thermo-elastic-gaas). This repo reproduces the **strain step** only.
