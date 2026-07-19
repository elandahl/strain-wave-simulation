# Strain Wave Simulation

Compute displacement and strain profiles from laser-excited film/substrate systems. This repository is split from the published [thermo-elastic-gaas](https://github.com/elandahl/thermo-elastic-gaas) package and is intended for extension to new strain models and material stacks.

## Current models

- **`ttm_dalembert_cr_gaas` (default)** — TTM near field + exact d'Alembert far-field propagation (dispersion-free; recovers the discrete pulse train of Sci. Rep. Fig. 3). Physically realistic choice for new work.
- **`ttm_dalembert_cr_si`**, **`ttm_dalembert_cr_ge`**,
  **`ttm_dalembert_cr_insb`** / **`ttm_dalembert_cr_substrate`** — same
  method for the four APS substrates registered in `strain_wave.materials`.
- **`ttm_fd_courant_cr_gaas`** — Courant-matched FD far field (acceptance-tested vs d'Alembert).
- **`ttm_cr_gaas`** — historical notebook leapfrog (dispersive; archival only).
- **`ttm_fd_courant_cr_gaas`** — same TTM/near field, with a boundary-driven
  GaAs finite-difference field at acoustic Courant number exactly one. It
  matches d'Alembert in the homogeneous source-free limit and is the validated
  FD foundation for future distributed carrier stress, reflections, and
  inhomogeneous-substrate physics.
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

# Sci. Rep. Fig. 2 (Cr/Si, Δt = 0.34 ns)
python scripts/run.py --preset paper_fig2_si --no-show

# Ge/InSb material smoke runs (replace default sample parameters with the
# actual Cr thickness/fluence/TBC before data analysis):
python scripts/run.py --model ttm_dalembert_cr_ge --no-show
python scripts/run.py --model ttm_dalembert_cr_insb --no-show

# Historical notebook-faithful leapfrog solver:
python scripts/run.py --preset paper_fig3_gaas --model ttm_cr_gaas --no-show

# Validated Courant-one FD substrate field:
python scripts/run.py --preset paper_fig3_gaas --model ttm_fd_courant_cr_gaas --no-show
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
| Carrier generation / transport / recombination and deformation-potential stress | extend the field/source interface in `models/courant_fd.py`; register a new material model |
| Phonon MFP strain models | `src/strain_wave/models/` + register with `register_model()` |
| Additional substrates | `src/strain_wave/materials/`; generalized d'Alembert solver already dispatches by material |
| Multilayers / variable thickness | extend `SimulationConfig` and add new model class |
| Different films | new material modules + model implementation |

Ge and InSb room-temperature bulk-property provenance, plus the distinction
between intrinsic material properties and sample-specific fit parameters, is
documented in
[`docs/MATERIALS_PROVENANCE.md`](docs/MATERIALS_PROVENANCE.md).

## Paper reproduction and validation

The combined pipeline for the published paper remains frozen at [thermo-elastic-gaas](https://github.com/elandahl/thermo-elastic-gaas) (tag `paper-v1.0`). This repo reproduces the **strain step** only.

- `--preset paper_fig3_gaas` mirrors the Sci. Rep. Fig. 3 parameter set (`src/strain_wave/presets.py`).
- `scripts/validation_matrix.py` compares the historical leapfrog and default
  d'Alembert strain models through both XRD instrument models and produces
  comparison figures (`docs/ACOUSTIC_MODELS.md`,
  `docs/images/matrix_strain_far.png`).
- `scripts/validate_fd_courant.py` requires the Courant-one FD field to match
  d'Alembert for the full 1.8 ns paper preset. The checked-in run record is
  `docs/fd_courant_acceptance.json`.
- `scripts/validate_physics.py` runs the **physics acceptance suite** (acoustic
  echo spacing/decay, thickness scaling, lumped thermalization, and the
  numerical-dispersion convergence diagnosis). Report: `docs/physics_acceptance.json`;
  full write-up: [`docs/VALIDATION.md`](docs/VALIDATION.md). Both validation
  layers are also enforced by `pytest`.

## Physics roadmap

For physics *beyond* the current Fourier + constant-G Cr/GaAs baseline —
including the APS 7ID / Pohang multi-thickness, four-substrate (GaAs, Si,
InSb, Ge) re-analysis campaign — see [`docs/PHYSICS_ROADMAP.md`](docs/PHYSICS_ROADMAP.md).
