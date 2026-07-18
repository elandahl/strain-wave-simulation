# Acoustic models

This repo supports switchable strain models (see `strain_wave.models.base`).
Two are currently registered; both use the *same* two-temperature (TTM)
thermal solve and differ only in how the elastic wave reaches the far field.

| model | far field | use when |
|---|---|---|
| `ttm_cr_gaas` | leapfrog FD everywhere | reproducing the notebook / thermo-elastic-gaas (tag `paper-v1.0`) bit-for-bit |
| `ttm_dalembert_cr_gaas` | exact d'Alembert translation beyond a monitor plane | physically meaningful far-field strain (e.g. Sci. Rep. Fig. 3 pulse train) |

Select with `--model` on the CLI, or `SimulationConfig(model=...)`.

```bash
python scripts/run.py --preset paper_fig3_gaas --model ttm_dalembert_cr_gaas --no-show
```

## Why a second model exists: numerical dispersion at tiny Courant number

The original solver uses a single time step for both heat diffusion and the
elastic wave, set by thermal stability:

- `dt = dz^2 / 2 / 0.002 ≈ 1.78 fs` (for `dz = 2.67 nm`)
- acoustic Courant number `C = v_GaAs * dt / dz ≈ 0.003`

The second-order leapfrog scheme has numerical phase velocity

    v_num / v = sin(k dz / 2) / (k dz / 2)   (limit C -> 0)

so short-wavelength components travel *slower* than `v`. The error is tiny
per cell but accumulates with propagation distance. For the Fig. 3 preset the
strain pulse travels ~8.5 um (~3200 cells); components with ~50 nm wavelength
lag the front by hundreds of nm, producing a chirped oscillatory wake behind
the wavefront. That wake is a numerical artifact: the published Fig. 3
(Jo et al., Sci. Rep. 12, 16606 (2022)) instead shows a discrete pulse train
with ~114 nm spacing (film round-trip time 2*80 nm / v_Cr times v_GaAs) and
amplitudes decaying by the interface reflection coefficient
`r = (Z_Cr - Z_GaAs)/(Z_Cr + Z_GaAs) ≈ 0.23` per round trip.

## How `ttm_dalembert_cr_gaas` works

Dispersion accumulates with distance, but everything near the interface —
pulse generation in the 80 nm film, film round trips, transmission into the
substrate — happens over tens of cells and is computed accurately by the
trusted FD solver. The substrate beyond the interface is homogeneous with an
absorbing far boundary, so the acoustic field there is purely right-going and
obeys d'Alembert's solution exactly: `eta(z, t) = f(t - z/v)`.

The model therefore:

1. runs the unmodified TTM + leapfrog solver, recording strain and substrate
   temperature time-histories at a monitor plane 10 cells (~27 nm) inside the
   substrate;
2. subtracts the quasi-static thermal strain `3*beta_GaAs*(T_s - T0)` from the
   record to isolate the right-going acoustic wave;
3. maps the acoustic record to depth exactly:
   `eta_ac(z, t_max) = eta_ac(z_monitor, t_max - (z - z_monitor)/v_GaAs)`;
4. adds back the local thermal strain and rebuilds displacement by
   integrating strain.

Film, interface, and near-field strain are bit-identical to `ttm_cr_gaas`
(enforced by `tests/test_smoke.py`); only `z > z_monitor` is replaced.

Remaining approximations: the sharpest edges of high-order film echoes still
carry some dispersion from their round trips inside the film itself
(n-th echo travels 2*n*80 nm of FD grid before transmission), but echo
amplitudes decay like `0.23^n`, so only the first few matter.

## Validation matrix result (paper_fig3_gaas preset)

`scripts/validation_matrix.py` runs both models and pushes both strain
profiles through xrd-strain-simulation with both instrument models.
Figures: `docs/images/matrix_strain_far.png` (here) and
`docs/images/matrix_rocking.png` (in xrd-strain-simulation).

![strain comparison](images/matrix_strain_far.png)

Key numbers from the 2026-07-18 run (`results/matrix/matrix_summary.json`):

- near field (0–800 nm): identical by construction.
- far field (7.9–8.5 um): leapfrog shows a dispersive wake with ~55 nm
  dominant wavelength; d'Alembert shows a sharp bipolar wavefront plus echoes
  spaced ~100–114 nm — the discrete pulse train seen in the published Fig. 3.
- peak strain at the wavefront nearly doubles (2.9e-3 -> 5.5e-3) once
  dispersion no longer smears the pulse.
- rocking curves: max |Δ log10 I| between the two strain models is ~0.6–0.8
  and confined to the weak shoulder ~+50 to +150 arcsec; the Bragg peak and
  near-peak structure are unchanged. This confirms the earlier diagnosis that
  the dispersion artifact barely affects the diffraction observable.
