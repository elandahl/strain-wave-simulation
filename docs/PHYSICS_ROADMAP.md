# Physics roadmap

This document records what physics the current strain models cover, what the
2022 Sci. Rep. paper invokes but does not simulate mechanistically, and what
we need next for the multi-material TRXD campaign. It is intentionally a
*roadmap*, not a commitment to implement everything at once: each item is
framed so it can become a switchable model in `strain_wave.models` (or a
crystal calculator in `xrd-strain-simulation`) without breaking the frozen
notebook-faithful reference.

## Project context: the two datasets

Beyond reproducing Sci. Rep. Fig. 3, this project will curate and re-analyze
two large time-resolved X-ray diffraction (TRXD) datasets of Cr films on
semiconductor substrates:

| Campaign | Facility / beamline | Strength | What it constrains |
|----------|---------------------|----------|--------------------|
| High angular resolution | APS 7ID-C (as in Jo et al., Sci. Rep. **12**, 16606 (2022)) | Fine rocking-curve shape, Brillouin fringes, instrument-aware angular resolution (~0.5 mdeg stated; ~22″ empirical effective) | Acoustic pulse-train morphology, G and k_e from fringe amplitude/location, near-interface thermal strain asymmetry |
| Dense long time base | Pohang Light Source (dense fill pattern) | Many delay points over a long temporal window | Thermal diffusion timescales, σ_TBC / k_s evolution, late-time thermal strain, possible carrier or recombination tails |

**Sample matrix (target coverage):**

- Film: Cr, **multiple thicknesses** (not only the 80 nm of Fig. 3)
- Substrates: **GaAs, Si, InSb, Ge**

Together these datasets demand a strain + XRD pipeline that is (i) multi-material,
(ii) thickness-parameterized, (iii) able to forward-model both early-time
acoustic Brillouin features and late-time thermal asymmetry, and
(iv) switchable so new physics can be added without invalidating the archival
Cr/GaAs notebook solver (`thermo-elastic-gaas` tag `paper-v1.0`).

## Current physics baseline (what we have)

Registered strain models in this repo (see `docs/ACOUSTIC_MODELS.md`):

| Model | Role |
|-------|------|
| `ttm_dalembert_cr_gaas` (**default**) | TTM + thermoelastic film/interface; exact d'Alembert far field. Physically preferred for opaque-film Cr/GaAs. |
| `ttm_fd_courant_cr_gaas` | Same near field; Courant-one FD substrate field. Validated against d'Alembert; foundation for future distributed sources / inhomogeneity. |
| `ttm_cr_gaas` | Historical notebook leapfrog (C ≈ 0.003). Archival / bit-for-bit with frozen repo only. |

**Physics already in the PDEs** (matches Jo et al. Eqs. 1–3, 8–9):

- Film two-temperature model: `C_e = γ T_e`, `k_e(T_e,T_l)`, constant G, laser source with exponential absorption
- Interface thermal boundary resistance (`R_ps = 1/σ_TBC`)
- Substrate Fourier heat diffusion with effective `k_s`
- Thermoelastic wave equation, impedance mismatch, free surface, absorbing far BC
- Acoustic pulse train from film round trips
- Downstream: GaAs (004) dynamical XRD + switchable instruments (`empirical`, `aps_7idc`)

**What this baseline can already do for the campaign:**

- Reproduce Sci. Rep. Fig. 3 class of Cr/GaAs results with paper-matched presets
- Vary Cr thickness (pulse-train spacing ∝ film thickness — already a paper observable)
- Fit effective G, k_e, σ_TBC, k_s against APS rocking curves and PLS delay scans
  *within* the Fourier + constant-G phenomenological framework

## Gaps relative to the Sci. Rep. paper

The paper's *forward* model is largely what we implement. The paper's
*interpretive* physics is often folded into effective fitted parameters rather
than simulated. Ranked by impact on the multi-substrate campaign:

### Priority A — needed to interpret the datasets mechanistically

**A1. Quasi-ballistic / non-Fourier substrate heat transport**

- Paper claim: Fourier's law breaks down; reduced k_s (Si ~4.3×, GaAs ~11× vs
  bulk) reflects quasi-ballistic phonon transport and MFP size effects near the
  interface; BTE / spectral phonon models are the proper description.
- Current code: single effective Fourier `k_s`.
- Needed: switchable substrate heat models, e.g.
  - two-channel (ballistic + diffusive),
  - hyperbolic / Cattaneo–Vernotte or Guyer–Krumhansl,
  - MFP-spectral / simplified BTE.
- Why it matters for the data: APS Brillouin features vs late PLS thermal
  asymmetry are the natural place to *partition* ballistic vs diffusive energy —
  exactly the paper's scientific punchline, now across four substrates and
  multiple film thicknesses.

**A2. Ballistic vs diffusive partition at the interface**

- Paper claim: only the diffusive component participates in thermal heating;
  ballistic phonons leave without thermalizing.
- Current code: single σ_TBC channel.
- Needed: interface model with (at least) two channels — ballistic (non-
  thermalizing) flux vs diffusive conductance — so σ_TBC is not forced to
  absorb both.

**A3. Multi-substrate materials library (GaAs, Si, InSb, Ge)**

- Paper: Cr/GaAs *and* Cr/Si.
- Campaign: also InSb and Ge.
- **Baseline complete:** `strain_wave.materials` now carries audited
  room-temperature bulk defaults for all four substrates; the generalized
  d'Alembert model has Cr/Si, Cr/Ge, and Cr/InSb aliases.
- **XRD baseline complete:** externally checked GaAs, Si, Ge, and InSb (004),
  10 keV production calculators.
- Remaining campaign work: attach actual per-sample Cr thickness, absorbed
  fluence, TBC/effective conductivity priors, and delay grids. These are not
  intrinsic material properties and must not be hidden in the registry.

### Priority B — refinements the paper invokes, currently lumped into fits

**B1. Temperature-dependent electron–phonon coupling G(T_e)**

- Paper Eq. (5): `G = G_RT [A_e/B_l (T_e + T_l) + 1]`; explains G ≈ 2× room-
  temperature value at high T_e.
- Current code: constant fitted G (while `k_e` *is* already T-dependent).
- Needed: optional `G_model = constant | kaganov` (or similar) in film TTM.
- Data impact: high-fluence APS runs and thickness series (different peak T_e).

**B2. Grain-boundary-derived k_e (Mayadas–Shatzkes)**

- Paper Eqs. (6–7): reduced k_e,eff from grain size (~30 nm) and reflection
  probability r.
- Current code: lumped `k_e_factor` (+ density factor).
- Needed: optional derivation of k_e,eff from (l_crystal, l_grain, r), with
  SEM-informed defaults per deposition batch.
- Data impact: comparing films of different thickness / grain morphology
  across the curated sample set.

**B3. Fluence-dependent thermalization localization**

- Paper: higher fluence → more localized energy → narrower, stronger pulse
  train (kinetic theory).
- Current code: linear TTM response at fixed G; won't capture nonlinear
  sharpening unless G(T) / nonlinear kinetics are on.
- Needed: mostly falls out of B1 + keeping nonlinear C_e(T_e); validate against
  multi-fluence subsets if present in the archives.

### Priority C — physics for *other* systems / companion papers (not this Sci. Rep.)

These are **out of scope for opaque thick Cr**, which absorbs all 800 nm light
in ~18.5 nm (≪ 80 nm film). They become central for thin / semi-transparent
films or bare semiconductors — and appear in the authors' companion literature
(carrier diffusion, Auger recombination, deformation-potential strain).

**C1. Photocarrier generation, transport, recombination**

- Generation over optical absorption depth in the semiconductor
- Ambipolar diffusion; Auger / radiative / SRH recombination heating
- Requires FD field solver with *sources* (use `ttm_fd_courant_*` foundation +
  implicit/subcycled diffusion; or generalized d'Alembert with Duhamel sources
  while the medium stays homogeneous)

**C2. Deformation-potential (electronic) strain**

- Stress from separated electron–hole pairs, distinct from thermal expansion
- Essential for interpreting non-opaque or semiconductor-only pump cases

**C3. Acoustic damping / anharmonicity**

- Frequency-dependent attenuation of the pulse train over long PLS delays
- Natural as a transfer function on the d'Alembert characteristic or a weak
  damping term in the FD propagator

**C4. Medium-changing acoustics** (when C=1 stops being magic)

- Depth-dependent v(z), buried interfaces, multilayers, nonlinear acoustics,
  2D/3D
- Requires conventional (refined) FD or spectral methods; use C=1 / d'Alembert
  as the homogeneous source-free acceptance test

## Acoustic numerics stance (already decided)

Documented in `docs/ACOUSTIC_MODELS.md` and the Jo/Soo follow-up email tag
`email-jo-soo-followup-c1-2026-07-18`:

- **Default for opaque-film work:** `ttm_dalembert_*` (exact continuum far field)
- **FD foundation for sources / future medium complexity:** `ttm_fd_courant_*`
  (C = 1; matches d'Alembert in the source-free homogeneous limit)
- **Historical only:** `ttm_cr_gaas` / frozen `thermo-elastic-gaas`

Coupled charge/heat sources remain compatible with C=1 via operator splitting
or implicit diffusion; physics that changes the *medium* is the true limit of
the magic step.

## Suggested implementation order (campaign-driven)

Aligned with curating APS + PLS data rather than with abstract completeness:

1. **Materials + XRD coverage (A3) — baseline done.** GaAs, Si, Ge, and InSb
   fixed-scope (004), 10 keV support is in place and externally checked.
2. **Presets / bookkeeping for the two campaigns — next.** Naming conventions for
   thickness × substrate × facility, fluence, delay grids; strain-profile
   metadata that records model id + fitted parameters for each fit.
3. **Non-Fourier / two-channel heat (A1–A2)** — highest scientific return for
   re-interpreting reduced k_s and the paper's ballistic/diffusive claim across
   four substrates and long PLS delays.
4. **G(T_e) and Mayadas–Shatzkes k_e (B1–B2)** — improve physical meaning of
   film parameters when comparing thicknesses and deposition batches.
5. **Carrier / deformation-potential stack (C1–C2)** — when thin-film or
   semiconductor-pumped subsets appear, or when companion-paper physics is
   revisited; build on `ttm_fd_courant` + source terms.
6. **Damping / medium complexity (C3–C4)** — as data demand (long delays,
   multilayers, high strain).

## Repository responsibilities

| Concern | Repo |
|---------|------|
| Archival notebook physics, Sci. Rep. Fig. 3 validation record | `thermo-elastic-gaas` (frozen; docs only) |
| Strain / heat / acoustics models, materials, presets, this roadmap | `strain-wave-simulation` |
| Dynamical diffraction, instruments, new crystal reflections | `xrd-strain-simulation` |
| Portable strain profiles between the two active repos | `docs/STRAIN_PROFILE.md` (v1 `.npz`) |

## Open questions for collaborators (Jo / Soo)

Useful to settle early in the data-curation phase:

1. Which Cr thicknesses and which of the four substrates are present in each
   archive (APS vs PLS), and at which fluences?
2. For the PLS dense-fill series: maximum delay and sampling cadence — does
   acoustic damping or late carrier physics become visible?
3. Should reduced k_s continue as a fitted effective Fourier parameter in v1
   re-analysis, with non-Fourier models as a parallel v2 science track?
4. Which Cr/Ge and Cr/InSb samples have the best independent thickness, XRR,
   fluence, and timing metadata for the first four-material fit?
5. Are SEM grain sizes / XRR thicknesses available per sample for Mayadas–
   Shatzkes and thickness-series modeling?

## References

- Jo et al., Sci. Rep. **12**, 16606 (2022), DOI 10.1038/s41598-022-20715-5
  (benchmark paper; Fig. 3 validation in `thermo-elastic-gaas`)
- Companion physics (carriers / Auger / electronic strain): Jo et al., Appl.
  Phys. Lett. **113**, 032107 (2018); Jo et al., Curr. Appl. Phys. **18**,
  1230 (2018); Lee et al., Appl. Sci. **9**, 4788 (2019)
- Acoustic model documentation: `docs/ACOUSTIC_MODELS.md`
- Courant / C=1 demonstration: `docs/images/courant_convergence.png`
