# TODO

## Invariants library — DONE ✓

The invariants library (`src/invariants/`) is stable and fully tested:
- All functions follow `calc_` prefix; derivatives use `d"num"_by_d"denom"` (e.g. `dJ2_by_dsig`)
- All argument and local variable names are lowercase snake_case
- `lode_angle` used everywhere (not `theta`)
- Full-form reference implementations in `test/mod_stress_invar_refs.f90` and `test/mod_strain_invar_refs.f90`
- All utilities in `src/utils/voigt_utils.f90` are pure functions with no dead arguments
- Tests use fixed stress states (no `random_number`)

## Repo-wide modernisation — partially done

Completed across the whole repo:
- All `double precision` replaced with `real(dp)` (`use stdlib_kinds, only: dp`)
- Dilation variable renamed from `Dp`/`D_p` → `dilation`/`dil_u` (avoids `dp` kind collision)
- `state_params_deriv`: `Get_dD_to_dI` → `calc_ddil_by_dI`, `Get_dD_to_dEpsP` → `calc_ddil_by_deps_p`

Still needed in `src/models/mohr_coulomb_sr/` (outside invariants scope):
- `yield_function.f90`: local vars `dpdsig`, `dqdSig`, `dJ3dsig`, `dlode_angle_dSig` → `_by_dsig` convention
- `plastic_potential.f90`: same
- General camelCase variable names in legacy MCSS and MCSR files

---

## New architecture — DONE ✓

All steps compiled and unit-tested (39 tests passing):

- **Step 1** `csm_model.f90` — abstract base type `csm_model_t` ✓
- **Step 2** `voigt_conventions.f90` — Voigt reordering and inflate/deflate ✓
- **Step 3** `euler_substep.f90` — modified Euler substepping integrator ✓
- **Step 4** Linear elastic model — full stack validated (base → integrator → UMAT) ✓
- **Step 5** MCSS model — functions, types, UMAT wiring complete ✓

---

## MCSS refactoring — in progress

Recent work on `src/models/mohr_coulomb_ss/`:

### Completed ✓
- Renamed all derivative functions to `d<num>_by_d<denom>` convention:
  - `calc_dF_dsig_abbo` → `calc_dF_by_dsig_abbo`
  - `calc_dK_dlode` → `calc_dK_by_dlode`
- Replaced inline dJ/dσ in `calc_dF_by_dsig_abbo` with `calc_dJ_by_dsig` from invariants library
- Added `calc_dJ_by_dsig` to `mod_stress_invar_deriv` (thin wrapper over dJ2/dσ ÷ 2J)
- Refactored `mcss_hardening_modulus` into orchestrator + two private helpers:
  - `calc_dF_by_dstate` — ∂F/∂c and ∂F/∂φ
  - `calc_dstate_by_deps_eq` — dc/dεq and dφ/dεq (exponential softening rates)
- Renamed `lode_t` → `lode_tr` throughout (`abbo_sloan_params_t` field + all uses)
- Power operations consistent throughout module (`*` not `**2`, avoids gfortran IPO instability)
- Moved `pm4sand` to `legacy/` (depended on Anura3D internals; was blocking fpm build)
- Fixed `umat_linear_elastic` to use `integrator_params_t` / `DEFAULT_INTEGRATOR_PARAMS`
- Added `abbo_sloan_preset(lode_tr_deg)` factory in `mod_abbo_sloan_presets`:
  - Returns pre-calibrated constants for 25.0° and 29.5° transition Lode angles
  - Exported `DEFAULT_SMOOTH_COEFF` from `mod_mcss_types`

### Still needed
- Integration tests: single-element driver for MCSS (explicitly deferred — architecture first)
- `umat_mcss.f90` wiring (UMAT wrapper not yet connected end-to-end)

---

## Hardening modulus — needs verification

`mcss_hardening_modulus` in `src/models/mohr_coulomb_ss/mcss_functions.f90` computes the
analytical H from ∂F/∂c and ∂F/∂φ (see comment block above the function). The smoke tests
(H=0 at residual, H>0 at peak, no NaN at J=0) all pass, but the derivatives of K and a
w.r.t. φ have not been verified against a finite-difference reference. Before production
use, add a finite-difference gradient check (similar to the existing `dF/dSigma` test) and
confirm the sign convention matches physical softening behaviour in a single-element driver.

---

## Longer term

- Pull NorSand into the new architecture (currently legacy only)
- Extract invariants library to its own fpm dependency once signatures are stable
