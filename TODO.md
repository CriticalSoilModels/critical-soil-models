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

## New architecture — pseudocode written, not yet compiled

The following files describe the target architecture but are marked `PSEUDOCODE — does not compile`:

- `src/core/csm_model.f90` — abstract base type `csm_model_t`
- `src/integrators/euler_substep.f90` — modified Euler substepping integrator
- `src/conventions/voigt_conventions.f90` — Voigt reordering and inflate/deflate
- `src/models/mohr_coulomb_ss/mcss_model.f90` — MCSS concrete model
- `src/models/mohr_coulomb_ss/umat_mcss.f90` — MCSS UMAT wrapper

Suggested order to make them compile:

### Step 1 — `csm_model.f90` (abstract base)
No dependencies on other new-arch files. Start here; it gates everything else.

### Step 2 — `voigt_conventions.f90`
Standalone utility module. Needed by all UMAT wrappers. Can be done in parallel with Step 1.

### Step 3 — `euler_substep.f90` (integrator)
Depends only on `csm_model_t`. Can compile against the abstract type once Step 1 is done.
Has no concrete model calls — purely polymorphic dispatch through deferred procedures.

### Step 4 — Linear elastic model (validate the scaffold)
Simplest possible concrete model: no state, no yield surface.
Write `linear_elastic_model.f90` (extends `csm_model_t`) and a UMAT wrapper.
This validates the full stack (base type → integrator → UMAT) with the least complexity.
A passing single-element driver test here proves the architecture works end-to-end.

### Step 5 — MCSS model
Once the scaffold is proven with linear elastic:
- Implement `mcss_model.f90` properly (the pseudocode is a detailed spec)
- Implement `mcss_functions.f90` (yield fn, flow rule, softening)
- Unit-test the constitutive functions in isolation
- Wire up `umat_mcss.f90`

---

## Hardening modulus — needs verification

`mcss_hardening_modulus` in `src/models/mohr_coulomb_ss/mcss_functions.f90` computes the
analytical H from ∂F/∂c and ∂F/∂φ (see comment block above the function). The smoke tests
(H=0 at residual, H>0 at peak, no NaN at J=0) all pass, but the derivatives of K and a
w.r.t. φ have not been verified against a finite-difference reference. Before production
use, add a finite-difference gradient check (similar to the existing `dF/dSigma` test) and
confirm the sign convention matches physical softening behaviour in a single-element driver.

---

## Voigt shear ordering — review needed

`calc_J3_inv` and `calc_voigt_square` hardcode the `[11,22,33,12,13,23]` shear pairing
in their arithmetic. Everything else (J2, p, q, elastic stiffness) is already agnostic.
Review and find a better approach to make J3 and the lode angle independent of shear
slot ordering.

---

## Longer term

- Pull NorSand into the new architecture (currently legacy only)
- Set up CI/CD (GitHub Actions with `fpm test`)
- Extract invariants library to its own fpm dependency once signatures are stable
