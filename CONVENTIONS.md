# Variable and Naming Conventions

Naming conventions for variables, types, and procedures in this library.
Follows the Modern Fortran style guide in `CLAUDE.md` with additional rules
specific to geotechnical and continuum mechanics notation.

---

## Derivative Notation

Two distinct patterns — never mix them:

### Increments (δ-type)

Prefix `d` directly on the variable name, no underscore:

```fortran
dsig        ! δσ  — stress increment
deps        ! δε  — total strain increment
deps_p      ! δεᵖ — plastic strain increment
deps_sub    ! δε  — substep strain increment
dsig_e      ! δσₑ — elastic stress increment
dsig_1      ! first Euler estimate of stress increment
dsig_2      ! second Euler estimate of stress increment
dlambda     ! δλ  — plastic multiplier increment
dt_sub      ! δt  — pseudo-time substep increment
```

### Partial Derivatives (∂-type)

Pattern `d"numerator"_by_d"denominator"` — the `_by_d` connector is unambiguous
and extends cleanly to second and cross derivatives:

```fortran
dF_by_dsig          ! ∂F/∂σ  — yield surface gradient w.r.t. stress
dG_by_dsig          ! ∂G/∂σ  — plastic potential gradient w.r.t. stress (flow direction)
dF_by_dkappa        ! ∂F/∂κ  — yield surface gradient w.r.t. hardening variable
dkappa_by_deps_p    ! ∂κ/∂εᵖ — hardening variable w.r.t. plastic strain

! Second and cross derivatives (consistent tangent)
d2F_by_dsig2        ! ∂²F/∂σ²
d2F_by_dsig_dkappa  ! ∂²F/(∂σ∂κ)
```

**Why keep the `d` in the denominator:** cross derivatives like `d2F_by_dsig_dkappa`
would be ambiguous without it — `d2F_by_sig_kappa` looks like a derivative w.r.t.
a single variable named `sig_kappa`.

---

## Standard Geotechnical Variables

### Stress and Strain

| Symbol | Fortran name | Notes |
|---|---|---|
| σ (Voigt vector) | `sig(6)` | Internal convention `[11,22,33,12,13,23]` |
| ε (total strain) | `eps(6)` | |
| εᵖ (plastic strain) | `eps_p(6)` | Accumulated |
| εᵉ (elastic strain) | `eps_e(6)` | |
| σ_tr (trial stress) | `sig_tr(6)` | Elastic predictor |
| σ₀ (stress at step start) | `sig_0(6)` | |

### Stress Invariants

| Symbol | Fortran name | Notes |
|---|---|---|
| p (mean stress) | `p` | Tension positive internally |
| q (deviatoric stress) | `q` | |
| θ (Lode angle) | `theta` | |

### Material Parameters (Elastic)

| Symbol | Fortran name | Notes |
|---|---|---|
| G (shear modulus) | `G` | [kPa] |
| ν (Poisson's ratio) | `nu` | [-] |
| K (bulk modulus) | `K` | [kPa] — derived where possible |
| λ, μ (Lamé constants) | `lame_1`, `lame_2` | Avoid `F1`, `F2` — conflicts with yield function |

### Material Parameters (Mohr-Coulomb)

| Symbol | Fortran name | Notes |
|---|---|---|
| c (cohesion) | `c` | [kPa] |
| φ (friction angle) | `phi` | Stored in radians internally |
| ψ (dilation angle) | `psi` | Stored in radians internally |
| c_p (peak cohesion) | `c_peak` | |
| c_r (residual cohesion) | `c_res` | |
| φ_p (peak friction) | `phi_peak` | |
| φ_r (residual friction) | `phi_res` | |

### Plasticity

| Symbol | Fortran name | Notes |
|---|---|---|
| F (yield function) | `F` | Universal in geotechnics literature — keep single letter |
| G (plastic potential) | `G_pot` | Distinguish from shear modulus `G` |
| λ (plastic multiplier) | `lambda_p` | `p` subscript = plastic |
| κ (hardening variable) | `kappa` | Model-specific meaning |
| εᵖ_eq (equiv. plastic strain) | `eps_p_eq` | |
| α (elastic-plastic split) | `alpha_ep` | In `[0,1]` |
| **D**ₑ (elastic stiffness) | `D_e(6,6)` | `_e` distinguishes from elastoplastic tangent |
| **D**ₑₚ (elastoplastic tangent) | `D_ep(6,6)` | Consistent tangent |

### Integrator Variables

| Symbol | Fortran name | Notes |
|---|---|---|
| t (pseudo-time) | `t_sub` | `[0,1]` over the step |
| Δt (pseudo-time increment) | `dt_sub` | |
| R (relative error) | `err_rel` | Substepping error estimate |
| α (elastic-plastic split) | `alpha_ep` | |

---

## Conversion Constants

```fortran
real(dp), parameter :: DEG_TO_RAD = acos(-1.0_dp) / 180.0_dp
real(dp), parameter :: RAD_TO_DEG = 180.0_dp / acos(-1.0_dp)
```

---

## Summary of Changes from Legacy Code

The following renames apply when migrating legacy code to the new architecture:

| Legacy name | New name | Reason |
|---|---|---|
| `stress` | `sig` | Matches σ notation |
| `dstrain` / `DSTRAN` | `deps` | Matches δε notation |
| `D` | `D_e` | Distinguishes from elastoplastic tangent |
| `F1`, `F2` | `lame_1`, `lame_2` | Avoids collision with yield function `F` |
| `T`, `DT` | `t_sub`, `dt_sub` | Fixes uppercase violation |
| `R` (error) | `err_rel` | Descriptive |
| `dSig1`, `dSig2` | `dsig_1`, `dsig_2` | Fixes mixed case |
| `df_dsig`, `dg_dsig` | `dF_by_dsig`, `dG_by_dsig` | Consistent derivative notation |
| `df`, `dg` (local) | `dF_by_dsig`, `dG_by_dsig` | Consistent with interface |
| `RAD` | `DEG_TO_RAD` | Clear direction |
| `dlambda` | `dlambda` | Already correct |
| `eps_p`, `deps_p` | `eps_p`, `deps_p` | Already correct |
| `p`, `q`, `theta` | `p`, `q`, `theta` | Standard — keep |
| `G`, `nu`, `c`, `phi`, `psi` | keep | Universal symbols |
