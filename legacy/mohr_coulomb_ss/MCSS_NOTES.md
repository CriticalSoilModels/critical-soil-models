# Legacy MCSS Implementation Notes

Review of `umat_mcss_legacy.f90` and `mcss_functions.f90` to inform the new architecture.

---

## Yield Function (Abbo & Sloan 1995)

The yield surface is a **smooth hyperbolic approximation** to Mohr-Coulomb:

```
F = p·sin(φ) + sqrt(J²·K² + a²·sin²(φ)) - c·cos(φ)
```

where:
- `p` = mean stress = (σ₁+σ₂+σ₃)/3
- `J` = sqrt(J₂) = deviatoric stress invariant
- `K(θ)` = Lode-angle-dependent shape function (see below)
- `a = 0.0005·c·cot(φ)` = tip smoothing parameter (small, makes F differentiable at apex)

This is **not** the simplified `F = q + η·p - c` form — it uses the proper Mohr-Coulomb strength envelope expressed in (p, J, θ) space. At LodeT (transition Lode angle), K switches between:

- `|θ| < LodeT`: `K = cos(θ) - (1/√3)·sin(φ)·sin(θ)` (exact Mohr-Coulomb)
- `|θ| ≥ LodeT`: `K = A - B·sin(3θ)` (polynomial rounding at corners)

The active constants use **LodeT = 29.5°** (A1, A2, B1, B2 hardcoded; LodeT=25° and LodeT=30° options commented out).

---

## Stress Invariants

```
p     = (1/3)(σ₁+σ₂+σ₃)
J₂    = (1/6)[(σ₁-σ₂)²+(σ₁-σ₃)²+(σ₂-σ₃)²] + τ₁₂²+τ₂₃²+τ₁₃²
J     = sqrt(J₂)
J₃    = det(s)   (deviatoric stress third invariant)
θ     = (1/3)·arcsin(-3J₃ / (2J³·(1/√3)))   [Lode angle, ∈ [-π/6, π/6]]
```

S3TA (sin(3θ)) is clamped to [-1, 1] before the arcsin to handle floating-point noise.

Voigt convention: `[σ₁₁, σ₂₂, σ₃₃, σ₁₂, σ₂₃, σ₁₃]` (internal; UMAT boundary reorders from Anura3D order).

---

## Plastic Potential (Non-Associated)

Same mathematical form as yield function but with dilation angle `ψ` replacing friction angle `φ`:

```
G = p·sin(ψ) + sqrt(J²·K_ψ² + a_ψ²·sin²(ψ)) - c·cos(ψ)
```

Special case near pure hydrostatic (J → 0): `ψ` is linearly interpolated toward `φ` over the range `J ∈ [0, J0]` (J0 = 0.001) to avoid singularity. If `φ == ψ`, DPPDSig = DFDSig (associated flow).

---

## Softening Law — Exponential Decay

Single scalar `factor` (shape parameter) controls the rate:

```fortran
c   = c_r + (c_p - c_r) * exp(-factor * eps_p_eq)
phi = phi_r + (phi_p - phi_r) * exp(-factor * eps_p_eq)
psi = psi_r + (psi_p - psi_r) * exp(-factor * eps_p_eq)
```

All three parameters soften with the **same rate** and the **same scalar equivalent plastic strain** `eps_p_eq`.

Derivatives with respect to `eps_p_eq`:
```fortran
dc/deps_p_eq   = -factor * (c_p - c_r) * exp(-factor * eps_p_eq)
dphi/deps_p_eq = -factor * (phi_p - phi_r) * exp(-factor * eps_p_eq)
dpsi/deps_p_eq = -factor * (psi_p - psi_r) * exp(-factor * eps_p_eq)
```

---

## Equivalent Plastic Strain

Deviatoric-only norm (consistent with J₂ plasticity convention):

```
eps_p_eq = sqrt( (2/3) * e_p : e_p )
```

where `e_p` is the deviatoric part of the plastic strain tensor. In Voigt form:

```
eps_p_M = (1/3)(eps_p(1)+eps_p(2)+eps_p(3))
e_p(i)  = eps_p(i) - eps_p_M   [i=1,2,3]

eps_p_eq = sqrt( (2/3)*(e_p(1)²+e_p(2)²+e_p(3)²)
               + (4/3)*(eps_p(4)²+eps_p(5)²+eps_p(6)²) )
```

Note: shear components get factor 4/3 = 2*(2/3) because engineering shear strain is 2*tensorial shear strain.

> **⚠ BUG — incorrect shear factor**
>
> The 4/3 shear factor above is wrong when the Voigt vector stores **engineering shear strains**
> (γ = 2ε_tensor), which is the convention used throughout this codebase.
>
> Correct derivation: ε:ε in full-tensor form includes a factor of 2 for each off-diagonal pair,
> so ε:ε = ε₁₁²+ε₂₂²+ε₃₃² + 2(ε₁₂²+ε₁₃²+ε₂₃²). Substituting γᵢⱼ = 2εᵢⱼ:
>
> ```
> ε:ε = ε₁₁²+ε₂₂²+ε₃₃² + γ₁₂²/2 + γ₁₃²/2 + γ₂₃²/2
> ```
>
> Therefore the correct formula with engineering shear strains is:
>
> ```
> eps_p_eq = sqrt( (2/3)*(e_p(1)²+e_p(2)²+e_p(3)²)
>                + (1/3)*(eps_p(4)²+eps_p(5)²+eps_p(6)²) )
> ```
>
> The 4/3 factor is only correct if the Voigt vector stores **tensor** (non-doubled) shear
> components — a different convention.
>
> This was confirmed by porting the reference Python implementation
> (`calc_dev_strain_invariant`) to Fortran and comparing against `calc_eps_q_inv`
> from the invariant library, which uses `calc_two_norm_tensor_strain` (applies 0.5
> to shear terms, equivalent to the 1/3 shear factor above). Both agreed. The 4/3
> formula gave values ~2× larger for states with significant shear plastic strain.
>
> **Fixed in**: `src/models/mohr_coulomb_ss/mcss_functions.f90` — hardening now uses
> `calc_eps_q_inv(calc_dev_strain(eps_p, calc_eps_vol_inv(eps_p)))` from the library.

Derivative `d(eps_p_eq)/d(eps_p)` — see `CalculateDerivativesEquivalentPlasticStrainRespectPlasticStrain`.

---

## Integration Schemes

Two schemes selectable via PROPS(10):

### 0: Modified Euler Substepping
1. Elastic predictor: `σ_trial = σ + D_e·Δε`
2. If `F(σ_trial) > FTOL`: find elastic-plastic split `α` via **Pegasus method** (secant iteration, max 1000 iters)
3. Substep the plastic integration with error control
4. Apply **end-of-step drift correction** (consistent tangent, fallback to normal correction)

### 1: Ortiz-Simo (Return Mapping)
Gradient descent from elastic predictor back to yield surface:
```
dλ = F / (n·D_e·m - H)
H  = (dF/dXs · dXs/deps_p_eq) · (deps_p_eq/deps_p · m)
σ  = σ - dλ · D_e · m
```
Iterates until `|F| < FTOL` or `counter > max_iterations`.

---

## State Variables

STATEV layout (9 entries):
- (1): c — current cohesion
- (2): φ — current friction angle (radians)
- (3): ψ — current dilation angle (radians)
- (4–9): eps_p(6) — accumulated plastic strain (Voigt)

Initialization: at KSTEP=1 & KINC=1, state is initialized from peak values.

---

## Design Questions for New Implementation

### 1. Yield Function Form
Use the full **Abbo & Sloan smooth hyperbolic** form — not simplified. The smoothing at the tip (parameter `a`) and the Lode-angle-dependent K function are both needed for robustness.

Keep the Lode-angle transition constants (LodeT, A1, A2, B1, B2) as named module parameters so they can be changed without hunting through the code.

### 2. Lode Angle Dependency
Yes, full K(θ) dependency — required to recover correct Mohr-Coulomb strength in triaxial extension vs. compression.

### 3. Softening Law
Exponential decay with a single `factor` (shape parameter). All three strength params (c, φ, ψ) soften with the same rate. This maps cleanly to `mcss_update_hardening`.

### 4. What to Rewrite vs. Reuse
**Rewrite** (new module, clean API):
- `calc_yield_fn_mcss(p, J, lode, s3ta, c, phi)` — pure function
- `calc_dF_by_dsig_mcss(sig, c, phi, psi)` — subroutine
- `calc_dg_plas_by_dsig_mcss(sig, c, phi, psi)` — subroutine
- `calc_eps_p_eq(eps_p)` — **do not port**; use `calc_eps_q_inv(calc_dev_strain(...))` from the invariant library (see shear factor bug note above)
- `calc_softening_params(eps_p_eq, ...)` — pure function
- `calc_dsoftening_by_deps_p_eq(...)` — pure function

**Keep logic, adapt interface**:
- `DetermineElasticProportionPegasusMethod` → rename, clean up intent/kind
- `EndOfStepCorrection` → adapt for new signature

**Skip (handled by generic integrator)**:
- The `DetermineDSigAndDEpsP` elastoplastic matrix computation — `euler_substep` does this generically.

### 5. Pending `mcss_yield_fn` Fix
The current `mcss_yield_fn` placeholder in `mcss_model.f90` uses `F = q - phi*p - c`, which is dimensionally wrong (phi is in radians, used as a multiplier on p directly). This needs to be replaced with the proper Abbo & Sloan form once `mcss_functions.f90` is ported.
