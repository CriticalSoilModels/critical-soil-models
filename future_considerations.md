# Future Considerations

Design questions that are out of scope for the current architecture prototype but need
to be resolved before the library can support the full range of geotechnical models.

---

## 1. Full tensor vs Voigt notation

The current abstract interface passes stress as `real(dp) :: stress(6)` (Voigt notation).
This works for isotropic models and models whose internal state is scalar or vector-valued.
It breaks for:

- **Fabric/structure tensor models** (SANISAND, anisotropic critical state models) where
  internal state includes a symmetric 3x3 tensor `α_ij`. Stored as 6 components in Voigt
  for STATEV, but the tensor needs to be rotated as a second-order tensor — not as a
  stress vector — when the material frame rotates. A flat `real(dp)` snapshot doesn't
  capture this.

- **Models that work in principal stress space** and need to reconstruct the full tensor
  to find eigenvectors. The Voigt vector is sufficient input but intermediate steps
  require the 3x3 form.

**Options to consider:**
- Add an optional full-tensor path through the abstract interface
- Handle tensor-valued state variables separately from scalar state in snapshot/restore
- Require models that need full tensor treatment to handle the conversion internally

---

## 2. Finite strain and objective stress rates

For large deformation applications (MPM, large-strain FEM), the Cauchy stress must be
updated using an objective stress rate to remain frame-indifferent. Common choices:

- **Jaumann rate** — uses spin tensor `W`
- **Green-Naghdi rate** — uses rotation `R` from polar decomposition of `F`

The UMAT interface already passes `DFGRD0`, `DFGRD1` (deformation gradients) and `DROT`
(rotation increment) for this purpose. The current architecture ignores them.

**What needs to change:**
- The abstract interface or the UMAT wrapper needs a finite-strain path
- Stress rotation should happen at the wrapper layer, not inside each model
- Models should declare whether they are small-strain or finite-strain capable

---

## 3. Variable stress/strain dimensions (2D vs 3D)

The current pseudocode hardcodes `stress(6)` and `dstrain(6)` everywhere. This breaks for:

- **Plane strain / plane stress**: 4 components `[11, 22, 33, 12]`
- **Axisymmetric**: 4 components `[11, 22, 33, 12]` (where 33 is hoop)
- **PM4Sand and other inherently 2D models**: only defined for plane strain

The UMAT passes `NDI` (normal components) and `NSHR` (shear components) to communicate
the dimensionality. `NTENS = NDI + NSHR`.

**What needs to change:**
- The abstract interface signatures need to accommodate variable `NTENS`
- Models need to declare which dimensionalities they support
- The integrators need to be dimension-aware or dispatch to dimension-specific versions

---

## 4. Time-dependent models

The current `update_hardening(deps_p)` interface has no time argument. This is
insufficient for:

- **Viscoplastic models** (Bingham): strain rate and `DTIME` drive the update
- **Rate-dependent models** (MC strain rate): the inertial coefficient depends on the
  physical time increment, not just the strain increment
- **Creep models**: viscous strain accumulates on clock time independently of load steps

**Options to consider:**
```fortran
! Option A: add dt to update_hardening interface
subroutine harden_iface(self, deps_p, dt)

! Option B: store dt inside the model before calling update_hardening
call model%set_time_increment(dt)
call model%update_hardening(deps_p)

! Option C: separate deferred procedure for time-dependent update
procedure(viscous_iface), deferred :: update_viscous
```

Option A is the simplest but adds a `dt` argument that purely rate-independent models
must accept and ignore.

---

## 5. Multiphase and unsaturated soils

Some models require inputs beyond stress and strain:

- **Unsaturated soils** (Barcelona Basic Model, others): require matric suction `s = ua - uw`
  and degree of saturation `Sr` to compute Bishop effective stress
  `σ' = σ - [ua - χ(ua - uw)]`
- **Coupled hydromechanical models**: pore pressure as a driving variable, not just a
  result

These are not stress, not strain, and not material state — they are additional
**environmental inputs** per time step. The UMAT passes some of these via `ADDITIONALVAR`
(Anura3D) but there is no slot for them in the current abstract interface.

**What needs to change:**
- Consider an optional `env_vars(:)` argument on the integrator or the update procedure
- Or a `set_environment(suction, saturation, ...)` call before integration
- Models that don't need this simply ignore it

---

---

## 6. Solver interface formats beyond UMAT

UMAT is a legacy Abaqus format with a large, clunky argument list (~30 arguments) most
of which are unused by any given model. The goal is for this library to eventually support
cleaner calling interfaces alongside UMAT.

**The problem with UMAT as the only interface:**
- Forces every model to accept ~30 arguments regardless of what it needs
- Flat `PROPS` and `STATEV` arrays mean the calling code has no type safety — wrong
  index is a silent bug
- Tightly couples the library to Abaqus/Anura3D conventions
- Makes it harder to call models from Python, C, or other Fortran programs directly

**What the architecture already does right:**
The current design already solves this internally — `mcss_from_props`, `mcss_load_state`,
and the typed `mcss_model_t` are the real interface. UMAT is just one adapter on top.
Adding a new calling format means writing a new thin wrapper, not changing any model code.

**Possible future interfaces:**

```fortran
! Option A: Minimal direct interface — named arguments, no unused variables
call csm_stress_update(model, stress, dstrain)
! Caller constructs the model type directly rather than passing PROPS

! Option B: Incremental driver interface — the existing fumat/incremental-driver format
! Already closer to what this library wants to be

! Option C: C-compatible interface for Python/C interop
! Would require ISO_C_BINDING wrappers that serialize model state to/from flat arrays
subroutine csm_umat_c(stress, statev, ddsdde, dstran, props, nprops, nstatev) &
   bind(C, name="csm_mcss_update")
```

**What needs to happen:**
- Document the internal `model_t` interface as the canonical API
- UMAT wrappers become one of several adapters, not the definition of the interface
- Python interop (Option C) is explicitly out of scope for now — modern Fortran stdlib
  covers the testing and calibration workflows that would otherwise require Python

---

## Summary table

| Consideration | Affects interface? | Affects integrator? | Priority |
|---|---|---|---|
| Full tensor state variables | Yes — snapshot/restore | No | Medium |
| Finite strain / stress rotation | Yes — wrapper layer | No | Medium |
| Variable NTENS (2D/3D) | Yes — inflate/deflate at wrapper | No | High |
| Time-dependent models | Yes — update_hardening | Possibly | Medium |
| Multiphase inputs | Yes — new argument needed | No | Low |
| Non-UMAT interfaces | No — new wrappers only | No | Medium |
