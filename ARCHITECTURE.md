# Architecture Design Notes

## Goals

- Decouple constitutive models from integration algorithms
- Enable unit testing of yield functions, flow rules, and hardening laws independently
- Share generic integrators (Euler substepping, Ortiz-Simo) across all models
- Set a single internal Voigt convention; translate at the solver boundary only

---

## Model Development Workflow

New constitutive models follow a five-step process that prioritises getting the physics
right before optimising for performance.

**Step 1 — Define POD structs**

Write `*_params_t` and `*_state_t` as plain-old-data types: no allocatables, no
procedure pointers, no OOP machinery. These are the primitive building blocks that
both the OOP interface and the GPU layer will share.

**Step 2 — Write the OOP interface with math inside**

Implement the five deferred procedures (`yield_fn`, `flow_rule`, `plastic_potential`,
`update_hardening`, `elastic_stiffness`) directly on the concrete model type. Put the
full constitutive math here. This keeps the physics readable and co-located with the
type definition during the development phase.

**Step 3 — Test and iterate**

Use the OOP interface and the generic `euler_substep` integrator against single-element
drivers. All unit tests operate at this level. Iterate on the physics until the model
behaviour is correct.

**Step 4 — Extract math into pure procedures**

Once the model is validated, move the constitutive math out of the OOP methods and into
module-level `pure` functions in `*_functions.f90`. These are annotated `!$acc routine seq`
(or equivalent) so they are GPU-callable. The functions take `params` and `state` as
explicit arguments — no `class()`, no dynamic dispatch.

**Step 5 — Make OOP methods thin wrappers**

Update each deferred procedure to delegate to the corresponding pure function:

```fortran
function mcss_yield_fn_method(self, sig) result(F)
   class(mcss_model_t), intent(in) :: self
   real(dp),            intent(in) :: sig(6)
   real(dp) :: F
   F = mcss_yield_fn(self%params, self%state, sig)  ! one line
end function
```

The OOP interface is behaviourally identical — all existing tests still pass. The pure
functions are now independently callable from GPU kernels without any polymorphism.

This workflow means the math is always visible during implementation (Step 2–3) and the
refactor to GPU-callable form (Steps 4–5) is mechanical, not a redesign.

---

## Layered Structure

```
┌─────────────────────────────────────────────┐
│  UMAT wrapper (thin)                        │  ← solver-facing (Abaqus, Anura3D, etc.)
│  - unpack PROPS/STATEV                      │
│  - translate Voigt convention               │
│  - call integrator                          │
│  - pack STATEV, fill DDSDDE                 │
├─────────────────────────────────────────────┤
│  OOP integrators (CPU / research path)      │  ← src/integrators/
│  - euler_substep (class(csm_model_t))       │
│  - ortiz_simo                               │
│  (know nothing about any specific model)    │
├─────────────────────────────────────────────┤
│  Abstract model type  (csm_model_t)         │  ← src/core/
│  - deferred: yield_fn, flow_rule,           │
│    plastic_potential, update_hardening,     │
│    elastic_stiffness, snapshot, restore     │
├─────────────────────────────────────────────┤
│  Concrete OOP models  (*_model.f90)         │  ← src/models/*/
│  - extend csm_model_t                       │
│  - hold params_t + state_t                  │
│  - deferred procedures are thin wrappers    │
│    around the pure function layer below     │
├─────────────────────────────────────────────┤
│  Pure function layer  (*_functions.f90)     │  ← src/models/*/
│  - no class(), no dynamic dispatch          │
│  - !$acc routine seq — GPU callable         │
│  - take (params_t, state_t, sig) explicitly │
│  - model-specific flat integrators live     │
│    here too (HPC / GPU path)                │
├─────────────────────────────────────────────┤
│  POD structs  (*_params_t, *_state_t)       │  ← src/models/*/
│  - no allocatables, no procedures           │
│  - shared by OOP and pure function layers   │
│  - SoA building blocks for MPM consumers   │
├─────────────────────────────────────────────┤
│  Shared infrastructure                      │  ← src/invariants/, src/conventions/
│  - stress/strain invariants                 │
│  - Voigt convention translation             │
└─────────────────────────────────────────────┘
```

**Two calling paths from the same library:**

- **Research / prototyping:** OOP integrator → abstract type → concrete model → pure functions
- **HPC / GPU:** flat integrator → pure functions directly (no polymorphism)

---

## Layer 1: Internal Convention

**Canonical internal Voigt order:** `[11, 22, 33, 12, 13, 23]`

This is the standard mechanics order and matches the incremental driver. All models and
integrators use this convention internally. Translation to/from solver conventions happens
only at the UMAT boundary.

```fortran
! src/conventions/voigt_conventions.f90
module mod_voigt_conventions

  integer, parameter :: ANURA3D_ORDER(6) = [1, 2, 3, 4, 6, 5]
  ! add others as needed: ABAQUS_ORDER, etc.

  function to_internal(v, from_order) result(v_out)
  function from_internal(v, to_order) result(v_out)

end module
```

---

## Layer 2: Abstract Model Type

The model type carries both **parameters** (fixed, set from PROPS) and **state** (evolving,
maps to/from STATEV). Stress is always passed separately — the integrator owns stress.

```fortran
! src/model_base/csm_model.f90
type, abstract :: csm_model_t
  ! No stress. No flat PROPS/STATEV arrays.
  ! Concrete types define their own named param + state fields.
contains
  procedure(yield_iface),     deferred :: yield_fn
    ! F = model%yield_fn(stress)

  procedure(flow_iface),      deferred :: flow_rule
    ! dF_by_dsig = model%flow_rule(sig)      [yield surface gradient]

  procedure(potential_iface), deferred :: plastic_potential
    ! dG_by_dsig = model%plastic_potential(sig)  [flow direction, non-associated]

  procedure(harden_iface),    deferred :: update_hardening
    ! call model%update_hardening(deps_p)    [updates internal state]

  procedure(stiff_iface),     deferred :: elastic_stiffness
    ! stiff_e = model%elastic_stiffness()    [6x6, depends only on current state]

end type
```

### Concrete model example (MCSS)

```fortran
type, extends(csm_model_t) :: mcss_model_t
  ! --- parameters (set once from PROPS) ---
  real(dp) :: G, nu
  real(dp) :: c_peak, c_res
  real(dp) :: phi_peak, phi_res
  real(dp) :: psi_peak, psi_res, shape_factor

  ! --- state (evolves, maps to/from STATEV) ---
  real(dp) :: c, phi, psi
  real(dp) :: eps_p(6)
contains
  procedure :: yield_fn          => mcss_yield_fn
  procedure :: flow_rule         => mcss_flow_rule
  procedure :: plastic_potential => mcss_plastic_potential
  procedure :: update_hardening  => mcss_update_hardening
  procedure :: elastic_stiffness => mcss_elastic_stiffness
end type
```

---

## Layer 3: Generic Integrators

Integrators live in `src/integrators/` and call only the deferred procedures — no
model-specific code.

```fortran
! src/integrators/euler_substep.f90
subroutine euler_substep(model, sig, deps, ftol, stol)
  class(csm_model_t), intent(inout) :: model
  real(dp), intent(inout) :: sig(6)
  real(dp), intent(in)    :: deps(6)
  real(dp), intent(in)    :: ftol, stol

  ! 1. elastic predictor
  stiff_e = model%elastic_stiffness()
  sig_tr  = sig + matmul(stiff_e, deps)

  ! 2. check yield
  F = model%yield_fn(sig_tr)
  if (F < ftol) then
    sig = sig_tr; return   ! elastic step
  end if

  ! 3. find elastic/plastic split (alpha_ep)
  ! ... pegasus / Newton on model%yield_fn ...

  ! 4. substep loop
  do while (t_sub < 1.0)
    dG_by_dsig = model%plastic_potential(sig)
    dF_by_dsig = model%flow_rule(sig)
    ! ... compute dsig_1, dsig_2, err_rel ...
    call model%update_hardening(deps_p)
    ! ... accept/reject step, adjust dt_sub ...
  end do

end subroutine
```

---

## Layer 4: UMAT Wrapper (thin)

The UMAT is the only place that knows about PROPS indexing, STATEV layout, and solver
conventions. Everything below this layer uses named fields and internal conventions.

```fortran
subroutine umat_mc_strain_softening(STRESS, STATEV, DDSDDE, ..., PROPS, ...)

  ! 1. Unpack — positional indexing lives only here
  model = mcss_from_props(PROPS)       ! params
  call mcss_load_state(model, STATEV)  ! state

  ! 2. Translate to internal convention
  sig  = to_internal(STRESS, ANURA3D_ORDER)
  deps = to_internal(DSTRAN, ANURA3D_ORDER)

  ! 3. Integrate (choice of integrator is a parameter)
  call euler_substep(model, sig, deps, ftol=PROPS(11), stol=1e-4_dp)
  ! or: call ortiz_simo(model, sig, dstran, ...)

  ! 4. Translate back
  STRESS = from_internal(sig, ANURA3D_ORDER)

  ! 5. Pack state
  call mcss_save_state(model, STATEV)

  ! 6. Tangent stiffness
  DDSDDE = deflate_stiffness(model%elastic_stiffness(), ptype)

end subroutine
```

---

## State rollback in the integrator

The modified Euler substep needs to temporarily advance state for the second estimate
then roll back. The two paths handle this differently:

**OOP integrator (CPU path):** uses the `snapshot` / `restore` deferred procedures.
The integrator holds an opaque `real(dp), allocatable :: saved(:)` — it never
interprets the contents. Index layout lives only inside each model's implementation.

```fortran
call model%snapshot(saved)
call model%update_hardening(deps_p1)
call euler_step(model, sig1, ...)
call model%restore(saved)
```

**Flat integrator (GPU path):** because `*_state_t` is a POD type with no allocatables,
snapshot/restore is a plain Fortran struct assignment — fixed-size, stack-allocated,
no subroutine call needed. This is GPU-safe and avoids any heap allocation inside a kernel.

```fortran
pure subroutine mcss_euler_substep(params, state, sig, deps, ftol, stol)
   type(mcss_state_t) :: saved_state   ! stack copy — compiler knows the size

   saved_state = state                 ! snapshot
   call mcss_update_hardening(params, state, deps_p1)
   ! ... compute error ...
   if (err > stol) state = saved_state ! restore
end subroutine
```

This is the payoff of requiring POD structs: snapshot/restore on GPU costs nothing
beyond a stack frame copy, with no abstraction overhead.

## Params vs state split

Both params and state are POD types (no allocatables, no procedure pointers).

- **Params** (`*_params_t`, fixed after construction from PROPS): named fields such as
  `G`, `phi_peak`. No wrapping overhead — params are read in every procedure and the
  extra indirection adds noise for no benefit. The OOP model holds `model%params`.

- **State** (`*_state_t`, evolves during integration): separate named type held as
  `model%state`. Accessed as `self%state%c`, `self%state%eps_p` etc. Separation means
  snapshot/restore is a clean struct copy, and the integrator never needs to know what
  fields exist. For the pure function layer, state is passed explicitly and returned by
  `intent(inout)` — the caller (MPM loop) owns the storage.

Both types are the SoA building blocks for MPM consumers: `real(dp) :: G(n_points)` in
the MPM code is just the SoA version of `params%G`.

**Hard rule:** `*_params_t` and `*_state_t` must never contain allocatable members,
pointers, or type-bound procedures. This is not a style preference — it is what makes
struct assignment a safe, fixed-size stack copy inside GPU kernels (snapshot/restore in
the flat integrator), and what allows MPM consumers to build SoA layouts from scalar
fields.

## Deferred procedure rules (Fortran)

- All arguments except `self` must match the interface exactly (type, kind, rank, intent).
- The `self` argument **can** be the concrete type (`class(mcss_model_t)`) even though
  the interface declares `class(csm_model_t)`. This is passed-object covariance.
- This means concrete implementations access named fields directly with no `SELECT TYPE`.

## Directory structure

```
src/
├── critical-soil-models.f90       ← top-level re-export
│
├── conventions/                   ← Voigt ordering, inflate/deflate, problem types
│   └── voigt_conventions.f90      ← kept at top level so finite strain wrappers
│                                     can use it independently of core/
├── core/                          ← framework only: abstract base type
│   └── csm_model.f90
│
├── integrators/                   ← generic integrators, no model-specific code
│   ├── euler_substep.f90
│   └── ortiz_simo.f90             ← future
│
├── invariants/                    ← future separate library; kept here until stable
│   ├── stress_invariants.f90
│   ├── stress_invariant_deriv.f90
│   ├── strain_invariants.f90
│   └── strain_invariant_deriv.f90
│
├── models/                        ← all concrete model implementations
│   ├── linear_elastic/
│   │   ├── linear_elastic_model.f90
│   │   └── umat_linear_elastic.f90
│   ├── bingham/
│   │   ├── bingham_model.f90
│   │   └── umat_bingham.f90
│   ├── mohr_coulomb_ss/
│   │   ├── mcss_model.f90
│   │   ├── mcss_functions.f90
│   │   └── umat_mcss.f90
│   ├── mohr_coulomb_sr/
│   │   ├── mcsr_model.f90
│   │   ├── mcsr_functions.f90
│   │   └── umat_mcsr.f90
│   └── norsand/
│       ├── norsand_model.f90
│       └── umat_norsand.f90
│
└── utils/                         ← audit needed; keep for now
    ├── array_helper.f90
    └── boolean_helper_funcs.f90
```

Each model directory contains exactly three files:
- `*_model.f90` — defines `*_params_t` and `*_state_t` POD structs; extends
  `csm_model_t`; deferred procedures are thin wrappers that delegate to `*_functions.f90`
- `*_functions.f90` — pure constitutive math (`!$acc routine seq`); takes params and
  state as explicit arguments; no `class()`; GPU-callable; independently unit-testable;
  also contains model-specific flat integrators for HPC use
- `umat_*.f90` — thin wrapper; only place PROPS/STATEV indexing and solver convention
  translation appear

---

## Invariants library

`src/invariants/` contains stress and strain invariant calculations (p, q, Lode angle,
derivatives). These are general continuum mechanics utilities with no dependency on soil
behaviour — long-term they belong in a separate library that this one references as an
fpm dependency.

They are kept here deliberately until:
- All invariant functions have passing tests
- Voigt convention is fully aligned with the internal `[11, 22, 33, 12, 13, 23]` standard
- Function signatures are stable enough that extraction won't require chasing changes
  across two repos

Do not add soil-model-specific logic to `src/invariants/`. It should remain pure
tensor/mechanics math so extraction is clean when the time comes.

## Type hierarchy (small strain vs finite strain)

The base type is named `csm_model_t` — neutral, not `csm_small_strain_t`. When finite
strain models are needed, the hierarchy expands by addition:

```
csm_model_t                  ← current base; interface is implicitly small-strain
                               (takes stress, not F)
└── csm_finite_strain_t      ← future; adds deferred procedures that take F, DROT etc.
```

Existing small-strain models extend `csm_model_t` directly and are untouched when
`csm_finite_strain_t` is added. Deferring this until a finite-strain model is actually
being implemented is the right call — it is a pure addition with no refactoring of
existing code.

## Public API — two levels

`src/critical-soil-models.f90` exports two levels of interface:

**Level 1 — Model types directly (primary API)**
Users construct a model, load state, and call integrators themselves. No UMAT involved.
Suitable for test drivers, calibration tools, and custom solvers.

```fortran
use critical_soil_models, only: mcss_model_t, mcss_from_props, &
                                 mcss_load_state, mcss_save_state, &
                                 euler_substep

type(mcss_model_t) :: model
model = mcss_from_props(props)
call mcss_load_state(model, statev)
call euler_substep(model, sig, deps, ftol=1e-8_dp, stol=1e-4_dp)
```

**Level 2 — UMAT wrappers (solver integration)**
Thin adapters for Abaqus, Anura3D, etc. Call the same model types internally.

```fortran
use critical_soil_models, only: umat_mc_strain_softening
```

`csm_model_t` and the integrators are also exported so users can write generic
code that works across models (e.g. a single-element driver that accepts any
`class(csm_model_t)`).

## Open Design Questions

- [ ] Which model to prototype first — linear elastic (validates scaffolding)
      or MCSS (real-world complexity)?
- [ ] Integrator tolerance parameters: pass as arguments (current) or store inside
      the model type alongside yield_tol / max_iters?
- [ ] Tangent stiffness (DDSDDE): elastic-only for now, consistent tangent later.
