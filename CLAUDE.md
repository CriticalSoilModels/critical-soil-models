# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

A Modern Fortran library of geotechnical constitutive models (soil mechanics). The library implements UMAT (User Material) subroutines compatible with FEA solvers like Abaqus, decoupled from any PDE solver for independent testing and validation.

## Build and Test Commands

```bash
# Build the library
fpm build

# Run all tests
fpm test

# Run a specific test suite
fpm test -- test_stress_invar_suite

# Run a specific test within a suite
fpm test -- test_stress_invar_suite mean_stress_check

# Run the main application (prints banner)
fpm run

# Generate documentation (requires ford)
ford fpm.toml

# Clean build artifacts (preserves dependencies — stdlib is NOT re-downloaded)
fpm clean

# Clean everything including dependencies (stdlib will be re-cloned and rebuilt)
fpm clean --all
```

### Environment Setup
```bash
conda env create -f environment.yml
conda activate fpm
```
Requires: gfortran, fpm (Fortran Package Manager), fortls, ford. `stdlib` is fetched
automatically as a git dependency on first build and cached in `build/dependencies/`.
Use `fpm clean` (not `fpm clean --all`) to avoid re-downloading it.

## Architecture

See `ARCHITECTURE.md` for full design notes. Summary:

- `src/core/csm_model.f90` — abstract base type `csm_model_t` with seven deferred procedures
- `src/conventions/` — Voigt convention translation, inflate/deflate for 2D/3D
- `src/integrators/` — generic stress integrators (Euler substepping, Ortiz-Simo)
- `src/invariants/` — stress/strain invariants; future separate library, kept here until stable
- `src/models/` — one subdirectory per model, each with `*_model.f90`, `*_functions.f90`, `umat_*.f90`

### Constitutive Models

| Model | Directory | Status |
|-------|-----------|--------|
| Linear Elastic | `src/models/linear_elastic/` | Stable |
| Bingham viscoplastic | `src/models/bingham/` | Stable (from Anura3D) |
| Mohr-Coulomb strain softening | `src/models/mohr_coulomb_ss/` | Legacy + new arch in progress |
| MC with strain rate | `src/models/mohr_coulomb_sr/` | Legacy + new arch in progress |
| NorSand | `src/models/norsand/` | In development, no unit tests yet |

### Key Files

- `fpm.toml` — build configuration; `stdlib` dependency path is `../../stdlib`
- `ARCHITECTURE.md` — design decisions and rationale
- `future_considerations.md` — deferred design questions

---

# Modern Fortran Style Guide

A practical style guide for modern Fortran (2008+), with emphasis on GPU-accelerated scientific computing using OpenACC and CUDA Fortran.

## Naming Conventions

### Files

- Module files: `module_name.f90` (lowercase with underscores)
- Preprocessed files: `module_name.F90` (uppercase extension for files needing preprocessing)
- Test files: `test_module_name.f90`
- One module per file (submodules may be separate)

### Modules

Use a short project prefix to namespace modules and avoid collisions:

```fortran
! Good — prefixed with project abbreviation
module csm_solver
module csm_boundary

! Bad
module ShallowWaterSolver   ! CamelCase
module solver               ! No prefix, collision risk
```

Choose a 2-4 letter prefix and use it consistently. For this project: `csm_`.

### Derived Types

All derived types use the `_t` suffix with snake_case:

```fortran
type :: simulation_state_t
type :: mcss_model_t
type :: mcss_state_t
```

### Variable Naming for Mathematical Quantities

See `CONVENTIONS.md` for the full variable naming standard including derivative
notation, standard geotechnical symbols, and the legacy rename table.
Key rules:
- Increments: `d"var"` — `dsig`, `deps`, `dlambda`
- Partial derivatives: `d"num"_by_d"denom"` — `dF_by_dsig`, `dG_by_dsig`
- Stress vector: `sig(6)`, strain increment: `deps(6)`, elastic stiffness: `D_e(6,6)`

### Variables and Procedures

Snake_case with descriptive names:

```fortran
integer :: num_cells
real(wp) :: total_energy
subroutine compute_flux(state, flux)
```

**Exception:** Single-letter variables are acceptable for loop indices (`i`, `j`, `k`) and mathematical formulas matching published notation.

### Function and Subroutine Naming

Verb + noun patterns:

```fortran
subroutine compute_flux(...)
subroutine apply_boundary_conditions(...)
function calculate_timestep(...) result(dt)

! Logical returns use is_/has_/can_ prefixes
function is_dry(depth) result(dry)
function has_converged(residual, tol) result(converged)
```

### Units in Comments

Document physical units in variable declarations:

```fortran
type :: mcss_state_t
   real(wp) :: c    !! Cohesion [kPa]
   real(wp) :: phi  !! Friction angle [rad]
   real(wp) :: psi  !! Dilation angle [rad]
end type

real(wp), parameter :: GRAVITY = 9.81_wp   ! [m/s^2]
```

### Constants

UPPERCASE with underscores:

```fortran
real(wp), parameter :: GRAVITY = 9.81_wp
real(wp), parameter :: DRY_TOLERANCE = 1.0e-6_wp
integer, parameter :: MAX_ITERATIONS = 1000
```

---

## Required Practices

### Use Statements with `only` Clause

```fortran
! Good
use iso_fortran_env, only: real64, int32
use mod_csm_model, only: csm_model_t

! Bad — pollutes namespace, hides dependencies
use iso_fortran_env
```

### Implicit None

Always in modules and programs:

```fortran
module csm_my_module
   implicit none
   private
end module
```

### Intent Declarations

Always declare intent for all procedure arguments:

```fortran
subroutine compute_flux(state, dt, flux)
   type(state_t), intent(in)    :: state
   real(wp),      intent(in)    :: dt
   real(wp),      intent(out)   :: flux(:,:)
```

### Functions Should Have No Side Effects

```fortran
! Good — pure computation
function kinetic_energy(mass, velocity) result(ke)
   real(wp), intent(in) :: mass, velocity
   real(wp) :: ke
   ke = 0.5_wp * mass * velocity**2
end function

! Bad — function mutates state (use a subroutine instead)
function update_and_return_energy(state) result(energy)
   type(state_t), intent(inout) :: state    ! Side effect hidden in function
```

### Private by Default

```fortran
module csm_solver
   implicit none
   private

   public :: solver_t
   public :: initialize_solver
```

### Limit Procedure Arguments

Public procedures should have **6 or fewer arguments**. Group related arguments into derived types:

```fortran
! Bad
subroutine run_simulation(nx, ny, dx, dy, dt, t_end, h, hu, hv, ...)

! Good
subroutine run_simulation(mesh, state, config, output)
```

**Performance note:** For performance-critical compute kernels called millions of times, explicit scalar/array arguments are acceptable over derived type indirection.

---

## Forbidden Practices

- **No `goto`** — use structured control flow
- **No arithmetic IF** — use `if-then-else` or `select case`
- **No COMMON blocks** — use module variables or derived types
- **No EQUIVALENCE**
- **No fixed-form source** — all files must be `.f90` / `.F90`
- **No assumed-size arrays** (`arr(*)`) — use assumed-shape (`arr(:)`)
- **No `external` statements** — use modules for explicit interfaces
- **No implicit save** — avoid module-level variables with initialisation; if needed use explicit `save`

---

## Recommended Practices

### Working Precision

```fortran
module csm_precision
   use iso_fortran_env, only: real64
   implicit none
   integer, parameter :: wp = real64
end module

! Never use literal kind numbers
real(8) :: x           ! Bad — non-portable
real(wp) :: x          ! Good
```

### Prefer Allocatable Over Pointer

```fortran
! Good
real(wp), allocatable :: values(:)

! Bad — manual cleanup required
real(wp), pointer :: values(:) => null()
```

Use pointers only when aliasing, linked structures, or polymorphic returns are needed.

### Pure and Elemental Procedures

```fortran
pure function kinetic_energy(mass, velocity) result(ke)
   real(wp), intent(in) :: mass, velocity
   real(wp) :: ke
   ke = 0.5_wp * mass * velocity**2
end function

elemental function degrees_to_radians(deg) result(rad)
   real(wp), intent(in) :: deg
   real(wp) :: rad
   rad = deg * PI / 180.0_wp
end function
```

### Associate for Readability

```fortran
associate(p => state%mean_stress, &
          q => state%dev_stress,  &
          c => model%state%c,     &
          phi => model%state%phi)
   F = q - phi * p - c
end associate
```

**Compiler caveat:** `associate` support varies. Flang tends to be most robust. Test in performance-critical code and fall back to explicit temporaries if issues arise with gfortran/nvfortran.

### No Magic Numbers

```fortran
real(wp), parameter :: YIELD_TOL = 1.0e-8_wp

! Bad
if (F < 1.0e-8) return

! Good
if (F < YIELD_TOL) return
```

### Avoid Deep Nesting

Maximum 3-4 levels. Use early `cycle` and `return`:

```fortran
! Good — flat structure
do i = 1, n
   if (.not. active(i)) cycle
   if (.not. valid(i)) cycle
   ! work here
end do
```

### Labeled Loops

When using `cycle`/`exit` with nested loops:

```fortran
outer: do i = 1, n
   inner: do j = 1, m
      if (found(i,j)) exit outer
      if (skip(j)) cycle inner
   end do inner
end do outer
```

### Documentation (FORD compatible)

```fortran
type :: mcss_model_t
   !! Mohr-Coulomb strain softening model
   !!
   !! Extends csm_model_t. Params fixed after construction; state evolves.

   real(wp) :: G     !! Shear modulus [kPa]
   real(wp) :: nu    !! Poisson's ratio [-]
```

---

## GPU Programming

### OpenACC Basics

```fortran
!$acc parallel loop collapse(2) default(present)
do j = 1, ny
   do i = 1, nx
      flux(i,j) = hu(i,j)**2 / h(i,j) + 0.5_wp * GRAVITY * h(i,j)**2
   end do
end do
```

### Data Management Strategy

Copy to GPU once, keep there, copy back at end:

```fortran
!$acc enter data copyin(state%stress, state%statev)

do while (t < t_end)
   call stress_update(state)   ! all GPU work
end do

!$acc exit data copyout(state%stress, state%statev)
```

### Memory Coalescing

In Fortran (column-major), inner loop index should be the first array index:

```fortran
! Good — contiguous access
!$acc parallel loop collapse(2)
do j = 1, ny
   do i = 1, nx
      a(i,j) = b(i,j) + c(i,j)
   end do
end do
```

### Do Concurrent (Use with Caution)

```fortran
do concurrent (i = 1:n)
   result(i) = input(i)**2
end do
```

No `exit`, `cycle`, `return`, or `goto` inside. Compiler support varies — prefer regular `do` when in doubt.

---

## Common AI/LLM Mistakes in Fortran

### Pi is Not a Built-in Constant

```fortran
real(wp), parameter :: PI = 4.0_wp * atan(1.0_wp)
```

### `random_number` is a Subroutine

```fortran
call random_number(x)   ! correct
x = random_number()     ! wrong
```

### No I/O in `pure` Procedures

```fortran
pure function compute(x) result(y)
   real(wp), intent(in) :: x
   real(wp) :: y
   ! print *, x   ! compiler error — I/O is a side effect
   y = x**2
end function
```

### Declarations Before Executable Code

```fortran
subroutine foo()
   real(wp) :: x   ! declarations first
   x = 1.0_wp      ! then executable statements
end subroutine
```

### Array Constructors

```fortran
integer :: arr(3) = [1, 2, 3]     ! modern (Fortran 2003+)
integer :: arr(3) = (/1, 2, 3/)   ! old style, still valid
integer :: arr(3) = (1, 2, 3)     ! wrong — compiler error
```

---

## Summary Table

| Category | Do | Don't |
|---|---|---|
| Types | `state_t`, `mcss_model_t` | `State`, `TModel` |
| Variables | `num_cells`, `yield_tol` | `nCells`, `yt` |
| Constants | `MAX_ITER`, `YIELD_TOL` | `maxIter`, `tol` |
| Imports | `use mod, only: x` | `use mod` |
| Arrays | `arr(:)` | `arr(*)` |
| Memory | `allocatable` | `pointer` (unless needed) |
| Control | `if/select/do` | `goto` |
