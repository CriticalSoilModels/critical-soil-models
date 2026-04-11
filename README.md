# critical-soil-models

[![Docs](https://github.com/CriticalSoilModels/critical-soil-models/actions/workflows/docs.yml/badge.svg)](https://criticalsoilmodels.github.io/critical-soil-models/)

A Modern Fortran library of geotechnical constitutive models, implemented as UMAT subroutines compatible with FEA solvers like Abaqus and decoupled from any PDE solver for independent testing and validation.

Geotechnical engineers spend too much time reimplementing and debugging constitutive relations that already exist. This library aims to collect, abstract, and accelerate constitutive model development by providing a framework that is independent of any one PDE solver. Decoupling the constitutive model from the solver makes it easier to isolate bugs, cross-validate implementations, and share models across projects.

## Compiler Support

| Compiler | Version(s) | Status |
|----------|-----------|--------|
| GFortran | 12, 13, 14 | [![GFortran](https://github.com/CriticalSoilModels/critical-soil-models/actions/workflows/ci-gcc.yml/badge.svg)](https://github.com/CriticalSoilModels/critical-soil-models/actions/workflows/ci-gcc.yml) |
| Intel ifx | 2024.2 | [![Intel ifx](https://github.com/CriticalSoilModels/critical-soil-models/actions/workflows/ci-ifx.yml/badge.svg)](https://github.com/CriticalSoilModels/critical-soil-models/actions/workflows/ci-ifx.yml) |
| LFortran | latest | [![LFortran](https://github.com/CriticalSoilModels/critical-soil-models/actions/workflows/ci-lfortran.yml/badge.svg)](https://github.com/CriticalSoilModels/critical-soil-models/actions/workflows/ci-lfortran.yml) |

> **Note:** GFortran 15 has a known runtime crash when running the MCSR integration tests.
> GFortran 12–14 are unaffected.

## Implemented Models

### Linear Elastic
Standard isotropic linear elasticity. Fully integrated into the `csm_model_t` architecture.
Unit tested and validated through the generic Euler substepping integrator.

### Bingham Viscoplastic
Ported from Anura3D. Stable implementation, not yet integrated into the new architecture or unit tested.

### Mohr-Coulomb Strain Softening (MCSS)
Associative or non-associative MC model with exponential softening of cohesion and friction
angle from peak to residual. Both legacy and new-architecture implementations exist.
Function-level and integration-level unit tests pass. Hardening modulus derivatives are
pending a finite-difference verification before production use.

### Mohr-Coulomb with Strain Rate (MCSR)
MC model with rate-dependent elastic moduli and yield parameters via an inertial coefficient.
Based on [Martinelli et al. (2022)](https://doi.org/10.1680/jgeot.21.00192). Both legacy and
new-architecture implementations exist. Integration tests cover elastic steps, plastic return
mapping, and agreement between the Euler and cutting-plane (CPRM) integrators.

### NorSand
Legacy implementation only. Not yet integrated into the new architecture or unit tested.

---

## Setup

### Prerequisites

- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or Anaconda
- Git

### Create the environment

```bash
conda env create --file environment.yml
conda activate csm
```

This installs gfortran, flang, ifx, lfortran, fpm, fortls, ford, and graphviz.

> **Linux only:** the `ifx_linux-64` and Intel runtime packages in `environment.yml` are
> Linux-specific. Comment them out on macOS or Windows.

## Building and Testing

This project uses [fpm](https://fpm.fortran-lang.org/) with a `fpm.rsp` response file for
named build configurations. Invoke a configuration with `fpm @<tag>`.

```bash
# Run tests (default flags)
fpm @ifx_test
fpm @flang_test
fpm @lfortran_test

# Run tests with release optimisations
fpm @ifx_release_test
fpm @flang_release_test
fpm @lfortran_release_test

# Run a specific test suite
fpm @ifx_test -- test_mcss_integration

# Run a specific test within a suite
fpm @ifx_test -- test_mcss_integration "plastic step: stress on yield surface"

# Generate documentation (requires ford)
ford fpm.toml

# Clean build artifacts (preserves cached dependencies)
fpm clean
```

### Response file tag reference

| Tag | Action | Flags |
|-----|--------|-------|
| `@ifx` | build | default |
| `@ifx_test` | test | default |
| `@ifx_debug` | test | `-g -traceback -check all -fpe0 -warn all` |
| `@ifx_release` | build | `-O3 -xHost` |
| `@ifx_release_test` | test | `-O3 -xHost` |
| `@flang` | build | default |
| `@flang_test` | test | default |
| `@flang_debug` | test | `-g -fdebug-info-for-profiling -fcheck=bounds` |
| `@flang_release` | build | `-O3 -march=native -funroll-loops` |
| `@flang_release_test` | test | `-O3 -march=native -funroll-loops` |
| `@lfortran` | build | default |
| `@lfortran_test` | test | default |
| `@lfortran_debug` | test | `-g` |
| `@lfortran_release` | build | `-O3` |
| `@lfortran_release_test` | test | `-O3` |

**Naming convention:**
- `@<compiler>` — build only, default flags
- `@<compiler>_test` — run tests, default flags
- `@<compiler>_debug` — run tests with debug/bounds-checking flags
- `@<compiler>_release` — build only with release optimisation flags
- `@<compiler>_release_test` — run tests with release optimisation flags

---

## License

Distributed under the [GNU General Public License v3.0](LICENSE).

## Citations

If you use this library in research or a published manuscript, please cite it and the paper
the model implementation was taken from. See `CITATION.cff` for citation details.
