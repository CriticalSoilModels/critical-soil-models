# critical-soil-models

[![Docs](https://github.com/CriticalSoilModels/critical-soil-models/actions/workflows/docs.yml/badge.svg)](https://criticalsoilmodels.github.io/critical-soil-models/)

A Modern Fortran library of geotechnical constitutive models, implemented as UMAT subroutines compatible with FEA solvers like Abaqus and decoupled from any PDE solver for independent testing and validation.

Geotechnical engineers spend too much time reimplementing and debugging constitutive relations that already exist. This library aims to collect, abstract, and accelerate constitutive model development by providing a framework that is independent of any one PDE solver. Decoupling the constitutive model from the solver makes it easier to isolate bugs, cross-validate implementations, and share models across projects.

---

## Implemented Models

### Linear Elastic
- Standard isotropic linear elasticity
- Fully integrated into the `csm_model_t` architecture
- Unit tested and validated through the generic Euler substepping integrator

### Bingham Viscoplastic
- Ported from Anura3D
- Stable implementation, not yet integrated into the new architecture or unit tested

### Mohr-Coulomb Strain Softening (MCSS)
- Associative or non-associative MC model with exponential softening of cohesion and friction angle from peak to residual
- Both legacy and new-architecture implementations exist
- Function-level and integration-level unit tests pass
- Hardening modulus derivatives pending finite-difference verification before production use

### Mohr-Coulomb with Strain Rate (MCSR)
- MC model with rate-dependent elastic moduli and yield parameters via an inertial coefficient
- Based on [Martinelli et al. (2022)](https://doi.org/10.1680/jgeot.21.00192)
- Both legacy and new-architecture implementations exist
- Integration tests cover elastic steps, plastic return mapping, and agreement between the Euler and cutting-plane integrators

### NorSand
- Legacy implementation only
- Not yet integrated into the new architecture or unit tested

---

## Setup

[pixi](https://pixi.sh) is a cross-platform package manager that installs all required compilers and tools automatically. No manual compiler installation needed.

**Install pixi** (one-time, then restart your terminal):

| Platform | Command |
|----------|---------|
| Linux / macOS | `curl -fsSL https://pixi.sh/install.sh \| bash` |
| Windows (PowerShell) | `iwr -useb https://pixi.sh/install.ps1 \| iex` |

**Install project dependencies:**

```bash
git clone https://github.com/CriticalSoilModels/critical-soil-models.git
cd critical-soil-models
pixi install
```

This installs the compilers, [fpm](https://fpm.fortran-lang.org/) (build tool), fortls (language server), and ford (documentation generator). Compiler availability by platform:

| Compiler | Linux | macOS | Windows |
|----------|-------|-------|---------|
| Intel ifx | ✓ | — | ✓ |
| LLVM flang | ✓ | — | — |
| lfortran | ✓ | ✓ | ✓ |
| gfortran | ✓ | ✓ | ✓ |

> **Note:** gfortran has a known runtime crash with the test suite on this project. Use ifx or lfortran instead.

---

## Building and Testing

Activate the environment first, then use `fpm` directly:

```bash
pixi shell       # enter the environment (compilers are now on your PATH)
fpm @ifx_test    # run the test suite
exit             # leave when done
```

Or prefix individual commands with `pixi run` without entering the shell:

```bash
pixi run fpm @ifx_test
```

### Choosing a compiler

| Platform | Use |
|----------|-----|
| Linux | `fpm @ifx_test` |
| Windows | `fpm @ifx_test` |
| macOS | `fpm @lfortran_test` |

### Common commands

```bash
# Run the full test suite
fpm @ifx_test

# Run with debug flags (bounds checking, tracebacks)
fpm @ifx_debug

# Run a specific test suite
fpm @ifx_test -- test_norsand_integration

# Run a specific test by name
fpm @ifx_test -- test_norsand_integration "elastic step: euler and cpa identical"

# Build without running tests
fpm @ifx

# Run examples (output written to example/outputs/)
fpm @ifx_example

# Generate documentation
ford fpm.toml

# Clean build artifacts (keeps downloaded dependencies)
fpm clean
```

### All fpm.rsp profiles

| Tag | Action | Compiler flags |
|-----|--------|----------------|
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

---

## CI Status

| Compiler | Platform | Status |
|----------|----------|--------|
| GFortran 12–14 | Linux | [![GFortran](https://github.com/CriticalSoilModels/critical-soil-models/actions/workflows/ci-gcc.yml/badge.svg)](https://github.com/CriticalSoilModels/critical-soil-models/actions/workflows/ci-gcc.yml) |
| Intel ifx 2024.2 | Linux | [![Intel ifx](https://github.com/CriticalSoilModels/critical-soil-models/actions/workflows/ci-ifx.yml/badge.svg)](https://github.com/CriticalSoilModels/critical-soil-models/actions/workflows/ci-ifx.yml) |
| Intel ifx 2024.2 | Windows | [![Intel ifx (Windows)](https://github.com/CriticalSoilModels/critical-soil-models/actions/workflows/ci-ifx-windows.yml/badge.svg)](https://github.com/CriticalSoilModels/critical-soil-models/actions/workflows/ci-ifx-windows.yml) |
| LFortran latest | Linux | [![LFortran](https://github.com/CriticalSoilModels/critical-soil-models/actions/workflows/ci-lfortran.yml/badge.svg)](https://github.com/CriticalSoilModels/critical-soil-models/actions/workflows/ci-lfortran.yml) |

---

## License

Distributed under the [GNU General Public License v3.0](LICENSE).

## Citations

If you use this library in research or a published manuscript, please cite it and the paper the model implementation was taken from. See `CITATION.cff` for citation details.
