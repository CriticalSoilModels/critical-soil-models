---
title: User Guide
---

# Critical Soil Models

A Modern Fortran library of geotechnical constitutive models, implemented as
UMAT subroutines compatible with FEA solvers such as Abaqus. The library is
decoupled from any PDE solver so models can be tested and validated independently.

## Getting started

### Building

```bash
conda activate fpm
fpm build
fpm test
```

### Using a model in your solver

All models are exposed through a single `umat` entry point in `src/umat.f90`.
Select a model by passing its name via the `cmname` argument, following the
Abaqus UMAT convention.

## Architecture

See the [API reference](../index.html) for full module and type documentation.
Key entry points:

- **`csm_model_t`** — abstract base type every model extends
- **`integrate_stress`** — generic stress integrator (Euler substepping, Ortiz-Simo)
- **`mod_stress_invariants`** — p, q, J2, J3, Lode angle and their derivatives
