---
# Model Checklist

Tracks implementation status, testing, and documentation for each constitutive model.

## Linear Elastic
- [x] Implemented
- [x] Unit tested
- [x] Integrated into new architecture (`csm_model_t`)
- [ ] Paper / reference documented

## Bingham viscoplastic
- [x] Implemented (ported from Anura3D)
- [ ] Unit tested
- [ ] Integrated into new architecture
- [ ] Paper / reference documented

## Mohr-Coulomb strain softening (MCSS)
- [x] Implemented (legacy + new architecture)
- [x] Unit tested (functions + integration)
- [x] Integrated into new architecture (`csm_model_t`)
- [ ] `umat_mcss.f90` wiring verified end-to-end
- [ ] Hardening modulus ∂F/∂c and ∂F/∂φ verified against finite-difference reference
- [ ] Paper / reference documented

## Mohr-Coulomb with strain rate (MCSR)
- [x] Implemented (legacy + new architecture)
- [x] Unit tested (integration tests: elastic, plastic, cprm vs euler)
- [x] Integrated into new architecture (`csm_model_t`)
- [ ] `umat_mcsr.f90` wiring verified end-to-end
- [ ] Paper: [Constitutive modelling of non-cohesive soils under high-strain rates](https://doi.org/10.1680/jgeot.21.00192)

## NorSand
- [x] Implemented (legacy only)
- [ ] Unit tested
- [ ] Integrated into new architecture
- [ ] Paper / reference documented

## PM4Sand (2D — plane strain)
- [ ] Implemented
- [ ] Unit tested
- [ ] Integrated into new architecture
- [ ] Paper / reference documented

## Modified Cam Clay
- [ ] Implemented
- [ ] Unit tested
- [ ] Integrated into new architecture
- [ ] Paper / reference documented

## Drucker-Prager
- [ ] Implemented
- [ ] Unit tested
- [ ] Integrated into new architecture
- [ ] Paper / reference documented

---

## Lower priority

- [ ] Original Cam-Clay
- [ ] Mohr-Coulomb with Cap
- [ ] Hypoplastic models

---

## General goals

- [ ] Expose API for users to plug in custom models via `csm_model_t`
- [ ] Interface from C / C++ / Python / Julia
- [ ] Extract invariants library to its own fpm dependency once signatures are stable
- [ ] Pull NorSand into new architecture
