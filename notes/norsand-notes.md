# NorSand Implementation Notes

## Source

Ported from a legacy implementation by Luis E. Zambrano-Cruzatty (MIT license).
Original code included strain-rate hardening extensions on top of the base NorSand model
(Jefferies 1993, Jefferies & Shuttle 2002).

## Deliberate Simplifications

### Strain-rate extensions removed

The original code includes rate-dependent corrections to G, χ (dilatancy), and p_i via
the Isotache approach (Mesri et al. 1995, Xu & Zhang 2015). These are controlled by:

- `alpha_G`, `alpha_chi`, `alpha_pi` — strain rate scaling factors
- `RefRate` — reference strain rate
- `switch_smooth`, `N_S` — moving-average strain rate smoothing
- `Dashpot_method` — toggle between dashpot and original consistency approaches

**Decision:** These are excluded from the base implementation. The base NorSand model is
rate-independent. A separate `NorSandRate` model can extend the base if rate effects are
needed in the future.

### Drainage parameter S removed

The original code includes a parameter S (0 = drained, 1 = undrained) that scales a stress
ratio correction term T_s in the p_i hardening law:

```
dp_i/dε_eq += -S · T_s · p_i
```

At S=0 (drained) this term vanishes entirely. **Decision:** S is excluded; the base
implementation is always drained. Undrained behavior should be enforced at the FEA level
via volumetric strain constraints, not inside the constitutive model.

### Debug STATEV entries removed

The original stored `Rel_error_max` and `Drift_yield_max` (STATEV 24–25) for integrator
diagnostics. These are not part of the constitutive state and are excluded.

---

## Base NorSand PROPS layout (13 entries)

| Index | Symbol | Description | Units |
|-------|--------|-------------|-------|
| 1 | G_0 | Reference shear modulus | kPa |
| 2 | p_ref | Reference mean stress | kPa |
| 3 | nG | Shear modulus exponent | — |
| 4 | nu | Poisson's ratio | — |
| 5 | e_o | Initial void ratio | — |
| 6 | Gamma | CSL altitude at p = 1 kPa | — |
| 7 | lambda_c | CSL slope (ln scale) | — |
| 8 | R | OCR (used to initialise p_i) | — |
| 9 | M_tc | Critical friction ratio (TC) | — |
| 10 | N | Nova volumetric coupling | — |
| 11 | chi_tc | Dilatancy coefficient | — |
| 12 | H_0 | Hardening modulus intercept | — |
| 13 | H_y | Hardening modulus slope | — |

## Base NorSand STATEV layout (15 entries)

| Index | Symbol | Description | Units |
|-------|--------|-------------|-------|
| 1 | G | Current shear modulus | kPa |
| 2 | K | Current bulk modulus | kPa |
| 3 | p | Current mean effective stress | kPa |
| 4 | e | Current void ratio | — |
| 5 | psi | State parameter ψ = e − e_c | — |
| 6 | chi_tce | Current dilatancy coefficient | — |
| 7 | p_i | Image mean stress | kPa |
| 8 | M_i | Image stress ratio | — |
| 9 | switch_yield | Yielding flag (0.0/1.0) | — |
| 10–15 | eps_p(6) | Accumulated plastic strain (Voigt) | — |

Note: total 15 entries (9 scalars + 6 plastic strain components).

Note: `p` is stored in STATEV so that `elastic_stiffness()` (which takes no stress
argument in the abstract interface) can compute the pressure-dependent G = G_0*(p/p_ref)^nG.
`p` is updated at the start of each stress integration step before the integrator is called.

---

## Known limitations

- Cap surface (inner yield surface F2 = p − p_max) is not implemented. In the original code
  the cap is rarely activated in practice. A TODO comment marks where it would be added.
- Lode-angle dependence of M_i follows the standard Sheng–Sloan–Yu (2000) interpolation.
- N (Nova coupling) modifies the yield surface shape; set N=0 to recover standard NorSand.
