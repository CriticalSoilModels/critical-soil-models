# Documentation Setup Progress

## Strategy

- **API docs**: Ford (generates HTML from `!!` docstrings in Fortran source)
- **Mathematical theory**: Typst (separate PDF artifacts, not embedded in Ford)
- **Environment**: `fpm` conda env (neither ford nor typst installed yet)

## Completed

### Typst theory document
- Created `docs/theory/invariants.typ` — full port of LaTeX invariants notes
  - Stress invariants: p, J1, J2, J3, Lode's angle
  - Derivatives: ∂p/∂σ, ∂J2/∂σ, ∂J/∂σ, ∂q/∂σ, ∂J3/∂σ, ∂θ/∂σ (full derivation)
  - Strain invariants: volumetric εᵥ, equivalent εq and their derivatives
  - Frobenius norm section
  - Notation macros defined at top of file (`#let sig`, `#let pd`, etc.)
  - Discrepancy with Potts eq. (VII.8) flagged in a callout block
- Created `docs/theory/refs.bib` — bibliography for all cited works

### Ford config
- Already configured in `fpm.toml` under `[extra.ford]`
  - `src_dir = ["src", "example"]`
  - `output_dir = "build/doc"`

## Remaining

- [ ] Install ford into `fpm` conda env
- [ ] Install typst into `fpm` conda env (or system)
- [ ] Verify ford builds successfully (`ford fpm.toml`)
- [ ] Verify typst compiles (`typst compile docs/theory/invariants.typ`)
- [ ] Write build script (`docs/build.sh` or `Makefile`) with:
  - `--all` — builds both ford and typst
  - `--ford` — ford only
  - `--typst` — typst only
- [ ] Confirm typst source lives in `docs/` (separate from `build/doc/`) so `fpm clean` cannot touch it
- [ ] Add `docs/` to `.gitignore` exclusions if needed (keep source, ignore compiled PDFs?)

## Risk: Will `fpm clean` delete Typst docs?

**No.** `fpm clean` only removes the `build/` directory. Ford outputs to
`build/doc/` (per `fpm.toml`). Typst source lives in `docs/theory/` and
compiled PDFs will go to `docs/theory/` as well — entirely outside `build/`.
No risk of collision.
