# Drautz 2020 — tesseral-harmonic cluster expansion

**Citation.** R. Drautz, "Atomic cluster expansion of scalar, vectorial, and
tensorial properties including magnetism and charge transfer," *Phys. Rev. B*
**102**, 024104 (2020). DOI: 10.1103/PhysRevB.102.024104.

## Why it matters here

Fixes the **real (tesseral) spherical-harmonic `Zₗₘ` convention** (normalization
and Condon–Shortley sign) and the angular-momentum coupling used to build
symmetry-adapted basis functions over clusters of unit spin directions. The
from-scratch `basis/Harmonics.jl` and `basis/AngularMomentum.jl` reimplement this
convention; correctness is anchored by:

- closed-form `Zₗₘ` for low `l` and on-sphere central-difference of `∂ᵢZlm`
  (convention-independent), and
- agreement with the Magesty.jl oracle (which follows the same convention).

## Key conventions to preserve

- Per-site harmonic normalization carries `(4π)^(−1/2)`; an N-body cluster product
  carries `(4π)^(−N/2)`, cancelled by a `(4π)^(N/2)` design-matrix factor.
- Time-reversal even channel only (`Σ l_s` even), via the `(−1)^{Σl_s}` term in
  the symmetry projector.

(PDF, if present, lives at `papers/drautz-2020-tesseral.pdf` — local only.)
