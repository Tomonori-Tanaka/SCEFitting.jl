# MagestyRebuild.jl

A clean, extensible, Julia-native rebuild of **Magesty.jl** — fitting
**spin-cluster expansion (SCE)** models to noncollinear DFT data.

> **Status: work in progress (v0 vertical slice).** The numerical core (tesseral
> spherical harmonics, Clebsch–Gordan coupling, symmetry-adapted basis, design
> matrices, regression) is reimplemented from scratch and validated against
> Magesty.jl as a pinned numerical oracle. This is an architectural exploration;
> for production use see Magesty.jl.

## What it does

Given a crystal structure and a set of noncollinear DFT spin configurations
(energies + per-atom moments / effective fields), fit

```
E({e_a}) = j0 + Σ_φ J_φ Φ_φ({e_a})
```

together with the associated per-atom torques, where the `Φ_φ` are
symmetry-adapted linear combinations of real tesseral spherical harmonics over
clusters of spin directions `e_a` (unit vectors).

## Design highlights

- **Pluggable seams** via Julia package extensions: DFT data sources
  (`AbstractDFTSource`), estimators (`AbstractEstimator` + StatsAPI), symmetry
  backends (`AbstractSymmetryBackend`; Spglib optional). The lightweight core
  loads with no heavy dependencies.
- **Generalized cutoff neighbor list** — no fixed image grid; correct for
  triclinic cells and cutoffs spanning many lattice translations.
- **Canonical `SALCKey` column addressing** — a design-matrix column always means
  a fixed interaction, independent of construction order.
- **Spin-spiral-ready hooks** — neighbor pairs and cluster orbits retain the
  inter-site lattice translation `R`, leaving a clean seam for generalized-Bloch
  `E(q)` training data.

## References

- R. Drautz & M. Fähnle, *Phys. Rev. B* **69**, 104404 (2004) — spin-cluster
  expansion.
- R. Drautz, *Phys. Rev. B* **102**, 024104 (2020) — tesseral-harmonic /
  cluster-expansion formalism.
- T. Tanaka & Y. Gohda, *Phys. Rev. Research* **8**, 023300 (2026).
