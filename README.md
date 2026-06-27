# MagestyRebuild.jl

A clean, extensible, Julia-native rebuild of **Magesty.jl** — fitting
**spin-cluster expansion (SCE)** models to noncollinear DFT data.

> **Status: work in progress (v0 vertical slice).** The numerical core (tesseral
> spherical harmonics, Clebsch–Gordan coupling, symmetry-adapted basis, design
> matrices, regression) is reimplemented from scratch and validated against
> Magesty.jl as a pinned numerical oracle. This is an architectural exploration;
> for production use see Magesty.jl. Both SCE observables — **energy** and per-atom
> **torque** — are fitted; see [Status](#status) for what is and isn't implemented.

## What it does

Given a crystal structure and a set of noncollinear spin configurations with their
energies, fit

```
E({e_a}) = j0 + Σ_φ J_φ Φ_φ({e_a})
```

where the spins `e_a` are unit vectors and the `Φ_φ` are symmetry-adapted,
time-reversal-even scalar invariants built from real tesseral spherical harmonics
over clusters of spins. Fitting recovers the cluster coefficients `J_φ`.

The same coefficients also fix the per-atom **torque** `τ_a = e_a × ∂E/∂e_a`, the
SCE's other DFT observable. Passing per-configuration torques to `SCEDataset` and a
`torque_weight ∈ (0, 1]` to `fit` runs an energy+torque co-fit that minimizes
`(1 − w)·MSE_energy + w·MSE_torque`.

## Usage

```julia
using MagestyRebuild
import Spglib                     # load it to activate the SpglibBackend extension
                                 # (`import`, not `using`, to avoid a `Lattice` name clash)
using LinearAlgebra

# A 4-atom chain of identical spins along z
lat = Lattice([8.0 0 0; 0 8.0 0; 0 0 10.0])
frac = [0 0 0 0; 0 0 0 0; 0.0 0.25 0.5 0.75]
chain = Crystal(lat, frac, [1, 1, 1, 1], ["Fe"])

# nearest-neighbor 2-body interaction, isotropic (Heisenberg) channel only
interaction = Interaction(; nbody = 2, pair_cutoff = 2.6, lmax = [1], isotropy = true)
basis = SCEBasis(chain, interaction; backend = SpglibBackend())

# synthetic Heisenberg data E = J Σ_⟨ij⟩ e_i·e_j, then fit
configs = [mapreduce(_ -> (v = randn(3); v / norm(v)), hcat, 1:4) for _ in 1:30]
heis = basis.salcs.salcs[1]
J = 0.0137
E = [J * 0.5 * sum(c[:, m.atoms[1]]' * c[:, m.atoms[2]] for m in heis.members) for c in configs]

f = fit(SCEFit, SCEDataset(basis, configs, E), OLS())
r2_energy(f)                      # ≈ 1.0
2 * sqrt(3) * coef(f)[1]          # ≈ J  (recovered coupling)
```

A standalone, runnable version is in [`examples/heisenberg_chain.jl`](examples/heisenberg_chain.jl).

## Design highlights

- **Pluggable seams** via multiple dispatch + Julia package extensions: symmetry
  backends (`AbstractSymmetryBackend`; `NoSymmetry` in-tree, `SpglibBackend` in an
  extension) and estimators (`AbstractEstimator`; `OLS`/`Ridge` in-tree). The
  lightweight core loads with no heavy dependencies.
- **Generalized cutoff neighbor list** — no fixed image grid; correct for
  triclinic cells and cutoffs spanning many lattice translations.
- **Real Wigner-D from the package's own `Zₗₘ`** by an exact least-squares fit —
  convention-consistent by construction, handles improper rotations natively.
- **Orbit–stabilizer SALC projection** with a deterministic gauge and canonical
  **`SALCKey`** column addressing — a design-matrix column always means a fixed
  interaction, independent of construction order.
- **Spin-spiral-ready hooks** — neighbor pairs and cluster orbits retain the
  inter-site lattice translation `R`, leaving a clean seam for generalized-Bloch
  `E(q)` training data.

See [`docs/design-notes.md`](docs/design-notes.md) for the rationale behind these
refinements over Magesty.jl, and [`SPEC.md`](SPEC.md) for the realized architecture.

## Status

Implemented and validated (v0): geometry → symmetry (pluggable backend) → cluster
orbits → SALC basis (isotropic and anisotropic channels) → **energy and torque**
design matrices → `OLS`/`Ridge` fit (energy-only or energy+torque co-fit) →
`predict_energy` / `predict_torque`.

Not yet implemented (follow-ups): 3-body clusters (the machinery is `N`-generic;
enumeration is capped at 2-body), extensions for GLMNet estimators / VASP I/O /
Sunny export, `Tables.jl` results, and basis persistence.

## References

- R. Drautz & M. Fähnle, *Phys. Rev. B* **69**, 104404 (2004) — spin-cluster
  expansion.
- R. Drautz, *Phys. Rev. B* **102**, 024104 (2020) — tesseral-harmonic /
  cluster-expansion formalism.
- T. Tanaka & Y. Gohda, *Phys. Rev. Research* **8**, 023300 (2026).
