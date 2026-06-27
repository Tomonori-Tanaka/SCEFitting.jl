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

Standalone runnable versions are in [`examples/heisenberg_chain.jl`](examples/heisenberg_chain.jl)
and [`examples/kagome_threebody.jl`](examples/kagome_threebody.jl) (3-body / multi-term SALCs).

### Persistence and input files

Save a fitted model (or just a basis) to a self-contained, human-readable **TOML**
document and reload it later. Coefficients re-pair to the basis by `SALCKey`, so a
reloaded model predicts identically:

```julia
MagestyRebuild.save("model.toml", SCEModel(f))     # or save("basis.toml", basis)
model = MagestyRebuild.load(SCEModel, "model.toml")
predict_energy(model, configs)
```

A basis can also be built from a human-authored `input.toml` (inline crystal +
interaction + optional symmetry) instead of constructing `Crystal` / `Interaction`
in Julia. Training data and the estimator stay in Julia:

```toml
# input.toml
[structure]
lattice = [[8.0, 0.0, 0.0], [0.0, 8.0, 0.0], [0.0, 0.0, 10.0]]  # each entry = one lattice vector
positions = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.25], [0.0, 0.0, 0.5], [0.0, 0.0, 0.75]]
species = [1, 1, 1, 1]
species_labels = ["Fe"]

[interaction]
nbody = 2
pair_cutoff = 2.6
lmax = [1]
isotropy = true

[symmetry]
backend = "spglib"   # or "none"
tol = 1.0e-5
```

```julia
basis = SCEBasis("input.toml")     # reads [symmetry] backend/tol from the file
```

See [`examples/persist_and_input.jl`](examples/persist_and_input.jl) for the full loop.

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

Implemented and validated (v0): geometry → symmetry (pluggable backend) →
**arbitrary-body-order** cluster orbits → SALC basis (isotropic and anisotropic
channels, including the `N ≥ 3` coupling-path / `l`-ordering mixing) → **energy and
torque** design matrices → `OLS`/`Ridge` fit (energy-only or energy+torque co-fit) →
`predict_energy` / `predict_torque`. Cross-validated against Magesty.jl through
3-body (invariant-subspace dimensions agree exactly).

Basis/model **persistence** and a human-authored **`input.toml`** are implemented
(self-contained, human-readable TOML; coefficients re-pair by `SALCKey`).

Not yet implemented (follow-ups): extensions for GLMNet estimators / VASP I/O /
Sunny export, and `Tables.jl` results.

## References

- R. Drautz & M. Fähnle, *Phys. Rev. B* **69**, 104404 (2004) — spin-cluster
  expansion.
- R. Drautz, *Phys. Rev. B* **102**, 024104 (2020) — tesseral-harmonic /
  cluster-expansion formalism.
- T. Tanaka & Y. Gohda, *Phys. Rev. Research* **8**, 023300 (2026).
