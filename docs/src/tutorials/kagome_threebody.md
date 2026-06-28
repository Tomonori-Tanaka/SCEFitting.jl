# Kagome three-body

```@meta
CurrentModule = MagestyRebuild
```

Beyond pairs, the new physics is that a space-group operation can *permute* the sites of a
cluster. On the kagome lattice the three sites of a triangle are symmetry-equivalent
(``C_{3v}`` permutes them), so a 3-body term with unequal site angular momenta — e.g.
``l = (1,1,2)`` — must combine the three ``l``-orderings into a single **multi-term SALC**.
This tutorial builds the arbitrary-body-order basis, exhibits that multi-term channel, and
recovers a synthetic 3-body model from energies and torques.

## The kagome basis

```@example kagome
using MagestyRebuild
import Spglib
using LinearAlgebra, Random

rng = MersenneTwister(2026)
randcfg(nat) = reduce(hcat, (v = randn(rng, 3); v / norm(v)) for _ = 1:nat)

a, c = 2.0, 8.0
lat  = Lattice([a -a/2 0.0; 0.0 a*sqrt(3)/2 0.0; 0.0 0.0 c])
frac = [0.0 0.5 0.5; 0.5 0.0 0.5; 0.0 0.0 0.0]
kagome = Crystal(lat, frac, [1, 1, 1], ["Fe"])

interaction = Interaction(; nbody = 3, pair_cutoff = 1.2, lmax = [2])   # up to 3-body, l ≤ 2
basis = SCEBasis(kagome, interaction; backend = SpglibBackend())

(space_group = basis.spacegroup.symbol, n_salc = nsalc(basis))
```

## A multi-term SALC

The ``l = (1,1,2)``, ``L_f = 0`` three-body channel combines the distinct ``l``-orderings
into one SALC, stored as several *terms* per member:

```@example kagome
s112 = first(s for s in basis.salcs.salcs
             if s.key.body == 3 && s.ls == [1, 1, 2] && s.Lf == 0)
length(s112.members[1].terms)        # number of l-orderings folded into this one SALC
```

A single-`l` or equal-`l` channel collapses to one term and matches the ordinary
pairwise construction; the multi-term machinery only switches on when permutation-equivalent
sites carry unequal `l`. The mechanism is described in
[Architecture](../theory/architecture.md#Combined-space-projection).

## Recover a synthetic 3-body model

We draw random "true" couplings, generate in-span energies and the corresponding torques,
and recover the couplings with an energy + torque co-fit.

```@example kagome
m       = nsalc(basis)
configs = [randcfg(3) for _ = 1:200]
skel    = SCEDataset(basis, configs, zeros(length(configs)))   # for its design matrix
J_true  = randn(rng, m)
j0_true = 0.5
energies = j0_true .+ skel.X_E * J_true

model0  = SCEModel(basis, j0_true, J_true, basis.salcs.keys)
torques = [predict_torque(model0, c) for c in configs]

f = fit(SCEFit, SCEDataset(basis, configs, energies, torques), OLS(); torque_weight = 0.3)
(r2_energy = round(r2_energy(f); digits = 12),
 r2_torque = round(r2_torque(f); digits = 12),
 max_dJ    = round(maximum(abs, coef(f) .- J_true); sigdigits = 3))
```

Both observables are reproduced and every one of the ``m`` couplings — across all body
orders, `l`-channels, and `Lf` — is recovered to numerical precision. The multi-term
channel is fit on exactly the same footing as the bilinear ones, because
[`evaluate`](@ref) and the torque-gradient kernel both loop over a SALC's terms.

## Notes

- The same `nbody = K` construction works for any `K`; the cost grows with `K` and `lmax`.
- For a real fit you would supply DFT energies and torques (see
  [Persistence and I/O](../guide/io.md)) rather than synthetic in-span data, and likely a
  regularized estimator (see [Data and fitting](../guide/fitting.md)).
