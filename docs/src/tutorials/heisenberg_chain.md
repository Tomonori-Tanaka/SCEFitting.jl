# Heisenberg chain

```@meta
CurrentModule = MagestyRebuild
```

This tutorial walks the whole pipeline on the simplest non-trivial system: a chain of
identical spins with a nearest-neighbor Heisenberg coupling. We build the basis, generate
synthetic DFT-like data from a known coupling, fit it, and recover the coupling — then
repeat with an energy + torque co-fit.

## Setup

```@example heis
using MagestyRebuild
import Spglib                       # activate the SpglibBackend extension
using LinearAlgebra, Random

rng = MersenneTwister(2026)
randspin()    = (v = randn(rng, 3); v / norm(v))
randcfg(nat)  = mapreduce(_ -> randspin(), hcat, 1:nat)
nothing # hide
```

## The crystal and the basis

A 4-atom chain along ``z`` (spaced so neighbors differ), with a nearest-neighbor, 2-body,
isotropic interaction — the Heisenberg channel only.

```@example heis
lat   = Lattice([8.0 0 0; 0 8.0 0; 0 0 10.0])
frac  = [0 0 0 0; 0 0 0 0; 0.0 0.25 0.5 0.75]
chain = Crystal(lat, frac, [1, 1, 1, 1], ["Fe"])

interaction = Interaction(; nbody = 2, pair_cutoff = 2.6, lmax = [1], isotropy = true)
basis       = SCEBasis(chain, interaction; backend = SpglibBackend())

(space_group = basis.spacegroup.symbol, n_salc = nsalc(basis))
```

There is a single SALC: the nearest-neighbor Heisenberg invariant
``\sum_{\langle ij\rangle}\hat{\boldsymbol e}_i\cdot\hat{\boldsymbol e}_j``.

## Synthetic data and the fit

We generate energies from a known coupling ``J_\text{true}`` over random spin
configurations, then fit with ordinary least squares.

```@example heis
heis    = basis.salcs.salcs[1]              # the Heisenberg SALC
J_true  = 0.0137
configs = [randcfg(4) for _ = 1:40]
E = [J_true * 0.5 * sum(dot(c[:, m.atoms[1]], c[:, m.atoms[2]]) for m in heis.members)
     for c in configs]

f = fit(SCEFit, SCEDataset(basis, configs, E), OLS())
J_recovered = 2 * sqrt(3) * coef(f)[1]

(r2 = round(r2_energy(f); digits = 12), J_true, J_recovered)
```

The fit is exact and recovers the coupling. The factor ``2\sqrt 3`` is the fixed
normalization between the physical coupling and the fitted SALC coefficient — a useful
sanity check: a Heisenberg model must satisfy ``J = 2\sqrt 3\,j_\varphi``.

## Adding the torque

The same coupling fixes the per-atom torque. For the Heisenberg energy,
``\partial E/\partial\hat{\boldsymbol e}_a`` is a sum of neighbor spins; the torque is its
cross product with ``\hat{\boldsymbol e}_a``.

```@example heis
function heis_torque(c, J)
    nat = size(c, 2); G = zeros(3, nat)
    for m in heis.members
        i, j = m.atoms[1], m.atoms[2]
        G[:, i] .+= (J * 0.5) .* c[:, j]
        G[:, j] .+= (J * 0.5) .* c[:, i]
    end
    reduce(hcat, cross(c[:, a], G[:, a]) for a = 1:nat)
end
torques = [heis_torque(c, J_true) for c in configs]

fc = fit(SCEFit, SCEDataset(basis, configs, E, torques), OLS(); torque_weight = 0.5)
(r2_energy = round(r2_energy(fc); digits = 12), r2_torque = round(r2_torque(fc); digits = 12))
```

Both observables are reproduced. Because [`predict_torque`](@ref) is the analytic
derivative of the surface [`predict_energy`](@ref) evaluates, the recovered model predicts
the torque on a fresh configuration too:

```@example heis
isapprox(predict_torque(fc, configs[1]), torques[1]; atol = 1e-8)
```

## Where this generalizes

- Drop `isotropy = true` to keep the anisotropic (`Lf > 0`) channels — the Dzyaloshinskii–
  Moriya and symmetric-Γ bilinears.
- Raise `lmax` to `[2]` for biquadratic / single-ion channels.
- Raise `nbody` to `3` for the [kagome three-body tutorial](kagome_threebody.md).
- Swap `OLS()` for `Ridge`, `Lasso`, or `ElasticNet` (see [Data and fitting](../guide/fitting.md)).
