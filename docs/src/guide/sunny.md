# Sunny export

```@meta
CurrentModule = MagestyRebuild
```

A fitted SCE model is an energy surface; for linear spin-wave theory you want it as a
[Sunny.jl](https://sunnysuite.github.io/Sunny/) `System`. [`to_sunny`](@ref) builds a real
`Sunny.System` from the Sunny-representable channels of a model. It is provided by an
extension that activates under `using Sunny`.

```julia
using Sunny
sys = to_sunny(model; spins = 3/2)        # a real Sunny.System
```

## What is exported

Sunny represents two channels, and those are exactly what the conversion keeps:

| SCE channel | Sunny term |
|-------------|------------|
| `ls = [1,1]` 2-body (bilinear) | a `3 × 3` exchange matrix (Heisenberg + Dzyaloshinskii–Moriya + symmetric Γ in one) |
| `ls = [2]` 1-body | a traceless-symmetric single-ion anisotropy tensor |
| `ls = [0…]` | folded into the constant `j0` |

Every other channel — three-body and higher, higher `l` — **cannot** be represented by a
Sunny `System` and is **reported as skipped**, never silently dropped:

```julia
sys = to_sunny(model; spins = 3/2)        # @warn lists any skipped channels
```

## Spin length and mode

Spins enter only at assembly. The SCE couplings are fit with *unit* directions, so the
exchange is rescaled `J = M / (Sₐ S_b)` to make Sunny's length-``S`` dipoles reproduce the
unit-vector energy. The single-ion term carries the classical (`mode = :dipole_uncorrected`,
factor ``1/s^2``) or quantum rank-2 (`mode = :dipole`, factor ``2/(s(2s-1))``, requires
``s > 1/2``) rescaling.

```julia
sys = to_sunny(model; spins = 1.5, g = 2, mode = :dipole)
sys = to_sunny(model; spins = Dict("Fe" => 3/2, "O" => 1))   # per-species
```

## Two placements: exact vs. unfolded

```julia
to_sunny(model; spins, placement = :explicit)    # exact, on the training supercell
to_sunny(model; spins, placement = :primitive)   # unfolded onto the chemical primitive cell
to_sunny(model; spins, placement = :auto)        # primitive if clean, else fall back (default)
```

- **`:explicit`** places the couplings directly on the training supercell. The classical
  energy reproduces `predict_energy − j0` exactly, but the magnon dispersion is folded
  into the supercell Brillouin zone.
- **`:primitive`** recovers the chemical primitive cell from the space group's pure
  translations, groups supercell atoms into sublattices, and places one Sunny bond per
  *primitive* bond (no multiplicity — Sunny's periodic replication restores it), giving the
  *unfolded* dispersion. This is exact only when the model genuinely lives on the primitive
  cell; a `clean` check detects when it does not (the interaction range reaches the
  supercell boundary) and falls back to the exact supercell route.
- **`:auto`** (the default) uses `:primitive` when clean, otherwise `:explicit`.

## How it is validated

The numerically load-bearing conversion math — the harmonic-to-Cartesian matrices, the
directed-member folding, the primitive unfold — lives in a **dependency-free core layer**
and is gated *without Sunny* by reconstructing the energy
(`_reconstruct_energy ≈ predict_energy − j0`). The thin extension only assembles the
`Sunny.System` from the core-computed matrices, and a separate test environment checks the
real `Sunny.System` energy against the SCE energy for both placements. The rationale is in
[Architecture](../theory/architecture.md#The-core/extension-split).
