# Building the basis

```@meta
CurrentModule = SCEFitting
```

The first half of the workflow turns a crystal and an interaction specification into a
symmetry-adapted SALC basis — the fixed set of invariants ``\Phi_\varphi`` whose
coefficients the fit will recover. This page covers the four ingredients:
[`Crystal`](@ref), [`BasisSpec`](@ref), the periodic image selection, and the symmetry
backend, all assembled by [`SCEBasis`](@ref).

## Geometry: `Lattice` and `Crystal`

A [`Lattice`](@ref) holds the cell vectors as columns; a [`Crystal`](@ref) adds
fractional atomic positions, integer species labels, and per-species names.

```@example basis
using SCEFitting
import Spglib

lat = Lattice([3.0 0 0; 0 3.0 0; 0 0 3.0])           # columns are the cell vectors aᵢ
cr  = Crystal(lat,
              [0.0 0.5; 0.0 0.5; 0.0 0.5],            # 3 × n_atoms fractional positions
              [1, 1],                                 # species index per atom
              ["Fe"])                                 # label per species
(n_atoms(cr), cr.species_labels)
```

Periodic directions default to all three axes; pass `pbc = (true, true, false)` to a
`Lattice` for a slab. Fractional positions are wrapped into `[0, 1)` — a precondition the
neighbor list relies on.

## The interaction specification

An [`BasisSpec`](@ref) fixes the body order, the pair cutoff, the per-species maximum
angular momentum `l`, and whether to keep only the isotropic channel:

```@example basis
interaction = BasisSpec(; nbody = 2, pair_cutoff = 2.7, lmax = [1], isotropy = true)
```

- `nbody` — maximum cluster size. `2` is pairwise; `nbody = 3` adds triplets, and so on
  (see [Body order beyond pairs](#Body-order-beyond-pairs)).
- `pair_cutoff` — the radius (Å) within which two atoms form an edge. It may be `Inf`,
  meaning *every resolvable pair* — the whole Wigner–Seitz cell of the supercell.
- `lmax` — the maximum harmonic degree per species. `[1]` is dipolar (the only `l` that
  builds the Heisenberg/DM/Γ bilinears); `[2]` adds quadrupolar single-ion and biquadratic
  channels.
- `isotropy` — `true` keeps only the rotation-invariant `Lf = 0` channel (e.g. Heisenberg
  `eᵢ·eⱼ`); `false` keeps the anisotropic channels too.

## Periodic resolvability: which pairs are physical

A plain-PBC supercell can only resolve interactions whose displacement lies in its
**Wigner–Seitz cell** — a farther periodic image of an atom carries the *same spin*, so its
interaction is an alias of the minimum-image one and cannot be fit independently. The
[`MinimumImage`](@ref) selection (the default) enumerates exactly this set, with
boundary ties kept as distinct members; [`AllImages`](@ref) keeps every image within the
cutoff and is reserved for the (future) spin-spiral / generalized-Bloch path.

```@example basis
nl = build_neighbor_list(cr, Inf, MinimumImage())     # full Wigner–Seitz cell
length([p for p in nl.pairs if (p.i, p.j) == (1, 2)])  # the 8-fold body-diagonal corner tie
```

This is a load-bearing distinction with subtleties at the cell boundary (and for `N ≥ 3`
clusters); it has its own chapter, [Periodic resolvability](../theory/resolvability.md).
You rarely call `build_neighbor_list` directly — `SCEBasis` threads the selection through
for you — but `pair_cutoff = Inf` in the `BasisSpec` is how you ask for the whole cell.

## Symmetry backends

Symmetry collapses the raw cluster instances into orbits and projects the harmonic
products onto the space-group-invariant SALCs. The backend is pluggable:

- [`NoSymmetry`](@ref) — the identity only (in-tree, no dependency). Every directed cluster
  is its own orbit, so the basis is larger but still correct.
- [`SpglibBackend`](@ref) — real space-group symmetry via [Spglib](https://github.com/singularitti/Spglib.jl).
  Provided by an extension: `import Spglib` activates it.

```@example basis
basis = SCEBasis(cr, interaction; backend = SpglibBackend())
(basis.spacegroup.symbol, basis.spacegroup.number, n_salcs(basis))
```

!!! tip "Loading Spglib"
    Use `import Spglib`, not `using Spglib` — Spglib exports `Lattice` / `Crystal`
    names that would clash with this package's. Importing only activates the extension.

## Body order beyond pairs

`nbody = K` enumerates clusters of up to `K` atoms as *pairwise-within-cutoff cliques* of
distinct atoms. At `N ≥ 3` something genuinely new happens: a space-group operation can
permute equivalent sites, mixing the coupling paths and — for unequal `l` — the
`l`-orderings. The SALC projection runs over the combined (ordering × coupling-path ×
`Lf`) space, so a single SALC can carry several `l`-orderings as *terms*. The
[kagome three-body tutorial](../tutorials/kagome_threebody.md) shows such a multi-term
channel; the mechanism is described in [Architecture](../theory/architecture.md).

## What `SCEBasis` produces

[`SCEBasis`](@ref) bundles the crystal, the space group, the SALC basis, and the
interaction. Each basis function is a [`SALC`](@ref) addressed by a canonical
[`SALCKey`](@ref) — a stable identity (body order, orbit, sorted `l`-multiset, `Lf`,
block) that names a fixed interaction independent of construction order. The design-matrix
columns, the persisted coefficients, and [`coeftable`](@ref) all key off it.

You can also build a basis from a human-authored `input.toml` instead of constructing the
objects in Julia — see [Persistence and I/O](io.md).

Next: [Data and fitting](fitting.md).
