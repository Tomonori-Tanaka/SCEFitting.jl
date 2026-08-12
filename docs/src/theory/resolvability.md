# Periodic resolvability

```@meta
CurrentModule = SCEFitting
```

A spherical cutoff is the wrong primitive for *what a finite supercell can resolve*. This
chapter explains why the default selection enumerates the **Wigner–Seitz cell**, what the
boundary ties are, and how the same logic extends to three- and four-body clusters.

## Only minimum-image interactions are resolvable

Under plain periodic boundary conditions (ordinary DFT, excluding the generalized-Bloch
spin-spiral case), every periodic image of an atom carries the **same spin**. So if atom
``A`` reaches both the minimum image of ``B`` (at distance ``d``) and a farther image of
``B`` (at ``d' > d``), the two "interactions" are the *same* ``\hat{\boldsymbol e}_A \cdot
\hat{\boldsymbol e}_B`` in every training configuration — perfectly collinear design-matrix
columns. A supercell can only independently resolve the **minimum-image** set: the
displacements lying inside the Wigner–Seitz cell of the (super)lattice.

## The Wigner–Seitz cell is a polyhedron, not a ball

For a cubic cell of side ``L`` the WS cell is the cube ``[-L/2, L/2]^3``. Its inscribed
sphere has radius ``L/2`` (the face centers), but the farthest resolvable displacement is
the body-diagonal corner ``(L/2, L/2, L/2)`` at ``\sqrt 3\,L/2 \approx 0.866\,L``. A
spherical cutoff therefore *cannot* express "all resolvable pairs": to reach the corners it
must exceed ``L/2`` along the faces — exactly where it starts sweeping in the aliased
(non-resolvable) images. This face-vs-corner mismatch is the classic source of
``L/2``-plane double-counting errors.

So `cutoff` is not a sphere radius applied blindly. The default [`MinimumImage`](@ref)
selection keeps, per atom pair, only the minimum-image displacement(s), trimmed by the
radial cutoff; `cutoff = Inf` keeps the whole WS cell (Magesty spells this
`cutoff = -1`).

## Boundary ties

On the WS boundary several images are exactly equidistant, and they are kept as
**distinct members** of one cluster orbit. In the infinite crystal they *are*
different bonds; on the finite reference cell they connect the **same** pair of
reference-cell atoms, so their design-matrix contributions are evaluated on the same
spin data and need not be independent:

| WS boundary feature | Multiplicity |
|---------------------|-------------:|
| face (``L/2``) | 2 |
| edge | 4 |
| corner ``(L/2,L/2,L/2)`` | 8 |

For a body-centered cubic cell, the cross pair sits at all eight corners at once:

```@example res
using SCEFitting
cr = Crystal(Lattice([3.0 0 0; 0 3.0 0; 0 0 3.0]), [0.0 0.5; 0.0 0.5; 0.0 0.5], [1, 1], ["Fe"])
nl = SCEFitting.build_neighbor_list(cr, Inf, MinimumImage())   # public-unexported: qualify
length([p for p in nl.pairs if (p.i, p.j) == (1, 2)])     # 8 equidistant corner images
```

Self-pairs (``i = i + R``) are dropped: both ends share ``\hat{\boldsymbol e}_i``, so the
term is a constant (``L_f = 0``) or a one-body alias (``L_f > 0``), never an independent
pair. Likewise an ``N``-body cluster must use distinct atoms.

### What a tie costs on a finite cell

Because the tied images connect the same reference-cell atom pair, the orbit sum over
a tie can make some of an orbit's SALC content **identically zero as a function of
cell-periodic spin data** — the column is not merely small, it vanishes for every
configuration the reference cell can express. The content is *unidentifiable from
this cell*, not absent from the physics: a fit sees a zero design column and the
coefficient is unconstrained. A tie is *necessary* for this (with a unique minimum
image every member has its own atom content and nothing can cancel), but which
content dies depends on how the space group relates the tied members, not on the tie
alone. The pattern to expect from the simplest two-fold tie (a pair at half a lattice
vector) is bond reversal ``\hat{\boldsymbol r} \to -\hat{\boldsymbol r}``, which an
invariant of bond rank ``L_f`` carries as ``(-1)^{L_f}`` while the site factors are
unchanged — the odd-``L_f`` content of such a channel goes, the antisymmetric
DMI-like ``L_f = 1`` one included. Treat that as a guide to *where to look*, not as a
rule: a larger tie can remove even-``L_f`` content too (measured upstream on the
eight-fold bcc corner tie with ``l \le 2`` spin factors: 4 of the orbit's 7 SALC
columns identically zero), and a smaller stabilizer can leave odd content standing.

The remedy is a reference cell that breaks the tie **in every direction whose
displacement carries a half-lattice-vector component** — doubling one axis alone
moves only the ties along that axis. (The joint-family SLCE.jl classifies and
freezes such columns mechanically — `unresolvable_columns` there; this package
documents the phenomenon and leaves the columns in place, so watch for exactly-zero
design columns on high-symmetry cells.)
[Adapted from SLCE.jl de79b92/3e68fc1.]

## The third edge: compact clusters at `N ≥ 3`

For pairs the WS boundary only multiplies tied images. For clusters of three or more atoms
it adds a genuinely new constraint. A triangle ``\{i, j, k\}`` is admissible only if **all
three edges sit at their atom-pair minimum image simultaneously** — the *compact-cluster*
criterion. Having each pair individually minimum-image-resolvable is **not** enough: the
images that make ``i\!-\!j`` and ``i\!-\!k`` minimal may force ``j\!-\!k`` onto a longer,
non-minimum image, which must reject the cluster.

A sharp illustration: three atoms equally spaced around a one-dimensional ring have every
pair minimum-image at the same distance, yet admit *zero* compact triangles — you cannot
realize all three minimal edges at once, just as an equilateral triangle does not embed on
a ring. The enumeration checks every pair of a candidate clique on its *actual chosen
images*, so it counts exactly the resolvable compact clusters, with the boundary ties
multiplying them as they do pairs.

This combinatorics is easy to get subtly wrong, so the full count — the candidate set for
``N = 2, 3, 4`` (including the whole WS cell) and the symmetry-orbit partition — is pinned
against an independent brute-force enumeration on cells deliberately seeded with face,
edge, and corner ties.

## The spin-spiral seam

[`AllImages`](@ref) keeps *every* image within the cutoff, each tagged with its lattice
translation ``R``. For plain-PBC fitting this over-counts (it admits aliases of shorter
bonds), so it is not the default — but it is the representation a **generalized-Bloch /
spin-spiral** extension needs, where the phase ``e^{i\boldsymbol q\cdot\boldsymbol R}``
distinguishes images that a single supercell cannot. This is why [`NeighborPair`](@ref) and
the cluster members retain ``R``. Below half the smallest perpendicular cell width,
``\text{cutoff} < \min_d d_i/2``, the two selections coincide (each in-cutoff image is
already the minimum).
