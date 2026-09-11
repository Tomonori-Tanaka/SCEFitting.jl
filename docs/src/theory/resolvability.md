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

```@example resolvability
using SCEFitting
crystal = Crystal(Lattice([3.0 0 0; 0 3.0 0; 0 0 3.0]),
                  [0.0 0.5; 0.0 0.5; 0.0 0.5], [1, 1], ["Fe"])
# public-unexported: qualify
neighbors = SCEFitting.build_neighbor_list(crystal, Inf, MinimumImage())
length([p for p in neighbors.pairs if (p.i, p.j) == (1, 2)])     # 8 equidistant corner images
```

Self-pairs (``i = i + R``) are dropped: both ends share ``\hat{\boldsymbol e}_i``, so the
term is a constant (``L_f = 0``) or a one-body alias (``L_f > 0``), never an independent
pair. Likewise an ``N``-body cluster must use distinct atoms.

### What a tie costs on a finite cell

Because the tied images connect the same reference-cell atom pair, the orbit sum over
a tie can make some of an orbit's SALC content **identically zero as a function of
cell-periodic spin data** — the column is not merely small, it vanishes for every
configuration the reference cell can express. The content is *unidentifiable from
this cell*, not absent from the physics: without intervention a fit would see a zero
(or linearly dependent) design column and the coefficient would be unconstrained. A
tie is *necessary* for this (with a unique minimum image every member has its own
atom content and nothing can cancel), but which content dies depends on how the
space group relates the tied members, not on the tie alone. The pattern to expect
from the simplest two-fold tie (a pair at half a lattice vector) is bond reversal
``\hat{\boldsymbol r} \to -\hat{\boldsymbol r}``, which an invariant of bond rank
``L_f`` carries as ``(-1)^{L_f}`` while the site factors are unchanged — the
odd-``L_f`` content of such a channel goes, the antisymmetric DMI-like ``L_f = 1``
one included. Treat that as a guide to *where to look*, not as a rule: a larger tie
can remove even-``L_f`` content too (measured upstream on the eight-fold bcc corner
tie with ``l \le 2`` spin factors: 4 of the orbit's 7 SALC columns identically
zero), and a smaller stabilizer can leave odd content standing.

The basis builder handles this **exactly, at construction**: each orbit's SALCs are
expanded into their aggregated (shift-blind) monomial coefficients, and combinations
that aggregate to zero — or become linearly dependent within the orbit — are dropped
with a warning naming the orbit, channel, and reason. Per orbit (and for members
with all-distinct atoms) the surviving functions are linearly independent by
construction, so the fit and every coefficient-level readout stay well-posed; the
model space is unchanged (the drop is exact). Dependence **across** orbits — tied
images the point group does not relate, sitting in separate orbits — is a different
face of the same tie and is handled separately: see
[When symmetry does not fuse the tie](@ref) below. Measured on bulk MnTe with SOC (3×3×3 supercell, ``P6_3/mmc``): 51 raw SALCs
of which 14 aggregate to zero, previously reaching OLS as a silently rank-deficient
design with `max|coef| ~ 1e7`. The remedy for the *physics* (recovering the dropped
content) is still a reference cell that breaks the tie **in every direction whose
displacement carries a half-lattice-vector component** — doubling one axis alone
moves only the ties along that axis. (The joint-family SLCE.jl classifies and
freezes such columns at readout — `unresolvable_columns` there; this package removes
them from the basis at build time.)
[Adapted from SLCE.jl de79b92/3e68fc1; reduction added with the MnTe fix.]

## When symmetry does not fuse the tie

Everything above needed the point group to permute the tied images, which is what
puts them in **one** orbit whose sum weights them equally. When no operation relates
the two equidistant images they sit in **different orbits** carrying independent
couplings. Nothing cancels and no column is zero. What collapses is the **span**:
under cell-periodic evaluation every image of an atom carries the same spin, so the
two orbits are the same function of anything this cell can express, and the data fix
only how much of that function is used in **total** — never how the total divides
between the two couplings.

This is not a low-symmetry curiosity. On the conventional Nd₂Fe₁₄B cell (68 atoms,
``P4_2/mnm``, 16 operations) with isotropic pairs enumerated out to the Wigner–Seitz
boundary, **ten** pairs of distinct pair orbits alias in exactly this way: each joins
the same atom pair at the same distance (6.2–8.7 Å) through two images that differ
by ``c`` (nine of them) or by ``a + b`` (one), and the group contains no operation
that swaps the "+" and "−" image at that atom pair. Their design columns coincide,
the torque design loses ten of its 179 ranks, and every penalized estimator returns
some split of each pair's sum without saying so.

### How this package handles it: one column, an equal split by convention

The basis builder detects these **cross-orbit alias groups** structurally, once, at
construction (no training data involved): within one channel
``(N, \{l\}, L_S, L_f)``, orbits whose members join the same atom sets are compared
through the same aggregated shift-blind function vectors the per-orbit reduction
uses, and orbits whose vectors are proportional (relative residual below
`alias_rtol`, default `1e-6`) form a group. The group is then

- **tied in the fit**: its SALC columns are summed into one design column
  ([`SCEFitting.n_columns`](@ref) counts those), with per-bond weights that are
  exactly ``\pm 1`` for transported copies of one representative tensor, so the
  data determine the group's sum and nothing else is asked of them;
- **read back with an equal per-bond split**: every orbit of the group receives the
  tied column's coefficient (times its weight). Summing the columns and assigning
  the tied coefficient to every member is what makes the split equal *per bond*
  even when the group's orbits have different multiplicities — a representative
  column with a ``J/k`` read-out would only be right for equal multiplicities;
- **disclosed everywhere the number can leave**: the [`coeftable`](@ref) columns
  `alias_group` and `split` (`:convention` on such coefficients), a `split` entry
  per coupling in the saved model (schema v7; older files load as `:legacy`), the
  build-time report naming the orbits, atom pairs and image offsets, and a
  one-time warning from [`multipole_terms`](@ref) / [`bilinear_terms`](@ref) /
  [`to_sunny`](@ref);
- **merged in the cost-aware selection**: [`SCEFitting.salc_groups`](@ref) gives the
  orbits of a group one label, so a group-sparse estimator can only keep or drop
  them together — a selection that kept one orbit and dropped the other would be
  deciding the split by itself.

The equal split is a **convention, not a measurement**: the training cell samples
``J(\boldsymbol q)`` only on its commensurate mesh, and the tied images alias exactly
there. Dividing the sum equally is the minimum-bandwidth interpolation between those
samples (the Nyquist term shared half-and-half between ``\pm\pi``, as phonopy
distributes a supercell force constant over equidistant images). What it leaves
**exactly invariant**: every training-cell energy and torque, every ``\boldsymbol q =
0`` quantity (``J(0)``, the mean-field ``T_c``, any configuration periodic in the
training cell), and the isotropic stiffness ``\sum_R J(R)\,|\boldsymbol d|^2``.
What **depends on it**: ``J(\boldsymbol q)`` off ``\Gamma`` — the magnon dispersion
along the tied direction is linear in ``J_+ - J_-`` at ``\boldsymbol q\cdot\Delta
\boldsymbol R = \pi`` —, the energy of a tiled supercell in a configuration that is
not periodic in the training cell (so the ordering temperature and the ordered state
of a Monte-Carlo run on the tiled cell), the off-diagonal stiffness, and any
published bond-resolved ``J(R)`` table.

The remedy that turns the convention into a measurement is the same as for the fused
tie: a training cell doubled along **every** axis the group's image offsets reach
(the build-time report lists them; `alias_groups(basis)[k].delta_shifts` carries
them), or — for isotropic channels only — spin-spiral training data at a wavevector
with ``\boldsymbol q\cdot\Delta\boldsymbol R \notin 2\pi\mathbb Z``, which this package
does not yet consume.

Two things the detection does **not** do. A class of orbits whose functions are
linearly dependent without being pairwise proportional — an anisotropic channel
whose member tensors carry their own bond geometry — has no equal split to define,
so it is reported (`kind == :span_collapsed`) and left untied, as is a proportional
set whose per-member tensor norms differ (`:unequal_norm`, the premise behind the
``\pm 1`` weights failing); the `OLS` rank warning still applies there and the
remedy is the doubled cell. A tie the point group fuses only *partially* (one orbit
of two images plus two single-image orbits, say) is handled: the tied column sums
every member bond, so the read-out is the bond-weighted mean. And nothing is merged or dropped in the basis itself: the orbits, their
`SALCKey`s and members stay distinct, so a datum type that breaks the periodicity
could separate them on the same basis object.

`alias_rtol = nothing` turns the detection off (every SALC its own column; the
design is then rank deficient on such a cell and `OLS` warns as before).

The joint-family SLCE.jl takes the opposite decision on the same face: it freezes
every column of every orbit that shares an atom set with another orbit
(`unresolvable_columns` there), discarding the determined sum as well so that a fit
to data containing that shell fails loudly. That is the right choice when the
deliverable is a dispersion the data never constrained; here the deliverable is a
model that has to be tiled and sampled, where a dropped shell is a missing coupling,
so the sum is kept and the split is recorded as what it is.

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
