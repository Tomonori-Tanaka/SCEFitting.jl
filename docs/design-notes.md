# Design notes

Why `MagestyRebuild` is built the way it is — the deliberate refinements over
`Magesty.jl`. Bit-for-bit agreement with Magesty is **not** a goal; these are
cleaner constructions that are validated to be *physically* correct (symmetry
invariance, finite differences, known-coupling recovery), and may differ from
Magesty by a gauge or by rounding.

## 1. Generalized cutoff neighbor list (vs. a fixed 27-cell grid)

Magesty enumerates neighbors over a hard-coded 3×3×3 image grid, which silently
caps any interaction range at one lattice translation and stores a dense
`O(n²)` `min_distance_pairs` matrix. The rebuild (`geometry/neighborlist.jl`)
drives the image range by the cutoff: along each periodic axis,
`N_d = ceil(cutoff · ‖b_d‖)` with `b_d` the reciprocal row, which the projection
bound `|n_d + Δf_d| ≤ cutoff/spacingᵈ` proves encloses every in-cutoff image —
correct for triclinic cells and arbitrary cutoffs, and a sparse `Vector` of pairs.
Each `NeighborPair` keeps the integer lattice translation `R` (`shift`), which the
orbit reduction propagates; this is the hook for future reciprocal-space /
generalized-Bloch `E(q)` training data. Validated against an *independent*
over-large-shell brute force (the production tight range is never checked against
itself), including a strongly sheared triclinic cell.

## 2. Real Wigner-D from the package's own `Zₗₘ` (vs. Euler angles + complex D + c2r)

Magesty builds the real Wigner-D as `conj(C)·wignerD(l,α,β,γ)·transpose(C)` —
Euler-angle decomposition, the complex Wigner-D (a dependency), and the complex→
real matrix `C`. The rebuild (`AngularMomentum.wignerD_real`) instead **fits** the
matrix directly: sample the sphere (Fibonacci points), require
`Zₗ,m(R·r̂) = Σ_{m'} Δ[m,m'] Zₗ,m'(r̂)`, and solve the (exact, oversampled) linear
system. Advantages: it is consistent with *our* `Zₗₘ` convention by construction
(no transpose/phase bookkeeping), uses real arithmetic only, needs no `WignerD`
dependency, and handles improper rotations (reflections, inversion) with no special
case. Validated by the functional identity on fresh points, orthogonality,
`l=0/1/inversion` special cases, and agreement with Magesty's `Δl`.

## 3. Orbit–stabilizer SALC projection (vs. a whole-supercell Reynolds average)

Magesty averages over the whole space group acting on explicit supercell atoms.
The rebuild (`basis/salcbasis.jl`) uses the orbit–stabilizer theorem: project the
representative cluster's coefficient tensor onto the trivial irrep of its **site
stabilizer** (small `(2Lf+1)×(2Lf+1)` matrices), gauge-fix, then transport to the
orbit members. The invariant count is identical, but the projectors are tiny and
the cluster model stays primitive-cell + integer `R`.

The subtle point — found by an adversarial review — is parity. The per-`Lf`
projector represents an operation on the final multiplet by `wignerD_real(Lf, R)`,
whose value for an *improper* `R` carries the genuine harmonic parity `(−1)^Lf`.
But the physical object is a product of `N` site harmonics with total parity
`(−1)^{Σl} = +1` in the time-reversal-even (even-`Σl`) sector. These agree for even
`Lf` and disagree for odd `Lf`, which broke invariance for anisotropic channels on
centrosymmetric cells. The fix: represent each op by its **proper part**,
`wignerD_real(Lf, is_proper ? R : −R)` (matching `det = +1`); this correctly drops
symmetry-forbidden odd-`Lf` channels (e.g. a Dzyaloshinskii–Moriya-like term on a
bond with an inversion center) and keeps the genuinely allowed ones. The
ground-truth gate is the invariance test `Φ(g·e) = Φ(e)` with *non-collinear*
spins for all `Lf` — note that collinear test spins make every odd-`Lf` invariant
vanish and silently hide such bugs.

## 4. Canonical `SALCKey` column addressing (vs. implicit construction order)

In Magesty the design-matrix column order is the SALC iteration order, an implicit
"load-bearing" coupling that silently invalidates saved models if reordered. The
rebuild gives every SALC a structural `SALCKey` `(body, orbit_id, sorted ls, Lf,
block)`; the basis is sorted by key and the keys travel with the coefficients
(`SALCBasis.keys`, a gauge-independent `fingerprint`). A model re-pairs `jphi` to a
freshly built basis *by key*, not by position. (`orbit_id` is a stable ordinal
within a fixed `(crystal, spacegroup, cutoff)`, protected by the fingerprint.)

## 5. Pluggable backends via package extensions

Symmetry analysis and (future) regularized estimators are behind abstract types
(`AbstractSymmetryBackend`, `AbstractEstimator`). The *types* live in the core
(`NoSymmetry`/`SpglibBackend`, `OLS`/`Ridge`), but heavy methods live in
extensions: `analyze_symmetry(::SpglibBackend, …)` is in `ext/` and lights up only
on `import Spglib`. The core thus loads with no heavy dependencies, and calling an
unloaded backend hits a friendly "load the backend" error on the abstract
supertype rather than a bare `MethodError`. This also removes Magesty's habit of
force-loading Spglib at `using` time.

## 6. Oracle methodology

`test/oracle/` is a separate environment that `dev`s a pinned `Magesty.jl`; the
core suite never depends on Magesty. The oracle compares only **convention-fixed
kernels** bit-for-bit (`Zₗₘ`, Clebsch–Gordan vs. `WignerSymbols`, real Wigner-D,
coupled tensors) and **gauge-invariant aggregates** (space-group/spacegroup
numbers, orbit counts, held-out predictions, the recovered Heisenberg `J`). Raw
SALC coefficients and individual design-matrix columns are *never* compared — they
are gauge-dependent. Absolute conventions (signs, normalization) are pinned
independently by closed forms and finite differences, because a self-consistent
wrong convention would pass every gauge-invariant comparison.
