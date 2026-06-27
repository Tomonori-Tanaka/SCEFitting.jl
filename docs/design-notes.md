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
representative cluster's coefficient onto the trivial irrep of its **site
stabilizer**, gauge-fix, then transport to the orbit members. The invariant count
is identical, but the projectors are small and the cluster model stays
primitive-cell + integer `R`.

For **arbitrary body order** the projection runs over the combined
`(ordering o, coupling path p, Mf)` space (the next section explains why). The
key construction choice: each stabilizer operation acts by rotating every **site**
axis by `wignerD_real(l_i, R)` (full `R`) and relabeling axes by the induced
permutation, and the action matrix is read off by **contracting against the
package's own orthonormal coupled tensors** — no 6j/9j recoupling, in the same
spirit as building `wignerD_real` itself. Because the orthonormal coupled tensors
make the action's matrix exact, `P = mean_g M_g` is a true projector and the
eigenvalues come out exactly `0`/`1` (an in-band idempotency assertion guards this).

A subtlety found earlier (an adversarial review) is **parity**: an *improper* `R`
carries the harmonic parity `(−1)^l` per site, while the time-reversal-even sector
has total product parity `(−1)^{Σl} = +1`. A first version that represented an op on
the final multiplet directly by `wignerD_real(Lf, R)` disagreed for odd `Lf` and
broke invariance on centrosymmetric cells (it needed an explicit proper-part fix).
The site-axis construction above **handles this automatically**: rotating each site
axis by the full `R` carries the genuine `(−1)^{l_i}` per site, and even `Σl` makes
the net improper action correct — so symmetry-forbidden odd-`Lf` channels (e.g. a
Dzyaloshinskii–Moriya-like term on a bond with an inversion center) drop out with no
special case. The ground-truth gate is the invariance test `Φ(g·e) = Φ(e)` with
*non-collinear* spins for all `Lf` and all body orders (collinear test spins make
every odd-`Lf` invariant vanish and would silently hide such bugs), plus a
linear-independence check (design-matrix rank = #SALC).

## 3b. Combined-space projection and multi-term SALCs (`N ≥ 3`)

At `N ≥ 3` a stabilizer operation can **permute symmetry-equivalent sites**. That
permutation does two things a per-multiplet projector cannot express:

1. It mixes **coupling paths** — different intermediate `L` sequences that share the
   same final `Lf` — so the invariant lives in a span of paths, folded into one
   tensor (this already bites equal-`l` channels like `(2,2,2)`).
2. When the permuted sites carry **unequal `l`** (e.g. `(1,1,2)` on an equilateral
   triangle), it mixes **`l`-orderings**: the invariant is a combination of
   `(1,1,2)`, `(1,2,1)`, `(2,1,1)` that no single ordering contains.

So a SALC is a sum over orderings, each with its own per-site `ls` — a
**multi-term** object (`SALCTerm`); `evaluate` and the torque gradient loop over
terms. When the site permutations are a *proper* subgroup of `Sₙ`, a degenerate
multiset can split into more than one ordering orbit (e.g. a mirror-only triangle
puts `l = 2` on the apex vs. on a base site — two distinct SALCs that share the
sorted label `[1,1,2]`); `SALCKey` keeps `block` running across them so it stays
injective (the by-key column contract, §4).

This was **cross-validated against Magesty**: for a kagome `P6/mmm` cell the number
of independent invariants per `(body, ls, Lf)` channel — a convention-free integer —
matches exactly through 3-body (all 42 channels, including the multi-term and
multi-path ones), and Magesty's own SALCs independently pass the invariance test. Two
from-scratch constructions agreeing on every invariant-subspace dimension is strong
mutual evidence that both are correct.

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

## 6. Torque as the energy surface's exact derivative

The torque `τ_a = e_a × ∂E/∂e_a` is the SCE's second observable. Rather than
treat the torque design matrix as an independently-derived object (a place a sign
or normalization can silently drift out of step with the energy kernel), the
rebuild builds the per-site gradient `accumulate_grad!` from the *same*
`μ = idx − ls − 1` mapping, `ls`, `folded`, and `(4π)^(N/2)` scale as the energy
kernel `evaluate` — only the innermost `Zₗᵢμᵢ` of each product is swapped for
`∇Zₗᵢμᵢ` (product rule, summed over the cluster's sites). The consequence is that
`predict_torque` is *by construction* the analytic gradient of the surface
`predict_energy` evaluates, and the gate is exactly that: an on-sphere
finite-difference check `predict_torque ≈ e × ∇E_FD`, plus the Heisenberg closed
form `τ_a = J e_a × Σ_{nn} e_j`. This is convention-independent ground truth — it
would catch a self-consistent wrong convention that every gauge-invariant
comparison misses. (Magesty's reverse-mode, cluster-major gradient accumulator is
an optimization the v0 deliberately forgoes: the simple `O(N²)`-per-multi-index
leave-one-out product is clearer and fast enough for 2-body, and shares its inner
loop with `evaluate` so the two cannot diverge.)

The co-fit follows Magesty: minimize `L = (1−w)·MSE_E + w·MSE_T` by whitening each
block (`√((1−w)/n_E)`, `√(w/n_T)`), with `j0` profiled out of the energy block
analytically — the torque block carries no intercept, so it never sees `j0`.

## 7. Persistence and input files (TOML, not XML/JSON)

Magesty persists the basis/model as **XML** (EzXML) and takes run setup from an
**`input.toml`**. The rebuild keeps the two-artifact split but unifies on one
zero-dependency format:

- **Input** (`sce/input.jl`): a human-authored `input.toml` — `[structure]` (inline
  crystal: `lattice` as a list of the three lattice vectors, fractional `positions`,
  per-atom `species`, `species_labels`), `[interaction]` (`nbody`, `pair_cutoff`,
  per-species `lmax`, `isotropy`), and optional `[symmetry]` (`backend`, `tol`).
  `SCEBasis("input.toml")` builds the basis; keyword arguments override the file's
  backend/tol. Like Magesty, training data and the estimator are kept **out** of this
  file (loaded / chosen in Julia) — `input.toml` specifies only what defines the basis.
- **Persistence** (`sce/persist.jl`): a **self-contained** TOML document — the crystal,
  the space-group ops, the interaction, and the *full* SALC basis (every member, term,
  and folded tensor), plus (for a model) `j0` and per-`SALCKey` coefficients. Reload
  reconstructs the basis **verbatim** (no re-projection), so a saved model keeps working
  even if the construction code's gauge later changes — the file is the ground truth.
  Coefficients re-pair to the basis **by key** (§4), not by position; a scrambled
  on-disk coefficient order still reloads correctly.

Why TOML for both (the format was deliberated): the persistence artifact is a deep
nested numeric structure, which first suggested JSON. But Julia's **stdlib** `TOML`
round-trips `Float64` *exactly* (verified, including `-0.0` and `nextfloat(1.0)`) and
expresses the full nested SALC document, so a single zero-dependency format serves
both the human-written input and the machine-written dump — the most Julia-native
outcome. The `struct ⇄ Dict` schema layer (`_to_doc` / `_from_doc`) is kept
**format-agnostic** and is unit-tested with no serializer at all, so a JSON (or other)
backend would be a thin future addition rather than a rewrite. Two robustness details:
the `UInt64` structural fingerprint is stored as a *string* (TOML/JSON numbers lose
precision past `2^53`) and is **recomputed** on load rather than trusted (`hash` is
Julia-version dependent); and `-0.0` is normalized to `+0.0` on write so two builds of
the same object serialize byte-identically.

## 8. Oracle methodology

`test/oracle/` is a separate environment that `dev`s a pinned `Magesty.jl`; the
core suite never depends on Magesty. The oracle compares only **convention-fixed
kernels** bit-for-bit (`Zₗₘ`, Clebsch–Gordan vs. `WignerSymbols`, real Wigner-D,
coupled tensors) and **gauge-invariant aggregates** (space-group/spacegroup
numbers, orbit counts, held-out predictions, the recovered Heisenberg `J`). Raw
SALC coefficients and individual design-matrix columns are *never* compared — they
are gauge-dependent. Absolute conventions (signs, normalization) are pinned
independently by closed forms and finite differences, because a self-consistent
wrong convention would pass every gauge-invariant comparison.
