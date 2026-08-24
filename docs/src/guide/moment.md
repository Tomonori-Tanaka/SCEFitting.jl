# Adiabatic site moments

```@meta
CurrentModule = SCEFitting
```

The energy expansion of the other guides treats the per-atom moment magnitude as
a fixed property of each site. Constrained-noncollinear DFT also returns the **bare
site moment** `M_a` of every configuration, and that quantity moves with the spin
environment: the adiabatic map `m_a(e)` — the moment a site carries once the
electrons relax around a given direction pattern `e` — is itself a smooth
function of the configuration. This page fits that map on a **pointed** (site-marked)
symmetry-adapted basis, so a Monte-Carlo or dynamics consumer can read a
configuration-dependent moment off the same kind of model it reads the energy from.

The objects are [`MomentSpec`](@ref) / [`MomentBasis`](@ref) (the pointed basis),
[`MomentDataset`](@ref) (rows, gate, resolvability), `fit(MomentFit, …)` with any
estimator, [`MomentModel`](@ref) / [`predict_moment`](@ref) (prediction), and a set
of diagnostics that ride on the same rows.

## The pointed expansion in one paragraph

The energy is a scalar invariant of the whole cluster; a site moment is a vector
attached to one site, so it is **covariant**: it rotates with the spins. The
expansion therefore marks one site `i` per term and reads the moment through its
projection on an evaluation axis `ê_i`,

```math
y_i(e) = \hat{\boldsymbol e}_i \cdot \boldsymbol m_i(e)
       = \sum_\varphi V_\varphi\, \Phi_\varphi(i;\, e, \hat{\boldsymbol e}_i),
```

where each `Φ_φ(i; ·)` is an ordinary SALC of the decor engine with the site `i`
**marked** — the mark is the displacement decor `SiteDecor(disp = (1, 0))`, so
projection, canonical members and evaluation are the unchanged machinery of
[Building the basis](basis.md) — and the marked site's own `ê` factor is the
decor's spin part. One coefficient vector `V` serves every symmetry-equivalent
site. A mark of rank `l = 0` gives the **per-orbit intercept** `μ₀` (the moment the
site carries in a featureless environment); higher mark ranks and environment
factors describe how the neighbours' directions inflate or quench it.

## Spec and basis

[`MomentSpec`](@ref) is the truncation. Two fields have no energy-side analogue and
are load-bearing:

- `sampled` — which species the downstream consumer actually samples. Every
  species whose spins enter the environment (`lmax_env > 0`) must be sampled,
  because a basis that reads an unsampled spin cannot be evaluated at run time.
  The spec refuses the inconsistent combination rather than fixing it up.
- `marked` — whose site moments are expanded (default: every species). A marked
  species may have `lmax_env = 0`: its own moment is modelled while its spins never
  act as environment for others.

The cutoffs are **mark-aware**: `cutoff_pair` bounds the mark–environment bond of a
2-body cluster, `cutoff_star` the mark bonds of a star (`nbody ≥ 3`), and the
environment–environment edge of a star is free. `isotropy = true` (the default)
keeps the `L_S = 0` blocks only — the adiabatic map is taken to be spin-rotation
covariant, exactly like an `isotropy = true` energy basis.

```@example moment
using SCEFitting, LinearAlgebra, Random
import Spglib                     # activates the SpglibBackend extension

# A 4-atom chain of Fe along z (8 Å of vacuum in x, y; nearest neighbours 2.5 Å apart)
lat = Lattice([8.0 0 0; 0 8.0 0; 0 0 10.0])
frac = [0 0 0 0; 0 0 0 0; 0.0 0.25 0.5 0.75]
chain = Crystal(lat, frac, [1, 1, 1, 1], ["Fe"])

spec = MomentSpec(; lmax_env = [1], sampled = [true], lmax_mark = 1,
                  nbody = 2, cutoff_pair = 2.6)
mb = MomentBasis(chain, spec; backend = SpglibBackend())
```

With `lmax_mark = 1` and `lmax_env = 1` the basis is the `l = 0` mark (the
intercept `μ₀`) plus the isotropic mark–neighbour pair `ê_i · e_j`. The marked
reference-cell atoms — the design's row atoms — are `mb.marked_atoms`.

The same truncation can be written as the `[moment]` section of the TOML setup file
that already holds the crystal and the energy basis, and built with
`MomentBasis("input.toml")` — see [the `[moment]` section](io.md#The-pointed-site-moment-basis)
in the I/O guide. Every key maps one-to-one onto the keywords above (`sampled` is
spelled as species labels there, `sampled = ["Fe"]`).

### Periodic resolvability of the pointed columns

A finite training cell cannot determine every pointed column: some cancel
identically on cell-periodic data, and some only appear in fixed combinations.
[`moment_resolvability`](@ref) is the structural (data-free) gate that names them:

```@example moment
res = moment_resolvability(mb)
(; rank = res.rank, vanishing = res.vanishing, n_dependent = length(res.null_combinations))
```

The dataset constructor runs this gate first, freezes the vanishing columns to
exact zero in every fit, and discloses the dependent combinations in its
construction warning. A full-rank result certifies that the **basis** is
resolvable on this cell — not that any particular training set identifies it.

## Data

A [`SpinDatum`](@ref) carries the moment channel's trio, each `nothing` when absent:

- `moments_bare` — the bare moment vectors `M_a` (`3 × n_atoms`, μ_B; VASP
  `M_int`), the projection target. Kept as a signed vector: it may pass through
  zero, which is exactly what makes the target analytic there.
- `constraint_mode` — how the source calculation constrained the spins: `4`
  (direction-pinning, `I_CONSTRAINED_M = 4`) or `1` (transverse penalty).
- `constraint_axes` — the penalty axes of a mode-1 calculation (`3 × n_atoms`,
  zero columns where a site was unconstrained).

The file readers of [Persistence and I/O](io.md) populate the trio from the
`mint` / `mconstr` columns (`read_extxyz`) or from an `M_int` file next to the
`MW_int` one (`read_embset_pair`), and re-verify the axes against the converged
moments through `SCEFitting.check_moment_gates`. Here the data are synthetic: a
planted map `m_a = (μ₀ + c Σ_{j∈nn(a)} e_a · e_j) e_a`, exactly decomposable along
the spin, read out in mode 4:

```@example moment
rng = MersenneTwister(1)
unit(v) = v / norm(v)
μ0, c = 2.2, 0.1
function planted_datum(rng)
    e = reduce(hcat, [unit(randn(rng, 3)) for _ in 1:4])
    M = similar(e)
    for a in 1:4
        s = dot(e[:, a], e[:, mod1(a - 1, 4)]) + dot(e[:, a], e[:, mod1(a + 1, 4)])
        M[:, a] = (μ0 + c * s) * e[:, a]
    end
    # `moments` (the smoothed MW) supplies directions and magnitudes; `moments_bare`
    # is the channel's target. No constraining field in this synthetic set.
    return SpinDatum(0.0, M, zeros(3, 4); moments_bare = M, constraint_mode = 4)
end
data = [planted_datum(rng) for _ in 1:40]
length(data)
```

## The dataset: rows, axis rule, gate

[`MomentDataset`](@ref) assembles one row per `(configuration, marked atom)` with
the target `y = ê · M`. The axis `ê` follows the **mode rule**: the spin direction
itself in mode 4, the constraint axis in mode 1 (a mode-1 atom with a zero axis
has no defined projection; its rows are excluded and reported). Modes may be mixed
in one dataset.

Rows pass a **decomposability gate**, `g = ‖M⊥‖²/|M| = |M| sin²θ ≤ gate_eps`: a
moment that is not parallel to its axis has a transverse remainder the projection
discards, and the tolerance is yours to state (`gate_eps` has no default). A
`coverage_floor` refuses a dataset in which too few rows of some marked orbit
survive. Two further doors: every `directions` column must be unit, and every atom
the basis references must carry a nonzero moment (`zero_moment_atol`; pass the
reader's own value if you changed it there).

```@example moment
ds = MomentDataset(mb, data; gate_eps = 1e-8)
```

The dataset keeps the per-row bookkeeping a consumer of the diagnostics needs:
`keep` / `gate` (the gate decision and statistic), `row_config` / `row_atom`,
`orbit_rep` and the per-orbit `orbit_report`, and `order` — the marked-sublattice
order parameter `|⟨e⟩|` of every configuration.

## Fitting and predicting

`fit(MomentFit, ds, estimator)` solves the gated rows and, for disclosure, the
ungated ones; both coefficient sets are stored. There is **no centering and no
global intercept**: the `l = 0` mark columns are the per-orbit intercepts, so the
design reaches the estimator exactly as built.

A regularized estimator built from the basis leaves those `μ₀` columns
**unpenalized**: [`penalty_metric`](@ref)`(mb)` gives them a scale of exactly `0`,
which every estimator reads as "do not penalize this column". Shrinking the reference
moment magnitude toward zero has no physical meaning, and the exemption is not a
special case in `fit` — it is the same per-column scale that makes the penalty
invariant under the basis's column conventions (see the
[fitting guide](fitting.md)). Pass `free_intercepts = false` to penalize them anyway,
or `metric = nothing` for the plain unweighted penalty.

```@example moment
f = fit(MomentFit, ds, OLS())
rmse_moment(f)
```

```@example moment
model = MomentModel(f)
e = data[1].directions
maximum(abs, predict_moment(model, e) .- ds.y[1:4])   # the first configuration's rows
```

[`predict_moment`](@ref) returns `ê_a · m_a(e)` for each marked atom; its
`axes = e` default is the mode-4 identity (the magnitude along each spin, what a
Monte-Carlo consumer asks for), and explicit axes reproduce mode-1 training rows.

For group-adaptive shrinkage, `salc_groups(mb)` labels the pointed columns at
**mark-class** granularity and `GroupAdaptiveRidge(mb; lambda)` builds the
matching estimator; `fit` reduces it to the active columns alongside the freeze.
`Ridge(mb; lambda)` and `AdaptiveRidge(mb; lambda)` are the ungrouped forms, and all
three carry the metric.

## Choosing λ

[`cross_validate`](@ref)`(ds, estimator)` is the honest criterion. Its folds are
grouped **by configuration**, so the rows of one configuration — one per marked atom,
all sharing its spin directions — never split across the train/holdout boundary:

```@example moment
est = Ridge(mb; lambda = 1e-4, metric_nconfig = 512)
cv = cross_validate(ds, est; nfolds = 3)
cv.pooled_rmse_moment
```

For a λ sweep, move that estimator along the path with
`SCEFitting.with_lambda(est, λ)`: it carries the penalty metric forward, where a hand
rebuild would drop it and a fresh `Ridge(mb; lambda = λ)` would rebuild the reference
ensemble at every point.

Each fold re-solves with the same frozen column set as the full dataset, and a fold
whose training rows miss a marked orbit entirely is refused by name — that orbit's
`μ₀` would be unidentified on the fold rather than merely noisy. Training and scoring
run on the gate-kept rows; `score_defined` reports the rejected rows separately, as
disclosure rather than as a criterion (their targets carry a transverse component no
coefficient can fit).

[`effective_dof`](@ref)`(f)` and [`gcv`](@ref)`(f)` are the fast reference for a
linear estimator, computed on the design the fit actually solved. GCV treats rows as
exchangeable, which the rows of one configuration are not, so it runs optimistic —
use it to scan, and cross-validate the shortlist.

## Diagnostics

All three ride on the dataset's rows, so they read the same gate and the same
marked atoms as the fit.

**Coverage-band residual profile.** Per-configuration mean residuals organized
along `|⟨e⟩|`: a systematic trend across bands on held-out data is the signature
of an insufficient basis, and is reported next to any σ.

```@example moment
prof = moment_band_profile(f; nbins = 4)
(; slope = prof.slope, r = prof.r, bands = [(b.lo, b.hi, b.n) for b in prof.bands])
```

**Local field and coverage.** The pair-consistent local field `h₁ = Σ_j e_j` over
the marked atom's `cutoff_pair` neighbours gives two coordinates per row, `‖h₁‖`
and the alignment `ê · ĥ₁`; [`moment_coverage`](@ref) compares a new set's
coordinates with the training set's so an extrapolation is named before it is
trusted.

```@example moment
lf_train = moment_local_field(mb, data)
lf_new = moment_local_field(mb, [planted_datum(rng) for _ in 1:10])
cov = moment_coverage(lf_train, lf_new)
(; frac_beyond = cov.frac_beyond, frac_anti = cov.frac_anti)
```

**Simple-feature floor.** On exactly the fit's kept rows, a trivial per-orbit
model (intercept + Legendre shell sums `Σ_j P_l(ê_i · ê_j)`) gives the
performance the SALC basis must beat. Whether the floor's features lie in the
SALC column span is **reported** (`inclusion`), not assumed; the nested bound
`sigma_model ≤ sigma_floor` applies to an unregularized fit only.

```@example moment
floor = moment_simple_floor(f, data; lmax = 1)
(; sigma_model = floor.sigma_model, sigma_floor = floor.sigma_floor,
   inclusion = round.(floor.inclusion; sigdigits = 3))
```

On this planted map both residuals are rounding noise and the features lie in the
SALC span (`inclusion ≈ 0`); on real data the floor is a number to beat, and a SALC
fit that does not beat it is mis-assembled or under-resolved.

## Scope and limits

- **Not persisted.** `save` refuses a moment basis, fit, or model by name; rebuild
  the basis from the crystal and spec.
- **One setup, one reference configuration per dataset.** A `SpinDatum` carries no
  provenance, so mixing computational setups biases `μ₀` invisibly — batch the
  training set upstream.
- **Reference geometry only.** This package carries no displacements, so every
  datum sits at the reference geometry by construction.
- The gate removes rows, never reweights them; a large gap between the gated and
  ungated coefficient sets (`f.coeffs` vs `f.coeffs_ungated`) means the gate is
  doing physics, not cleanup, and belongs in the report.
