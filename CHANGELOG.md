# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/); this package predates a tagged
release, so everything lives under *Unreleased*.

## [Unreleased]

### Added — thread-parallel design-matrix assembly and batch prediction

- The two design-matrix builders (`_design_energy`, `_design_torque`) and the vector
  `predict_energy` / `predict_torque` forms now parallelize over independent columns /
  configurations with `Threads.@threads`. Each task owns whole columns (or output slots),
  so writes are disjoint and the result is **identical at any thread count** (the per-atom
  gradient buffer `G` in the torque builder is now task-local). Speedup scales with
  `JULIA_NUM_THREADS` / `julia -t`; serial (1 thread) is unchanged. This is the cost that
  grows with dataset size — large active-learning batches and bigger supercells.
- New `test/unit/test_threading.jl` pins the threaded output to a race-free serial
  reference and to the independent scalar predict path (a data race would fail it under
  `julia -t N>1`).

### Added — `to_sunny` spin scaling routes (`:moment` / `:coupling`)

- `to_sunny` gains a `scaling` keyword so it can export a **non-half-integer** effective spin
  (e.g. the itinerant ``S_{\text{eff}} = m/(g\mu_B) \approx 1.1`` of bcc Fe), which Sunny's
  `Moment` (half-integer only) cannot carry directly:
  - `:moment` — put `S_eff` into the `Moment` and rescale exchange by `1/(SₐS_b)`; static
    energy *and* dispersion exact, but `S_eff` must be a half-integer (the previous, only
    behavior).
  - `:coupling` — keep `Moment` at a placeholder `s₀ = 1` and fold `S_eff` into the couplings
    (`J = M/(s₀√(SᵢSⱼ))`, single-ion `1/(s₀ Sᵢ)`); works for **any** positive `S_eff`. Only
    the magnon *dispersion* is physical (invariant under the overall spin scale `sᵢ→c sᵢ,
    J→J/c`); the static energy is rescaled.
  - `:auto` (default) picks `:moment` for half-integer spins, else `:coupling`. `mode` also
    gains `:auto` (`:dipole` for half-integer, `:dipole_uncorrected` otherwise). A `:coupling`
    placeholder cannot carry the quantum quadrupole, so single-ion + `:dipole` + `:coupling`
    is rejected. Modeled on the Magesty.jl Sunny export.
  - Validated: the `:moment` and `:coupling` dispersions agree to ``<10^{-7}`` for a
    half-integer spin where both apply, the `:coupling` static energy equals
    `(predict_energy − j0)/S`, and a non-half-integer `S_eff` now builds and disperses.

### Added — bcc Fe worked-example tutorial

- New tutorial [`docs/src/tutorials/case1_bcc_fe.md`](docs/src/tutorials/case1_bcc_fe.md):
  a real (non-synthetic) end-to-end fit for body-centered cubic iron — a 128-atom
  ``4\times4\times4`` supercell, isotropic two-body basis (``Im\bar3m``, 13 SALCs), fit from
  noncollinear spin-DFT energies and torques, with in-sample validation (parity plots),
  the isotropic ``J_{ij}`` read back out of [`bilinear_terms`](@ref) against pair distance,
  and a **live** [`to_sunny`](@ref) magnon dispersion (``\Gamma\text{–}H\text{–}N\text{–}\Gamma\text{–}P\text{–}H``
  and a ``\Gamma\text{–}N`` comparison to neutron data) using the new `:coupling` route at
  ``S_{\text{eff}} = 1.1``. The fit uses `torque_weight = 1.0` (a torque fit, constraining the
  energy gradient that sets the dispersion), reproducing the reference Magesty.jl tutorial's
  ``J_{ij}`` and magnon spectrum. The whole pipeline runs **live** in the docs (the ~10× `SCEBasis`
  speedup below makes the 128-atom build a few seconds). Ships the reference data — `POSCAR`,
  a 50-configuration `EMBSET`, and the experimental `febcc_spinwave.csv` — under
  `docs/src/tutorials/case1_inputs/`; the upstream VASP I/O and mean-field sampling that
  produced it live in the companion `SCETools.jl`. Adds `Sunny` to the docs environment.

### Performance — `SCEBasis` build (clusters + SALC)

- The symmetry-orbit and SALC construction were the dominant cost of building a basis on
  a large supercell (minutes for the 128-atom bcc Fe tutorial cell). Four numerics-
  preserving changes cut it to ~1.5 s (≈10×), and removed the GiB-scale transient
  allocations:
  - `build_clusters` reduces orbits by **growing each orbit from one representative**
    (`O(n_orbits · n_ops)`) instead of computing a canonical key per candidate
    (`O(n_candidates · n_ops)`); the canonical-key inner loop is now allocation-free
    (`Val(N)` static site tuples).
  - `build_salc_basis` connects all orbit members to the representative in **one
    `O(n_ops)` sweep** (was an `O(n_ops)` scan per member), and **memoizes the real
    Wigner-D matrices** `D^l(R_g)` by `(l, g)` across the build.
- Output is **byte-identical** (same orbits, same SALC keys/coefficients): verified by the
  full unit suite and the gauge-invariant Magesty oracle. Benchmarks and before/after
  numbers live in [`bench/`](bench/) and `.claude/bench_log.md`.

### Changed — public API naming consistency (BREAKING)

- A naming-and-usability pass renamed several exported symbols for consistency. All are
  **clean renames** (no deprecation aliases); update call sites accordingly:
  - **Count accessors** unified to an `n_*` prefix: `num_atoms` → `n_atoms`,
    `nsalc` → `n_salcs`. (`nobs` / `dof` are StatsAPI names and are unchanged.)
  - **SALC access** de-stuttered: the `SCEBasis` field `.salcs` (a `SALCBasis`) is renamed
    to `.salc_basis`, and a new exported accessor `salcs(basis)::Vector{SALC}` returns the
    basis functions directly — write `salcs(basis)[k]` instead of the old
    `basis.salcs.salcs[k]`.
  - **Predictor type** `SCEModel` → **`SCEPredictor`** (the lightweight, persistable
    predictor): the type, its constructor `SCEPredictor(fit)`, and `load(SCEPredictor, …)`.
  - **Setup reader** `read_input` → `read_setup` (returns the parsed crystal + interaction +
    symmetry setup; pairs with `read_configs` for training data).
  - **SALC kernel** `evaluate` → `evaluate_salc` (a less generic, collision-resistant name
    for the exported invariant evaluator).
- The persisted-model/basis **TOML schema is unaffected** (the on-disk tags `sce-basis` /
  `sce-model` and all doc keys are unchanged), so model and basis files written by an
  earlier build still load.
- Docs additions: the `SCEFit` (heavyweight, data-bearing) vs `SCEPredictor` (lightweight,
  persistable) roles are now contrasted in their docstrings, the README, and Getting
  Started; the `coeftable` columns (`body` / `orbit_id` / `ls` / `Lf` / `block` / `J`) gained
  a legend in the I/O guide.

### Added — lattice figures in the tutorials

- The Heisenberg-chain and kagome three-body tutorials now open with a generated lattice
  figure (CairoMakie) so the system is shown, not just described. Each figure draws the
  **unit (calculation) cell** and the periodic connectivity from the **actual** computed
  geometry — sites from `cartesian_positions`, bonds from `build_neighbor_list`, and the
  highlighted kagome triangle read back from `build_clusters` (its sites and periodic
  `shifts`). The chain figure shows the four sites closing into a ring (the dashed
  ``1\text{–}4`` bond across the cell boundary); the kagome figure shows the cell's three
  sites and a 3-body cluster that borrows two of its sites as periodic images from the
  neighbouring cell. The pictures cannot drift from what the basis is built on. `CairoMakie`
  is a `docs/` dependency only; the package itself gains no plotting dep.

### Changed — package renamed `MagestyRebuild` → `SCEFitting`

- The package, its module, and the repository directory were renamed from
  `MagestyRebuild` (`Magesty_rebuild.jl`) to **`SCEFitting`** (`SCEFitting.jl`), unifying the
  naming with the companion `SCETools.jl` under a shared `SCE*` family. The UUID is unchanged,
  so the package identity is preserved. The package extensions are now
  `SCEFittingGLMNetExt`, `SCEFittingSpglibExt`, and `SCEFittingSunnyExt`. The persisted-model
  TOML schema tags changed from `magesty-rebuild/sce-{basis,model}` to
  `scefitting/sce-{basis,model}` (`schema_version` is still `1`); any model TOML written by an
  earlier build must be re-saved. Downstream code updates `using MagestyRebuild` to
  `using SCEFitting`. The legacy `Magesty.jl` package (the design-reference original) is
  unaffected and keeps its name.

### Changed — VASP I/O moved to `SCETools.jl`

- The concrete VASP adapter (`SCEFitting.VASP`: `read_poscar`, `write_poscar`, `Oszicar`)
  has been **moved to the `SCETools.jl` package** (`SCETools.VASP`), joining the INCAR writer
  so all VASP I/O lives in one place. The core now keeps only the **code-agnostic DFT-data
  seam** — `AbstractDFTSource`, `AbstractTrainingDatum`, `SpinDatum`, `read_configs`, and
  `SCEDataset(basis, src)` — so the SCE pipeline stays code-agnostic. To read VASP training
  data, `using SCETools` and `SCETools.VASP.read_poscar` / `Oszicar` (the `SCEDataset(basis,
  src)` seam is unchanged). The `test/unit/test_vaspio.jl` unit tests, the
  `examples/vasp_dft_source.jl` example, and the VASP-vs-Magesty oracle cross-check moved to
  SCETools with it.

### Changed — sampling extracted into `SCETools.jl`; fitted-model introspection added

- The mean-field spin-configuration **sampler** (the P0–P4 work documented below) has been
  **moved out of this package** into the new auxiliary package `SCETools.jl`, which depends
  on `SCEFitting`. This package is now focused on building and fitting SCE models;
  generating spin configurations (and, later, active learning) lives in `SCETools.jl`. The
  removed exports are `AbstractSampler`, `MFASampler`, `MFASample`, `ExchangeModel`,
  `MultipoleField`, `sample`, `mfa_temperature_scale`, `mfa_sublattice_m`,
  `thermal_averaged_m`, and `tau_from_magnetization` (now exported by `SCETools`).
- **Added a code-neutral fitted-model introspection surface** so a downstream consumer reads
  the fitted Hamiltonian without reaching into the SALC-basis internals:
  - `multipole_terms(model)` returns a flat `Vector{MultipoleTerm}` — one record per cluster
    member / `l`-ordering of every SALC with a nonzero coefficient, carrying the raw `jϕ`
    coefficient (the per-N scale `(4π)^(body/2)` left for the consumer), the `body`, member
    `atoms` / `shifts`, per-site `ls`, and the `folded` tensor;
  - `bilinear_terms(model)` returns `(; pairs, onsites, skipped)` — the bilinear (`ls=[1,1]`)
    and single-ion (`ls=[2]`) channels as Cartesian `3×3` matrices, reusing the validated
    Sunny export conversion;
  - `num_atoms(model::SCEModel)` is a new method of the exported `num_atoms`.
  The gate is energy reconstruction (`test/unit/test_introspect.jl`): summing the per-term
  tesseral contraction reproduces `predict_energy − j0`.
- The tesseral spherical-harmonic submodule `SCEFitting.Harmonics` (`Zlm`, `lm_index`) is
  documented as a stable surface for downstream packages.

### Added — mean-field spin-configuration sampling: P4 (full multipole / many-body)

> The P0–P4 sampler entries below are retained as history; the code now lives in
> `SCETools.jl` (see the *Changed* entry above).

- **`MultipoleField` and `MFASampler(model::SCEModel; reference)`**: the full multipole
  mean-field sampler over **all** SCE clusters and harmonic orders (higher-order /
  many-body). The mean-field decoupling of any cluster term factorizes
  (`⟨∏ Z⟩ → ∏⟨Z⟩`), giving the generalized molecular field
  `h_a^{lm} = Σ_φ jφ·(4π)^(N/2)·folded · ∏_{b≠a} ⟨Z_{l_b}^{m_b}(e_b)⟩` — built by contracting
  each folded coefficient tensor against the *other* sites' multipole averages (the
  `accumulate_grad!` leave-one-out structure with the site-`a` harmonic left symbolic). The
  order parameters generalize from the magnetization to the **full per-atom multipole
  averages `⟨Z_lm⟩_a`** (`l ≤ lmax`), iterated to self-consistency by sphere quadrature; the
  `l=1` (bilinear) Perron eigenvalue sets `T_MF = ρ/3`, and the single-site Bingham /
  higher-multipole distribution is drawn with the Metropolis engine.
- `MFASampler(model::SCEModel; reference)` keeps every channel (bilinear, single-ion, and
  higher-order); `MFASampler(ExchangeModel(model); reference)` remains the bilinear-only
  truncation. Validated by the exact reduction to the single-global Langevin curve for a
  pure-bilinear model, by coupling-scale invariance, and by matching the single-site
  potential to the conditional mean SCE energy `⟨E | e_a⟩` of a biquadratic model (the
  many-body factorization check).

### Added — mean-field spin-configuration sampling: P3 (tensorial + single-ion)

- **Tensorial `ExchangeModel`**: `ExchangeModel` now carries the full bilinear tensor
  `bilinear[a,b] = S_ab` (Heisenberg + DM + anisotropic exchange) and single-ion anisotropy
  `onsite[a] = A_a` (the `ls=[2]` channel), with an `isotropic` flag for the P2 fast path.
  `ExchangeModel(model::SCEModel)` now extracts **all** of these (only the higher-order /
  higher-`l` SALCs are dropped); `ExchangeModel(Jiso; onsite)` and
  `ExchangeModel(bilinear; onsite)` are the raw constructors.
- **`MFASampler(exch; reference)` — tensorial path**: with DMI / anisotropic exchange or
  single-ion anisotropy the single-site potential `V_a(e) = β(e·g_a + e' A_a e)` (molecular
  field `g_a = Σ_b S_ab m_b ê_b`, `β = 3/(ρτ)`) gains an `l=2` Bingham factor that the
  closed-form vMF cannot represent, so the magnetizations are solved as `m_a = ⟨e·ê_a⟩` by
  sphere quadrature and the configurations are drawn with the **Metropolis** engine
  (per-atom chains, proposal scaled to the peak sharpness). The longitudinal molecular-field
  matrix `A[a,b] = −ê_a' S_ab ê_b` (so `T_MF = ρ/3`) generalizes the P2 isotropic case,
  which still takes the closed-form vMF path. **Noncollinear references** are supported
  (rigid cone axis, decision D2; a warning flags a non-stationary reference, e.g. DMI
  canting a collinear state). Because the single-ion enters with the same `β` as the
  exchange, the anisotropy persists above `T_MF` (an anisotropy-weighted paramagnet).
- An easy-axis single-ion sharpens the cone (and keeps order above the exchange `T_MF`); an
  easy-plane single-ion gives a girdle; the sampled `⟨Z_2m⟩` match the quadrature
  self-consistency. The shared sphere quadrature now caps its auto-sized node count (the
  large molecular fields of the `τ → 0` limit would otherwise blow it up).

### Added — mean-field spin-configuration sampling: P2 (multi-sublattice, isotropic)

- **`ExchangeModel`** (`docs/specs/mfa-sampling.md`): the neutral carrier of the isotropic
  bilinear exchange the mean-field sampler needs. `ExchangeModel(model::SCEModel)` extracts
  the symmetric `Jiso[a,b] = Σ_R J_iso(a,b,R)` from a fitted SCE by reusing the Sunny
  bilinear extraction and keeping the Heisenberg part `tr(M)/3` of each bond; the DMI /
  anisotropic, single-ion, and higher-order channels are dropped (P3/P4) and reported via
  `@warn`. `ExchangeModel(Jiso)` takes a raw symmetric matrix (an external-`Jij` reader is
  a P5 target).
- **`MFASampler(exch::ExchangeModel; reference)`**: the multi-sublattice sampler. The
  per-atom magnetizations `m_a(τ)` are solved from the coupled mean-field self-consistency
  `m_a = L(3(Ā m)_a/τ)`, where the molecular-field matrix `A[a,b] = −Jiso[a,b](ê_a·ê_b)`
  folds the reference directions in (so ferro / antiferro / ferri order all become
  ferromagnetic in the magnitude variables) and `Ā = A/ρ` is normalized by the Perron
  eigenvalue `ρ`, with `T_MF = ρ/3`. Distinct sublattices disorder at distinct rates under
  a single `T_MF`; because `Ā` is scale-free, only the coupling *ratios* matter — `m_a(τ)`
  is invariant under an overall coupling rescaling (decision D4). Each spin is then drawn
  from `vMF(ê_a, κ_a)` with the self-consistent per-atom concentration `κ_a = 3(Ā m)_a/τ`.
- **`mfa_sublattice_m(sampler, τ)`** returns the per-atom `m_a(τ)`. The coupled solve uses
  depth-1 Anderson acceleration (plain iteration suffers critical slowing as `τ → 1⁻`), and
  the construction verifies the reference is a stationary, sign-definite ordered state
  (warns on a non-stationary noncollinear or frustrated reference — the rigid-axis MFA is
  exact only for collinear references, decision D2). `MFASample.m` is now a per-atom vector
  per config.

### Added — mean-field spin-configuration sampling: P1 (single global, isotropic)

- **`MFASampler(reference)`** and the **`sample`** verb (`docs/specs/mfa-sampling.md`):
  the single global, isotropic mean-field sampler, building on the P0 single-site engine.
  Each spin is drawn from a von Mises–Fisher distribution `vMF(ê_a, κ)` about its reference
  direction, with one global concentration `κ = 3m/τ` fixed by the classical-Heisenberg
  mean-field self-consistency `m = L(3m/τ)` (`L` = Langevin function) in the reduced
  temperature `τ = T/T_MF`. Numerically equivalent to Magesty's `MfaSampling`, but with an
  explicit seeded `rng::AbstractRNG` and no Roots.jl dependency (a self-written bisection
  solves the monotone self-consistency).
- **`sample(sampler, n; tau | m, …)`** draws `n` configs at one control value;
  **`sample(sampler; tau | m, nsamples, …)`** sweeps a collection. Both return an
  **`MFASample`** (decision D1): `.configs::Vector{Matrix{Float64}}` plus parallel labels
  `.tau` and `.m`, iterable/indexable as its configs. Keywords `fixed` / `uniform` /
  `randomize` carry over from Magesty's `mfa_sweep`.
- **`AbstractSampler`** is the dispatch seam for the later model-backed samplers (P2+);
  **`thermal_averaged_m`** / **`tau_from_magnetization`** expose the self-consistency and
  its inverse; **`mfa_temperature_scale`** returns `T_MF` (decision D4 — `1.0`, reduced
  units, for the coupling-free global sampler).

### Added — `refit` and the regression-diagnostic accessors

- **`refit(f, estimator = OLS(); threshold = 0.0)`**: re-solve on the **support** of an
  existing fit — the de-biasing step that follows a sparse fit. A column survives when its
  scaled-magnitude contribution `|coef(f)[j]|·‖X[:, j]‖` exceeds `threshold` (`0` keeps the
  nonzero support exactly); coefficients off the support are zeroed and `j0` is recovered
  analytically. An empty support returns an all-zero `jϕ` (with a warning) and `j0 =
  mean(y_E)`. A `PrecomputedPilot` (or an `AdaptiveLasso` carrying one) is rejected — its
  fixed coefficient vector has the original column count, not the refit support length.
- **Diagnostic accessors**: `dof` (`length(coef(f)) + 1`), `rss_energy` / `rss_torque`
  (residual sums of squares), and `residuals_energy` / `residuals_torque` (the raw residual
  vectors). The energy and torque blocks are reported separately throughout, matching the
  existing `r2_*` / `rmse_*` split; `r2_energy` / `rmse_energy` / `r2_torque` / `rmse_torque`
  are refactored to build on the new `rss_*` so each metric has a single source of truth.
- Internal: the `(X, y)` centering / whitening / `groups` assembly is factored out of `fit`
  into `_assemble_problem`, shared by `fit` and `refit` so the two build identical designs
  (the oracle confirms `fit` is numerically unchanged).

### Added — adaptive / L0-approximating estimators

- **`AdaptiveRidge`** (in-tree, no extension): iterative reweighted ridge
  (Frommlet & Nuel 2016) that approximates an L0 penalty. It refits the analytic weighted
  ridge `(X'X + λ·Diagonal(w)) \ X'y` with `wⱼ = 1/(βⱼ² + ε)` until the relative ∞-norm
  change drops below `tol` (or `max_iter` steps), so large coefficients keep a light
  penalty and small ones are driven toward zero. `lambda = 0` reduces to `OLS`;
  `islinear ⇒ true` (a linear smoother in the converged-weight sense).
- **`AdaptiveLasso`** (type in core, GLMNet solve in `ext/SCEFittingGLMNetExt`): the
  one-shot Adaptive Lasso (Zou 2006). A `pilot` estimator (default `OLS`, any estimator
  allowed) supplies `β̂`, then a weighted Lasso is solved with per-column penalty factor
  `wⱼ = 1/max(|β̂ⱼ|, ε)^γ`. `gamma = 0` reduces exactly to a plain `Lasso`. It shares
  `ElasticNet`'s `lambda` / `standardize` / grouped-CV behavior (`lambda = nothing` selects
  λ by configuration-grouped CV with the adaptive weights held fixed).
- **`PrecomputedPilot`** (in-tree adapter): returns a fixed coefficient vector from
  `solve_coefficients` (length-checked against `size(X, 2)`), so an `AdaptiveLasso` can
  reuse a prior fit's `coef(f)` as its pilot instead of re-running a pilot regression.
- All three honor the column-centered `(X, y)` / analytic-`j0` contract. The GLMNet
  plumbing is shared between `ElasticNet` and `AdaptiveLasso` (an optional `penalty_factor`).
  Validated in `test/glmnet/` (support recovery, `γ = 0` ≡ plain Lasso, pilot
  pluggability, grouped-CV co-fit) and `test/unit/test_fit.jl` (core `AdaptiveRidge` /
  `PrecomputedPilot` solves, construction / validation, deferred-backend error).

### Changed — torque sign convention → Landau–Lifshitz / physical torque

- The per-atom torque is now `τ_a = −e_a × ∂E/∂e_a` (the physical / Landau–Lifshitz
  torque `m_a × B_eff,a`), matching the published *General spin models* paper
  (Phys. Rev. Research **8**, 023300 (2026)). Previously it was the energy-rotation-gradient
  `+e_a × ∂E/∂e_a` (the methods-paper Eq. 15 convention), the opposite sign.
- The change flips **both** sides together, so the fit is unchanged: `predict_torque` /
  `_design_torque` now compute `∂Φ/∂e × e = −e × ∂Φ/∂e`, and the DFT training target from a
  constrained OSZICAR is `τ_a = m_a × B_a` (was `−m_a × B_a`). Because the torque design
  matrix `X_T` and the target `y_T` both negate, the co-fit objective `‖X_T J − y_T‖²` is
  invariant — **`j0` and every coefficient `J` are identical**; only the *reported* torque
  sign changes. The energy fit, energies, and `predict_energy` are untouched.
- `predict_torque(model, config)` and the torque returned by `read_configs` / `SpinDatum`
  flip sign for downstream consumers. The finite-difference self-consistency gate
  (`test_torque.jl`, `test_nbody.jl`) and the oracle's analytic Heisenberg torque now check
  the `−e × ∇E` convention. Docs, examples, and the design-matrix convention notes updated.

### Added — Documenter.jl documentation site

- A full browser-viewable documentation site under `docs/` (Documenter.jl): home, a
  getting-started page, a four-part guide (building the basis, data and fitting,
  persistence and I/O, Sunny export), two narrated tutorials (Heisenberg chain, kagome
  three-body) with executed `@example` blocks, a three-part theory section (the SCE
  formalism, periodic resolvability, the architecture), and a complete API reference.
  Build with `make -C docs serve` (or `build` / `open`). Local-only for now
  (`remotes = nothing`, no `deploydocs` until a remote exists).
- **`docs` fixed several exported public APIs that silently had no runtime docstring** —
  a comment or a sibling definition sat between the docstring and the documented binding,
  so Julia never attached it. `coeftable`, `load`, `intercept`, `nobs`, `nsalc`,
  `rmse_energy`, `rmse_torque`, and the `AbstractDFTSource` / `AbstractTrainingDatum`
  abstract types now carry their own docstrings; every exported binding is documented
  (the API reference builds with `checkdocs = :exports`). No behavior change.

### Tested — N-body Wigner–Seitz cluster counting (`N ≥ 3`)

- Added `test/unit/test_ws_nbody.jl`, pinning the count of 3- and 4-body clusters on
  the Wigner–Seitz boundary against an **independent brute-force enumeration** of all
  compact clusters. The candidate set (for `N = 2, 3, 4`, including `pair_cutoff = Inf`)
  and the symmetry-orbit partition are checked to match exactly on cells deliberately
  seeded with face / edge / corner ties (cubic face-atoms, fcc, skewed hexagonal), and
  every emitted member is re-verified to have all `C(N,2)` edges at the minimum image.
- Documents and guards the **compact-cluster / third-edge** criterion: a cluster is
  admitted only when *all* its pairwise edges sit at their atom-pair minimum image
  simultaneously — individually minimum-image-resolvable pairs are not enough (the images
  minimizing `i–j` and `i–k` may push `j–k` onto a longer image). The regression includes
  an equal-spaced 1-D ring where every pair is minimum-image yet no compact triangle
  exists, and an over-/under-merge guard on the orbit reduction under a cubic point group.
  No behavior change — the enumeration was already correct; this makes it guaranteed.

### Added — GLMNet estimators (Lasso / elastic-net)

- **Estimator types in core** (`fitting/estimators.jl`): `ElasticNet(; alpha, lambda,
  standardize, nfolds, select, seed, nlambda)` and the `Lasso(; …)` convenience
  (`alpha = 1`). Like the Spglib/Sunny seam, the types live in the core package — named,
  validated, and dispatched on without the heavy dependency — while the actual solve
  lights up only under `using GLMNet`.
- **`solve_coefficients(::ElasticNet, X, y; groups)`** in the new
  `SCEFittingGLMNetExt` extension. GLMNet minimizes
  `(1/2n)·‖y − Xβ‖² + λ·[(1−α)/2·‖β‖₂² + α·‖β‖₁]` on the column-centered `(X, y)` the
  fit hands it, with `intercept = false` (so `j0` is still recovered analytically) and
  column `standardize` (the penalty acts per-column at `λ·std`, returning β on the
  original scale). `lambda = nothing` selects the penalty by K-fold cross-validation over
  GLMNet's automatic λ path — `select = :lambda_min` (lowest CV error) or `:lambda_1se`
  (sparsest within one SE); a numeric `lambda` fits at exactly that penalty.
- **Configuration-grouped, reproducible CV.** `fit` now passes per-row `groups` labels
  to `solve_coefficients`; for an energy+torque co-fit a configuration's energy row and
  all its torque-component rows share a label, so CV folds never split one configuration
  (which would leak within-configuration structure and bias λ selection). Folds are
  assigned deterministically by a seeded `hash` ranking — reproducible without taking on
  a `Random` dependency in the extension.
- The `AbstractEstimator` contract gains an optional `groups` keyword on
  `solve_coefficients`; `OLS`/`Ridge` accept and ignore it.
- Validated in a separate `test/glmnet/` environment (heavy Fortran-backed dependency,
  mirroring `test/sunny/` and `test/oracle/`): tiny-λ ≈ OLS, the analytic-`j0` centering
  invariant, CV support recovery and sparsity, `:lambda_1se` shrinkage, seeded
  reproducibility, ElasticNet ≠ Lasso, and an energy+torque grouped-CV co-fit. Core-only
  construction/validation and the deferred-backend error are covered in the main suite
  (`test/unit/test_fit.jl`).

### Added — Sunny.jl export (supercell + primitive-cell routes)

- **Conversion core** (`sce/sunny.jl`): turns a fitted `SCEModel` into the
  Sunny-representable channels — `ls=[1,1]` 2-body → a 3×3 exchange matrix
  (`M = (3/4π)·folded`, carrying Heisenberg / Dzyaloshinskii–Moriya / symmetric-Γ
  in one matrix) and `ls=[2]` 1-body → a traceless-symmetric single-ion tensor.
  `ls=[0…]` channels fold into `j0`; every other SALC (3-body+, higher `l`) is
  **skipped and reported**, since Sunny cannot represent it. The directed cluster
  members fold into one matrix per undirected supercell bond
  (`_sunny_supercell_terms`), and the whole conversion is gated **without Sunny** by
  reconstructing the energy (`_reconstruct_energy ≈ predict_energy − j0`).
- **`to_sunny(model; spins, g, mode, placement)`** in the `SCEFittingSunnyExt`
  extension (loaded by `using Sunny`): builds a real `Sunny.System` on the training
  supercell (P1, inhomogeneous), placing `set_exchange_at!` per bond (rescaled
  `J = M/(SₐS_b)`) and `set_onsite_coupling_at!` per atom. The system's classical
  energy reproduces `predict_energy − j0` exactly; exchange is independent of the
  spin length and Sunny mode, the single-ion term carries the
  classical/`:dipole`-quantum rescaling. Skipped channels are surfaced via `@warn`.
- **Primitive-cell unfold** (`placement = :primitive`/`:auto`, `_sunny_primitive`):
  recovers the chemical primitive cell from the space group's pure translations, groups
  supercell atoms into sublattices, and folds the supercell bonds onto one Sunny bond
  per **primitive** bond `(i, j, n)` (no multiplicity — Sunny's periodic replication
  restores it), for *unfolded* spin-wave dispersion. A `clean` flag detects when the
  model does not live on the primitive cell (interaction range reaching the supercell
  boundary), and the export falls back to the exact supercell route.
- Conversion math is a **core dependency-free** layer; only the `Sunny.System`
  assembly lives in the extension. Validated by `test/unit/test_sunny.jl` (Sunny-free:
  the `Z₁`/`Z₂` contractions, classification, energy reconstruction, the primitive fold
  and its clean detection, skip reporting) and the separate `test/sunny/` environment
  (the real `Sunny.System` energy vs the SCE energy across modes/spins for both routes
  — the primitive system reshaped back to the supercell reproduces the SCE energy,
  confirming the unfolded bonds and offsets — plus the skip warning and per-species
  spins).

### Changed — minimum-image periodic resolvability (Wigner–Seitz cell)

- **Periodic-image selection** (`geometry/neighborlist.jl`): `AbstractImageSelection`
  with `MinimumImage` (new default) and `AllImages`. Previously the neighbor list kept
  every periodic image within a *spherical* cutoff, which over-counts beyond `L/2`: a
  farther image of an atom carries the **same spin** as its minimum image, so the two
  interactions are not independently resolvable from a finite supercell (their design
  columns are collinear). `MinimumImage` keeps only the minimum-image, Wigner–Seitz-cell
  pairs — the physically resolvable set, which reaches the body-diagonal corner
  `(L/2,L/2,L/2)` at `√3·L/2`, **not** a sphere of radius `L/2` — with WS-boundary ties
  (faces 2-fold, edges 4-fold, corners 8-fold) kept as distinct members and `i==j`
  self-pairs dropped (same spin ⇒ a constant or a 1-body alias, never an independent
  pair). The image-box search is adaptive (provably sufficient on skewed / non-reduced
  cells). `AllImages` retains the old every-image behavior as the **generalized-Bloch /
  spin-spiral seam**, where `e^{iq·R}` resolves what one supercell cannot.
- **`pair_cutoff = Inf`** (`Interaction`) now means "every resolvable pair" — the whole
  WS cell — under `MinimumImage` (cf. Magesty's `-1` sentinel); `≤ 0` / `NaN` are
  rejected. `SCEBasis(...; images = MinimumImage())` and the `input.toml`
  `[interaction].images` / `pair_cutoff = inf` keys thread the choice through.
- **N-body min-image consistency** (`clusters/enumerate.jl`): a `MinimumImage` clique
  requires every edge at its atom-pair minimum-image distance and all atoms distinct, so
  a cluster is never built from an aliased bond or a reused atom image. The `AllImages`
  edge cutoff now uses the same relative tolerance as the tie band (was an absolute
  `+1e-9`).
- For a cutoff below half the smallest perpendicular cell width, `MinimumImage` and
  `AllImages` coincide (each in-cutoff image is already the minimum), so the existing
  physics is unchanged: the Heisenberg `J = 2√3·jϕ` recovery, the kagome 3-body co-fit,
  and the oracle all use cutoffs in this regime. Bit-for-bit Magesty agreement is not a
  goal; a single-atom cell's self-pair "bond" (1 orbit under Magesty/`AllImages`) is now
  0 orbits under the default `MinimumImage` — a deliberate refinement.
- Validated (`test/unit/test_imageselection.jl`): cutoff `< L/2` equivalence, the cubic
  8-fold corner and `L/2` 2-fold face ties, alias rejection above `L/2`, `AllImages`
  rejects `Inf`, hexagonal / skewed-cell search-box sufficiency, no-self-pair and
  distinct-atom clusters, a full-WS design matrix with full column rank under symmetry,
  `input.toml` `inf` + `images`, and a `pair_cutoff = Inf` persistence round-trip.

### Added — VASP I/O and a code-agnostic DFT-source seam

- **DFT-source boundary** (`io/dftsource.jl`): `AbstractDFTSource` +
  `read_configs(src) -> Vector{SpinDatum}`, the `SpinDatum` training datum (energy,
  `3×n` unit spin directions, moment magnitudes, constraining field, and the derived
  torque target `τ_a = −m_a × B_a`, eV), and `SCEDataset(basis, src | data; use_torque)`.
  This is the *only* thing the SCE pipeline consumes — the originating DFT code is
  irrelevant once you hold a `SpinDatum`/`SCEDataset`.
- **VASP adapter** (`io/vasp.jl`, `module SCEFitting.VASP`): `read_poscar` /
  `write_poscar` (POSCAR/CONTCAR ↔ `Crystal`; scaling incl. negative-volume,
  Direct/Cartesian, Selective dynamics, VASP4/5) and `Oszicar`, an `AbstractDFTSource`
  over constrained-noncollinear OSZICARs (energy `F=`/`E0`, `MW_int`/`M_int` moments,
  `lambda*MW_perp` field, SAXIS `Rz(α)·Ry(β)` rotation). Code-specific I/O is a
  **namespaced submodule** kept out of the core: the core and its export list do not grow
  as DFT codes are added; only the code-agnostic boundary is exported.
- A torque-carrying `SCEDataset` rejects all-zero torque targets (no constraining field
  found) so unconstrained data cannot be silently fit as "torque = 0".
- Validated (`test/unit/test_vaspio.jl`): POSCAR (Direct/Cartesian/VASP4/negative-volume/
  selective-dynamics/round-trip), OSZICAR (energy kinds, `mint`, SAXIS, multi-step,
  missing field, errors), the torque formula, and source → `SCEDataset`. The from-scratch
  POSCAR/OSZICAR readers are cross-checked **bit-for-bit against Magesty.jl**'s parsers in
  the oracle (synthetic files, since real VASP outputs are not vendored).

### Added — tabular coefficient output (Tables.jl)

- **`coeftable(f)` / `coeftable(model)` → `SCECoefficients`** (`sce/coeftable.jl`): a
  Tables.jl source — one row per SALC (`body`, `orbit_id`, `ls` as a comma string,
  `Lf`, `block`, `J`) — so the fitted coefficients drop into `DataFrame` /
  `CSV.write` / `Arrow.write`. The library owns the internal-storage → labeled-row
  mapping (the `J` column pairs with the basis keys positionally); the caller brings
  the table/IO package. `j0` is the intercept (`intercept(c)`), not a row. Tables.jl
  is a lightweight core dependency, the same seam that would later open tabular
  training-data ingestion.

### Added — persistence + TOML input files

- **Persistence** (`sce/persist.jl`): `SCEFitting.save(path, x)` and
  `SCEFitting.load(SCEBasis | SCEModel, path)` serialize a self-contained,
  human-readable **TOML** document — the crystal, the space-group ops, the
  interaction, and the *full* SALC basis (every member / term / folded tensor); a
  model adds `j0` and per-`SALCKey` coefficients. Reload rebuilds the basis verbatim
  (no re-projection) and re-pairs coefficients to the basis **by key**, not by
  position (a `SCEModel` saved by one build reloads correctly into a structurally
  identical basis). The `struct ⇄ Dict` schema (`_to_doc` / `_from_doc`) is
  format-agnostic and unit-tested without any serializer; TOML is the stdlib (no
  dependency) and round-trips `Float64` exactly.
- **TOML input** (`sce/input.jl`): `read_input(path)` and the new
  `SCEBasis(path::AbstractString; backend, tol)` constructor build a basis from a
  human-authored `input.toml` (`[structure]` inline crystal, `[interaction]`,
  optional `[symmetry]`); keyword arguments override the file's backend/tol.
  Training data and the estimator are kept out of the file (loaded/chosen in Julia),
  mirroring the basis/data separation. `read_input` is exported.
- Considered JSON vs TOML deliberately: stdlib TOML round-trips `Float64` exactly and
  expresses the full nested SALC document, so the persistence artifact and the input
  file share one zero-dependency format. The schema layer stays format-agnostic, so a
  JSON (or other) backend remains a thin future addition.
- Validated (`test/unit/test_persist.jl`, `test/unit/test_input.jl`): basis / model /
  fit round-trips (predictions bit-identical, coefficients re-paired by key under a
  scrambled on-disk order, multi-op space-group ops, empty basis), input parsing +
  defaults + keyword overrides + error paths.

### Added — arbitrary body order (N-body clusters)

- **Cluster enumeration** (`clusters/enumerate.jl`): `candidate_clusters` generalized
  to `N`-body pairwise-within-cutoff cliques (`N = 2` stays the directed neighbor pairs).
- **SALC projection** (`basis/salcbasis.jl`): general site stabilizer with induced
  permutations; projection over the combined (ordering × coupling-path × `Mf`) space
  with the action matrix built by contracting against the package's own orthonormal
  coupled tensors (no 6j/9j). Handles, at `N ≥ 3`, coupling-path mixing and — for
  unequal `l` on symmetry-equivalent sites (e.g. `l=(1,1,2)` on a triangle) —
  `l`-ordering mixing. Improper-op parity is automatic (no proper-part special case).
- **Multi-term SALCs** (`basis/salc.jl`): a SALC carries one `SALCTerm` (own `ls` +
  `folded`) per `l`-ordering; `evaluate` and `accumulate_grad!` loop over terms.
- **`SALCKey`** stays injective when a proper-subgroup stabilizer splits a degenerate
  multiset into several ordering orbits (`block` runs across them).
- Validated by ground-truth invariance / time-reversal / linear-independence tests at
  `N = 3` (incl. the multi-term channel) and an energy+torque 3-body recovery;
  cross-checked against Magesty.jl — per-`(body, ls, Lf)` invariant-subspace dimensions
  agree exactly through 3-body, and Magesty's own SALCs independently pass invariance.

### Added — torque observable (energy + torque co-fit)

- **Gradient kernel** (`basis/salc.jl`): `accumulate_grad!` sums `jϕ·∂Φ/∂e_a` per
  site (product rule over cluster members), sharing the energy kernel's `μ`-mapping
  and `(4π)^(N/2)` scale.
- **Torque design matrix** `X_T` and per-config torque targets via the four-argument
  `SCEDataset(basis, configs, energies, torques)`; entry `(4π)^(N/2)·(e_a × ∂Φ/∂e_a)`,
  rows flattened config/atom/`xyz`.
- **Energy+torque co-fit**: `fit(SCEFit, dataset, est; torque_weight = w)` minimizes
  `(1−w)·MSE_E + w·MSE_T` by per-block whitening; `j0` stays analytic from the
  energy block.
- **Prediction & metrics**: `predict_torque`, `r2_torque`, `rmse_torque`,
  `has_torque`. `predict_torque = e × ∇(predict_energy)` by construction (validated
  by on-sphere finite differences and the Heisenberg closed form).

### Added — v0 vertical slice (energy fitting end-to-end)

- **Geometry**: `Lattice`, `Crystal`, and a generalized cutoff `NeighborList`
  (`build_neighbor_list`) that replaces a fixed 27-cell image grid; `NeighborPair`
  retains the inter-site lattice translation `R`.
- **Harmonics** (`basis/Harmonics.jl`): real tesseral `Zₗₘ` + tangent-projected
  gradient (Drautz convention).
- **Angular momentum** (`basis/AngularMomentum.jl`): Clebsch–Gordan (Racah),
  `wignerD_real` (least-squares from the package's own `Zₗₘ`), complex→real coupled
  tensors; `CoupledBasis{R}`.
- **Symmetry**: `AbstractSymmetryBackend` with in-tree `NoSymmetry` and a
  `SpglibBackend` whose method lives in `ext/SCEFittingSpglibExt`.
- **Clusters**: orbit reduction with `R`-carrying members (`build_clusters`).
- **SALC basis**: orbit–stabilizer projector (isotropic and anisotropic channels),
  deterministic gauge, canonical `SALCKey` column addressing (`build_salc_basis`).
- **Fitting / API**: `Interaction`, `SCEBasis`, `SCEDataset` (energy design
  matrix), `SCEModel`/`SCEFit`, `fit`, `predict_energy`, `OLS`/`Ridge`,
  `coef`/`intercept`/`nobs`/`r2_energy`/`rmse_energy`.
- A runnable Heisenberg-chain example (`examples/heisenberg_chain.jl`) that
  recovers `J`; an oracle test environment (`test/oracle/`) validating against a
  pinned Magesty.jl.

### Notes

- Bit-for-bit agreement with Magesty.jl is explicitly not a goal — see
  `docs/design-notes.md`.
