# The pointed (site-marked) SALC basis for the adiabatic site-moment expansion
# m_i(e) — the moment channel's counterpart of `SCEBasis`. Design record:
# _brain_storming/adiabatic-moment-sce (D1–D10 + the D8 addendum, M3-1/M3-2).
# Ported from SLCE.jl at `3d54abb` = `6960276` (the basis) + `bd03517` (the
# mark→term index fast design path, asserted value-identical to the full
# evaluation) + `28dfa24` (the resolvability cache AND the fix for the null
# report of a WIDE signature block, which the economy SVD under-enumerated — the
# QR completion below). The one spelling change is the screen keyword
# (`isotropy = true` keeps L_S = 0 only — this package's name for upstream's
# `soc = false`; the two are the same screen under opposite names).
#
# The MARK is the displacement decor `SiteDecor(disp = (1, 0))` — the polar l = 0,
# k = 1 radial factor |u|² — evaluated on the synthetic indicator field u = x̂ at the
# marked atom: the factor is 1 on the marked site and 0 elsewhere, so the existing
# decor SALC engine, evaluation kernels, and canonical-member machinery carry the
# pointed basis without modification. The marked site's own ê factor is the decor's
# SPIN part (`SiteDecor(spin = l, disp = (1, 0))`); at evaluation the spin matrix's
# marked COLUMN is substituted by the evaluation axis ê (the D8 addendum's
# marked-column substitution — identity for a mode-4 datum, the constraint axis for
# mode 1), which is exact because every pointed label carries exactly one mark, so a
# member not marked at the row's atom contributes zero before it can read the
# substituted column.
#
# Three mark-aware rules the generic engine cannot express, all confined here:
#   * the marked site's ê factor is allowed for ANY species (an E-inactive species'
#     induced moment is exactly what the channel predicts — a literal per-species
#     l-cap would delete it), while ENVIRONMENT spin factors are allowed only for
#     species the consumer samples (M3-1 decision A; slaving makes the mediated
#     physics expressible through sampled-species clusters);
#   * the 3-body cutoff is MARK–ENVIRONMENT-BOND based, not all-edge: the triangle
#     is pinned by the two mark bonds (each minimum-image within its radius), and
#     the environment–environment edge is free — an all-edge cut at the same radius
#     keeps 3 of the 15 nn star pairs on FeGe and loses 20–32 % in σ (M2-5);
#   * time reversal: y = ê·M is TR-even, so only even-Σl labels exist (the mark's
#     spin rank counts) — enforced by the engine's existing screen.

# ── refusal vocabulary shared with the resolvability gate ──────────────────────────

"""
    UnclassifiableBasis <: Exception

Thrown when a gate cannot classify a basis because a member carries two spin
factors on ONE reference-cell atom (two periodic images of one neighbor), so two
harmonics of the SAME unit vector multiply and the monomial key the gate is built on
no longer names a function uniquely. The honest answer is "unknown", never "none";
the gate refuses rather than reports. Two gates raise it:

- [`moment_resolvability`](@ref) — two ENVIRONMENT spin factors on one atom (the
  mark factor is exempt: it reads the independent axis variable), where the signature
  expansion would overcount the rank.
- [`SCEDataset`](@ref) — an energy-side basis with a repeated-atom (`AllImages`
  self-image) member, whose columns are redundant on the reference cell and which the
  per-orbit function-space reduction cannot certify. Only the fitting door refuses;
  building such a basis as a tiling template stays legal.

(Upstream SLCE.jl defines the same exception for its energy-side gate.)

!!! note "What is and is not public here"
    The **type** and the doors that raise it are public API: catch it by type. The
    **`reason` text is not** — never classify a failure by matching the message.
"""
struct UnclassifiableBasis <: Exception
    reason::String
end
Base.showerror(io::IO, e::UnclassifiableBasis) =
    print(io, "UnclassifiableBasis: ", e.reason)

# Relative cut on (assembled value) / (gross accumulated mass): "did this cancel?",
# never "is this small?". Cancellation to roundoff lands at ~N·eps of the gross
# mass; a column that merely has few terms still has a ratio of order 1. A cut
# taken against the matrix's own maximum cannot see the case where EVERYTHING
# cancels, which is exactly what a Wigner–Seitz-tied cell produces. (Same value
# and reading as upstream's `_CANCELLATION_RTOL`.)
const _CANCELLATION_RTOL = 1e-10

# ── MomentSpec ─────────────────────────────────────────────────────────────────────

"""
    MomentSpec(; lmax_env, sampled, lmax_mark = 2, marked = nothing, nbody = 3,
               cutoff_pair, cutoff_star = nothing, lsum = nothing, isotropy = true)

The truncation spec of a pointed moment basis ([`MomentBasis`](@ref)). Species are
indexed like `Crystal.species`; `lmax_env` fixes the species count. The same spec can
be written as the `[moment]` section of a TOML setup file ([`read_setup`](@ref),
`MomentBasis(path)`).

- `lmax_env::Vector{Int}` — per-species cap on **environment** spin factors
  (`0` = this species' spins never enter the environment).
- `sampled::Vector{Bool}` — which species the downstream consumer (MC / LLG)
  actually samples. **Required, and load-bearing**: every species with
  `lmax_env > 0` must be sampled (a basis reading unsampled spins cannot be
  evaluated at run time) — refused loudly, never fixed up.
- `lmax_mark::Int` — cap on the marked site's own ê factor (any species; the mark
  rank counts toward the even-Σl time-reversal screen, so odd ranks pair with odd
  environment content).
- `marked::Vector{Bool}` — which species' site moments the basis expands
  (default: every species).
- `nbody` — 1, 2, or 3 (1-body = the per-orbit intercepts μ₀ plus even-rank
  single-site ê invariants).
- `cutoff_pair` — mark–environment bond radius (Å) for 2-body clusters: a scalar
  or a symmetric per-species-pair matrix.
- `cutoff_star` — mark–environment bond radius for 3-body stars (default:
  `cutoff_pair`). Only the two mark bonds are constrained; the
  environment–environment edge is free.
- `lsum` — optional cap on the total spin rank of a label (`nothing` = uncapped).
- `isotropy` — keep only the `L_S = 0` blocks (default `true`: the adiabatic map
  is treated as spin-rotation covariant, exactly like an `isotropy = true` energy
  basis). `isotropy = false` keeps every `L_S` — upstream SLCE.jl spells this
  `soc = true`; the two are the same screen under opposite names.
"""
struct MomentSpec
    nbody::Int
    lmax_mark::Int
    lmax_env::Vector{Int}
    sampled::Vector{Bool}
    marked::Vector{Bool}
    cutoff_pair::Matrix{Float64}
    cutoff_star::Matrix{Float64}
    lsum::Int
    isotropy::Bool
end

function MomentSpec(; lmax_env::AbstractVector{<:Integer},
                    sampled::AbstractVector{Bool},
                    lmax_mark::Integer = 2,
                    marked::Union{Nothing,AbstractVector{Bool}} = nothing,
                    nbody::Integer = 3,
                    cutoff_pair::Union{Real,AbstractMatrix{<:Real}},
                    cutoff_star::Union{Nothing,Real,AbstractMatrix{<:Real}} = nothing,
                    lsum::Union{Nothing,Integer} = nothing,
                    isotropy::Bool = true)::MomentSpec
    nkd = length(lmax_env)
    nkd >= 1 || throw(ArgumentError("lmax_env must name at least one species"))
    all(l -> l >= 0, lmax_env) ||
        throw(ArgumentError("lmax_env entries must be ≥ 0; got $lmax_env"))
    lmax_mark >= 0 || throw(ArgumentError("lmax_mark must be ≥ 0; got $lmax_mark"))
    1 <= nbody <= 3 || throw(ArgumentError("nbody must be 1, 2, or 3; got $nbody"))
    length(sampled) == nkd ||
        throw(ArgumentError("sampled has $(length(sampled)) entries for $nkd species"))
    for s = 1:nkd
        lmax_env[s] > 0 && !sampled[s] &&
            throw(ArgumentError("species $s carries environment spin factors " *
                                "(lmax_env = $(lmax_env[s])) but is not sampled by " *
                                "the consumer — a basis reading unsampled spins " *
                                "cannot be evaluated at run time. Either sample the " *
                                "species or set its lmax_env to 0 (M3-1: induced " *
                                "physics is expressible through sampled-species " *
                                "clusters)"))
    end
    mk = marked === nothing ? fill(true, nkd) : collect(Bool, marked)
    length(mk) == nkd ||
        throw(ArgumentError("marked has $(length(mk)) entries for $nkd species"))
    any(mk) || throw(ArgumentError("no species is marked — nothing to expand"))
    _cut(c) = c isa Real ? fill(Float64(c), nkd, nkd) : Matrix{Float64}(c)
    cp = _cut(cutoff_pair)
    cs = cutoff_star === nothing ? copy(cp) : _cut(cutoff_star)
    for (name, m) in (("cutoff_pair", cp), ("cutoff_star", cs))
        size(m) == (nkd, nkd) ||
            throw(ArgumentError("$name is $(size(m)) for $nkd species"))
        m == m' || throw(ArgumentError("$name must be symmetric"))
        all(v -> !isnan(v) && v >= 0, m) ||
            throw(ArgumentError("$name entries must be ≥ 0 Å"))
    end
    ls = lsum === nothing ? typemax(Int) : Int(lsum)
    ls >= 0 || throw(ArgumentError("lsum must be ≥ 0; got $lsum"))
    return MomentSpec(Int(nbody), Int(lmax_mark), collect(Int, lmax_env),
                      collect(Bool, sampled), mk, cp, cs, ls, isotropy)
end

_mark_decor(l::Int)::SiteDecor =
    l == 0 ? SiteDecor(disp = (1, 0)) : SiteDecor(spin = l, disp = (1, 0))

# In a pointed label the disp factor IS the mark.
is_marked(d::SiteDecor)::Bool = has_disp(d)

# All sorted decor multisets for body order N: exactly one marked slot (spin rank
# 0…lmax_mark), every environment slot a pure spin factor (rank ≥ 1), total spin
# rank even (TR) and ≤ lsum. Which SITE may carry which decor is the per-assignment
# `admit` closure's business, not the label's.
function _moment_labels(spec::MomentSpec, N::Int)::Vector{Vector{SiteDecor}}
    lem = maximum(spec.lmax_env; init = 0)
    labs = Vector{Vector{SiteDecor}}()
    if N == 1
        for lm = 0:spec.lmax_mark
            (iseven(lm) && lm <= spec.lsum) || continue
            push!(labs, [_mark_decor(lm)])
        end
    elseif N == 2
        for lm = 0:spec.lmax_mark, le = 1:lem
            (iseven(lm + le) && lm + le <= spec.lsum) || continue
            push!(labs, sort([_mark_decor(lm), SiteDecor(spin = le)]))
        end
    elseif N == 3
        for lm = 0:spec.lmax_mark, l1 = 1:lem, l2 = l1:lem
            (iseven(lm + l1 + l2) && lm + l1 + l2 <= spec.lsum) || continue
            push!(labs, sort([_mark_decor(lm), SiteDecor(spin = l1),
                              SiteDecor(spin = l2)]))
        end
    end
    sort!(labs)
    unique!(labs)
    return labs
end

# ── pointed 3-body star candidates ─────────────────────────────────────────────────

# Star clusters {mark, env₁, env₂}: both mark–environment bonds are minimum-image
# neighbor pairs within the mark–env star radius for their species pair; the
# environment–environment edge is FREE (the triangle is pinned by the two mark
# bonds, so no periodic alias hides there — M2-5). The candidate set is closed
# under the space group by construction (species and minimum-image distances are
# symmetry invariants), which `_orbits_from_members`' closure assertion re-checks.
#
# MULTIPLICITY CONVENTION: `candidate_clusters` lists every physical instance once
# per SITE ORDERING (3! = 6 anchored variants for a 3-body) — the space the SALC
# projection needs, and the multiplicity `_canonicalize_members` folds into the
# member weights. A pointed candidate set with fewer orderings per instance yields
# SALCs scaled down by the missing factor (measured upstream: 3 of 6 orderings
# halved the closed-star column against the prototype's 6.0 geometric oracle),
# silently breaking cross-orbit coefficient comparability. So the enumeration here
# first finds each triangle once (per translation class), then expands EVERY class
# to all 3! re-anchored orderings — identical members to what `candidate_clusters`
# would emit for the clusters it also admits.
const _PERMS3 = ((1, 2, 3), (1, 3, 2), (2, 1, 3), (2, 3, 1), (3, 1, 2), (3, 2, 1))

function _pointed_star_candidates(crystal::Crystal, nl::NeighborList,
                                  spec::MomentSpec)::Vector{ClusterMember}
    nat = n_atoms(crystal)
    sp = crystal.species
    fac = 1.0 + nl.tol
    nbrs = [Tuple{Int,SVector{3,Int}}[] for _ = 1:nat]
    for p in nl.pairs
        spec.marked[sp[p.i]] || continue
        p.distance <= spec.cutoff_star[sp[p.i], sp[p.j]] * fac || continue
        push!(nbrs[p.i], (p.j, p.shift))
    end
    z = SVector{3,Int}(0, 0, 0)
    classes = Dict{Any,ClusterMember}()
    for i = 1:nat
        ns = nbrs[i]
        for x = 1:length(ns), y = (x + 1):length(ns)
            (j, Rj) = ns[x]
            (k, Rk) = ns[y]
            # Only an exact (atom, shift) repeat is skipped. Two DIFFERENT minimum
            # images of one neighbor (j == k, Rj != Rk — a small cell's tie) are
            # kept, as upstream keeps them, so the pointed columns agree with
            # SLCE.jl on the same cell. Under plain PBC both environment factors
            # then read the one spin e_j and the member reduces to a lower-body
            # function; `moment_resolvability` refuses such a basis as
            # `UnclassifiableBasis`, and `MomentDataset` runs that gate at its
            # door — a hard refusal, never a silent overcount. (The energy side's
            # `candidate_clusters` drops these members instead; the pointed
            # enumeration does not, deliberately.)
            (j, Rj) == (k, Rk) && continue
            m = ClusterMember([i, j, k], [z, Rj, Rk])
            get!(classes, _member_sig(m), m)
        end
    end
    out = ClusterMember[]
    for sig in sort!(collect(keys(classes)))     # deterministic emission order
        m = classes[sig]
        for p in _PERMS3
            atoms2 = [m.atoms[p[1]], m.atoms[p[2]], m.atoms[p[3]]]
            s1 = m.shifts[p[1]]
            shifts2 = [m.shifts[p[1]] - s1, m.shifts[p[2]] - s1, m.shifts[p[3]] - s1]
            push!(out, ClusterMember(atoms2, shifts2))
        end
    end
    return out
end

# ── MomentBasis ────────────────────────────────────────────────────────────────────

"""
    MomentBasis(crystal, spec; backend = NoSymmetry(), tol = 1e-5,
                tie_tol = same-distance default)

The pointed SALC basis of the adiabatic site-moment expansion `m_i(e)`: per marked
reference-cell atom `i`, the design row is the SALC vector evaluated with the mark
on `i` (`X[(config, i), φ] = Φ_φ(i; e, ê_i)`), and one shared coefficient vector
serves every symmetry-equivalent site. Built on the ordinary decor SALC engine —
full space group, orbits split under the site group through the SALC projection —
with the mark-aware admission rules of the design record (see the file header).

Fields: `crystal`, `spacegroup`, `spec`, `salc_basis` (keys sorted, addressable),
`records` (one NamedTuple per SALC: representative geometry for reporting —
`body`, `species`, `edges`, `nmem`), and `marked_atoms` (the design's row atoms).
`tie_tol` records the boundary-tie band the basis was built with. The trailing
`resolvability` field is [`moment_resolvability`](@ref)'s default-`rtol` cache (a
`Ref`, `nothing` until the first call) — derived data, never set by hand.
"""
struct MomentBasis
    crystal::Crystal
    spacegroup::SpaceGroup
    spec::MomentSpec
    salc_basis::SALCBasis
    records::Vector{NamedTuple}
    marked_atoms::Vector{Int}
    tie_tol::Float64
    resolvability::Base.RefValue{Union{Nothing,NamedTuple}}
end

n_salcs(mb::MomentBasis)::Int = length(mb.salc_basis.salcs)
salcs(mb::MomentBasis)::Vector{SALC} = mb.salc_basis.salcs

function Base.show(io::IO, mb::MomentBasis)
    print(io, "MomentBasis(", length(mb.salc_basis.salcs), " SALCs, ",
          length(mb.marked_atoms), " marked atoms, nbody = ", mb.spec.nbody, ")")
end

function MomentBasis(crystal::Crystal, spec::MomentSpec;
                     backend::AbstractSymmetryBackend = NoSymmetry(),
                     tol::Real = 1e-5,
                     tie_tol::Real = _SAME_DIST_RTOL)::MomentBasis
    nkd = length(crystal.species_labels)
    length(spec.lmax_env) == nkd ||
        throw(ArgumentError("spec names $(length(spec.lmax_env)) species, crystal " *
                            "has $nkd"))
    (isfinite(tie_tol) && 0 <= tie_tol < _TIE_TOL_MAX) ||
        throw(ArgumentError("tie_tol must be in [0, $( _TIE_TOL_MAX)); got $tie_tol"))
    sg = analyze_symmetry(backend, crystal; tol = tol)
    sp = crystal.species
    nat = n_atoms(crystal)
    marked_atoms = [a for a = 1:nat if spec.marked[sp[a]]]
    isempty(marked_atoms) &&
        throw(ArgumentError("no atom of a marked species in the reference cell"))

    # clusters: bodies 1–2 from the ordinary enumeration at the pair radii (the
    # single edge IS the mark–env bond), 3-body stars from the pointed enumeration
    nl2 = build_neighbor_list(crystal, spec.cutoff_pair, MinimumImage();
                              tol = tie_tol)
    cs = build_clusters(crystal, nl2, sg; nbody = min(spec.nbody, 2))
    orbits = Tuple{Int,Int,ClusterOrbit}[]           # (body, per-body id, orbit)
    for body = 1:min(spec.nbody, 2)
        for (k, O) in enumerate(get(cs.by_body, body, ClusterOrbit[]))
            push!(orbits, (body, k, O))
        end
    end
    dmin2_star = Matrix{Float64}(undef, 0, 0)
    if spec.nbody >= 3
        nl3 = build_neighbor_list(crystal, spec.cutoff_star, MinimumImage();
                                  tol = tie_tol)
        stars = _pointed_star_candidates(crystal, nl3, spec)
        for (k, O) in enumerate(_orbits_from_members(crystal, sg, stars, 3))
            push!(orbits, (3, k, O))
        end
        dmin2_star = _dmin2_matrix(nl3, nat)
    end

    maxl = max(spec.lmax_mark, maximum(spec.lmax_env; init = 0))
    wcache = _build_wig_cache(sg, maxl)
    cart = cartesian_positions(crystal)
    A = Matrix(crystal.lattice.vectors)
    fac = 1.0 + tie_tol

    out = SALC[]
    recs = NamedTuple[]
    for (body, oid, O) in orbits
        labels = _moment_labels(spec, body)
        isempty(labels) && continue
        rep = O.representative
        # representative-site cartesian positions (image shifts included) → edges
        pos = [SVector{3,Float64}(cart[:, rep.atoms[s]]) +
               SVector{3,Float64}(A * Float64.(rep.shifts[s])) for s = 1:body]
        edges = [norm(pos[s] - pos[t]) for s = 1:body, t = 1:body]
        admit = function (t::Vector{SiteDecor})
            for s in eachindex(t)
                d = t[s]
                if is_marked(d)
                    # the mark: any marked species, its own ê-rank cap, and — for
                    # stars — every mark–env bond minimum-image within the radius
                    spec.marked[O.species[s]] || return false
                    d.spin_l <= spec.lmax_mark || return false
                    if body == 3
                        for u in eachindex(t)
                            u == s && continue
                            r = spec.cutoff_star[O.species[s], O.species[u]]
                            edges[s, u] <= r * fac || return false
                            edges[s, u]^2 <=
                                dmin2_star[rep.atoms[s], rep.atoms[u]] * fac^2 ||
                                return false
                        end
                    end
                else
                    # environment spins: sampled species only, per-species cap
                    d.spin_l >= 1 || return false
                    d.spin_l <= spec.lmax_env[O.species[s]] || return false
                end
            end
            return true
        end
        got = _orbit_salcs_decors(crystal, sg, body, oid, O, labels, wcache;
                                  isotropy = spec.isotropy, admit = admit)
        for s in got
            push!(out, s)
            push!(recs, (; body = body, species = Tuple(O.species),
                          edges = Tuple(sort([round(edges[a, b]; digits = 6)
                                              for a = 1:body for b = (a + 1):body])),
                          nmem = length(O.members)))
        end
    end
    isempty(out) &&
        throw(ArgumentError("the moment basis is empty: no pointed SALC survives " *
                            "the spec — check the cutoffs, lmax_mark / lmax_env, " *
                            "and that a marked species has admissible neighbors"))
    perm = sortperm(out; by = s -> s.key)
    keyvec = [out[j].key for j in perm]
    allunique(keyvec) || error("duplicate pointed SALC keys — enumeration bug")
    return MomentBasis(crystal, sg, spec, SALCBasis(out[perm], keyvec), recs[perm],
                       marked_atoms, Float64(tie_tol),
                       Ref{Union{Nothing,NamedTuple}}(nothing))
end

# ── evaluation (marked-column substitution) ────────────────────────────────────────

# One design entry: Φ_φ(a; e, ê) with the spin matrix's marked column substituted by
# the evaluation axis and the mark realized by the indicator field u = x̂ at `a`
# (the mark factor |u|²R₀₀ is direction-independent). Exact, not approximate: every
# pointed label carries exactly one mark, so a member not marked at `a` is killed by
# its |u_b|² = 0 factor before any substituted value could reach it.
function _eval_moment_entry(salc::SALC, esub::Matrix{Float64}, u::Matrix{Float64},
                            scratch::SALCScratch)::Float64
    return evaluate_salc(salc, esub, u, scratch)
end

"""
    _design_moment(mb, configs, axes; member_index = true) -> Matrix{Float64}

The moment channel's design matrix: rows are (configuration-major, marked-atom
minor) pairs `(c, a)` over `mb.marked_atoms`, columns the pointed SALCs. `axes[c]`
is the per-atom evaluation-axis matrix of configuration `c` (the D8 addendum's
mode rule resolves it at the dataset layer: mode 4 → the datum's `directions`,
mode 1 → its `constraint_axes`); only the marked column of `e` is substituted per
row, environment columns stay configuration coordinates in both modes.
`member_index = true` (default) evaluates each row through the mark→term index —
value-identical to the full per-SALC evaluation (`member_index = false`, the
in-tree oracle path); see `_mark_term_index`. Precondition, not checked here: the
columns of every `configs[c]` and the MARKED columns of every `axes[c]` are unit
vectors (the harmonic kernels assume it); the public doors enforce it —
`MomentDataset` validates every datum's `directions` (and the `SpinDatum` ctor its
`constraint_axes`), `predict_moment` validates `e` and the marked `axes` columns.
A `SpinDatum` built field-by-field carries NO direction check of its own, so a
new caller of this function is a new door and must validate first.
"""
function _design_moment(mb::MomentBasis, configs::Vector{Matrix{Float64}},
                        axes::Vector{Matrix{Float64}};
                        member_index::Bool = true)::Matrix{Float64}
    length(configs) == length(axes) ||
        throw(ArgumentError("$(length(configs)) configs, $(length(axes)) axes"))
    sal = salcs(mb)
    atoms = mb.marked_atoms
    nat = n_atoms(mb.crystal)
    nrow = length(configs) * length(atoms)
    X = Matrix{Float64}(undef, nrow, length(sal))
    idx = member_index ? _mark_term_index(sal, atoms) : nothing
    # Shape checks once, serially: a throw from inside the threaded loop surfaces
    # as a TaskFailedException wrapping the ArgumentError.
    for (ci, e) in enumerate(configs)
        size(e) == (3, nat) ||
            throw(ArgumentError("config $ci is $(size(e)), expected (3, $nat)"))
        size(axes[ci]) == (3, nat) ||
            throw(ArgumentError("axes $ci is $(size(axes[ci])), expected (3, $nat)"))
    end
    Threads.@threads for j = 1:length(sal)
        scratch = SALCScratch()
        esub = Matrix{Float64}(undef, 3, nat)
        u = zeros(3, nat)
        for (ci, e) in enumerate(configs)
            for (ai, a) in enumerate(atoms)
                copyto!(esub, e)
                esub[1, a] = axes[ci][1, a]
                esub[2, a] = axes[ci][2, a]
                esub[3, a] = axes[ci][3, a]
                u[1, a] = 1.0
                X[(ci - 1) * length(atoms) + ai, j] = idx === nothing ?
                    _eval_moment_entry(sal[j], esub, u, scratch) :
                    _eval_moment_members(sal[j], idx[j][ai], esub, u, scratch)
                u[1, a] = 0.0
            end
        end
    end
    return X
end

# ── fast design path: mark → term index ────────────────────────────────────────────

# For each pointed SALC and each marked atom, the (member, term) pairs whose mark
# sits on that atom: the design row (c, a) reads exactly these — every other term
# dies on its |u|² = 0 factor before contributing. The granularity is the TERM,
# not the member: a member's terms enumerate distinct decor→site assignments, so
# one member can carry its mark on different sites in different terms. Skipping
# the dead terms removes only exact-zero additions from each accumulated sum
# (live terms keep their relative order), so the fast path is VALUE-IDENTICAL to
# the full evaluation up to the sign of a zero; the equality gates in
# test_momentbasis.jl hold both paths to elementwise ==.
function _mark_term_index(sal::Vector{SALC},
                          atoms::Vector{Int})::Vector{Vector{Vector{Tuple{Int,Int}}}}
    pos = Dict(a => i for (i, a) in enumerate(atoms))
    idx = [[Tuple{Int,Int}[] for _ in atoms] for _ in sal]
    for (j, s) in enumerate(sal), (mi, m) in enumerate(s.members),
        (ti, t) in enumerate(m.terms)

        # One mark per pointed label (the marked-column-substitution invariant):
        # the FIRST DISP slot is therefore THE mark, the same slot every other
        # consumer (e.g. moment_resolvability) sees. The exact-zero skip argument
        # additionally needs the mark's radial power k ≥ 1 (0.0^0 == 1.0 would
        # make skipped terms nonzero) — pinned loudly, not assumed.
        site = 0
        nmarks = 0
        for sl in t.slots
            if sl.factor.channel == DISP
                sl.factor.k >= 1 ||
                    error("pointed SALC $j member $mi term $ti: mark radial " *
                          "power k = $(sl.factor.k) < 1 breaks the dead-term skip")
                nmarks += 1
                site = sl.site
            end
        end
        site == 0 && error("pointed SALC $j member $mi term $ti carries no mark")
        # The marked-column substitution is exact only with ONE mark per label;
        # every consumer of "the mark slot" in this file rests on it — asserted,
        # not assumed.
        nmarks == 1 || error("pointed SALC $j member $mi term $ti carries $nmarks " *
                             "marks; the marked-column substitution needs exactly one")
        a = m.atoms[site]
        # A mark on an atom outside `atoms` contributes to no design row in either
        # path (no row ever raises its |u|²); it is simply absent from the index.
        haskey(pos, a) && push!(idx[j][pos[a]], (mi, ti))
    end
    return idx
end

# Subset-sum evaluation of a pointed SALC over the terms marked at the row's
# atom — same scale and same relative accumulation order as
# `evaluate_salc(salc, e, u, scratch)` restricted to the live terms.
function _eval_moment_members(salc::SALC, tids::Vector{Tuple{Int,Int}},
                              esub::Matrix{Float64}, u::Matrix{Float64},
                              scratch::SALCScratch)::Float64
    n_spin = count(has_spin, salc.decors)
    scale = (4π)^(n_spin / 2)
    total = 0.0
    @inbounds for (mi, ti) in tids
        m = salc.members[mi]
        t = m.terms[ti]
        total += _eval_term_mixed(t.folded, t.slots, m.atoms, esub, u, scratch)
    end
    return scale * total
end

# ── pointed resolvability gate ─────────────────────────────────────────────────────

# Row identity of the pointed signature expansion WITHIN one marked atom's block:
# the mark factor's (l, μ) in the INDEPENDENT axis variable ê_a and the sorted
# environment factors (reference-cell atom, l, μ) in the configuration variable e.
# The marked reference-cell atom itself is the block index, not part of the key —
# the blocks of different marked atoms share no row. Periodic images fold onto
# reference-cell atoms — exactly the folding where a Wigner–Seitz tie's aliasing
# becomes visible.
struct _MomentRowKey
    mark_lm::Tuple{Int,Int}
    env::Vector{Tuple{Int,Int,Int}}
end
Base.hash(k::_MomentRowKey, h::UInt) = hash(k.env, hash(k.mark_lm, hash(:mrk, h)))
Base.:(==)(a::_MomentRowKey, b::_MomentRowKey) =
    a.mark_lm == b.mark_lm && a.env == b.env

"""
    moment_resolvability(mb; rtol = nothing) -> NamedTuple

The pointed periodic-resolvability gate (design record D9′): what THIS reference
cell can and cannot determine about the pointed columns. Returns

- `vanishing::Vector{Int}` — columns whose signature expansion cancels identically
  on cell-periodic data (judged per column against its own gross accumulation, the
  `_CANCELLATION_RTOL` convention);
- `rank::Int` and `null_combinations::Vector{Vector{Tuple{Int,Float64}}}` — the
  numerical rank of the kept signature block and, per flat direction, the column
  indices with their weights (the gate names the dependent columns, not just the
  count);
- `kept::Vector{Int}` — the columns that entered the rank computation (the
  complement of `vanishing`);
- `census::Vector` — per cluster orbit, `n_mark_atoms` = the number of DISTINCT
  reference-cell atoms that carry the mark somewhere in the orbit (what the code
  counts; not the number of stabilizer-inequivalent mark placements of one
  representative — with one species unmarked an Fe–Ge orbit still reports 4). It
  is bookkeeping for the face-(b) hazard (same unmarked cluster, different marks,
  whose columns can collapse to a determined sum under a boundary tie), not a
  decision criterion by itself.

Purely structural (the symbolic expansion, never sampled data); a full-rank result
certifies resolvability of the basis on this cell, not identifiability from any
particular training set.

The default-`rtol` (`rtol = nothing` → `1e-10`) result is cached on the basis (the
answer is a pure function of the basis, and the moment dataset layer runs this gate
at every construction — a train/held-out pair would otherwise pay the symbolic
expansion twice). Single-threaded, cached calls return the SAME object (`===`);
treat it as read-only. Passing an explicit `rtol` — the default value included —
always recomputes and is never cached. Concurrent first calls may compute twice and
race the store — the value is deterministic, so both compute equal (not `===`)
results and the race is benign. An [`UnclassifiableBasis`](@ref) refusal is
deliberately not cached: the gate re-throws loudly on every call.
"""
function moment_resolvability(mb::MomentBasis; rtol::Union{Nothing,Real} = nothing)
    if rtol === nothing
        mb.resolvability[] === nothing &&
            (mb.resolvability[] = _moment_resolvability(mb, 1e-10))
        return mb.resolvability[]::NamedTuple
    end
    return _moment_resolvability(mb, Float64(rtol))
end

function _moment_resolvability(mb::MomentBasis, rtol::Float64)
    sal = salcs(mb)
    # A member carrying two ENVIRONMENT spin slots on one reference-cell atom (two
    # periodic images of one neighbor) multiplies two harmonics of the SAME unit
    # vector, and the signature expansion — which keys factors as independent
    # monomials — then OVERCOUNTS the rank (the product reduces by the sphere's
    # Clebsch–Gordan relations; upstream measured 108 symbolic vs 98 actual on
    # the FeGe primitive cell). The mark factor is exempt: it reads the
    # independent axis variable ê, never the same sphere as an environment factor.
    for (j, s) in enumerate(sal), mem in s.members, t in mem.terms
        env_atoms = Int[]
        mark_site = 0
        for sl in t.slots
            sl.factor.channel == DISP && (mark_site = sl.site)
        end
        for sl in t.slots
            (sl.factor.channel == SPIN && sl.site != mark_site) || continue
            push!(env_atoms, mem.atoms[sl.site])
        end
        allunique(env_atoms) ||
            throw(UnclassifiableBasis("pointed SALC $j (key $(s.key)) has a member " *
                                      "with two environment spin factors on one " *
                                      "reference-cell atom (two periodic images of " *
                                      "one neighbor): the symbolic signature cannot " *
                                      "classify harmonic products on a single " *
                                      "sphere, so the gate refuses rather than " *
                                      "overcounting the rank. Use a reference cell " *
                                      "in which the images are distinct atoms"))
    end
    # The signature block S (rows = distinct signature keys, columns = pointed
    # SALCs) is never formed densely: on a supercell its row count scales as
    # (atoms × neighbours × harmonic components) — tens of millions of rows on a
    # 36×36 torus, i.e. tens of GB dense — while everything read off it (column
    # norms, rank, right singular vectors) is a function of SᵀS alone. Rows are
    # grouped by the marked atom (blocks of different marked atoms share no key),
    # so each block B_a (rows marked at atom a × all columns) is assembled densely,
    # its column norms accumulated, and it is folded into an upper-triangular R by
    # R ← qr([R; B_a]).R. Stacked QR is exact: R = Qᵀ S for an orthogonal Q, so
    # svd(R[:, kept]) has the singular values and right singular vectors of
    # S[:, kept]. Duplicate (key, column) entries within a block are SUMMED before
    # the norm is taken (the cancellation the vanishing test looks for), exactly as
    # the dense S[r, j] += w did. [Found on Miyazaki's 36×36 KLM torus, 2026-08-23:
    # the dense block was SIGKILLed at 8–16 GB on every basis past 201 columns.]
    nat = n_atoms(mb.crystal)
    c = length(sal)
    term_index = _mark_term_index(sal, collect(1:nat))    # (member, term) marked at each atom
    scales = [(4π)^(count(has_spin, s.decors) / 2) for s in sal]
    gross = zeros(c)
    colnorm2 = zeros(c)
    nrows_total = 0
    R = zeros(0, c)
    rows = Dict{_MomentRowKey,Int}()
    entries = Vector{Vector{Tuple{Int,Float64}}}()
    for a = 1:nat
        empty!(rows); empty!(entries)
        for j = 1:c
            s = sal[j]
            for (mi, ti) in term_index[j][a]
                mem = s.members[mi]
                t = mem.terms[ti]
                slots = t.slots
                mslot = findfirst(sl -> sl.factor.channel == DISP, slots)
                mslot === nothing && error("pointed SALC without a mark slot")
                mark_site = slots[mslot].site
                for index in CartesianIndices(t.folded)
                    w = t.folded[index] * scales[j]
                    w == 0.0 && continue
                    mark_lm = (0, 0)
                    env = Tuple{Int,Int,Int}[]
                    for i in eachindex(slots)
                        sl = slots[i]
                        if sl.factor.channel == DISP
                            continue                     # the mark's |u|²R₀₀ ≡ const
                        elseif sl.site == mark_site
                            mark_lm = (sl.factor.l, index[i] - sl.factor.l - 1)
                        else
                            push!(env, (mem.atoms[sl.site], sl.factor.l,
                                        index[i] - sl.factor.l - 1))
                        end
                    end
                    sort!(env)
                    key = _MomentRowKey(mark_lm, env)
                    r = get!(rows, key) do
                        push!(entries, Tuple{Int,Float64}[])
                        length(entries)
                    end
                    push!(entries[r], (j, w))
                    gross[j] += abs(w)
                end
            end
        end
        isempty(entries) && continue
        B = zeros(length(entries), c)
        for (r, es) in enumerate(entries), (j, w) in es
            B[r, j] += w
        end
        for j = 1:c
            colnorm2[j] += sum(abs2, @view B[:, j])
        end
        nrows_total += size(B, 1)
        W = vcat(R, B)                  # one copy; qr! factors it in place
        R = Matrix(qr!(W).R)
    end
    colnorm = sqrt.(colnorm2)
    vanishing = [j for j = 1:c if gross[j] > 0.0 && colnorm[j] <= _CANCELLATION_RTOL * gross[j]]
    kept = setdiff(1:c, vanishing)
    rank = 0
    null_combinations = Vector{Vector{Tuple{Int,Float64}}}()
    if !isempty(kept)
        F = svd(R[:, kept])
        cut = maximum(F.S; init = 0.0) * max(Float64(rtol),
                                             min(nrows_total, c) * eps(Float64))
        rank = count(>(cut), F.S)
        ck = length(kept)
        _push_comb!(v) = push!(null_combinations,
                               [(kept[t], v[t]) for t in eachindex(kept)
                                if abs(v[t]) > 1e-8])
        for q = (rank + 1):size(F.V, 2)
            _push_comb!(F.V[:, q])
        end
        # A WIDE S (more kept columns than signature rows) has flat directions
        # the economy SVD cannot list: its V spans only min(r, c) directions, and
        # the orthogonal complement of span(V) is null too. Found upstream by the
        # P1 face-(b) control (rank 20 of 74 kept, null_combinations empty — the
        # dataset door's dependency disclosure silently missed all 54). The
        # complement is read off a QR completion of V (never a full SVD: the
        # row side can be huge and its full U is never needed).
        if size(F.V, 2) < ck
            Qfull = qr(F.V).Q * Matrix{Float64}(I, ck, ck)
            for q = (size(F.V, 2) + 1):ck
                _push_comb!(Qfull[:, q])
            end
        end
    end
    # census: per cluster orbit (body, orbit_id), the stabilizer-inequivalent
    # admissible mark placements of the representative
    seen = Set{Tuple{Int,Int}}()
    census = NamedTuple[]
    for s in sal
        key = (s.key.body, s.key.orbit_id)
        key in seen && continue
        push!(seen, key)
        # reconstruct the orbit's mark classes from the SALCs sharing the orbit:
        # each member's marked atom, folded under "same reference-cell atom"
        marks = Set{Int}()
        for s2 in sal
            (s2.key.body, s2.key.orbit_id) == key || continue
            for mem in s2.members, t in mem.terms
                for sl in t.slots
                    sl.factor.channel == DISP || continue
                    push!(marks, mem.atoms[sl.site])
                end
            end
        end
        push!(census, (; body = key[1], orbit_id = key[2],
                        n_mark_atoms = length(marks)))
    end
    return (; vanishing, rank, kept, null_combinations, census)
end
