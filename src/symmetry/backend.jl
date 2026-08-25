"""
    AbstractSymmetryBackend

Strategy for obtaining a crystal's space-group operations. Implement a backend by
subtyping this and defining

    analyze_symmetry(::MyBackend, crystal::Crystal; tol) -> SpaceGroup

Backends only need to supply raw `(rotations, translations, symbol, number)`; a
shared in-tree assembler (`_assemble_spacegroup`) derives Cartesian rotations,
properness, and the atom permutation table `map_sym`, so every backend shares one
convention.
"""
abstract type AbstractSymmetryBackend end

"""
    NoSymmetry()

In-tree fallback backend: returns the trivial P1 group (identity only). Lets the
whole pipeline run with no external dependency; the SALC projector still averages
over `{E} × {time reversal}`.
"""
struct NoSymmetry <: AbstractSymmetryBackend end

"""
    SpglibBackend()

Space-group analysis via Spglib. The type lives in the core package (so it can be
named without Spglib loaded), but `analyze_symmetry(::SpglibBackend, …)` is
provided only when `Spglib` is loaded (`using Spglib`).
"""
struct SpglibBackend <: AbstractSymmetryBackend end

"""
    analyze_symmetry(backend, crystal; tol = 1e-5) -> SpaceGroup

Compute the space group of `crystal` using `backend`.

A backend analyzes the cell as a fully periodic 3D crystal — Spglib is not told about
`Lattice(...; pbc)` and has no way to be. When the crystal declares an aperiodic axis,
the assembled group is intersected with what that declaration allows: operations that
close only through the periodicity along the aperiodic axis are dropped (with a warning
naming how many), and the translations of the survivors are re-seated on the
representative the finite structure actually has, so no cell shift is ever produced
along an axis that has no cells. The result is a subgroup, so the basis built from it is
larger than the fully periodic one, never short. A fully periodic crystal is untouched.
"""
function analyze_symmetry(backend::AbstractSymmetryBackend, crystal::Crystal;
                          tol::Real = 1e-5)::SpaceGroup
    error("analyze_symmetry has no method for $(typeof(backend)); load the backend " *
          "package first (e.g. `using Spglib` for SpglibBackend).")
end

function analyze_symmetry(::NoSymmetry, crystal::Crystal; tol::Real = 1e-5)::SpaceGroup
    return _assemble_spacegroup(crystal, [SMatrix{3,3,Float64}(I)],
                                [SVector{3,Float64}(0, 0, 0)], "P1", 1; tol = tol)
end

"""
    _assemble_spacegroup(crystal, rotations, translations, symbol, number; tol) -> SpaceGroup

Shared assembler: turn raw fractional `rotations`/`translations` into `SymOp`s
(deriving Cartesian rotations, properness, pure-translation flags) and derive
`map_sym` in-tree. Backends call this so they all share one convention.

The assembled set is validated as a group (see `_validate_ops`) before it is
returned.
"""
function _assemble_spacegroup(crystal::Crystal,
                              rotations::AbstractVector,
                              translations::AbstractVector,
                              symbol::AbstractString, number::Integer;
                              tol::Real)::SpaceGroup
    length(rotations) == length(translations) ||
        throw(ArgumentError("rotations and translations must have equal length"))
    A = crystal.lattice.vectors
    Ainv = crystal.lattice.reciprocal
    ops = Vector{SymOp}(undef, length(rotations))
    @inbounds for i in eachindex(rotations)
        W = SMatrix{3,3,Float64}(rotations[i])
        Rcart = A * W * Ainv
        t = translations[i]
        tfrac = SVector{3,Float64}(ntuple(k -> abs(t[k]) >= tol ? Float64(t[k]) : 0.0, 3))
        ops[i] = SymOp(W, Rcart, tfrac, det(Rcart) > 0, isapprox(W, I; atol = tol))
    end
    _validate_ops(ops, tol)
    ops, symbol = _restrict_to_pbc(crystal, ops, String(symbol), tol)
    map_sym = _build_map_sym(crystal, ops, tol)
    translation_ops = [i for i in eachindex(ops) if ops[i].is_translation]
    return SpaceGroup(symbol, Int(number), ops, map_sym, translation_ops, Float64(tol))
end

"""
    _restrict_to_pbc(crystal, ops, symbol, tol) -> (ops, symbol)

Keep only the operations that are symmetries of the crystal under its DECLARED
periodicity, returning them with a symbol that says so.

Backends analyse the cell as a fully periodic 3D crystal — Spglib is not told about
`Lattice(...; pbc)` and has no way to be. For a slab (`pbc = (true, true, false)`) that
list can therefore contain operations that only close through the artificial
periodicity along the aperiodic axis. They are not harmless: `build_clusters` demands
that the candidate set be closed under the group, while the neighbour list refuses to
emit an image along an aperiodic axis, so such an operation shows up as a closure
failure blamed on image selection and the tie tolerance.

The kept set is a **subgroup**: the aperiodic axes are matched without folding, so the
condition composes (`g` maps `x_i` to `x_j` exactly, `h` maps `x_j` to `x_l` exactly,
hence `g∘h` maps `x_i` to `x_l` exactly), the identity satisfies it, and an operation
that permutes the atoms bijectively has an inverse that does too. Using a subgroup
never over-reduces: the basis it builds is larger than necessary, never short.

A fully periodic crystal keeps every operation and never reaches the filter.
"""
function _restrict_to_pbc(crystal::Crystal, ops::Vector{SymOp}, symbol::String,
                          tol::Real)::Tuple{Vector{SymOp},String}
    pbc = crystal.lattice.pbc
    all(pbc) && return ops, symbol
    kept = SymOp[]
    ndrop = 0
    for (i, op) in enumerate(ops)
        r = _op_permutation(crystal, op, tol, i, false)
        if r === nothing
            ndrop += 1
            continue
        end
        # Re-seat the translation on the representative that is the physical
        # operation, so every later consumer — `_build_map_sym`, and `_site_image`,
        # whose cell shift is `round(W·x_a + t − x_b)` — sees an exact match and a
        # zero shift along the aperiodic axes.
        push!(kept, SymOp(op.rotation_frac, op.rotation_cart,
                          op.translation_frac - r[2], op.is_proper, op.is_translation))
    end
    # `kept` is returned even when nothing was dropped: the re-seated translations are
    # the point on their own, and leaving them at the backend's mod-1 representative
    # would keep handing `_site_image` a cell shift along an axis that has no cells.
    ndrop == 0 && return kept, symbol
    axes = join((("a", "b", "c")[k] for k = 1:3 if !pbc[k]), ", ")
    @warn "symmetry: $ndrop of $(length(ops)) operations " *
          "reported for the fully periodic cell close only through the periodicity " *
          "along $axes, which this lattice declares aperiodic — dropping them. The " *
          "basis is built from the remaining subgroup, so it is larger than the " *
          "fully periodic one, never short. If the cell is a vacuum-padded slab and " *
          "you meant the 3D group of that padded cell, pass the default " *
          "`pbc = (true, true, true)`: the vacuum and the cutoff already keep the " *
          "images apart." group = symbol kept = length(kept)
    _validate_ops(kept, tol)
    return kept, symbol * " (pbc subgroup)"
end

# The largest operation count any conventional crystallographic setting can reach:
# 48 point operations × 4 centering translations (F). Anything above it is a
# supercell, whose extra operations are pure translations.
const _MAX_CLOSURE_OPS = 192

# Are two fractional translations equal modulo a lattice vector?
_tclose(a::SVector{3,Float64}, b::SVector{3,Float64}, tol::Real) =
    all(k -> (d = abs(a[k] - b[k]) % 1.0; min(d, 1.0 - d) <= tol), 1:3)

# Index of the operation `(W|t)` in `ops` (0 if absent). `bykey` groups operation
# indices by their integer rotation, so only the same-rotation coset is scanned.
function _find_op(ops::Vector{SymOp}, bykey::Dict{SMatrix{3,3,Int,9},Vector{Int}},
                  W::SMatrix{3,3,Int,9}, t::SVector{3,Float64}, tol::Real)::Int
    idx = get(bykey, W, nothing)
    idx === nothing && return 0
    for i in idx
        _tclose(ops[i].translation_frac, t, tol) && return i
    end
    return 0
end

"""
    _validate_ops(ops, tol)

Check that an assembled operation list really is a space group of the lattice.
[Backported from SLCE.jl f8a529d.]

Backends are trusted for *which* group a crystal has, not for handing back a
well-formed one — and downstream code consumes a `SpaceGroup` as a group, not as
a list: `build_clusters` reduces orbits by stabilizer counting and
`build_salc_basis` projects with `(1/|G|) Σ_g`. A set that is merely a plausible
list of matrices therefore does not fail; it silently yields wrong orbit
multiplicities and a projector that is not idempotent. (Measured on a
deliberately non-closed set: `sum(multiplicity) = 8` over `6` candidate
clusters.)

Checked, in increasing cost: each `W` integral with `|det W| = 1` and a Cartesian
image that is orthogonal (i.e. a symmetry of *this* lattice's metric, which
integrality alone does not imply — a shear is integral too); `(I|0)` present; no
two operations equal modulo a lattice translation; the rotation parts closed
under multiplication (always cheap — at most 48 distinct rotations); every
operation's inverse present. The full pairwise `(W|t)` closure runs for up to
`_MAX_CLOSURE_OPS` operations, which covers every conventional setting; beyond
that (supercells) the preceding checks stand in for it.
"""
function _validate_ops(ops::Vector{SymOp}, tol::Real)
    n = length(ops)
    n == 0 && throw(ArgumentError("the operation list is empty; a space group " *
                                  "contains at least the identity"))
    otol = max(Float64(tol), 1e-8)
    Wint = Vector{SMatrix{3,3,Int,9}}(undef, n)
    for i = 1:n
        W = ops[i].rotation_frac
        isapprox(W, round.(W); atol = tol) || throw(ArgumentError(
            "symmetry op $i: the fractional rotation is not an integer matrix " *
            "($(W)) — a space-group operation maps the lattice onto itself"))
        Wi = SMatrix{3,3,Int,9}(round.(Int, W))
        dW = round(Int, det(Wi))             # exact: Wi is integral
        abs(dW) == 1 || throw(ArgumentError(
            "symmetry op $i: |det W| = $(abs(dW)) ≠ 1 — the operation is not a " *
            "bijection of the lattice"))
        R = ops[i].rotation_cart
        maximum(abs, R' * R - I) <= otol || throw(ArgumentError(
            "symmetry op $i: the Cartesian rotation is not orthogonal to $otol — " *
            "the fractional matrix is integral but does not preserve this " *
            "lattice's metric, so it is not a symmetry of the crystal"))
        Wint[i] = Wi
    end

    bykey = Dict{SMatrix{3,3,Int,9},Vector{Int}}()
    for i = 1:n
        push!(get!(Vector{Int}, bykey, Wint[i]), i)
    end
    for (_, idx) in bykey, a in eachindex(idx), b in (a + 1):length(idx)
        _tclose(ops[idx[a]].translation_frac, ops[idx[b]].translation_frac, tol) &&
            throw(ArgumentError(
                "symmetry ops $(idx[a]) and $(idx[b]) are the same operation modulo " *
                "a lattice translation — duplicates inflate every orbit multiplicity"))
    end

    id = SMatrix{3,3,Int,9}(I)
    _find_op(ops, bykey, id, zero(SVector{3,Float64}), tol) == 0 &&
        throw(ArgumentError("the identity (I|0) is missing from the operation list"))

    rots = Set(Wint)
    for a in rots, b in rots
        a * b in rots || throw(ArgumentError(
            "the rotation parts are not closed under multiplication ($(a) ∘ $(b) is " *
            "absent) — the operation list is not a group"))
    end

    for i = 1:n
        Wi = SMatrix{3,3,Int,9}(round.(Int, inv(ops[i].rotation_frac)))
        ti = -(SMatrix{3,3,Float64}(Wi) * ops[i].translation_frac)
        _find_op(ops, bykey, Wi, ti, tol) == 0 && throw(ArgumentError(
            "symmetry op $i has no inverse in the operation list — the list is not " *
            "a group"))
    end

    n <= _MAX_CLOSURE_OPS || return nothing
    for i = 1:n, j = 1:n
        W = Wint[i] * Wint[j]
        t = ops[i].rotation_frac * ops[j].translation_frac + ops[i].translation_frac
        _find_op(ops, bykey, W, t, tol) == 0 && throw(ArgumentError(
            "the composition of symmetry ops $i ∘ $j is not in the operation list — " *
            "the list is not closed, so orbit multiplicities and the SALC projector " *
            "are both wrong"))
    end
    return nothing
end

# Fractional separation along axis `k`, folded to the minimum image ONLY when that
# axis is periodic. On an APERIODIC axis there is no period to fold with: `f = ε` and
# `f = 1 − ε` are `|a_k|` apart, and folding them together is precisely how a symmetry
# that exists only by virtue of the cell's artificial periodicity gets mistaken for a
# symmetry of the declared structure.
@inline function _frac_sep(d::Float64, periodic::Bool)::Float64
    periodic || return abs(d)
    v = abs(d) % 1.0
    return min(v, 1.0 - v)
end

# Does `W` keep the periodic and aperiodic axes in separate blocks? An operation that
# sends an in-plane lattice translation into the aperiodic direction (or the reverse)
# does not map the *declared* translation lattice onto itself — there are no lattice
# vectors along an aperiodic axis to receive it — so it cannot be a symmetry however
# well it happens to permute the atoms.
@inline function _axis_blocks_ok(W::SMatrix{3,3,Float64,9}, pbc::SVector{3,Bool},
                                 tol::Real)::Bool
    @inbounds for k = 1:3, m = 1:3
        pbc[k] == pbc[m] && continue
        abs(W[k, m]) <= tol || return false
    end
    return true
end

# The atom permutation induced by `op`, matched by fractional distance within `tol`,
# together with the integer translation the APERIODIC axes have to be shifted by for
# the match to hold exactly. Returns `nothing` when the operation is not a symmetry
# under the crystal's declared periodicity.
#
# Backends report `t` modulo a lattice translation, and along a periodic axis that is
# the whole truth. Along an aperiodic axis it is not: the coset has no lattice to be
# taken modulo, so ONE representative is the physical operation and the others are
# not. A slab centred at `z = 1/2` has a mirror there, and Spglib is entitled to hand
# it back as the mirror at `z = 0` with `t_z = 0` — the same operation modulo `c`, and
# the only one of the two that the finite slab does not have. So the test is not
# "does `t` match exactly" but "is there an integer shift of `t` along the aperiodic
# axes that makes every atom match exactly" — a single shift for the whole operation,
# which is what makes the surviving set a group and forces `_site_image` to produce
# zero cell shifts along those axes.
#
# `strict` (a fully periodic crystal) short-circuits all of this: nothing may be
# dropped there, so a backend handing back an operation the crystal does not have is a
# fault, and "no image" stays the error it has always been.
function _op_permutation(crystal::Crystal, op::SymOp, tol::Real, o::Integer,
                         strict::Bool)
    nat = n_atoms(crystal)
    x = crystal.frac_positions
    species = crystal.species
    pbc = crystal.lattice.pbc
    tol2 = tol * tol
    strict || _axis_blocks_ok(op.rotation_frac, pbc, tol) || return nothing
    col = zeros(Int, nat)
    hit = falses(nat)                      # each column must be a *permutation*
    dt = zero(MVector{3,Float64})          # the aperiodic-axis representative shift
    W = op.rotation_frac
    t = op.translation_frac
    @inbounds for iat = 1:nat
        xi = SVector{3,Float64}(x[1, iat], x[2, iat], x[3, iat])
        xn = W * xi + t
        found = 0
        for jat = 1:nat
            species[jat] == species[iat] || continue
            d2 = 0.0
            for k = 1:3
                v = abs(xn[k] - x[k, jat]) % 1.0
                d2 += min(v, 1.0 - v)^2
            end
            if d2 < tol2
                found = jat
                break
            end
        end
        if found == 0
            strict || return nothing
            error("map_sym: atom $iat has no image under operation $o (tol=$tol)")
        end
        # Two atoms sharing an image means the operation does not map the crystal
        # onto itself at this tolerance. `SpaceGroup` documents these columns as
        # permutations and the orbit code inverts them, so catch it here — that is a
        # malformed operation either way, never something to filter away quietly.
        if hit[found]
            prev = findfirst(==(found), view(col, 1:(iat - 1)))
            error("map_sym: atoms $prev and $iat share the image $found under " *
                  "operation $o (tol=$tol) — the column is not a permutation")
        end
        hit[found] = true
        col[iat] = found
        strict && continue
        # One shift for the whole operation: the first atom fixes it, the rest must
        # agree. Disagreement means the atoms are only each other's images by hopping
        # different numbers of cells along an axis that has none.
        for k = 1:3
            pbc[k] && continue
            n = round(xn[k] - x[k, found])
            iat == 1 ? (dt[k] = n) : (dt[k] == n || return nothing)
        end
    end
    return col, SVector{3,Float64}(dt)
end

# map_sym[iat, op] = the atom that iat maps to under (W|t). Every operation reaching
# here has already survived `_restrict_to_pbc`, so `_op_permutation` cannot come back
# empty; `strict = true` states that — a missing image at this point is a fault.
function _build_map_sym(crystal::Crystal, ops::Vector{SymOp}, tol::Real)::Matrix{Int}
    map_sym = Matrix{Int}(undef, n_atoms(crystal), length(ops))
    @inbounds for (o, op) in enumerate(ops)
        map_sym[:, o] = first(_op_permutation(crystal, op, tol, o, true))
    end
    return map_sym
end
