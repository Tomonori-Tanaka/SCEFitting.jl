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
    map_sym = _build_map_sym(crystal, ops, tol)
    translation_ops = [i for i in eachindex(ops) if ops[i].is_translation]
    return SpaceGroup(String(symbol), Int(number), ops, map_sym, translation_ops, Float64(tol))
end

# map_sym[iat, op] = the atom that iat maps to under (W|t), matched by minimum-image
# fractional distance within `tol` and the same species.
function _build_map_sym(crystal::Crystal, ops::Vector{SymOp}, tol::Real)::Matrix{Int}
    nat = n_atoms(crystal)
    x = crystal.frac_positions
    species = crystal.species
    tol2 = tol * tol
    map_sym = zeros(Int, nat, length(ops))
    @inbounds for (o, op) in enumerate(ops)
        W = op.rotation_frac
        t = op.translation_frac
        for iat = 1:nat
            xi = SVector{3,Float64}(x[1, iat], x[2, iat], x[3, iat])
            xn = W * xi + t
            found = 0
            for jat = 1:nat
                species[jat] == species[iat] || continue
                d2 = 0.0
                for k = 1:3
                    v = abs(xn[k] - x[k, jat]) % 1.0
                    v = min(v, 1.0 - v)
                    d2 += v * v
                end
                if d2 < tol2
                    found = jat
                    break
                end
            end
            found == 0 &&
                error("map_sym: atom $iat has no image under operation $o (tol=$tol)")
            map_sym[iat, o] = found
        end
    end
    return map_sym
end
