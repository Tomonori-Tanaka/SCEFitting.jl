"""
    MagestyRebuild.VASP

VASP file adapters — a code-specific I/O namespace kept **outside** the code-agnostic
core. Nothing in the SCE pipeline depends on this module; it only *produces* the
core boundary objects: a [`Crystal`](@ref) (from a POSCAR) and
[`SpinDatum`](@ref)s (from constrained-noncollinear OSZICARs, via the
[`AbstractDFTSource`](@ref) / [`read_configs`](@ref) seam). Once you hold those, the
originating DFT code is irrelevant. Use as `MagestyRebuild.VASP.read_poscar(...)`
or `using MagestyRebuild.VASP: read_poscar, Oszicar`. A second DFT code is one more
sibling submodule — the core and its exports do not change.

Conventions follow VASP and Magesty.jl: the POSCAR scaling factor is a length scale
(or, if negative, a target cell volume); the OSZICAR moment is the `MW_int` column
(or `M_int` with `mint = true`) and the constraining field is `lambda*MW_perp`; both
are rotated from the `SAXIS` quantization frame into Cartesian by `Rz(α)·Ry(β)`. The
torque target is then `τ_a = −m_a × B_a` (eV).
"""
module VASP

using StaticArrays
using LinearAlgebra: norm, det
using ..MagestyRebuild: Crystal, Lattice, AbstractDFTSource, SpinDatum
import ..MagestyRebuild: read_configs

export read_poscar, write_poscar, Oszicar

# ── POSCAR ────────────────────────────────────────────────────────────────────

"""
    read_poscar(path) -> Crystal

Read a VASP POSCAR/CONTCAR file into a [`Crystal`](@ref). Handles the scaling factor
(a length scale, or a target volume when negative), `Direct`/`Cartesian` coordinates,
the optional `Selective dynamics` line, and POSCAR files with or without the element-
symbol line (synthesizing `"X1", "X2", …` when it is absent). Each of the three
lattice lines is one lattice vector (a column of `Lattice.vectors`).
"""
function read_poscar(path::AbstractString)::Crystal
    isfile(path) || throw(ArgumentError("POSCAR not found: $path"))
    lines = readlines(path)
    length(lines) >= 8 || throw(ArgumentError("POSCAR is too short to be valid"))

    scale = tryparse(Float64, strip(lines[2]))
    scale === nothing && throw(ArgumentError("invalid POSCAR scaling factor: $(lines[2])"))

    A_raw = MMatrix{3,3,Float64}(undef)
    for i = 1:3
        toks = split(strip(lines[2 + i]))
        length(toks) >= 3 || throw(ArgumentError("bad lattice vector on line $(2 + i)"))
        for j = 1:3
            A_raw[j, i] = parse(Float64, toks[j])   # column i = i-th lattice vector
        end
    end
    # Negative scale = target volume (VASP convention); otherwise a plain length scale.
    s = scale >= 0 ? scale : cbrt(abs(scale) / abs(det(SMatrix(A_raw))))
    A = SMatrix(A_raw) .* s

    toks6 = split(strip(lines[6]))
    counts = map(t -> tryparse(Int, t), toks6)
    local labels::Vector{String}, numbers::Vector{Int}, coordline::Int
    if all(!isnothing, counts) && all(c -> c > 0, counts)
        numbers = Int[c for c in counts]                # VASP4: no symbols line
        labels = String["X$i" for i = 1:length(numbers)]
        coordline = 7
    else
        labels = String.(toks6)                         # VASP5: symbols then counts
        cnt = map(t -> tryparse(Int, t), split(strip(lines[7])))
        all(!isnothing, cnt) && all(c -> c > 0, cnt) ||
            throw(ArgumentError("bad atom counts on POSCAR line 7: $(lines[7])"))
        numbers = Int[c for c in cnt]
        coordline = 8
    end
    length(labels) == length(numbers) ||
        throw(ArgumentError("POSCAR: $(length(labels)) species symbols vs $(length(numbers)) counts"))

    mode = lowercase(strip(lines[coordline]))
    if startswith(mode, "s")                            # optional Selective dynamics line
        coordline += 1
        coordline <= length(lines) ||
            throw(ArgumentError("POSCAR ends after 'Selective dynamics' (no coordinate-mode line)"))
        mode = lowercase(strip(lines[coordline]))
    end
    cartesian = startswith(mode, "c") || startswith(mode, "k")
    startswith(mode, "d") || cartesian ||
        throw(ArgumentError("invalid POSCAR coordinate mode: $(lines[coordline])"))

    nat = sum(numbers)
    pos = Matrix{Float64}(undef, 3, nat)
    first = coordline + 1
    for a = 1:nat
        li = first + a - 1
        li <= length(lines) || throw(ArgumentError("POSCAR ends before all $nat positions"))
        toks = split(strip(lines[li]))
        length(toks) >= 3 || throw(ArgumentError("bad position on POSCAR line $li"))
        for j = 1:3
            pos[j, a] = parse(Float64, toks[j])
        end
    end
    # Cartesian coords carry the same scaling; converting to fractional cancels it.
    frac = cartesian ? (A \ (s .* pos)) : pos

    species = Vector{Int}(undef, nat)
    a = 0
    for (si, c) in enumerate(numbers), _ = 1:c
        species[a += 1] = si
    end
    return Crystal(Lattice(A), frac, species, labels)
end

"""
    write_poscar(path, crystal; comment = "…", cartesian = false)

Write `crystal` to a VASP POSCAR file (scaling factor `1.0`, atoms grouped by species
in `species_labels` order, `Direct` coordinates unless `cartesian = true`). Species
with no atoms are omitted from the symbol/count lines. Note: atoms are **reordered** into
species groups, so any external per-atom array (forces, charges, …) indexed by the
original atom order must be permuted to match.
"""
function write_poscar(path::AbstractString, crystal::Crystal;
                      comment::AbstractString = "generated by MagestyRebuild",
                      cartesian::Bool = false)
    A = crystal.lattice.vectors
    nsp = length(crystal.species_labels)
    groups = [findall(==(s), crystal.species) for s = 1:nsp]
    present = [s for s = 1:nsp if !isempty(groups[s])]
    row(v) = join((string(v[k]) for k = 1:3), "  ")
    open(path, "w") do io
        println(io, comment)
        println(io, "1.0")
        for i = 1:3
            println(io, "  ", row(@view A[:, i]))       # line i = i-th lattice vector
        end
        println(io, join((crystal.species_labels[s] for s in present), " "))
        println(io, join((length(groups[s]) for s in present), " "))
        println(io, cartesian ? "Cartesian" : "Direct")
        for s in present, a in groups[s]
            f = SVector{3,Float64}(crystal.frac_positions[1, a], crystal.frac_positions[2, a],
                                   crystal.frac_positions[3, a])
            println(io, "  ", row(cartesian ? A * f : f))
        end
    end
    return path
end

# ── OSZICAR (constrained-noncollinear training data) ──────────────────────────

"""
    Oszicar(paths; saxis = [0, 0, 1], energy_kind = :free, mint = false)

An [`AbstractDFTSource`](@ref) over one or more VASP OSZICAR files — each contributes
one [`SpinDatum`](@ref) (in the given order). `paths` may be a single path or a vector.

# Keyword arguments
- `saxis`: the `SAXIS` quantization axis; moments and fields are rotated from this
  frame into Cartesian coordinates by `Rz(α)·Ry(β)` (default `[0,0,1]` = identity).
- `energy_kind`: `:free` (the `F=` free energy) or `:sigma0` (`E0`, energy σ→0).
- `mint`: read the moment from the `M_int` columns instead of `MW_int`. The constraining
  field `lambda*MW_perp` is referenced to the `MW_int` moment, so the default
  `mint = false` keeps the moment and field (hence the torque) mutually consistent.
"""
struct Oszicar <: AbstractDFTSource
    paths::Vector{String}
    saxis::SVector{3,Float64}
    energy_kind::Symbol
    mint::Bool
end

function Oszicar(paths::AbstractVector{<:AbstractString};
                saxis = SVector{3,Float64}(0, 0, 1),
                energy_kind::Symbol = :free, mint::Bool = false)
    energy_kind in (:free, :sigma0) ||
        throw(ArgumentError("energy_kind must be :free or :sigma0; got $(repr(energy_kind))"))
    return Oszicar(collect(String, paths), SVector{3,Float64}(saxis), energy_kind, mint)
end
Oszicar(path::AbstractString; kwargs...) = Oszicar([path]; kwargs...)

# SAXIS quantization-axis rotation Rz(α)·Ry(β), α/β the azimuth/polar of `saxis`.
function _saxis_rotation(saxis::SVector{3,Float64})::SMatrix{3,3,Float64,9}
    n = norm(saxis)
    (isfinite(n) && n > 0) ||
        throw(ArgumentError("saxis must be a nonzero, finite vector; got $saxis"))
    α = atan(saxis[2], saxis[1])
    β = atan(sqrt(saxis[1]^2 + saxis[2]^2), saxis[3])
    Rz = @SMatrix [cos(α) -sin(α) 0.0; sin(α) cos(α) 0.0; 0.0 0.0 1.0]
    Ry = @SMatrix [cos(β) 0.0 sin(β); 0.0 1.0 0.0; -sin(β) 0.0 cos(β)]
    return Rz * Ry
end

# Final-step energy and per-atom moment vectors (3 × n_atoms) from one OSZICAR.
# The moment block is the `ion … MW_int … M_int` table (7 columns: idx + MW_int xyz +
# M_int xyz); it ends at a `:` line or any non-7-column line (e.g. the `… F= … E0= …`
# summary). The last committed block / last `F=` line (final ionic step) wins.
function _oszicar_energy_moments(path::AbstractString, energy_kind::Symbol, mint::Bool)
    energy = 0.0
    found_e = false
    moments = SVector{3,Float64}[]
    tmp = SVector{3,Float64}[]
    collecting = false
    cols = mint ? (5, 6, 7) : (2, 3, 4)
    for line in eachline(path)
        if occursin("M_int", line)
            collecting = true
            empty!(tmp)
            continue
        end
        if collecting && (occursin(":", line) || length(split(line)) != 7)
            collecting = false
            moments = copy(tmp)
        end
        if collecting
            p = split(line)
            push!(tmp, SVector{3,Float64}(parse(Float64, p[cols[1]]),
                                          parse(Float64, p[cols[2]]),
                                          parse(Float64, p[cols[3]])))
        end
        if occursin("F=", line)
            # take the token right after the keyword (robust to other fields shifting),
            # rather than a fixed column index.
            p = split(line)
            key = energy_kind === :free ? "F=" : "E0="
            k = findfirst(==(key), p)
            if k !== nothing && k < length(p)
                v = tryparse(Float64, p[k + 1])
                v === nothing || (energy = v; found_e = true)
            end
        end
    end
    collecting && (moments = copy(tmp))                 # block running at EOF
    found_e || throw(ArgumentError("no energy (F= line) found in $path"))
    isempty(moments) &&
        throw(ArgumentError("no magnetic-moment (M_int) block found in $path"))
    return energy, reduce(hcat, moments)                # 3 × n_atoms
end

# Per-atom constraining field (3 × n_atoms) from the `lambda*MW_perp` block; rows are
# `idx Bx By Bz`, only for constrained atoms (others stay zero). The last block wins;
# a missing block (unconstrained run) yields a zero field.
function _oszicar_field(path::AbstractString, nat::Int)::Matrix{Float64}
    field = zeros(3, nat)
    tmp = zeros(3, nat)
    in_block = false
    got = false
    for line in eachline(path)
        if occursin("lambda*MW_perp", line)
            in_block = true
            fill!(tmp, 0.0)
            got = false
            continue
        end
        if in_block
            p = split(line)
            idx = length(p) == 4 ? tryparse(Int, p[1]) : nothing
            if idx !== nothing && 1 <= idx <= nat
                tmp[1, idx] = parse(Float64, p[2])
                tmp[2, idx] = parse(Float64, p[3])
                tmp[3, idx] = parse(Float64, p[4])
                got = true
                continue
            else
                in_block = false
                got && (field .= tmp)
            end
        end
    end
    in_block && got && (field .= tmp)                   # block running at EOF
    return field
end

function read_configs(src::Oszicar)::Vector{SpinDatum}
    R = _saxis_rotation(src.saxis)
    data = Vector{SpinDatum}(undef, length(src.paths))
    for (i, path) in enumerate(src.paths)
        isfile(path) || throw(ArgumentError("OSZICAR not found: $path"))
        energy, moments = _oszicar_energy_moments(path, src.energy_kind, src.mint)
        field = _oszicar_field(path, size(moments, 2))
        data[i] = SpinDatum(energy, R * moments, R * field)   # SAXIS → Cartesian frame
    end
    return data
end

end # module VASP
