"""
    Crystal(lattice, frac_positions, species, species_labels) -> Crystal

A crystal structure: a [`Lattice`](@ref) plus an atomic basis.

# Arguments
- `lattice::Lattice`
- `frac_positions::AbstractMatrix{<:Real}`: 3 × n_atoms fractional coordinates.
- `species::AbstractVector{<:Integer}`: per-atom species index (1-based, indexing
  `species_labels`).
- `species_labels::AbstractVector{<:AbstractString}`: chemical label per species
  index, e.g. `["Fe", "Co"]`.

# Notes
Fractional coordinates are wrapped into `[0, 1)` along periodic directions so the
neighbor-list image bound `|Δfᵢ| < 1` holds (see [`build_neighbor_list`](@ref)).
The wrapping/validation lives in the inner constructor so that even an exactly
field-typed call (`Matrix{Float64}` / `Vector{Int}` / `Vector{String}`) is routed
through it rather than bypassing to a default constructor.
"""
struct Crystal
    lattice::Lattice
    frac_positions::Matrix{Float64}
    species::Vector{Int}
    species_labels::Vector{String}

    function Crystal(lattice::Lattice, frac_positions::AbstractMatrix{<:Real},
                    species::AbstractVector{<:Integer},
                    species_labels::AbstractVector{<:AbstractString})
        size(frac_positions, 1) == 3 ||
            throw(ArgumentError("`frac_positions` must be 3 × n_atoms"))
        nat = size(frac_positions, 2)
        length(species) == nat ||
            throw(ArgumentError("`species` length $(length(species)) ≠ n_atoms $nat"))
        smax = isempty(species) ? 0 : maximum(species)
        smax <= length(species_labels) ||
            throw(ArgumentError("species index $smax exceeds $(length(species_labels)) labels"))
        fr = Matrix{Float64}(frac_positions)
        @inbounds for a = 1:nat, d = 1:3
            if lattice.pbc[d]
                fr[d, a] = mod(fr[d, a], 1.0)
            end
        end
        return new(lattice, fr, collect(Int, species), collect(String, species_labels))
    end
end

"""
    n_atoms(crystal) -> Int

Number of atoms in the crystal basis.
"""
n_atoms(crystal::Crystal)::Int = size(crystal.frac_positions, 2)

"""
    cartesian_positions(crystal) -> Matrix{Float64}

Cartesian coordinates (3 × n_atoms, Å) of the atoms in the home cell.
"""
cartesian_positions(crystal::Crystal)::Matrix{Float64} =
    crystal.lattice.vectors * crystal.frac_positions
