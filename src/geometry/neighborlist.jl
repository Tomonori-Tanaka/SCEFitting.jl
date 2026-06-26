"""
    NeighborPair(i, j, shift, offset, distance)

A directed neighbor pair: atom `i` in the home cell and atom `j` in the periodic
image translated by integer lattice `shift` (cartesian `offset =
lattice.vectors * shift`), separated by `distance` (Å).

The integer `shift` (the inter-site lattice translation `R`) is retained so that
reciprocal-space / spin-spiral design rows — which need `R` to form the phase
factor `e^{i q·R}` — remain constructible later without re-deriving the geometry.
"""
struct NeighborPair
    i::Int
    j::Int
    shift::SVector{3,Int}
    offset::SVector{3,Float64}
    distance::Float64
end

"""
    NeighborList(cutoff, pairs)

All directed neighbor pairs within `cutoff` (Å). Built by
[`build_neighbor_list`](@ref).
"""
struct NeighborList
    cutoff::Float64
    pairs::Vector{NeighborPair}
end

"""
    build_neighbor_list(crystal, cutoff) -> NeighborList

Enumerate every directed atom pair `(i, j, shift)` with interatomic distance
`≤ cutoff` (Å), generalizing a fixed image grid to a cutoff-driven one.

# Notes
Along each periodic direction `d` the image range is `N_d = ceil(cutoff / spacingᵈ)`
with `spacingᵈ = 1/‖row_d(reciprocal)‖` the interplanar spacing; non-periodic
directions use no images. With fractional coordinates wrapped to `[0, 1)`, the
projection bound `|n_d + Δf_d| ≤ cutoff/spacingᵈ` guarantees this box encloses
every in-cutoff image (valid for triclinic cells too). Self-pairs (`i == j` at zero
shift) are excluded; both `(i,j,R)` and `(j,i,−R)` appear (directed list).
"""
function build_neighbor_list(crystal::Crystal, cutoff::Real)::NeighborList
    cutoff > 0 || throw(ArgumentError("`cutoff` must be positive"))
    cut = Float64(cutoff)
    lat = crystal.lattice
    A = lat.vectors
    nat = num_atoms(crystal)
    # N_d = ceil(cut / spacingᵈ) = ceil(cut * ‖b_d‖); zero on non-periodic axes.
    nrange = ntuple(d -> lat.pbc[d] ? ceil(Int, cut * norm(lat.reciprocal[d, :])) : 0, 3)
    cart = cartesian_positions(crystal)
    cut2 = cut^2
    pairs = NeighborPair[]
    @inbounds for i = 1:nat
        ri = SVector{3,Float64}(cart[1, i], cart[2, i], cart[3, i])
        for j = 1:nat
            rj0 = SVector{3,Float64}(cart[1, j], cart[2, j], cart[3, j])
            for n1 = -nrange[1]:nrange[1], n2 = -nrange[2]:nrange[2], n3 = -nrange[3]:nrange[3]
                (i == j && n1 == 0 && n2 == 0 && n3 == 0) && continue
                offset = A * SVector{3,Float64}(n1, n2, n3)
                d2 = sum(abs2, rj0 + offset - ri)
                if d2 <= cut2
                    push!(pairs, NeighborPair(i, j, SVector{3,Int}(n1, n2, n3),
                                              offset, sqrt(d2)))
                end
            end
        end
    end
    return NeighborList(cut, pairs)
end
