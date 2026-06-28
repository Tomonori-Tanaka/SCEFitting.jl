"""
    ClusterOrbit

A symmetry-equivalence class of cluster instances: `members` are all candidate
clusters related by the space group, `representative` is the canonical one,
`multiplicity = length(members)`, and `species` is the per-site species of the
representative.
"""
struct ClusterOrbit
    body::Int
    representative::ClusterMember
    members::Vector{ClusterMember}
    multiplicity::Int
    species::Vector{Int}
end

"""
    ClusterSet

Cluster orbits grouped by body order. See [`build_clusters`](@ref).
"""
struct ClusterSet
    by_body::Dict{Int,Vector{ClusterOrbit}}
end

# Image of site (atom, R) under operation `g`: the image atom `b = map_sym[atom,g]`
# carries lattice translation `τ_{g,atom} + W·R`, where τ closes the home-cell map.
@inline function _site_image(crystal::Crystal, sg::SpaceGroup, g::Int,
                             atom::Int, R::SVector{3,Int})
    op = sg.ops[g]
    b = sg.map_sym[atom, g]
    xa = SVector{3,Float64}(crystal.frac_positions[1, atom],
                            crystal.frac_positions[2, atom],
                            crystal.frac_positions[3, atom])
    xb = SVector{3,Float64}(crystal.frac_positions[1, b],
                            crystal.frac_positions[2, b],
                            crystal.frac_positions[3, b])
    τ = round.(Int, op.rotation_frac * xa + op.translation_frac - xb)
    Wi = round.(Int, op.rotation_frac)
    return b, SVector{3,Int}(τ + Wi * R)
end

# Canonical key of a cluster: the lexicographic minimum, over all symmetry images
# and all choices of anchor site, of the site list sorted as (atom, Rx, Ry, Rz).
# Two clusters share a key iff they lie in the same orbit.
function _canonical_key(crystal::Crystal, sg::SpaceGroup, m::ClusterMember)
    N = length(m.atoms)
    best = nothing
    for g = 1:n_ops(sg)
        imgs = ntuple(k -> _site_image(crystal, sg, g, m.atoms[k], m.shifts[k]), N)
        for s = 1:N
            R0 = imgs[s][2]
            sites = sort!([(imgs[k][1], imgs[k][2][1] - R0[1], imgs[k][2][2] - R0[2],
                            imgs[k][2][3] - R0[3]) for k = 1:N])
            key = Tuple(sites)
            if best === nothing || key < best
                best = key
            end
        end
    end
    return best
end

"""
    build_clusters(crystal, neighbors, spacegroup; nbody = 2, selection = MinimumImage())
        -> ClusterSet

Enumerate candidate clusters and reduce them to symmetry orbits under
`spacegroup`. Orbits are ordered deterministically by their canonical key.

# Arguments
- `crystal::Crystal`, `neighbors::NeighborList`, `spacegroup::SpaceGroup`.

# Keyword arguments
- `nbody::Integer = 2`: maximum body order.
- `selection::AbstractImageSelection = MinimumImage()`: the image-admissibility rule
  for the cluster edges; must match the one that built `neighbors` (the same rule the
  internal `candidate_clusters` enumeration applies).

# Returns
- `ClusterSet`: orbits grouped by body order.
"""
function build_clusters(crystal::Crystal, neighbors::NeighborList, spacegroup::SpaceGroup;
                        nbody::Integer = 2,
                        selection::AbstractImageSelection = MinimumImage())::ClusterSet
    cand = candidate_clusters(crystal, neighbors, nbody; selection = selection)
    by_body = Dict{Int,Vector{ClusterOrbit}}()
    for body in sort(collect(keys(cand)))
        groups = Dict{Any,Vector{ClusterMember}}()
        keyorder = Any[]
        for m in cand[body]
            key = _canonical_key(crystal, spacegroup, m)
            haskey(groups, key) || (groups[key] = ClusterMember[]; push!(keyorder, key))
            push!(groups[key], m)
        end
        orbits = ClusterOrbit[]
        for key in sort(keyorder)
            ms = groups[key]
            rep = ms[1]
            species = [crystal.species[a] for a in rep.atoms]
            push!(orbits, ClusterOrbit(body, rep, ms, length(ms), species))
        end
        by_body[body] = orbits
    end
    return ClusterSet(by_body)
end
