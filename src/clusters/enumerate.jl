"""
    ClusterMember

One concrete cluster instance: site `atoms[k]` lives in the periodic image
translated by integer lattice vector `shifts[k]`. By convention candidate clusters
are anchored with the first site in the home cell (`shifts[1] = 0`). The integer
`shifts` (the inter-site lattice translations `R`) are retained for later
reciprocal-space / spin-spiral design rows.
"""
struct ClusterMember
    atoms::Vector{Int}
    shifts::Vector{SVector{3,Int}}
end

body_order(m::ClusterMember)::Int = length(m.atoms)

"""
    candidate_clusters(crystal, neighbors, nbody) -> Dict{Int,Vector{ClusterMember}}

Enumerate candidate clusters per body order (`1:nbody`), anchored with the first
site in the home cell. 1-body clusters are the individual atoms; 2-body clusters
are the directed neighbor pairs of `neighbors`. (Body orders ≥ 3 are deferred.)
"""
function candidate_clusters(crystal::Crystal, neighbors::NeighborList,
                            nbody::Integer)::Dict{Int,Vector{ClusterMember}}
    nbody >= 1 || throw(ArgumentError("nbody must be ≥ 1"))
    nbody <= 2 ||
        throw(ArgumentError("only 1- and 2-body clusters are implemented (got nbody=$nbody)"))
    out = Dict{Int,Vector{ClusterMember}}()
    z = SVector{3,Int}(0, 0, 0)
    out[1] = [ClusterMember([a], [z]) for a = 1:num_atoms(crystal)]
    if nbody >= 2
        out[2] = [ClusterMember([p.i, p.j], [z, p.shift]) for p in neighbors.pairs]
    end
    return out
end
