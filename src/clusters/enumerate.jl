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

# All ordered `k`-tuples of distinct indices drawn from `1:n` (permutations of
# `k`-subsets). Small `k` (= body order − 1) only; cheap.
function _ordered_subsets(n::Int, k::Int)::Vector{Vector{Int}}
    k == 0 && return [Int[]]
    out = Vector{Int}[]
    for sub in _ordered_subsets(n, k - 1)
        for i = 1:n
            i in sub && continue
            push!(out, vcat(sub, i))
        end
    end
    return out
end

# Per-atom-pair minimum-image squared distance, read off a neighbor list (Inf where
# the pair has no in-list image — i.e. its minimum image is beyond the cutoff).
function _dmin2_matrix(neighbors::NeighborList, nat::Int)::Matrix{Float64}
    dmin2 = fill(Inf, nat, nat)
    @inbounds for p in neighbors.pairs
        d2 = p.distance^2
        d2 < dmin2[p.i, p.j] && (dmin2[p.i, p.j] = d2)
    end
    return dmin2
end

"""
    candidate_clusters(crystal, neighbors, nbody; selection = MinimumImage())
        -> Dict{Int,Vector{ClusterMember}}

Enumerate candidate clusters per body order (`1:nbody`), anchored with the first
site in the home cell. 1-body clusters are the individual atoms; an `N`-body
cluster (`N ≥ 2`) is `N` distinct sites whose every pairwise **edge is admissible**,
generated as ordered tuples `[a, …]` anchored at a home-cell atom `a` with the
remaining sites drawn from `a`'s neighbors. For `N = 2` this is exactly the
directed neighbor pairs.

An edge `(ax, ay)` is admissible under [`MinimumImage`](@ref) iff it sits at the
*minimum-image distance* of atom pair `(ax, ay)` (within relative tolerance) — so a
clique cannot be built from an aliased, non-minimum-image bond — and under
[`AllImages`](@ref) iff it is within the radial cutoff. The `selection` should match
the one that built `neighbors`. For a cutoff below half the smallest perpendicular
cell width the two rules coincide (each in-cutoff image is already the minimum one).
"""
function candidate_clusters(crystal::Crystal, neighbors::NeighborList,
                            nbody::Integer;
                            selection::AbstractImageSelection =
                                MinimumImage())::Dict{Int,Vector{ClusterMember}}
    nbody >= 1 || throw(ArgumentError("nbody must be ≥ 1"))
    nat = n_atoms(crystal)
    z = SVector{3,Int}(0, 0, 0)
    out = Dict{Int,Vector{ClusterMember}}()
    out[1] = [ClusterMember([a], [z]) for a = 1:nat]
    nbody == 1 && return out

    # Per-home-atom neighbor sites (the directed neighbor list already has `i` in the
    # home cell), and a cartesian-position helper for the pairwise edge check.
    A = crystal.lattice.vectors
    cart = cartesian_positions(crystal)
    neigh = [Tuple{Int,SVector{3,Int}}[] for _ = 1:nat]
    for p in neighbors.pairs
        push!(neigh[p.i], (p.j, p.shift))
    end
    sitepos(b::Int, R::SVector{3,Int}) =
        SVector{3,Float64}(cart[1, b], cart[2, b], cart[3, b]) +
        A * SVector{3,Float64}(R[1], R[2], R[3])

    # Edge admissibility (see the docstring). MinimumImage compares each edge to its
    # atom-pair minimum-image distance; AllImages to the radial cutoff. Tolerances are
    # relative so a degenerate shell at the boundary is never split.
    fac = (1 + _SAME_DIST_RTOL)^2
    use_minimage = selection isa MinimumImage
    dmin2 = use_minimage ? _dmin2_matrix(neighbors, nat) : Matrix{Float64}(undef, 0, 0)
    # AllImages edge cutoff with the same relative tolerance as the tie band (so a
    # degenerate shell at the cutoff is admitted whole, not split by round-off).
    cut2 = neighbors.cutoff^2 * fac

    for N = 2:nbody
        members = ClusterMember[]
        for a = 1:nat
            sites = neigh[a]
            ns = length(sites)
            ns >= N - 1 || continue
            anchor = (a, z)
            for combo in _ordered_subsets(ns, N - 1)
                full = Tuple{Int,SVector{3,Int}}[anchor]
                for k in combo
                    push!(full, sites[k])
                end
                # MinimumImage clusters must use distinct atoms — two images of the
                # same atom share its spin in plain PBC, so reusing one aliases a
                # lower-body term. AllImages (spin-spiral) only forbids an exact
                # (atom, shift) repeat (distinct images of an atom are independent).
                (use_minimage ? _distinct_atoms(full) : _all_distinct(full)) || continue
                # every edge admissible (anchor↔chosen holds by construction for the
                # matching selection; this re-checks all pairs, incl. chosen↔chosen)
                ok = true
                for x = 1:N
                    ax = full[x][1]
                    px = sitepos(full[x]...)
                    for y = (x + 1):N
                        ay = full[y][1]
                        d2 = sum(abs2, sitepos(full[y]...) - px)
                        admit = use_minimage ?
                            (isfinite(dmin2[ax, ay]) && d2 <= dmin2[ax, ay] * fac) :
                            (d2 <= cut2)
                        if !admit
                            ok = false
                            break
                        end
                    end
                    ok || break
                end
                ok || continue
                push!(members, ClusterMember([s[1] for s in full],
                                             SVector{3,Int}[s[2] for s in full]))
            end
        end
        out[N] = members
    end
    return out
end

@inline function _all_distinct(sites::Vector{Tuple{Int,SVector{3,Int}}})::Bool
    n = length(sites)
    @inbounds for x = 1:n, y = (x + 1):n
        sites[x] == sites[y] && return false
    end
    return true
end

# Distinct atom *indices* — the MinimumImage requirement (see the call site).
@inline function _distinct_atoms(sites::Vector{Tuple{Int,SVector{3,Int}}})::Bool
    n = length(sites)
    @inbounds for x = 1:n, y = (x + 1):n
        sites[x][1] == sites[y][1] && return false
    end
    return true
end
