# SALC construction: project each cluster orbit's representative tensor onto the
# trivial irrep of its site stabilizer (orbit–stabilizer theorem), gauge-fix a
# deterministic orthonormal basis, then transport to all orbit members. v0 covers
# 1- and 2-body clusters with per-site l ≤ 2, lₛ ≥ 1, Σl even (so every surviving
# channel has l₁ = l₂ and the l₁≠l₂ site-swap subtlety never arises).

# Per-site l-assignments to enumerate for an orbit of the given species.
function _enumerate_ls(N::Int, species::Vector{Int}, lmax::Vector{Int})
    if N == 1
        s = species[1]
        return [[l] for l = 1:lmax[s] if iseven(l)]
    elseif N == 2
        s1, s2 = species[1], species[2]
        out = Vector{Int}[]
        for l1 = 1:lmax[s1], l2 = 1:lmax[s2]
            l1 <= l2 || continue            # unordered (sorted) multiset
            iseven(l1 + l2) || continue     # time-reversal even
            push!(out, [l1, l2])
        end
        return out
    else
        throw(ArgumentError("only 1- and 2-body clusters are implemented"))
    end
end

# Stabilizer of a 1-body representative: ops fixing the atom. (no swap)
function _stabilizer1(crystal::Crystal, sg::SpaceGroup, rep::ClusterMember)
    a = rep.atoms[1]
    return [(g, false) for g = 1:n_ops(sg) if sg.map_sym[a, g] == a]
end

# Stabilizer of a 2-body representative: ops mapping the (anchored) pair to itself,
# possibly swapping the two sites. Returns (op_index, swap::Bool).
function _stabilizer2(crystal::Crystal, sg::SpaceGroup, rep::ClusterMember, ls::Vector{Int})
    s1 = (rep.atoms[1], rep.shifts[1])
    s2 = (rep.atoms[2], rep.shifts[2])
    out = Tuple{Int,Bool}[]
    for g = 1:n_ops(sg)
        img1 = _site_image(crystal, sg, g, rep.atoms[1], rep.shifts[1])
        img2 = _site_image(crystal, sg, g, rep.atoms[2], rep.shifts[2])
        for Δ in (img1[2], img2[2])
            t1 = (img1[1], img1[2] - Δ)
            t2 = (img2[1], img2[2] - Δ)
            if t1 == s1 && t2 == s2
                push!(out, (g, false))
                break
            elseif ls[1] == ls[2] && t1 == s2 && t2 == s1
                push!(out, (g, true))
                break
            end
        end
    end
    return out
end

# Op connecting the representative to an orbit member (rep ↦ member), with its swap
# flag, so the member's transported tensor can be built.
function _connect(crystal::Crystal, sg::SpaceGroup, rep::ClusterMember,
                  member::ClusterMember, ls::Vector{Int})
    N = length(rep.atoms)
    if N == 1
        a = rep.atoms[1]
        for g = 1:n_ops(sg)
            sg.map_sym[a, g] == member.atoms[1] && return (g, false)
        end
    else
        m1 = (member.atoms[1], member.shifts[1])
        m2 = (member.atoms[2], member.shifts[2])
        for g = 1:n_ops(sg)
            img1 = _site_image(crystal, sg, g, rep.atoms[1], rep.shifts[1])
            img2 = _site_image(crystal, sg, g, rep.atoms[2], rep.shifts[2])
            for Δ in (img1[2], img2[2])
                t1 = (img1[1], img1[2] - Δ)
                t2 = (img2[1], img2[2] - Δ)
                if t1 == m1 && t2 == m2
                    return (g, false)
                elseif ls[1] == ls[2] && t1 == m2 && t2 == m1
                    return (g, true)
                end
            end
        end
    end
    error("no symmetry operation connects the representative to a member")
end

# Projector onto the stabilizer-invariant subspace of the (2Lf+1)-dim coefficient
# space, returning an orthonormal basis (columns) of the eigenvalue-1 subspace.
function _project_invariants(sg::SpaceGroup, stab::Vector{Tuple{Int,Bool}}, Lf::Int)
    d = 2Lf + 1
    P = zeros(Float64, d, d)
    for (g, swap) in stab
        ε = swap ? (iseven(Lf) ? 1.0 : -1.0) : 1.0   # CG exchange (−1)^Lf for swaps
        P .+= ε .* AngularMomentum.wignerD_real(Lf, sg.ops[g].rotation_cart)
    end
    P ./= length(stab)
    P .= (P .+ P') ./ 2
    F = eigen(Symmetric(P))
    all(λ -> abs(λ) < 1e-6 || abs(λ - 1) < 1e-6, F.values) ||
        error("SALC projector is not idempotent (eigenvalues $(F.values)); convention bug")
    cols = findall(λ -> abs(λ - 1) < 1e-8, F.values)
    return F.vectors[:, cols]
end

# Deterministic, BLAS-independent gauge: axis-pivoted modified Gram–Schmidt on the
# block projector, with a first-significant-component sign fix.
function _canonical_basis(V1::AbstractMatrix{Float64})
    n, d = size(V1)
    d == 0 && return zeros(Float64, n, 0)
    Q = V1 * V1'
    W = zeros(Float64, n, d)
    k = 0
    for j = 1:n
        u = Q[:, j]
        for i = 1:k
            u .-= dot(view(W, :, i), u) .* view(W, :, i)
        end
        nu = norm(u)
        if nu > 1e-8
            k += 1
            W[:, k] = u ./ nu
            _sign_canon!(view(W, :, k))
        end
        k == d && break
    end
    k == d || error("canonical gauge failed to span the invariant subspace")
    return W
end

function _sign_canon!(v)
    for x in v
        if abs(x) > 1e-8
            x < 0 && (v .*= -1)
            return
        end
    end
end

# Fold the Mf axis of a coupled tensor against a coefficient vector → rank-N tensor.
function _fold_mf(tensor::AbstractArray, c::AbstractVector)
    N = ndims(tensor) - 1
    sz = ntuple(i -> size(tensor, i), N)
    out = zeros(Float64, sz...)
    for idx in CartesianIndices(sz)
        s = 0.0
        for q in eachindex(c)
            s += c[q] * tensor[Tuple(idx)..., q]
        end
        out[idx] = s
    end
    return out
end

# Transport the representative folded tensor to a member: rotate each site axis by
# the connecting op, then reorder axes for a swap.
function _transport(folded_rep::Array{Float64}, ls::Vector{Int}, Rcart, swap::Bool)
    T = folded_rep
    for i in eachindex(ls)
        T = AngularMomentum.nmode_mul(T, AngularMomentum.wignerD_real(ls[i], Rcart), i)
    end
    if swap && length(ls) == 2
        T = permutedims(T, (2, 1))
    end
    return Array{Float64}(T)
end

"""
    build_salc_basis(crystal, spacegroup, clusters; lmax_by_species, isotropy = false) -> SALCBasis

Construct the symmetry-adapted (and time-reversal-even) SCE basis for every
cluster orbit. For each `(orbit, ls, Lf)` the stabilizer-invariant coefficient
subspace is found, gauge-fixed, and transported to all orbit members.

# Keyword arguments
- `lmax_by_species::AbstractVector{<:Integer}`: per-species maximum `l`.
- `isotropy::Bool = false`: keep only the scalar `Lf == 0` channel if `true`.

# Status
The isotropic (`Lf == 0`) scalar channel — Heisenberg and biquadratic invariants —
is validated (space-group invariance + time-reversal evenness, the ground-truth
tests). Anisotropic (`Lf > 0`) channels are constructed but have a known
member-transport site-axis-alignment issue on cells whose pair stabilizers contain
a site swap (the invariance test flags it); use `isotropy = true` for the
validated v0 path until that is resolved.
"""
function build_salc_basis(crystal::Crystal, spacegroup::SpaceGroup, clusters::ClusterSet;
                          lmax_by_species::AbstractVector{<:Integer},
                          isotropy::Bool = false)::SALCBasis
    lmax = collect(Int, lmax_by_species)
    salcs = SALC[]
    for N in sort(collect(keys(clusters.by_body)))
        for (orbit_id, O) in enumerate(clusters.by_body[N])
            rep = O.representative
            for ls in _enumerate_ls(N, O.species, lmax)
                stab = N == 1 ? _stabilizer1(crystal, spacegroup, rep) :
                       _stabilizer2(crystal, spacegroup, rep, ls)
                conns = [_connect(crystal, spacegroup, rep, m, ls) for m in O.members]
                for cb in coupled_bases(ls; isotropy = isotropy)
                    V1 = _project_invariants(spacegroup, stab, cb.Lf)
                    size(V1, 2) == 0 && continue
                    W = _canonical_basis(V1)
                    for block = 1:size(W, 2)
                        folded_rep = _fold_mf(cb.tensor, W[:, block])
                        members = SALCMember[]
                        anynz = false
                        for (m, (g, swap)) in zip(O.members, conns)
                            fm = _transport(folded_rep, ls,
                                            spacegroup.ops[g].rotation_cart, swap)
                            @inbounds for i in eachindex(fm)
                                abs(fm[i]) < 1e-10 && (fm[i] = 0.0)
                            end
                            norm(fm) > 1e-10 && (anynz = true)
                            push!(members, SALCMember(m.atoms, m.shifts, fm))
                        end
                        anynz || continue
                        key = SALCKey(N, orbit_id, collect(ls), cb.Lf, block)
                        push!(salcs, SALC(key, N, collect(ls), cb.Lf, members))
                    end
                end
            end
        end
    end
    sort!(salcs; by = s -> s.key)
    keyvec = SALCKey[s.key for s in salcs]
    return SALCBasis(salcs, keyvec, hash(keyvec))
end
