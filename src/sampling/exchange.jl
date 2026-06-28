# Mean-field sampler — P2: multi-sublattice isotropic exchange (see
# `docs/specs/mfa-sampling.md`).
#
# `ExchangeModel` is the neutral carrier of the bilinear couplings the mean-field sampler
# needs. P2 uses the **isotropic** (Heisenberg) part only: the total isotropic exchange
# `Jiso[a,b] = Σ_R J_iso(a,b,R)` between atom `a` and every periodic image of atom `b`
# (a symmetric `n × n` matrix). Tensorial (DMI / anisotropic) and single-ion channels are
# a P3 extension.
#
# Given a reference state `{ê_a}`, the mean-field decoupling of the isotropic energy
# `E = Σ_{bonds} J e_a·e_b` gives a single-site field on atom `a`,
#   h_a = −Σ_b Jiso[a,b] ⟨e_b⟩ = −Σ_b Jiso[a,b] m_b ê_b,
# whose projection on the (stationary) reference is `λ_a = ê_a·h_a = Σ_b A[a,b] m_b` with
#   A[a,b] = −Jiso[a,b] (ê_a·ê_b).
# Folding the reference directions into `A` turns any collinear order (ferro / antiferro /
# ferri) ferromagnetic in the magnitude variables `m_a`, so its leading (Perron)
# eigenvalue `ρ` with a positive eigenvector is the ordering mode and `T_MF = ρ/3`. In the
# reduced temperature `τ = T/T_MF` the self-consistency is
#   m_a = L(3 (Ā m)_a / τ),   Ā = A/ρ   (spectral radius 1),
# and the per-atom vMF concentration is `κ_a = 3 (Ā m)_a / τ`. Scaling all couplings
# scales `A` and `ρ` together, so `Ā` — and hence every `m_a(τ)` — is scale-free (D4):
# only the coupling *ratios* matter, not the absolute `T_MF`.

"""
    ExchangeModel(Jiso)
    ExchangeModel(model::SCEModel)

Neutral carrier of the isotropic bilinear exchange the mean-field sampler needs (P2).
`Jiso` is a symmetric `n × n` matrix whose entry `[a,b]` is the total isotropic exchange
`Σ_R J_iso(a,b,R)` between atom `a` and all periodic images of atom `b` (the sign
convention is that of the SCE energy `E = Σ_{bonds} J e_a·e_b`, with the diagonal the
sum over an atom's own image bonds).

`ExchangeModel(model)` extracts `Jiso` from a fitted [`SCEModel`](@ref): it reuses the
Sunny bilinear extraction (the per-bond `3×3` matrix `M`) and keeps the isotropic part
`tr(M)/3` of each bond. The DMI / anisotropic parts of `M`, single-ion (`ls=[2]`)
anisotropy, and higher-order channels are **dropped** (a P3/P4 extension) and reported via
`@warn`. Pass the matrix directly to supply an external `Jij` (a P5 reader target).
"""
struct ExchangeModel
    natoms::Int
    Jiso::Matrix{Float64}
    function ExchangeModel(Jiso::AbstractMatrix{<:Real})
        n = size(Jiso, 1)
        size(Jiso, 2) == n || throw(ArgumentError("Jiso must be square; got $(size(Jiso))"))
        J = Matrix{Float64}(Jiso)
        maximum(abs.(J - J'); init = 0.0) <= 1.0e-10 * (1 + maximum(abs.(J); init = 0.0)) ||
            throw(ArgumentError("Jiso must be symmetric (Jiso[a,b] = Jiso[b,a])"))
        J = 0.5 * (J + J')                           # symmetrize exactly
        return new(n, J)
    end
end

Base.show(io::IO, m::ExchangeModel) = print(io, "ExchangeModel(", m.natoms, " atoms)")

function ExchangeModel(model::SCEModel)
    terms = _sunny_supercell_terms(model)
    n = num_atoms(model.basis.crystal)
    J = zeros(Float64, n, n)
    aniso = 0.0
    scale = 0.0
    for ((a, b, _), M) in terms.pairs
        j = (M[1, 1] + M[2, 2] + M[3, 3]) / 3        # isotropic (Heisenberg) part tr(M)/3
        J[a, b] += j
        a != b && (J[b, a] += j)
        # track the dropped non-isotropic (DMI + anisotropic) magnitude for the warning
        aniso = max(aniso, maximum(abs.(M - j * I)))   # M is a 3×3 SMatrix (never empty)
        scale = max(scale, maximum(abs.(M)))
    end
    dropped = String[]
    aniso > 1.0e-8 * (1 + scale) && push!(dropped,
        "non-isotropic bilinear (DMI / anisotropic exchange) parts of $(length(terms.pairs)) " *
        "bond(s) — only the Heisenberg part tr(M)/3 is kept (P2)")
    isempty(terms.onsites) ||
        push!(dropped, "$(length(terms.onsites)) single-ion (ls=[2]) term(s)")
    isempty(terms.skipped) ||
        push!(dropped, "$(length(terms.skipped)) higher-order / higher-l SALC(s)")
    isempty(dropped) || @warn "ExchangeModel keeps the isotropic bilinear exchange only; " *
        "dropped: " * join(dropped, "; ") * " (tensorial / single-ion / many-body are P3+)."
    return ExchangeModel(J)
end

# The symmetric molecular-field matrix A[a,b] = −Jiso[a,b] (ê_a·ê_b) for a reference `ref`.
function _mfa_matrix(exch::ExchangeModel, ref::Matrix{Float64})::Matrix{Float64}
    n = exch.natoms
    A = zeros(Float64, n, n)
    @inbounds for b = 1:n
        eb = SVector{3,Float64}(ref[1, b], ref[2, b], ref[3, b])
        for a = 1:n
            ea = SVector{3,Float64}(ref[1, a], ref[2, a], ref[3, a])
            A[a, b] = -exch.Jiso[a, b] * dot(ea, eb)
        end
    end
    return A
end

# Largest eigenvalue ρ of the symmetric A and its eigenvector; the second return flags
# whether the leading eigenvector is sign-definite (the reference is a clean ordering mode).
function _perron(A::Matrix{Float64})
    F = eigen(Symmetric(A))
    ρ = F.values[end]
    v = F.vectors[:, end]
    tol = 1.0e-8 * maximum(abs.(v); init = 1.0)
    pos = count(>(tol), v)
    neg = count(<(-tol), v)
    return ρ, (pos == 0 || neg == 0)             # all entries one sign (ignoring ~0)
end

# Verify the reference is a stationary state of the isotropic model: the molecular field
# h_a = −Σ_b Jiso[a,b] ê_b must be parallel to ê_a (transverse part ≈ 0) and aligned
# (λ_a = ê_a·h_a > 0). Exact for any collinear reference; warns otherwise (D2).
function _check_reference_stationary(exch::ExchangeModel, ref::Matrix{Float64})
    n = exch.natoms
    worst_t = 0.0
    scale = 0.0
    antialigned = false
    @inbounds for a = 1:n
        ea = SVector{3,Float64}(ref[1, a], ref[2, a], ref[3, a])
        h = zero(SVector{3,Float64})
        for b = 1:n
            eb = SVector{3,Float64}(ref[1, b], ref[2, b], ref[3, b])
            h -= exch.Jiso[a, b] * eb
        end
        λ = dot(ea, h)
        trans = norm(h - λ * ea)
        worst_t = max(worst_t, trans)
        scale = max(scale, norm(h))
        λ < -1.0e-8 * (1 + norm(h)) && (antialigned = true)
    end
    worst_t > 1.0e-6 * (1 + scale) && @warn "the reference is not a stationary state of " *
        "the exchange model (molecular field not ∥ ê_a; max transverse $(worst_t)); the " *
        "rigid-direction MFA is approximate here — a rotating cone axis ê_a(τ) is a P3 extension."
    antialigned && @warn "the molecular field is antiparallel to the reference on some " *
        "atom; the reference is not a stable ordered state of the exchange model."
    return nothing
end

function MFASampler(exch::ExchangeModel; reference::AbstractMatrix{<:Real})
    ref = _normalize_reference(reference)
    size(ref, 2) == exch.natoms || throw(DimensionMismatch(
        "reference has $(size(ref, 2)) atoms but the ExchangeModel has $(exch.natoms)"))
    A = _mfa_matrix(exch, ref)
    ρ, signdef = _perron(A)
    ρ > 1.0e-12 * (1 + maximum(abs.(A); init = 0.0)) || throw(ArgumentError(
        "the exchange model has no ordering instability about this reference " *
        "(spectral radius ≈ 0); use the single global `MFASampler(reference)` instead"))
    signdef || @warn "the leading molecular-field eigenvector is not sign-definite; the " *
        "reference may be frustrated or not the model's ordered ground state — the " *
        "self-consistency still runs but T_MF/m_a(τ) may not reflect the intended order."
    _check_reference_stationary(exch, ref)
    return MFASampler(ref, A ./ ρ, ρ / 3)
end

# Is atom `a` coupled (a non-zero row in the normalized molecular-field matrix)? A free
# spin (all-zero row) feels no molecular field and disorders for every τ > 0.
function _is_coupled(Abar::Matrix{Float64}, a::Int)::Bool
    @inbounds for b = 1:size(Abar, 2)
        abs(Abar[a, b]) > 1.0e-12 && return true
    end
    return false
end

# The coupled per-atom self-consistency m_a = L(3 (Ā m)_a / τ) by damped fixed-point
# iteration from the fully ordered state m = 1 (selects the stable ordered branch over the
# trivial m = 0). Returns (ordered, κ::Vector, m::Vector).
function _coupled_state(Abar::Matrix{Float64}, τ::Float64, n::Int)
    if τ < _MFA_MIN_TAU
        # Fully ordered limit, per atom: coupled spins saturate (m=1, κ=Inf ⇒ the vMF draw
        # returns the reference), but a free spin feels no field and stays disordered
        # (m=0, κ=0 ⇒ isotropic). `ordered=false` so `_draw_config` honors the per-atom κ.
        m = Vector{Float64}(undef, n)
        κ = Vector{Float64}(undef, n)
        @inbounds for a = 1:n
            c = _is_coupled(Abar, a)
            m[a] = c ? 1.0 : 0.0
            κ[a] = c ? Inf : 0.0
        end
        return (false, κ, m)
    elseif τ > _MFA_MAX_TAU
        return (false, fill(_MFA_KAPPA_UNIFORM, n), zeros(n))
    end
    m, Δ = _coupled_solve(Abar, τ, n)
    # Tolerance 1e-6: deep in the critical regime (τ within ~0.01% of T_MF, where m ~ 1% and
    # the state is essentially disordered) depth-1 Anderson stagnates around a few ×1e-7,
    # accurate enough for the draw; a genuinely non-convergent (frustrated) reference leaves
    # a residual O(10⁻²) and still warns.
    Δ > 1.0e-6 && @warn "the coupled mean-field self-consistency did not converge at " *
        "τ = $τ (residual $Δ); m_a may be inaccurate (critical slowing near T_MF, or a " *
        "frustrated reference)."
    κ = Vector{Float64}(undef, n)
    f = Abar * m
    @inbounds for a = 1:n
        κ[a] = max(0.0, 3 * f[a] / τ)                 # see the clamp note in _scf_map!
    end
    return (false, κ, m)
end

# One self-consistency map step g(m)_a = L(3 (Ā m)_a / τ), clamped to [0,1]. The clamp is
# load-bearing only off the intended (sign-definite Ā) path: a frustrated reference (sign-
# indefinite Ā, @warn-flagged at construction) can drive a field negative, which is floored.
function _scf_map!(out::Vector{Float64}, Abar::Matrix{Float64}, m::Vector{Float64}, τ::Float64)
    mul!(out, Abar, m)
    @inbounds for a in eachindex(out)
        out[a] = clamp(_langevin(3 * out[a] / τ), 0.0, 1.0)
    end
    return out
end

# Solve m = g(m) by depth-1 Anderson acceleration from the ordered state m = 1 (selects the
# stable ordered branch over the trivial m = 0). Anderson is used instead of damped Picard
# because the contraction rate → 1 as τ → 1⁻ (critical slowing) makes plain iteration far
# too slow near T_MF. Returns (m, residual ‖g(m) − m‖∞).
function _coupled_solve(Abar::Matrix{Float64}, τ::Float64, n::Int)
    m = ones(n)
    g = similar(m)
    f = similar(m)              # residual g − m
    gp = zeros(n)               # previous g
    fp = zeros(n)               # previous residual
    df = similar(m)
    _scf_map!(g, Abar, m, τ)
    @. f = g - m
    Δ = maximum(abs, f)
    for k = 1:2000
        Δ <= 1.0e-13 && break
        if k == 1
            @. m = g                                  # plain Picard first step
        else
            @. df = f - fp
            d2 = dot(df, df)
            β = d2 > 1.0e-30 ? dot(f, df) / d2 : 0.0  # type-II Anderson mixing coefficient
            @. m = clamp(g - β * (g - gp), 0.0, 1.0)
        end
        @. gp = g
        @. fp = f
        _scf_map!(g, Abar, m, τ)
        @. f = g - m
        Δ = maximum(abs, f)
    end
    return m, Δ
end
