# Mean-field sampler — P1: single global magnetization, isotropic vMF draw
# (see `docs/specs/mfa-sampling.md`).
#
# Builds on the P0 single-site engine: each spin `a` is drawn from a von Mises–Fisher
# distribution `vMF(ê_a, κ)` about its reference direction, with one global
# concentration `κ = 3m/τ` fixed by the classical-Heisenberg mean-field
# self-consistency `m = L(3m/τ)` (`L` = Langevin function). This is the
# Magesty-equivalent isotropic sampler: it has no couplings, works purely in the
# reduced temperature `τ = T/T_MF`, and validates the engine end to end.
#
# `AbstractSampler` is the dispatch seam the later, model-backed samplers (per-
# sublattice / tensorial / higher-order MFA, P2+) slot into; `sample` is the one verb.

"""
    AbstractSampler

Dispatch seam for spin-configuration samplers. [`MFASampler`](@ref) (mean-field) is the
first; future Metropolis-MC / spin-spiral samplers can slot in behind the `sample` verb.
"""
abstract type AbstractSampler end

# Numerical guards on the reduced temperature τ = T/T_MF, mirroring the reference
# sampler so the ordered / near-uniform boundaries reproduce it exactly.
const _MFA_MIN_TAU = 1.0e-5        # below this: fully ordered, return the reference
const _MFA_MAX_TAU = 0.99999       # above this: use a vanishing concentration
const _MFA_KAPPA_UNIFORM = 1.0e-6  # the near-uniform-limit concentration
# Bracket for the self-consistent magnetization root m ∈ (0, 1). A distinct physical
# quantity from the temperature guards above, so it carries its own name. The lower end
# stays at 1e-5 (not smaller): L(κ) = coth κ − 1/κ cancels catastrophically as κ → 0, so
# evaluating it at κ = 3m/τ for m ≪ 1e-5 returns a garbage sign and breaks the bracket.
# The upper end sits at 1⁻ʳ so the near-saturated root for τ just above _MFA_MIN_TAU
# stays bracketed (large κ there, no cancellation) rather than clamping short.
const _MFA_M_MIN = 1.0e-5
const _MFA_M_MAX = 1.0 - 1.0e-9

# The Langevin function L(κ) = coth κ − 1/κ = ⟨cosθ⟩ for a vMF field of concentration κ.
_langevin(κ::Real)::Float64 = coth(κ) - 1 / κ

# Bisection root of a function that is negative at `lo` and positive at `hi` (both
# assumed bracketing). Used for the monotone mean-field self-consistency; avoids a
# Roots.jl dependency.
function _bisect(f, lo::Float64, hi::Float64; tol::Float64 = 1.0e-12, maxit::Int = 200)::Float64
    a, b = lo, hi
    fa = f(a)
    for _ = 1:maxit
        m = 0.5 * (a + b)
        fm = f(m)
        (abs(fm) <= tol || (b - a) <= tol) && return m
        if (fa < 0) == (fm < 0)
            a, fa = m, fm
        else
            b = m
        end
    end
    return 0.5 * (a + b)
end

"""
    thermal_averaged_m(τ) -> Float64

Solve the classical-Heisenberg mean-field self-consistency `m = L(3m/τ)` (equivalently
`m = coth(3m/τ) − τ/3m`) for the thermally averaged magnetization `m` at reduced
temperature `τ = T/T_MF`. Returns `1.0` for `τ < $(_MFA_MIN_TAU)` (fully ordered) and
`0.0` for `τ > $(_MFA_MAX_TAU)` (fully disordered).
"""
function thermal_averaged_m(τ::Real)::Float64
    τf = Float64(τ)
    τf < _MFA_MIN_TAU && return 1.0
    τf > _MFA_MAX_TAU && return 0.0
    # f(m) = m − L(3m/τ) is increasing on (0,1): f(m_min) = m(1−1/τ) < 0 (for τ<1) and
    # f(m_max) = 1⁻ − L(3/τ) > 0 across the whole interior τ range, so the root is
    # bracketed in [m_min, m_max].
    f(m) = m - _langevin(3m / τf)
    return _bisect(f, _MFA_M_MIN, _MFA_M_MAX)
end

"""
    tau_from_magnetization(m) -> Float64

Invert the mean-field self-consistency: the reduced temperature `τ = T/T_MF` whose
thermally averaged magnetization is `m`. Returns `1.0` for `m ≤ 0` and `0.0` for
`m ≥ 1`; magnetizations whose temperature falls outside `[$(_MFA_MIN_TAU),
$(_MFA_MAX_TAU)]` clamp to the disordered (`1.0`) or ordered (`0.0`) limit.
"""
function tau_from_magnetization(m::Real)::Float64
    mf = Float64(m)
    mf <= 0.0 && return 1.0
    mf >= 1.0 && return 0.0
    g(τ) = mf - _langevin(3mf / τ)               # increasing in τ
    g_lo = g(_MFA_MIN_TAU)
    g_hi = g(_MFA_MAX_TAU)
    if (g_lo < 0) == (g_hi < 0)                   # root outside the bracket: clamp
        return abs(g_lo) <= abs(g_hi) ? 0.0 : 1.0
    end
    return _bisect(g, _MFA_MIN_TAU, _MFA_MAX_TAU)
end

# Resolve a reduced-temperature value into the draw state (ordered flag, concentration
# κ, realized magnetization m). `κ` is unused when `ordered`.
function _mfa_state(τ::Float64)
    if τ < _MFA_MIN_TAU
        return (true, Inf, 1.0)
    elseif τ > _MFA_MAX_TAU
        return (false, _MFA_KAPPA_UNIFORM, 0.0)
    end
    m = thermal_averaged_m(τ)
    return (false, 3m / τ, m)
end

# Uniform random rotation in SO(3) via Shoemake's unit-quaternion construction. Used by
# the optional `randomize` global-frame randomization (a no-op on isotropic statistics,
# but it diversifies the realized frame for later anisotropic/SOC training).
function _random_rotation(rng::AbstractRNG)::SMatrix{3,3,Float64}
    u1, u2, u3 = rand(rng), rand(rng), rand(rng)
    x = sqrt(1 - u1) * sin(2π * u2)
    y = sqrt(1 - u1) * cos(2π * u2)
    z = sqrt(u1) * sin(2π * u3)
    w = sqrt(u1) * cos(2π * u3)
    return @SMatrix [
        1-2(y*y+z*z)  2(x*y-z*w)    2(x*z+y*w)
        2(x*y+z*w)    1-2(x*x+z*z)  2(y*z-x*w)
        2(x*z-y*w)    2(y*z+x*w)    1-2(x*x+y*y)
    ]
end

"""
    MFASampler(reference) <: AbstractSampler

Single global, isotropic mean-field sampler. `reference` is a `3 × n_atoms` matrix of
seed spin directions (the ordered reference state); its columns are normalized to unit
vectors on construction (per-atom magnitudes live separately in the dataset, so they are
discarded here). Draws follow [`sample`](@ref): every spin is sampled from `vMF(ê_a, κ)`
with one global concentration set by the reduced temperature `τ` (or magnetization `m`).
"""
struct MFASampler <: AbstractSampler
    reference::Matrix{Float64}     # 3 × n_atoms, unit columns
    function MFASampler(reference::AbstractMatrix{<:Real})
        size(reference, 1) == 3 ||
            throw(ArgumentError("reference must be 3 × n_atoms; got $(size(reference))"))
        n = size(reference, 2)
        n >= 1 || throw(ArgumentError("reference must have ≥ 1 atom"))
        ref = Matrix{Float64}(undef, 3, n)
        for a = 1:n
            v = SVector{3,Float64}(reference[1, a], reference[2, a], reference[3, a])
            nv = norm(v)
            nv > 1.0e-10 || throw(ArgumentError(
                "reference column $a has ~zero norm; cannot define a direction"))
            ref[:, a] = v / nv
        end
        return new(ref)
    end
end

Base.show(io::IO, s::MFASampler) =
    print(io, "MFASampler(", size(s.reference, 2), " atoms)")

"""
    mfa_temperature_scale(sampler) -> Float64

The mean-field ordering temperature `T_MF` in the sampler's energy units, so a caller
can convert `T = τ·T_MF` (decision D4). The single global [`MFASampler`](@ref) carries no
couplings, so it works in reduced units and returns `T_MF = 1.0` (then `τ` is itself the
temperature); model-backed samplers (P2+) return the linearized mean-field `T_MF`.
"""
mfa_temperature_scale(::MFASampler)::Float64 = 1.0

"""
    MFASample

Labeled result of [`sample`](@ref) (decision D1). `configs` is the bare
`Vector{Matrix{Float64}}` (each `3 × n_atoms` unit directions); the parallel `tau` and
`m` hold each config's reduced temperature and realized magnetization. The object is
iterable and indexable as its `configs`.
"""
struct MFASample
    configs::Vector{Matrix{Float64}}
    tau::Vector{Float64}
    m::Vector{Float64}
end

Base.length(s::MFASample) = length(s.configs)
Base.getindex(s::MFASample, i) = s.configs[i]
Base.eltype(::Type{MFASample}) = Matrix{Float64}
Base.iterate(s::MFASample, st::Int = 1) =
    st > length(s.configs) ? nothing : (s.configs[st], st + 1)
Base.show(io::IO, s::MFASample) =
    print(io, "MFASample(", length(s.configs), " configs)")

# Validate 1-based atom indices against the atom count.
function _check_atom_indices(idx::AbstractVector{<:Integer}, n::Int, name::AbstractString)
    for i in idx
        (1 <= i <= n) ||
            throw(ArgumentError("$name entry $i is outside 1:$n"))
    end
    return nothing
end

# Draw one configuration at concentration κ (ordered ⇒ copy the reference). `uniform`
# columns are redrawn isotropically (overriding the vMF draw); a global rotation and the
# `fixed` columns are applied last so `fixed` takes precedence over everything.
function _draw_config(rng::AbstractRNG, ref::Matrix{Float64}, ordered::Bool, κ::Float64,
                      randomize::Bool, fixed::AbstractVector{<:Integer},
                      uniform::AbstractVector{<:Integer})::Matrix{Float64}
    n = size(ref, 2)
    out = Matrix{Float64}(undef, 3, n)
    @inbounds for a = 1:n
        ea = SVector{3,Float64}(ref[1, a], ref[2, a], ref[3, a])
        v = ordered ? ea : sample_vmf(rng, ea, κ)
        out[1, a], out[2, a], out[3, a] = v[1], v[2], v[3]
    end
    for a in uniform
        u = _random_unit(rng)
        out[1, a], out[2, a], out[3, a] = u[1], u[2], u[3]
    end
    if randomize
        R = _random_rotation(rng)
        out = Matrix(R * out)
        for a in fixed
            out[:, a] = R * SVector{3,Float64}(ref[1, a], ref[2, a], ref[3, a])
        end
    else
        for a in fixed
            out[:, a] = SVector{3,Float64}(ref[1, a], ref[2, a], ref[3, a])
        end
    end
    return out
end

# Resolve exactly one of `tau` / `m` (scalars or collections) into a reduced-temperature
# iterator. Returns a vector of τ values.
function _resolve_taus(tau, m)::Vector{Float64}
    (tau === nothing) == (m === nothing) &&
        throw(ArgumentError("provide exactly one of `tau` or `m`"))
    if tau !== nothing
        return tau isa Real ? [Float64(tau)] : Float64[Float64(t) for t in tau]
    end
    return m isa Real ? [tau_from_magnetization(m)] :
           Float64[tau_from_magnetization(mi) for mi in m]
end

"""
    sample(sampler::MFASampler, n; tau, m, rng, randomize, fixed, uniform) -> MFASample
    sample(sampler::MFASampler; tau, m, nsamples = 1, rng, randomize, fixed, uniform) -> MFASample

Draw spin configurations from the mean-field sampler. Provide **exactly one** control
variable: the reduced temperature `tau = T/T_MF` or the magnetization `m`. Each spin is
drawn from `vMF(ê_a, κ)` about its reference direction with the global concentration
`κ = 3m/τ` set by the self-consistency [`thermal_averaged_m`](@ref).

The first form draws `n` configurations at a single control value. The second form sweeps
a **collection** `tau` (or `m`) and draws `nsamples` configurations per value, ordered
value-outer / sample-inner.

# Keyword arguments
- `rng::AbstractRNG = default_rng()`: explicit, seeded RNG (reproducible draws).
- `randomize::Bool = false`: apply one uniform random global rotation per configuration.
- `fixed::AbstractVector{<:Integer} = Int[]`: atoms held at the reference direction
  (rotated with the frame when `randomize`); takes precedence over `uniform`.
- `uniform::AbstractVector{<:Integer} = Int[]`: atoms redrawn isotropically (the
  disordered limit) regardless of the control value.

Returns an [`MFASample`](@ref): `.configs` (each `3 × n_atoms` unit directions) with
parallel labels `.tau` and `.m`.
"""
function sample(sampler::MFASampler, n::Integer; tau = nothing, m = nothing,
                rng::AbstractRNG = default_rng(), randomize::Bool = false,
                fixed::AbstractVector{<:Integer} = Int[],
                uniform::AbstractVector{<:Integer} = Int[])::MFASample
    n >= 0 || throw(ArgumentError("n must be ≥ 0; got $n"))
    taus = _resolve_taus(tau, m)
    length(taus) == 1 ||
        throw(ArgumentError("the positional `n` form takes a scalar `tau`/`m`; " *
                            "pass a collection without `n` to sweep"))
    return _sample_sweep(sampler, taus, n, rng, randomize, fixed, uniform)
end

function sample(sampler::MFASampler; tau = nothing, m = nothing, nsamples::Integer = 1,
                rng::AbstractRNG = default_rng(), randomize::Bool = false,
                fixed::AbstractVector{<:Integer} = Int[],
                uniform::AbstractVector{<:Integer} = Int[])::MFASample
    nsamples >= 0 || throw(ArgumentError("nsamples must be ≥ 0; got $nsamples"))
    taus = _resolve_taus(tau, m)
    return _sample_sweep(sampler, taus, nsamples, rng, randomize, fixed, uniform)
end

# Core sweep: for each τ draw `per` configurations, value-outer / sample-inner.
function _sample_sweep(sampler::MFASampler, taus::Vector{Float64}, per::Integer,
                       rng::AbstractRNG, randomize::Bool,
                       fixed::AbstractVector{<:Integer},
                       uniform::AbstractVector{<:Integer})::MFASample
    ref = sampler.reference
    natoms = size(ref, 2)
    _check_atom_indices(fixed, natoms, "fixed")
    _check_atom_indices(uniform, natoms, "uniform")
    total = length(taus) * per
    configs = Vector{Matrix{Float64}}(undef, total)
    tau_lab = Vector{Float64}(undef, total)
    m_lab = Vector{Float64}(undef, total)
    k = 0
    for τ in taus
        ordered, κ, mval = _mfa_state(τ)
        for _ = 1:per
            k += 1
            configs[k] = _draw_config(rng, ref, ordered, κ, randomize, fixed, uniform)
            tau_lab[k] = τ
            m_lab[k] = mval
        end
    end
    return MFASample(configs, tau_lab, m_lab)
end
