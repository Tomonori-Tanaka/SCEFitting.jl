"""
    Harmonics

Real (tesseral) spherical harmonics `Zₗₘ` and their analytic, tangent-projected
Cartesian gradients, following Drautz, *Phys. Rev. B* **102**, 024104 (2020).

Per-site normalization carries `(4π)^(−1/2)`. The polar part uses the Drautz
normalized associated Legendre `P̄ₗₘ = (−1)^m √((2l+1)/(4π)·(l−m)!/(l+m)!)
dᵐPₗ/dzᵐ` (the `(1−z²)^{m/2}` factor is absorbed into the azimuthal `(x+iy)^m`),
so on the unit sphere these reduce to the standard real solid harmonics
(`Z₁₁ ∝ x`, `Z₂₂ ∝ x²−y²`, …).

The gradient is the genuine on-sphere (tangent) gradient
`∇Z = ∂Z − r̂ (r̂·∂Z)`; it is what the torque design matrix needs.
"""
module Harmonics

using LegendrePolynomials: dnPl
using LinearAlgebra: norm
using StaticArrays

# Not exported: callers reach these via `Harmonics.Zlm` etc. The `_unsafe`
# variants skip validation for hot paths.

"""
    lm_index(l, m) -> Int

1-based linear index of `(l, m)` in the contiguous `Zₗₘ` ordering
`(0,0), (1,−1), (1,0), (1,1), (2,−2), …`. Inverse-friendly: `l² + l + m + 1`.
"""
@inline lm_index(l::Integer, m::Integer)::Int = l * l + l + m + 1

"""
    num_lm(lmax) -> Int

Number of `(l, m)` pairs with `l ≤ lmax`, i.e. `(lmax + 1)²`.
"""
@inline num_lm(lmax::Integer)::Int = (lmax + 1)^2

@inline _parity(n::Integer)::Int = isodd(n) ? -1 : 1

# √((2l+1)/(4π) · (l−m)!/(l+m)!) for m ≥ 0, formed without large factorials.
@inline function _plm_norm(l::Int, m::Int)::Float64
    acc = (2 * l + 1) / (4π)
    @inbounds for i = (l - m + 1):(l + m)
        acc /= i
    end
    return sqrt(acc)
end

# Drautz normalized associated Legendre P̄ₗₘ(z) and its z-derivative (m via |m|).
@inline _barP(l::Integer, m::Integer, z::Real)::Float64 =
    (n = abs(m); _parity(n) * _plm_norm(Int(l), n) * dnPl(z, l, n))
@inline _dbarP(l::Integer, m::Integer, z::Real)::Float64 =
    (n = abs(m); _parity(n) * _plm_norm(Int(l), n) * dnPl(z, l, n + 1))

@inline function _validate_lm(l::Integer, m::Integer)
    l >= 0 || throw(ArgumentError("need l ≥ 0 (got l=$l)"))
    -l <= m <= l || throw(ArgumentError("need -l ≤ m ≤ l (got l=$l, m=$m)"))
    return nothing
end

@inline function _validate_unit(u::AbstractVector{<:Real}; atol::Real = 1e-8)
    length(u) == 3 || throw(ArgumentError("direction must have length 3 (got $(length(u)))"))
    abs(norm(u) - 1) <= atol ||
        throw(ArgumentError("direction must be a unit vector (‖u‖ = $(norm(u)))"))
    return nothing
end

"""
    Zlm_unsafe(l, m, u) -> Float64

Tesseral harmonic `Zₗₘ(u)` for a unit vector `u`, without input validation.
"""
@inline function Zlm_unsafe(l::Integer, m::Integer, u::AbstractVector{<:Real})::Float64
    n = abs(m)
    plm = _barP(l, n, u[3])
    n == 0 && return plm
    c = _parity(n) * sqrt(2.0) * plm
    zpow = ComplexF64(u[1], u[2])^n
    return m > 0 ? c * real(zpow) : c * imag(zpow)
end

"""
    Zlm(l, m, u) -> Float64

Tesseral harmonic `Zₗₘ` at the unit direction `u` (3-vector).

# Arguments
- `l::Integer`: angular momentum (`l ≥ 0`).
- `m::Integer`: order (`-l ≤ m ≤ l`).
- `u::AbstractVector{<:Real}`: unit direction `[x, y, z]`.

# Returns
- `Float64`: `Zₗₘ(u)`.
"""
function Zlm(l::Integer, m::Integer, u::AbstractVector{<:Real})::Float64
    _validate_lm(l, m)
    _validate_unit(u)
    return Zlm_unsafe(l, m, u)
end

"""
    grad_Zlm_unsafe(l, m, u) -> SVector{3,Float64}

Tangent-projected Cartesian gradient `∇Zₗₘ(u) = ∂Z − u (u·∂Z)`, without
validation.
"""
@inline function grad_Zlm_unsafe(l::Integer, m::Integer,
                                 u::AbstractVector{<:Real})::SVector{3,Float64}
    x, y, z = u[1], u[2], u[3]
    n = abs(m)
    plm = _barP(l, n, z)
    dplm = _dbarP(l, n, z)
    if m == 0
        zz = z * dplm
        return SVector{3,Float64}(-x * zz, -y * zz, dplm - z * zz)
    end
    c = _parity(n) * sqrt(2.0)
    zxy = ComplexF64(x, y)
    zpn = zxy^n
    zpn1 = zxy^(n - 1)
    rn, iN = real(zpn), imag(zpn)
    rn1, in1 = real(zpn1), imag(zpn1)
    dZx, dZy, dZz = if m > 0
        (c * n * plm * rn1, -c * n * plm * in1, c * dplm * rn)
    else
        (c * n * plm * in1, c * n * plm * rn1, c * dplm * iN)
    end
    zz = x * dZx + y * dZy + z * dZz   # u · ∂Z, the radial part to remove
    return SVector{3,Float64}(dZx - x * zz, dZy - y * zz, dZz - z * zz)
end

"""
    grad_Zlm(l, m, u) -> SVector{3,Float64}

On-sphere (tangent-projected) Cartesian gradient `∇Zₗₘ(u)`.

# Arguments
- `l::Integer`: angular momentum (`l ≥ 0`).
- `m::Integer`: order (`-l ≤ m ≤ l`).
- `u::AbstractVector{<:Real}`: unit direction `[x, y, z]`.

# Returns
- `SVector{3,Float64}`: `(∂ₓ, ∂_y, ∂_z) Zₗₘ` with the radial component removed, so
  `u·∇Zₗₘ = 0`.
"""
function grad_Zlm(l::Integer, m::Integer, u::AbstractVector{<:Real})::SVector{3,Float64}
    _validate_lm(l, m)
    _validate_unit(u)
    return grad_Zlm_unsafe(l, m, u)
end

end # module Harmonics
