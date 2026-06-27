"""
    AngularMomentum

Angular-momentum coupling primitives, from scratch:

- `clebsch_gordan(j1, m1, j2, m2, J, M)` — Clebsch–Gordan coefficient via the
  Racah single-sum formula (integer momenta).
- `wignerD_real(l, R)` — the real Wigner-D matrix of a 3×3 rotation `R` in the
  tesseral-harmonic basis. It is built **directly from this package's own**
  `Harmonics.Zₗₘ` by solving the (exact, oversampled) linear system
  `Zₗ,m(R·r̂) = Σ_{m'} Δ[m,m'] Zₗ,m'(r̂)`. This guarantees the rotation matrix is
  consistent with the harmonic convention, uses real arithmetic only, and works
  for improper rotations (reflections, inversion) without special cases.
"""
module AngularMomentum

using LinearAlgebra: norm, I
using StaticArrays

using ..Harmonics: Zlm_unsafe

export clebsch_gordan, wignerD_real

# n! as Float64 (no overflow for the l-range used here); negative → 0.
@inline function _fact(n::Integer)::Float64
    n < 0 && return 0.0
    f = 1.0
    @inbounds for i = 2:n
        f *= i
    end
    return f
end

"""
    clebsch_gordan(j1, m1, j2, m2, J, M) -> Float64

Clebsch–Gordan coefficient `⟨j1 m1 j2 m2 | J M⟩` (integer angular momenta), via the
Racah single-sum formula. Returns `0.0` outside the selection rules.

# Arguments
- `j1, m1, j2, m2, J, M::Integer`: the coupling quantum numbers.

# Returns
- `Float64`: the coefficient.
"""
function clebsch_gordan(j1::Integer, m1::Integer, j2::Integer, m2::Integer,
                        J::Integer, M::Integer)::Float64
    (m1 + m2 == M) || return 0.0
    (j1 >= 0 && j2 >= 0 && J >= 0) || return 0.0
    (abs(m1) <= j1 && abs(m2) <= j2 && abs(M) <= J) || return 0.0
    (abs(j1 - j2) <= J <= j1 + j2) || return 0.0

    tri = _fact(j1 + j2 - J) * _fact(j1 - j2 + J) * _fact(-j1 + j2 + J) /
          _fact(j1 + j2 + J + 1)
    pref = sqrt((2J + 1) * tri) *
           sqrt(_fact(j1 + m1) * _fact(j1 - m1) * _fact(j2 + m2) *
                _fact(j2 - m2) * _fact(J + M) * _fact(J - M))

    kmin = max(0, j2 - m1 - J, j1 + m2 - J)
    kmax = min(j1 + j2 - J, j1 - m1, j2 + m2)
    s = 0.0
    @inbounds for k = kmin:kmax
        denom = _fact(k) * _fact(j1 + j2 - J - k) * _fact(j1 - m1 - k) *
                _fact(j2 + m2 - k) * _fact(J - j2 + m1 + k) * _fact(J - j1 - m2 + k)
        s += (isodd(k) ? -1.0 : 1.0) / denom
    end
    return pref * s
end

# Deterministic, well-spread sample directions on the unit sphere (Fibonacci
# sphere), used to fit the rotation matrix. No RNG so results are reproducible.
@inline function _fib_dir(k::Int, K::Int)::SVector{3,Float64}
    golden = π * (3.0 - sqrt(5.0))
    z = 1.0 - 2.0 * (k - 0.5) / K
    r = sqrt(max(0.0, 1.0 - z * z))
    θ = golden * k
    return SVector{3,Float64}(r * cos(θ), r * sin(θ), z)
end

"""
    wignerD_real(l, R) -> Matrix{Float64}

Real Wigner-D matrix `Δ` (size `(2l+1)×(2l+1)`) of the rotation `R` (3×3, proper or
improper) in the tesseral-harmonic basis, indexed by `m = -l:l`
(`row/col = m + l + 1`), satisfying

```
Zₗ,m(R·r̂) = Σ_{m'} Δ[m, m'] · Zₗ,m'(r̂)   for all r̂.
```

Built from the package's own `Zₗₘ` by an exact oversampled least-squares fit, so it
is consistent with the harmonic convention by construction.

# Arguments
- `l::Integer`: angular momentum (`l ≥ 0`).
- `R::AbstractMatrix{<:Real}`: a 3×3 orthogonal matrix.

# Returns
- `Matrix{Float64}`: the `(2l+1)×(2l+1)` real Wigner-D matrix.
"""
function wignerD_real(l::Integer, R::AbstractMatrix{<:Real})::Matrix{Float64}
    l >= 0 || throw(ArgumentError("need l ≥ 0 (got $l)"))
    size(R) == (3, 3) || throw(ArgumentError("R must be 3×3"))
    l == 0 && return fill(1.0, 1, 1)
    d = 2l + 1
    K = 3d                      # oversample for conditioning; fit is exact
    Rs = SMatrix{3,3,Float64}(R)
    A = Matrix{Float64}(undef, K, d)
    B = Matrix{Float64}(undef, K, d)
    @inbounds for k = 1:K
        u = _fib_dir(k, K)
        ru = Rs * u
        ru = ru / norm(ru)
        for (col, mp) in enumerate(-l:l)
            A[k, col] = Zlm_unsafe(l, mp, u)
            B[k, col] = Zlm_unsafe(l, mp, ru)
        end
    end
    # Z(R·) = A·Δᵀ  ⇒  Δᵀ = A \ B  (least squares; residual ≈ 0).
    return permutedims(A \ B)
end

end # module AngularMomentum
