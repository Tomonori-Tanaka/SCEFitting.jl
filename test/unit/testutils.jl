# Shared unit-test helpers, included exactly once by `runtests.jl` before the unit
# files. They used to be redefined per file in the shared test module — textually
# identical copies that overwrote one another (method-overwrite warnings on every
# run), where an edit to one copy would silently win for every later-included file.
# The bodies are kept byte-identical to the originals so every seeded fixture draws
# the same RNG stream.

using StaticArrays: SVector, SMatrix, @SMatrix
using LinearAlgebra: norm, I
using Random: AbstractRNG

# A uniform random unit direction.
rand_unit(rng::AbstractRNG) =
    (v = SVector{3,Float64}(randn(rng), randn(rng), randn(rng)); v / norm(v))

# Proper rotation from a random axis-angle (Rodrigues).
function rand_rotation(rng::AbstractRNG)
    a = SVector{3,Float64}(randn(rng), randn(rng), randn(rng))
    a = a / norm(a)
    θ = 2π * rand(rng)
    Kc = @SMatrix [0.0 -a[3] a[2]; a[3] 0.0 -a[1]; -a[2] a[1] 0.0]
    return SMatrix{3,3,Float64}(I) + sin(θ) * Kc + (1 - cos(θ)) * (Kc * Kc)
end

# A random spin configuration: `3 × nat`, unit columns (fresh draw per column, so
# neighboring spins differ and odd-Lf / anisotropic channels do not vanish).
function randcfg(rng::AbstractRNG, nat::Integer)
    M = Matrix{Float64}(undef, 3, nat)
    for a = 1:nat
        v = randn(rng, 3)
        M[:, a] = v / norm(v)
    end
    return M
end
