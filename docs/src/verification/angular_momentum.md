# Verification: angular momentum

```@meta
CurrentModule = SCEFitting
```

The SALC projection rests on three pieces of angular-momentum machinery, all
implemented from scratch in the public (unexported) [`AngularMomentum`](@ref)
submodule: the Clebsch–Gordan coefficients (Racah single-sum formula,
Condon–Shortley phase), the real Wigner-``D`` matrices in the tesseral basis, and
the complex→real unitary ``U^{(l)}``. This page renders the checks of
`test/unit/test_angmom.jl` in human-readable form — every table and summary below
is **computed when the documentation is built**, and each section ends in an
assertion, so a regression fails the (strict) docs build just as it fails the test
suite. The unit tests remain the machine gate; this page is the same evidence,
made readable.

## Clebsch–Gordan coefficients: known values

Closed-form values from the standard tables (Varshalovich, *Quantum Theory of
Angular Momentum*), including three selection-rule zeros. The notation
`⟨j₁ m₁; j₂ m₂ → J M⟩` stands for ``\langle j_1 m_1\, j_2 m_2 | J M \rangle``:

```@example am
using SCEFitting
using Printf, Markdown

const CG = SCEFitting.AngularMomentum.clebsch_gordan

known = [
    (1, 0, 1, 0, 0, 0,   "−1/√3",              -1 / sqrt(3)),
    (1, 1, 1, -1, 0, 0,  "1/√3",               1 / sqrt(3)),
    (1, 1, 1, -1, 1, 0,  "1/√2",               1 / sqrt(2)),
    (1, 1, 1, -1, 2, 0,  "1/√6",               1 / sqrt(6)),
    (2, 0, 2, 0, 0, 0,   "1/√5",               1 / sqrt(5)),
    (1, 1, 1, 1, 2, 2,   "1 (stretched)",      1.0),
    (1, 1, 1, 1, 0, 0,   "0 (M conservation)", 0.0),
    (1, 0, 1, 0, 3, 0,   "0 (triangle rule)",  0.0),
    (1, 0, 1, 0, 1, 0,   "0 (odd-J symmetry)", 0.0),
]

function known_value_table(rows; atol = 1e-14)
    io = IOBuffer()
    println(io, "| coefficient | exact | computed | abs. error | pass |")
    println(io, "|:---|:---|---:|---:|:---:|")
    worst = 0.0
    for (j1, m1, j2, m2, J, M, expr, exact) in rows
        val = CG(j1, m1, j2, m2, J, M)
        err = abs(val - exact)
        worst = max(worst, err)
        @printf(io, "| ⟨%d %d; %d %d → %d %d⟩ | %s | %.15f | %.1e | %s |\n",
                j1, m1, j2, m2, J, M, expr, val, err, err <= atol ? "✔" : "✗")
    end
    worst <= atol || error("Clebsch–Gordan known-value check failed (worst $worst)")
    return Markdown.parse(String(take!(io)))
end

known_value_table(known)
```

## Clebsch–Gordan coefficients: orthogonality sweeps

Individual values can be right by accident; the two orthogonality relations tie
the *whole* coefficient family together. First,

```math
\sum_{m_1 m_2} \langle j_1 m_1\, j_2 m_2 | J M \rangle
              \langle j_1 m_1\, j_2 m_2 | J' M \rangle = \delta_{J J'},
```

swept over every ``(j_1, j_2, J, J', M)`` with ``j_1, j_2 \le 3``:

```@example am
function cg_orthonormality(jmax)
    checked = 0
    maxdev = 0.0
    for j1 = 0:jmax, j2 = 0:jmax
        for J = abs(j1 - j2):(j1 + j2), Jp = abs(j1 - j2):(j1 + j2), M = -J:J
            s = 0.0
            for m1 = -j1:j1, m2 = -j2:j2
                s += CG(j1, m1, j2, m2, J, M) * CG(j1, m1, j2, m2, Jp, M)
            end
            maxdev = max(maxdev, abs(s - (J == Jp ? 1.0 : 0.0)))
            checked += 1
        end
    end
    return (sums_checked = checked, max_deviation = maxdev)
end

ortho = cg_orthonormality(3)
```

and the completeness relation in the other direction,

```math
\sum_{J M} \langle j_1 m_1\, j_2 m_2 | J M \rangle
           \langle j_1 m_1'\, j_2 m_2' | J M \rangle
         = \delta_{m_1 m_1'} \delta_{m_2 m_2'},
```

swept over every ``(m_1, m_2, m_1', m_2')`` for the same ``j_1, j_2 \le 3``:

```@example am
function cg_completeness(jmax)
    checked = 0
    maxdev = 0.0
    for j1 = 0:jmax, j2 = 0:jmax
        for m1 = -j1:j1, m2 = -j2:j2, m1p = -j1:j1, m2p = -j2:j2
            s = 0.0
            for J = abs(j1 - j2):(j1 + j2)
                s += CG(j1, m1, j2, m2, J, m1 + m2) *
                     CG(j1, m1p, j2, m2p, J, m1 + m2)
            end
            expected = (m1 == m1p && m2 == m2p) ? 1.0 : 0.0
            maxdev = max(maxdev, abs(s - expected))
            checked += 1
        end
    end
    return (sums_checked = checked, max_deviation = maxdev)
end

complete = cg_completeness(3)
```

```@example am
max(ortho.max_deviation, complete.max_deviation) <= 1e-12 ||
    error("Clebsch–Gordan orthogonality sweep failed")
println("both sweeps pass at the 1e-12 tolerance")
```

## Real Wigner-D matrices

[`wignerD_real`](@ref SCEFitting.AngularMomentum.wignerD_real) must satisfy the
functional identity ``Z_{l m}(R\hat{\boldsymbol u}) = \sum_{m'} \Delta^{(l)}_{m m'}(R)\,
Z_{l m'}(\hat{\boldsymbol u})`` for *every* direction — the SALC projector applies it to
arbitrary spin configurations. The sweep below checks that identity on fresh
random directions (never used in the internal least-squares construction), for
proper **and** improper operations (the second half of the rotations is composed
with the inversion), plus orthogonality ``\Delta \Delta^{\mathsf T} = I`` and the
inversion parity ``\Delta(-I) = (-1)^l I``:

```@example am
using LinearAlgebra, Random
using SCEFitting: Harmonics

const Dreal = SCEFitting.AngularMomentum.wignerD_real

rand_unit(rng) = normalize(randn(rng, 3))
function rand_rotation(rng)
    Q = Matrix(qr(randn(rng, 3, 3)).Q)
    Q[:, 1] .*= sign(det(Q))            # make it proper; improper = -Q below
    return Q
end

function wignerD_deviations(; lmax = 4, nrot = 20, npts = 5, seed = 7)
    rng = MersenneTwister(seed)
    dev_orth = 0.0
    dev_fun = 0.0
    for k = 1:2nrot
        R = rand_rotation(rng)
        k > nrot && (R = -R)            # improper half: inversion ∘ rotation
        for l = 0:lmax
            Δ = Dreal(l, R)
            dev_orth = max(dev_orth, opnorm(Δ * Δ' - I))
            for _ = 1:npts
                u = rand_unit(rng)
                lhs = [Harmonics.Zlm(l, m, R * u) for m = -l:l]
                rhs = Δ * [Harmonics.Zlm(l, m, u) for m = -l:l]
                dev_fun = max(dev_fun, maximum(abs.(lhs - rhs)))
            end
        end
    end
    parity = maximum(0:lmax) do l
        opnorm(Dreal(l, -Matrix(1.0I, 3, 3)) - (-1.0)^l * I)
    end
    return (orthogonality = dev_orth, functional_identity = dev_fun,
            inversion_parity = parity)
end

wd = wignerD_deviations()
```

```@example am
maximum(wd) <= 1e-9 || error("Wigner-D verification failed")
println("all Wigner-D checks pass at the 1e-9 tolerance ",
        "(200 operations × l ≤ 4 × 5 fresh directions)")
```

## Complex→real unitary

[`c2r_matrix`](@ref SCEFitting.AngularMomentum.c2r_matrix) is compared against the
closed form of the tesseral transformation (the paper's Eq. B2) — a
convention-independent anchor written out digit by digit, plus unitarity:

```@example am
const c2r = SCEFitting.AngularMomentum.c2r_matrix

function c2r_deviations(lmax)
    dev_form = 0.0
    dev_unit = 0.0
    for l = 0:lmax
        U = c2r(l)
        ix(m) = m + l + 1
        ref = zeros(ComplexF64, 2l + 1, 2l + 1)
        ref[ix(0), ix(0)] = 1
        for m = 1:l
            s = 1 / sqrt(2.0)
            ref[ix(m), ix(m)] = (-1.0)^m * s
            ref[ix(m), ix(-m)] = s
            ref[ix(-m), ix(-m)] = im * s
            ref[ix(-m), ix(m)] = -im * (-1.0)^m * s
        end
        dev_form = max(dev_form, maximum(abs.(U - ref)))
        dev_unit = max(dev_unit, opnorm(U' * U - I))
    end
    return (closed_form = dev_form, unitarity = dev_unit)
end

u = c2r_deviations(3)
```

```@example am
maximum(u) <= 1e-12 || error("complex→real unitary verification failed")
println("closed form and unitarity pass at the 1e-12 tolerance (l ≤ 3)")
```

## What this page does and does not gate

- Every number above is recomputed on each docs build; `warnonly = false` means a
  failed assertion breaks the build, so the published page can never show a stale
  green.
- The sweeps here mirror (and slightly extend — the completeness relation, ``j
  \le 3``) the unit tests in `test/unit/test_angmom.jl`; the test suite is still
  the primary, faster gate run on every change.
- Downstream consistency — that these pieces assemble into symmetry-invariant
  SALCs — is checked elsewhere: the invariance and idempotency gates of
  `test/unit/test_salc.jl` / `test_nbody.jl` and the pinned-oracle comparison in
  `test/oracle/`.
