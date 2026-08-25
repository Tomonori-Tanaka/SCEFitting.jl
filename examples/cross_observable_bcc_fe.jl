# Cross-observable consistency on REAL DFT data: do the constrained-moment torques
# and the total energies of the same calculations imply the SAME couplings?
#
# Every other gate in this repository closes a loop the code itself draws. The
# analytic Heisenberg torque in `heisenberg_chain.jl` pins `predict_torque` and the
# torque design matrix against `−e × ∇E`, and the unit suite pins `predict_torque`
# against finite differences of `predict_energy`. None of them can see the ONE step
# that leaves the package: VASP writes a constraining field `B_a`, `SpinDatum` turns
# it into the torque target `τ_a = m_a × B_a`, and nothing internal can tell whether
# that is the same torque as `−e_a × ∂E/∂e_a` of the energies in the same file. A
# sign or a factor there would fit happily and be wrong by that factor.
#
# The check: the energies and the fields are INDEPENDENT outputs of the same DFT
# runs, so fit each block alone and predict the other. Both directions have to work,
# and the two coefficient vectors have to be the same vector — same sign, same size.
#
# Data: `docs/src/tutorials/case1_inputs`, the bcc Fe 4×4×4 supercell (128 atoms,
# 50 constrained configurations) of the case-1 tutorial.
#
# Run:  julia --project=examples examples/cross_observable_bcc_fe.jl

using SCEFitting
import Spglib           # `import` only: activate the SpglibBackend extension
using LinearAlgebra

inputs = joinpath(pkgdir(SCEFitting), "docs", "src", "tutorials", "case1_inputs")

function read_poscar(path)
    lines   = readlines(path)
    scale   = parse(Float64, strip(lines[2]))
    A       = reduce(hcat, [parse.(Float64, split(strip(lines[i]))) for i = 3:5]) .* scale
    species = String.(split(strip(lines[6])))
    counts  = parse.(Int, split(strip(lines[7])))
    nat     = sum(counts)
    frac    = reduce(hcat, [parse.(Float64, split(strip(lines[8 + k]))[1:3]) for k = 1:nat])
    kinds   = reduce(vcat, [fill(s, c) for (s, c) in enumerate(counts)])
    return (; lattice = A, frac, kinds, species, nat)
end

poscar  = read_poscar(joinpath(inputs, "POSCAR"))
data    = read_embset(joinpath(inputs, "EMBSET"))     # derives ê and τ = m × B
crystal = Crystal(Lattice(poscar.lattice), poscar.frac, poscar.kinds, poscar.species)

interaction = BasisSpec(; nbody = 2, cutoff = Inf, lmax = [1], isotropy = true)
basis   = SCEBasis(crystal, interaction; backend = SpglibBackend())
dataset = SCEDataset(basis, data; use_torque = true)
println("space group : ", basis.spacegroup.symbol)
println("atoms       : ", poscar.nat, "   configurations: ", length(data),
        "   SALCs: ", n_salcs(basis))

fit_E = fit(SCEFit, dataset, OLS(); torque_weight = 0.0)   # ENERGIES only
fit_T = fit(SCEFit, dataset, OLS(); torque_weight = 1.0)   # TORQUES only

r2_T_from_E = r2_torque(fit_E)      # torques the energy fit never saw
r2_E_from_T = r2_energy(fit_T)      # energies the torque fit never saw
println("\nenergy-only fit → R² on the torques it never saw : ",
        round(r2_T_from_E; digits = 5))
println("torque-only fit → R² on the energies it never saw : ",
        round(r2_E_from_T; digits = 5))

# The two coefficient vectors, compared as vectors: the least-squares scale c in
# coef_T ≈ c·coef_E. c ≈ 1 is agreement; c < 0 is a sign convention flipped
# somewhere between `τ = m × B` and the torque design; c ≈ 2 or ≈ 0.5 is a factor
# slipped into one of them. (c is not exactly 1 — the truncated isotropic pair
# basis is not the true energy surface, so the two observables weight its error
# differently. The 2026-08-25 measurement on this fixture: 0.8729.)
cE, cT = coef(fit_E), coef(fit_T)
c_scale = dot(cE, cT) / dot(cE, cE)
println("least-squares scale in coef_T ≈ c·coef_E        : ", round(c_scale; digits = 5))

# Gates. The bands are set from the mutations they have to resolve, measured on
# this fixture (2026-08-25) by rescaling the torque targets — the equivalent of a
# slip in `τ = m × B` or in the torque design matrix:
#
#     targets ×      R²(τ | E-fit)   R²(E | τ-fit)        c
#     ---------      -------------   -------------   --------
#      1  (as is)         0.9805          0.9864       0.8729
#     −1  (sign flip)    −3.5023         −2.5710      −0.8729
#      2                  0.8055          0.3883       1.7458
#      0.5               −0.5608          0.6913       0.4365
#
# Every mutation trips at least two of the three gates; the `c` band trips on all
# three, which is why it is here alongside the two R².
@assert r2_T_from_E > 0.9 "energy-only fit does not predict the DFT torques (R² = $r2_T_from_E)"
@assert r2_E_from_T > 0.9 "torque-only fit does not predict the DFT energies (R² = $r2_E_from_T)"
@assert 0.6 < c_scale < 1.5 "the two observables imply different couplings (c = $c_scale)"
println("\n✓ energies and constrained-field torques imply the same couplings")
