# Heisenberg chain end-to-end: build the SCE basis for a spin chain, fit synthetic
# energies E = J Σ_⟨ij⟩ e_i·e_j, and recover the coupling J.
#
# Run:  julia --project=examples examples/heisenberg_chain.jl

using MagestyRebuild
import Spglib           # `import` (not `using`): just load it to activate the
                        # SpglibBackend extension, without importing names that
                        # clash with MagestyRebuild's `Lattice`/`Crystal`.
using LinearAlgebra
using Random

rng = MersenneTwister(2026)
randspin() = (v = randn(rng, 3); v / norm(v))
randcfg(nat) = mapreduce(_ -> randspin(), hcat, 1:nat)

# A 4-atom chain of identical spins along z (so neighboring spins differ).
lat = Lattice([8.0 0 0; 0 8.0 0; 0 0 10.0])
frac = [0 0 0 0; 0 0 0 0; 0.0 0.25 0.5 0.75]
chain = Crystal(lat, frac, [1, 1, 1, 1], ["Fe"])

# Nearest-neighbor 2-body interaction, isotropic (Heisenberg) channel only.
interaction = Interaction(; nbody = 2, pair_cutoff = 2.6, lmax = [1], isotropy = true)
basis = SCEBasis(chain, interaction; backend = SpglibBackend())
println("space group : ", basis.spacegroup.symbol, " (#", basis.spacegroup.number, ")")
println("# SALCs     : ", nsalc(basis))

# Synthetic Heisenberg training data.
heis = basis.salcs.salcs[1]                       # the nearest-neighbor Heisenberg SALC
J_true = 0.0137
configs = [randcfg(4) for _ = 1:40]
E = [J_true * 0.5 * sum(dot(c[:, m.atoms[1]], c[:, m.atoms[2]]) for m in heis.members)
     for c in configs]

# Fit and recover J (the SALC normalization gives J = 2√3 · jphi).
f = fit(SCEFit, SCEDataset(basis, configs, E), OLS())
J_recovered = 2 * sqrt(3) * coef(f)[1]

println("R²          : ", round(r2_energy(f); digits = 12))
println("J (true)    : ", J_true)
println("J (recovered): ", J_recovered)
@assert isapprox(J_recovered, J_true; rtol = 1e-6)
println("✓ recovered the Heisenberg coupling")
