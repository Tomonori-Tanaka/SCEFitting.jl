# The pinned fixtures: five real crystals, carried as RAW DATA so this file
# depends on no package type.  They are the magnetic rows of the integration
# roster (`test/integration/roster.jl`) minus the 2x2x2 supercell, which is the
# same content on 8x the cell and would multiply the pin's size for nothing.
#
#   bcc_Fe    Im-3m, one species, 3-body        -- the high-symmetry baseline
#   B2_FeRh   Pm-3m, two species, 3-body        -- two magnetic sublattices
#   hcp_Co    P6_3/mmc, non-orthogonal cell     -- two like sites, hexagonal
#   wz_GaN    P6_3mc, NO inversion centre       -- four atoms, no centre
#   rs_MnO    Fm-3m, a species with lmax = 0    -- the spin-free-species path
#
# Specs are pure spin (dense, sector-less): that is the surface both packages
# share, and the only one a pin can be written against once.

_pin_diag(a) = [a 0.0 0.0; 0.0 a 0.0; 0.0 0.0 a]

const PIN_FIXTURES = [
    (id = "bcc_Fe", L = _pin_diag(2.87), frac = [0.0 0.5; 0.0 0.5; 0.0 0.5],
     species = [1, 1], labels = ["Fe"], nbody = 3, lmax = [2], cutoff = 2.6),
    (id = "B2_FeRh", L = _pin_diag(2.99), frac = [0.0 0.5; 0.0 0.5; 0.0 0.5],
     species = [1, 2], labels = ["Fe", "Rh"], nbody = 3, lmax = [2, 2], cutoff = 2.7),
    (id = "hcp_Co",
     L = [2.507 -2.507/2 0.0; 0.0 2.507*sqrt(3)/2 0.0; 0.0 0.0 4.069],
     frac = [1/3 2/3; 2/3 1/3; 1/4 3/4], species = [1, 1], labels = ["Co"],
     nbody = 2, lmax = [2], cutoff = 2.6),
    (id = "wz_GaN",
     L = [3.189 -3.189/2 0.0; 0.0 3.189*sqrt(3)/2 0.0; 0.0 0.0 5.185],
     frac = [1/3 2/3 1/3 2/3; 2/3 1/3 2/3 1/3; 0.0 0.5 0.377 0.877],
     species = [1, 1, 2, 2], labels = ["Ga", "N"], nbody = 2, lmax = [1, 1],
     cutoff = 2.0),
    (id = "rs_MnO", L = _pin_diag(4.446),
     frac = [0.0 0.5 0.5 0.0 0.5 0.5 0.0 0.0;
             0.0 0.5 0.0 0.5 0.5 0.0 0.5 0.0;
             0.0 0.0 0.5 0.5 0.5 0.0 0.0 0.5],
     species = [1, 1, 1, 1, 2, 2, 2, 2], labels = ["Mn", "O"], nbody = 2,
     lmax = [2, 0], cutoff = 3.2),
]

# The pointed (site-marked) moment fixtures.  A SEPARATE roster: the moment
# channel has its own basis type, its own design rows (one per marked atom per
# configuration) and its own fit, so nothing above can stand in for it.
#
#   B2_FeRh_pointed   Pm-3m 2x2x2 supercell, Fe marked / Rh not
#
# The cell is a SUPERCELL on purpose.  On the 2-atom primitive cell every
# neighbour is reachable through several minimum images, the star enumeration
# keeps those tied instances, and `moment_resolvability` refuses the basis
# outright (`UnclassifiableBasis`) — measured on the 2-atom bcc Fe, B2 FeRh and
# 8-atom rock-salt MnO cells.  2x2x2 is the smallest cell here with untied
# nearest-neighbour stars.
#
# `cutoff_star` admits only the Fe-Rh nearest-neighbour bonds, which is what
# keeps the payload to the size of the pure-spin pins: the pointed member count
# carries a factor N! per instance, so a wider star multiplies the `folded`
# dump by an order of magnitude for no extra structure.

_pin_super(a, base, sp, L) = begin
    pos = Float64[]
    species = Int[]
    for i = 0:L-1, j = 0:L-1, k = 0:L-1, (b, s) in zip(base, sp)
        append!(pos, [(i + b[1]) / L, (j + b[2]) / L, (k + b[3]) / L])
        push!(species, s)
    end
    (_pin_diag(Float64(L) * a), reshape(pos, 3, :), species)
end

const _PIN_FERH_L, _PIN_FERH_FRAC, _PIN_FERH_SP =
    _pin_super(2.99, [(0.0, 0.0, 0.0), (0.5, 0.5, 0.5)], [1, 2], 2)

const PIN_MOMENT_FIXTURES = [
    (id = "B2_FeRh_pointed", L = _PIN_FERH_L, frac = _PIN_FERH_FRAC,
     species = _PIN_FERH_SP, labels = ["Fe", "Rh"], nbody = 3, lmax_env = [1, 1],
     lmax_mark = 1, marked = [true, false], cutoff_pair = 2.7, cutoff_star = 2.7,
     lsum = 2),
]
