# Cross-orbit alias groups (src/basis/aliases.jl): distinct orbits that join the same
# reference-cell atoms through periodic images the space group does not relate. On a
# cell-periodic training set their design columns coincide, so the data fix only the
# sum of their couplings. The package ties them into one design column and reads the
# coefficient back with equal per-bond weight, recorded as a convention.
#
# Every fixture uses `NoSymmetry()` so that each cluster instance is its own orbit —
# the situation a low-symmetry (or multi-site, as Nd2Fe14B P4_2/mnm) cell produces.
# Oracles: hand geometry (which images tie), the documented isotropic-pair SALC
# closed form `Φ = 2√3·(e_a·e_b)` (the v0 done-line), and the closed-form
# consequence `J_tied = (J₊ + J₋)/2` of tying two identical columns. The 2×1×1
# recovery uses a HAND-WRITTEN bond-sum energy, never the package's predictor.

using Test
using SCEFitting
using SCEFitting: alias_groups, n_columns, salcs, evaluate_salc, salc_groups,
    group_costs, penalty_metric, multipole_terms
using LinearAlgebra
using StaticArrays
using Random

const MA = SCEFitting

# P1 cell: atoms 1 and 2 differ by exactly a/2 along x, so their pair has two
# equidistant minimum images (±x) that no operation relates. Atom 3 sits off every
# tie. Hand check: d₊ = (1.5, 0.6, 0.3), d₋ = (−1.5, 0.6, 0.3), |d| = 1.6432 Å;
# pairs 1–3 (1.72 Å) and 2–3 (1.25 Å) are unique below the 2.6 Å cutoff and their
# next images lie beyond it (2.73 Å).
const _P1_FRAC = [0.0 0.5 0.25; 0.0 0.2 0.3; 0.0 0.1 0.42]
_p1_crystal() = Crystal(Lattice(Matrix(3.0 * I(3))), _P1_FRAC, [1, 1, 1], ["Fe"])
_p1_spec() = BasisSpec(; nbody = 2, cutoff = 2.6, lmax = [1], isotropy = true)

# 2×1×1 supercell of the same crystal: the two images become the distinct bonds
# 1–2 / 1'–2' (d₊) and 1'–2 / 1–2' (d₋). Atom order: 1, 2, 3, 1', 2', 3'.
function _p1_supercell()
    fr = hcat(_P1_FRAC ./ [2.0, 1.0, 1.0], (_P1_FRAC .+ [1.0, 0.0, 0.0]) ./ [2.0, 1.0, 1.0])
    return Crystal(Lattice([6.0 0 0; 0 3.0 0; 0 0 3.0]), fr, fill(1, 6), ["Fe"])
end

# Cartesian separation of a pair member (atom₂ + R₂) − (atom₁ + R₁).
function _pair_vector(cr::Crystal, m)
    A = cr.lattice.vectors
    fr = cr.frac_positions
    return A * (SVector{3,Float64}(fr[:, m.atoms[2]]) + m.shifts[2] -
                SVector{3,Float64}(fr[:, m.atoms[1]]) - m.shifts[1])
end

@testset "cross-orbit alias groups" begin
    rng = MersenneTwister(2026)
    cr = _p1_crystal()
    spec = _p1_spec()

    @testset "P1 a/2 tie: detection, closed form, tie, equal split" begin
        basis = @test_logs (:info, r"cross-orbit alias groups: 1 group") match_mode = :any SCEBasis(cr, spec)
        @test n_salcs(basis) == 4
        @test n_columns(basis) == 3
        groups = alias_groups(basis)
        @test length(groups) == 1
        g = groups[1]
        @test g.kind === :proportional
        @test g.body == 2
        @test length(g.salcs) == 2 && length(unique(g.orbit_ids)) == 2
        @test g.weights == [1.0, 1.0]                      # transported copies: exact ±1
        @test g.distance ≈ norm([1.5, 0.6, 0.3]) rtol = 1e-12
        # the two orbits' matching members differ by one lattice vector along a
        @test g.delta_shifts[1] == [zero(SVector{3,Int}), zero(SVector{3,Int})]
        @test abs(g.delta_shifts[2][2][1]) == 1 && g.delta_shifts[2][2][2:3] == [0, 0]
        # both tied members are the 1–2 pair, one per image (hand geometry)
        for j in g.salcs
            m = salcs(basis)[j].members
            @test length(m) == 1 && m[1].atoms == [1, 2]
            @test norm(_pair_vector(cr, m[1])) ≈ g.distance rtol = 1e-12
        end
        @test Set(sign(_pair_vector(cr, salcs(basis)[j].members[1])[1]) for j in g.salcs) ==
              Set([1.0, -1.0])
        # closed form: both aliased SALCs ARE 2√3·(e₁·e₂) on cell-periodic spins
        for _ = 1:5
            e = randcfg(rng, 3)
            for j in g.salcs
                @test evaluate_salc(salcs(basis)[j], e) ≈ 2sqrt(3) * dot(e[:, 1], e[:, 2]) rtol = 1e-12
            end
        end
        # the design carries the tied column once: bitwise the sum of the two SALC
        # columns of the untied build (the fold is an exact `x + x` for ±1 weights),
        # at whatever thread count the suite runs
        cfgs = [randcfg(rng, 3) for _ = 1:30]
        ds = SCEDataset(basis, cfgs, zeros(30))
        @test size(ds.X_E, 2) == 3
        @test ds.X_E[:, 1] ≈ [2 * 2sqrt(3) * dot(c[:, 1], c[:, 2]) for c in cfgs] rtol = 1e-12
        untied = SCEDataset(SCEBasis(cr, spec; alias_rtol = nothing), cfgs, zeros(30),
                            [zeros(3, 3) for _ = 1:30])
        tied = SCEDataset(basis, cfgs, zeros(30), [zeros(3, 3) for _ = 1:30])
        @test tied.X_E[:, 1] == untied.X_E[:, 1] .+ untied.X_E[:, 2]
        @test tied.X_E[:, 2:3] == untied.X_E[:, 3:4]
        @test tied.X_T[:, 1] == untied.X_T[:, 1] .+ untied.X_T[:, 2]
        @test tied.X_T[:, 2:3] == untied.X_T[:, 3:4]

        # synthetic truth with J₊ ≠ J₋: the fit returns (J₊ + J₋)/2 to BOTH orbits
        # (closed form: E = (J₊ + J₋)·Φ on periodic spins, the tied column is 2Φ)
        Jp, Jm = 0.7, 0.3
        jtrue = zeros(n_salcs(basis))
        jtrue[g.salcs[1]] = Jp
        jtrue[g.salcs[2]] = Jm
        jtrue[3] = -0.2
        jtrue[4] = 0.45
        gen = SCEPredictor(basis, 0.1, jtrue)          # hand-set model: split = :free
        @test gen.split == fill(:free, 4)
        E = predict_energy(gen, cfgs)
        T = [predict_torque(gen, c) for c in cfgs]
        f = fit(SCEFit, SCEDataset(basis, cfgs, E, T), OLS(); torque_weight = 0.5)
        @test length(coef(f)) == 3                      # column space
        @test coef(f)[1] ≈ (Jp + Jm) / 2 rtol = 1e-10
        m = SCEPredictor(f)
        @test length(m.jphi) == 4                       # SALC space
        @test m.jphi[g.salcs[1]] == m.jphi[g.salcs[2]]  # bit-equal: one number
        @test m.jphi[g.salcs[1]] ≈ (Jp + Jm) / 2 rtol = 1e-10
        @test m.jphi[3] ≈ -0.2 rtol = 1e-10
        @test m.jphi[4] ≈ 0.45 rtol = 1e-10
        @test m.split == [:convention, :convention, :free, :free]
        # the training cell is reproduced exactly — energies AND torques
        @test predict_energy(m, cfgs) ≈ E rtol = 1e-12
        for (c, t) in zip(cfgs, T)
            @test predict_torque(m, c) ≈ t atol = 1e-12
        end
        @test intercept(f) ≈ 0.1 rtol = 1e-10
        # coeftable discloses the group and the convention
        ct = coeftable(f)
        @test length(ct) == 4
        @test [r.alias_group for r in ct] == [1, 1, 0, 0]
        @test [r.split for r in ct] == [:convention, :convention, :free, :free]
        @test [r.J for r in ct] == m.jphi
        # read-out doors warn once
        @test_logs (:warn, r"multipole_terms") match_mode = :any multipole_terms(m)
        @test_logs (:warn, r"bilinear_terms") match_mode = :any bilinear_terms(m)
        @test_logs multipole_terms(gen)                  # hand-set: nothing to warn about
    end

    @testset "tetragonal P cell, Δz = c/2" begin
        # a = b = 3.4, c = 6: atom 2 at (0.3, 0.1, 0.5) has two images ±c/2 apart from
        # atom 1, d = (1.02, 0.34, ±3), |d| = 3.187 Å; the next image (x − a) is 3.85 Å
        # and the self pairs sit at 3.4 Å, both above the 3.3 Å cutoff.
        crt = Crystal(Lattice([3.4 0 0; 0 3.4 0; 0 0 6.0]), [0.0 0.3; 0.0 0.1; 0.0 0.5],
                      [1, 1], ["Fe"])
        bt = SCEBasis(crt, BasisSpec(; nbody = 2, cutoff = 3.3, lmax = [1], isotropy = true))
        @test n_salcs(bt) == 2 && n_columns(bt) == 1
        gt = alias_groups(bt)
        @test length(gt) == 1 && gt[1].kind === :proportional
        @test gt[1].distance ≈ norm([1.02, 0.34, 3.0]) rtol = 1e-12
        @test abs(gt[1].delta_shifts[2][2][3]) == 1 && gt[1].delta_shifts[2][2][1:2] == [0, 0]
        @test SCEFitting.salc_groups(bt) == [1]
    end

    @testset "the split is decided by the cell: 2×1×1 recovers J₊ ≠ J₋, 1×1×1 their mean" begin
        Jp, Jm, J13, J23, j0 = 0.7, 0.3, -0.2, 0.45, 0.1
        cr2 = _p1_supercell()
        b2 = SCEBasis(cr2, spec)
        @test isempty(alias_groups(b2))                     # the doubling breaks the tie
        # hand count of the admitted pairs (cutoff 2.6 Å): d₊ 1–2, 1'–2'; d₋ 1'–2, 1–2';
        # 1–3, 1'–3' (1.72 Å); 2–3, 2'–3' (1.25 Å); and 2–3', 2'–3 (2.46 Å — the image
        # of 2–3 that the 1×1×1 cell drops as non-minimum, a distinct pair here)
        @test n_salcs(b2) == 10 && n_columns(b2) == 10
        # hand-written bond-sum energy on the doubled cell (atoms 1,2,3,1',2',3'):
        #   d₊ bonds 1–2, 1'–2'; d₋ bonds 1'–2, 1–2'; 1–3, 1'–3'; 2–3, 2'–3'; the two
        #   2.46 Å bonds carry no coupling in the truth (and must be recovered as 0)
        Φ(e, a, b) = 2sqrt(3) * dot(e[:, a], e[:, b])
        E2(e) = j0 + Jp * (Φ(e, 1, 2) + Φ(e, 4, 5)) + Jm * (Φ(e, 4, 2) + Φ(e, 1, 5)) +
                J13 * (Φ(e, 1, 3) + Φ(e, 4, 6)) + J23 * (Φ(e, 2, 3) + Φ(e, 5, 6))
        cfg2 = [randcfg(rng, 6) for _ = 1:60]
        f2 = fit(SCEFit, SCEDataset(b2, cfg2, E2.(cfg2)), OLS())
        m2 = SCEPredictor(f2)
        @test intercept(f2) ≈ j0 rtol = 1e-8
        # classify every SALC of b2 by the bond it carries and compare with the truth
        for (k, s) in enumerate(salcs(b2))
            mem = s.members
            @test length(mem) == 1
            a, b = mem[1].atoms
            d = _pair_vector(cr2, mem[1])
            mod1(a, 3) > mod1(b, 3) && (d = -d)           # orient from the 1-type to the 2-type atom
            pair = Set([mod1(a, 3), mod1(b, 3)])          # which primitive-cell pair
            expected = pair == Set([1, 2]) ? (d[1] > 0 ? Jp : Jm) :
                       pair == Set([1, 3]) ? J13 :
                       norm(d) < 2.0 ? J23 : 0.0
            if expected == 0.0
                @test abs(m2.jphi[k]) < 1e-8
            else
                @test m2.jphi[k] ≈ expected rtol = 1e-8
            end
            @test m2.split[k] === :free
        end
        # fold the same truth onto the 1×1×1 cell: periodic configurations have
        # e' = e, so E₂ = 2·E₁ and the tied coefficient must be the per-bond mean
        b1 = SCEBasis(cr, spec)
        g = alias_groups(b1)[1]
        cfg1 = [randcfg(rng, 3) for _ = 1:30]
        E1 = [E2(hcat(c, c)) / 2 for c in cfg1]
        f1 = fit(SCEFit, SCEDataset(b1, cfg1, E1), OLS())
        m1 = SCEPredictor(f1)
        for j in g.salcs
            @test m1.jphi[j] ≈ (Jp + Jm) / 2 rtol = 1e-8
            @test m1.split[j] === :convention
        end
        @test intercept(f1) ≈ j0 / 2 rtol = 1e-8
    end

    @testset "mutation: alias_rtol = nothing disables the tie and the rank warning returns" begin
        b0 = SCEBasis(cr, spec; alias_rtol = nothing)
        @test isempty(alias_groups(b0))
        @test n_columns(b0) == n_salcs(b0) == 4
        @test b0.ties.trivial
        cfgs = [randcfg(rng, 3) for _ = 1:30]
        gen = SCEPredictor(b0, 0.0, [0.7, 0.3, -0.2, 0.45])
        E = predict_energy(gen, cfgs)
        @test_logs (:warn, r"rank deficient") match_mode = :any fit(SCEFit, SCEDataset(b0, cfgs, E), OLS())
        @test_throws ArgumentError SCEBasis(cr, spec; alias_rtol = -1.0)
        @test_throws ArgumentError SCEBasis(cr, spec; alias_rtol = Inf)
        # `alias_rtol = 0` keeps an EXACT alias (the P1 copies are bit-identical), which
        # is why the mutation above uses `nothing` rather than a zero threshold
        @test length(alias_groups(SCEBasis(cr, spec; alias_rtol = 0.0))) == 1
    end

    @testset "threshold band (change detector, not correctness)" begin
        # Pins the residual scale of a genuine alias on this fixture far below
        # `_ALIAS_RTOL`; recapture only if the SALC gauge fixing changes. The upper end
        # of the band (non-alias same-channel residual ≥ 1e-2) is measured on MnTe
        # upstream and is not reproducible on a fixture this small.
        basis = SCEBasis(cr, spec)
        g = alias_groups(basis)[1]
        va, ra = MA._function_vector(salcs(basis)[g.salcs[1]])
        vb, rb = MA._function_vector(salcs(basis)[g.salcs[2]])
        res, c = MA._proportional_residual(va, MA._dictnorm(va), vb, MA._dictnorm(vb))
        @test res <= 1e-12
        @test res < MA._ALIAS_RTOL
        @test c ≈ 1.0 rtol = 1e-12
        @test MA._ALIAS_RTOL == 1e-6
    end

    @testset "column-space helpers: groups, costs, metric, group estimator" begin
        basis = SCEBasis(cr, spec)
        b0 = SCEBasis(cr, spec; alias_rtol = nothing)
        g = alias_groups(basis)[1]
        lab = salc_groups(basis)
        @test length(lab) == n_columns(basis) == 3
        @test lab == [1, 2, 3]
        lab0 = salc_groups(b0)
        @test lab0 == [1, 2, 3, 4]
        # the tied column's Monte-Carlo cost is the sum of the two orbits' costs — both
        # bonds are priced, because after tiling they are distinct bonds
        c = group_costs(basis, lab)
        c0 = group_costs(b0, lab0)
        @test c[1] == c0[1] + c0[2]
        @test c[2:3] == c0[3:4]
        # identical columns tied by summation: Var[2Φ] = 4·Var[Φ] on the same ensemble
        mtr = penalty_metric(basis)
        mtr0 = penalty_metric(b0)
        @test length(mtr) == 3
        @test mtr[1] ≈ 4 * mtr0[1] rtol = 1e-12
        @test mtr[2:3] == mtr0[3:4]
        mt = penalty_metric(basis; torque_weight = 0.5)
        mt0 = penalty_metric(b0; torque_weight = 0.5)
        @test mt[1] ≈ 4 * mt0[1] rtol = 1e-12
        # the group estimator and the λ-path run in column space
        cfgs = [randcfg(rng, 3) for _ = 1:40]
        gen = SCEPredictor(basis, 0.0, [0.7, 0.3, -0.2, 0.45])
        E = predict_energy(gen, cfgs)
        T = [predict_torque(gen, c) for c in cfgs]
        ds = SCEDataset(basis, cfgs, E, T)
        est = GroupAdaptiveRidge(basis; lambda = 1e-6, torque_weight = 0.5)
        @test length(est.column_groups) == 3
        fg = fit(SCEFit, ds, est; torque_weight = 0.5)
        @test length(coef(fg)) == 3
        mg = SCEPredictor(fg)
        @test mg.jphi[1] == mg.jphi[2]
        @test mg.jphi[1] ≈ 0.5 rtol = 1e-4
        fr = refit(fg)
        @test length(coef(fr)) == 3
        @test SCEPredictor(fr).jphi[1] ≈ 0.5 rtol = 1e-10
    end

    @testset "persistence: v7 carries the split, pre-v7 loads as :legacy" begin
        basis = SCEBasis(cr, spec)
        cfgs = [randcfg(rng, 3) for _ = 1:30]
        gen = SCEPredictor(basis, 0.1, [0.7, 0.3, -0.2, 0.45])
        E = predict_energy(gen, cfgs)
        m = SCEPredictor(fit(SCEFit, SCEDataset(basis, cfgs, E), OLS()))
        doc = MA._to_doc(m)
        @test Int(doc["schema_version"]) == 7
        @test [c["split"] for c in doc["couplings"]] == ["convention", "convention", "free", "free"]
        m7 = MA._model_from_doc(doc)
        @test m7.split == m.split
        @test m7.jphi == m.jphi
        @test length(alias_groups(m7.basis)) == 1               # recomputed on load
        @test predict_energy(m7, cfgs) == predict_energy(m, cfgs)
        # a v6 document has no split: unknown provenance inside the group
        doc6 = MA._to_doc(m)
        doc6["schema_version"] = 6
        for c in doc6["couplings"]
            delete!(c, "split")
        end
        m6 = MA._model_from_doc(doc6)
        @test m6.split == [:legacy, :legacy, :free, :free]
        @test m6.jphi == m.jphi
        # a v7 coupling without its split, or with an unknown one, is refused
        bad = MA._to_doc(m)
        delete!(bad["couplings"][1], "split")
        @test_throws ArgumentError MA._model_from_doc(bad)
        bad2 = MA._to_doc(m)
        bad2["couplings"][1]["split"] = "measured"
        @test_throws ArgumentError MA._model_from_doc(bad2)
        # the hand-set model round-trips as :free
        @test MA._model_from_doc(MA._to_doc(gen)).split == fill(:free, 4)
    end

    @testset "alias-free bases are untouched (identity fast path)" begin
        crf = Crystal(Lattice(Matrix(3.0 * I(3))), [0.2 -0.2; 0.0 0.0; 0.0 0.0], [1, 1], ["Fe"])
        bf = SCEBasis(crf, BasisSpec(; nbody = 2, cutoff = 1.5, lmax = [2], isotropy = false))
        bf0 = SCEBasis(crf, BasisSpec(; nbody = 2, cutoff = 1.5, lmax = [2], isotropy = false);
                       alias_rtol = nothing)
        @test isempty(alias_groups(bf)) && bf.ties.trivial
        @test n_columns(bf) == n_salcs(bf)
        cfgs = [randcfg(rng, 2) for _ = 1:12]
        tq = [randn(rng, 3, 2) for _ = 1:12]
        d1 = SCEDataset(bf, cfgs, collect(1.0:12), tq)
        d0 = SCEDataset(bf0, cfgs, collect(1.0:12), tq)
        @test d1.X_E == d0.X_E && d1.X_T == d0.X_T             # bitwise
        @test salc_groups(bf) == salc_groups(bf0)
        @test penalty_metric(bf) == penalty_metric(bf0)
        f = fit(SCEFit, d1, OLS(); torque_weight = 0.5)
        @test coef(f) == SCEPredictor(f).jphi
        @test SCEPredictor(f).split == fill(:free, n_salcs(bf))
        @test [r.alias_group for r in coeftable(f)] == zeros(Int, n_salcs(bf))
        @test_logs multipole_terms(SCEPredictor(f))
    end
end

# A tie the point group fuses only PARTIALLY: a=b=4, c=5, atoms at (0,0,0),
# (0.5,0.5,0.25), (0.3,0.3,0.5) with the diagonal mirror x↔y (space group Cm, assembled
# by hand). The 1–2 pair has four equidistant corner images (±2, ±2, 1.25), |d| = 3.09 Å
# < 3.1 Å; the mirror fuses (2,−2) with (−2,2) into one orbit of two members and leaves
# (2,2) and (−2,−2) as one-member orbits — member counts [1, 2, 1]. Hand oracle: with
# per-orbit truths J_i the cell determines Σ n_i J_i only, and the tied read-out must be
# the bond-weighted mean Σ n_i J_i / Σ n_i. This is the case that separates "per bond"
# from "per orbit": a J/k read-out would return Σ J_i / 3.
@testset "partially fused corner tie (Cm): equal split is per bond, not per orbit" begin
    lat = Lattice(SMatrix{3,3,Float64}([4.0 0 0; 0 4.0 0; 0 0 5.0]))
    crystal = Crystal(lat, [0.0 0.5 0.3; 0.0 0.5 0.3; 0.0 0.25 0.5], [1, 1, 1], ["Fe"])
    mirror = SMatrix{3,3,Float64}([0 1 0; 1 0 0; 0 0 1])
    sg = SCEFitting._assemble_spacegroup(crystal, [SMatrix{3,3,Float64}(I), mirror],
                                         [zero(SVector{3,Float64}), zero(SVector{3,Float64})],
                                         "Cm(manual)", 8; tol = 1e-6)
    spec = BasisSpec(; nbody = 2, cutoff = 3.1, lmax = [1], isotropy = true)
    nl = build_neighbor_list(crystal, SCEFitting._superset_cutoff(spec), MinimumImage())
    clusters = build_clusters(crystal, nl, sg; nbody = 2, selection = MinimumImage(),
                              cutoff = spec.cutoff)
    sb = build_salc_basis(crystal, sg, clusters; lmax_by_species = spec.lmax,
                          lsum_by_body = spec.lsum, isotropy = true)
    basis = SCEBasis(crystal, sg, sb, spec)
    groups = alias_groups(basis)
    # two groups: the 1–2 corner tie under test, and the 1–3 pair's ±c/2 tie (atom 3
    # sits at z = 0.5: images (1.2, 1.2, ±2.5), 3.02 Å, two one-member orbits)
    @test length(groups) == 2
    @test sort([g.atoms for g in groups]) == [[1, 2], [1, 3]]
    g = groups[findfirst(g -> g.atoms == [1, 2], groups)]
    @test g.kind === :proportional
    nmem = [length(salcs(basis)[j].members) for j in g.salcs]
    @test sort(nmem) == [1, 1, 2]
    @test sum(nmem) == 4
    @test g.weights == ones(3)                      # equal per-member norms → exact +1
    @test n_columns(basis) == n_salcs(basis) - 3   # 3 → 1 and 2 → 1 columns
    # every member of the group is a 1–2 image at the corner distance
    for j in g.salcs, m in salcs(basis)[j].members
        @test m.atoms == [1, 2]
        @test norm(_pair_vector(crystal, m)) ≈ norm([2.0, 2.0, 1.25]) rtol = 1e-12
    end
    rng = MersenneTwister(7)
    jtrue = zeros(n_salcs(basis))
    Jset = [0.9, 0.2, 0.1]                          # per-bond mean ≠ per-orbit mean
    for (k, j) in enumerate(g.salcs)
        jtrue[j] = Jset[k]
    end
    gen = SCEPredictor(basis, 0.0, jtrue)
    cfgs = [randcfg(rng, 3) for _ = 1:40]
    E = predict_energy(gen, cfgs)
    T = [predict_torque(gen, c) for c in cfgs]
    m = SCEPredictor(fit(SCEFit, SCEDataset(basis, cfgs, E, T), OLS(); torque_weight = 0.5))
    expected = sum(nmem .* Jset) / sum(nmem)        # the bond-weighted mean
    for j in g.salcs
        @test m.jphi[j] ≈ expected rtol = 1e-10
        @test m.split[j] === :convention
    end
    @test !(expected ≈ sum(Jset) / 3)               # the per-orbit mean is a different number
    @test predict_energy(m, cfgs) ≈ E rtol = 1e-12
    # the hard cap on alias_rtol, as for tie_tol
    @test_throws ArgumentError SCEBasis(crystal, sg, sb, spec; alias_rtol = 1e-2)
    @test_throws ArgumentError SCEBasis(crystal, sg, sb, spec; alias_rtol = 0.5)
end
