# The extended-XYZ training-data container (src/io/extxyz.jl) and the moment
# channel's axis-consistency gates (check_moment_gates, src/io/dftsource.jl), plus
# the legacy EMBSET pair reader (read_embset_pair, src/io/embset.jl). The headline
# gates: (1) every stored value survives the file bit-exactly (shortest-round-trip
# printing); (2) spin-only vs joint is MEASURED from positions — and this pure-spin
# package refuses a joint file by name instead of flattening it; (3) the axis gates
# fire on generation AND on load — archived constraint axes are re-verified, never
# believed.

using Test
using SCEFitting
using LinearAlgebra
using Random

@testset "extxyz container + moment gates + EMBSET pair" begin
    tmp = mktempdir()
    xt = Crystal(Lattice(3.0 .* [1.0 0 0; 0 1 0; 0 0 1]),
                 [0.0 0.5; 0.0 0.5; 0.0 0.5], [1, 2], ["Fe", "Ge"])
    nat = 2
    rng = MersenneTwister(0x5CE2)
    mkdirs() = (m = randn(rng, 3, nat);
                for a = 1:nat; m[:, a] ./= norm(m[:, a]); end; m)
    function mkdat(i; mode = 4)
        dirs = mkdirs()
        mags = [1.2, 0.1] .+ 0.01 .* rand(rng, nat)
        M = 1.5 .* (mags' .* dirs) .+ 0.001 .* randn(rng, 3, nat)
        SpinDatum(-1.0 - 0.1i, mags' .* dirs, 0.01 .* randn(rng, 3, nat);
                  moments_bare = M, constraint_axes = dirs, constraint_mode = mode)
    end

    @testset "spin-only round trip: stored values survive bit-exactly" begin
        data = [mkdat(i) for i = 1:4]
        f = joinpath(tmp, "spin.extxyz")
        write_extxyz(f, data, xt; field_sign = "vasp:+B", source = "test")
        back = read_extxyz(f)
        @test length(back) == 4
        for i = 1:4
            @test back[i].energy == data[i].energy                    # bitwise
            @test back[i].field == data[i].field                      # bitwise
            @test back[i].moments_bare == data[i].moments_bare        # bitwise
            @test back[i].constraint_axes == data[i].constraint_axes  # bitwise
            @test back[i].constraint_mode == 4
            # torques re-derive from the written mw column (= ‖m‖·(m/‖m‖), an ulp
            # off the original m), so they are exact to rounding, not bitwise
            @test maximum(abs, back[i].torques .- data[i].torques) < 1e-15
            # directions/magmoms are re-derived from the written moment vectors —
            # exact up to the unit normalization, deliberately not bitwise
            @test maximum(abs, back[i].directions .- data[i].directions) < 1e-14
            @test maximum(abs, back[i].magmoms .- data[i].magmoms) < 1e-14
        end
        # with the reference crystal: the stored structure must match it bit for bit
        @test length(read_extxyz(f; reference = xt)) == 4
        # a second round trip is stable: the stored channels stay bitwise (the mw
        # column re-derives from magmoms .* directions, so IT may move by an ulp —
        # which is exactly why the bare channels are stored as their own columns)
        f2 = joinpath(tmp, "spin2.extxyz")
        write_extxyz(f2, back, xt; field_sign = "vasp:+B", source = "test")
        back2 = read_extxyz(f2)
        for i = 1:4
            @test back2[i].energy == back[i].energy
            @test back2[i].field == back[i].field
            @test back2[i].moments_bare == back[i].moments_bare
            @test back2[i].constraint_axes == back[i].constraint_axes
            @test maximum(abs, back2[i].directions .- back[i].directions) < 1e-14
        end
        # the file is the shared dialect: species pos(3) mw(3) bcon(3) mint(3)
        # mconstr(3), config_type=spin-only, units_field=eV/muB
        txt = read(f, String)
        @test occursin("Properties=species:S:1:pos:R:3:mw:R:3:bcon:R:3:mint:R:3:" *
                       "mconstr:R:3", txt)
        @test occursin("config_type=spin-only", txt) && occursin("units_field=eV/muB", txt)
        @test occursin("Lattice=\"3.0 0.0 0.0 0.0 3.0 0.0 0.0 0.0 3.0\"", txt)
    end

    @testset "joint data are refused by name, never flattened" begin
        data = [mkdat(i) for i = 1:2]
        f = joinpath(tmp, "base.extxyz")
        write_extxyz(f, data, xt)
        txt = read(f, String)
        lines = split(txt, "\n")
        # (a) positions that differ across frames (a displaced frame)
        tok = split(lines[7])                     # frame 2, atom 1
        tok[2] = "0.01"
        lines2 = copy(lines); lines2[7] = join(tok, " ")
        disp = joinpath(tmp, "disp.extxyz")
        write(disp, join(lines2, "\n"))
        err = try; read_extxyz(disp); nothing; catch e; e; end
        @test err isa ArgumentError && occursin("SLCE.jl", err.msg)
        @test occursin("positions differ", err.msg)
        # (b) a forces column
        forces = joinpath(tmp, "forces.extxyz")
        write(forces, replace(txt, ":mconstr:R:3" => ":mconstr:R:3:forces:R:3",
                              r"(?m)^(Fe|Ge)(.*)$" => s"\1\2 0.0 0.0 0.0"))
        err = try; read_extxyz(forces); nothing; catch e; e; end
        @test err isa ArgumentError && occursin("forces", err.msg) &&
              occursin("SLCE.jl", err.msg)
        # (c) a config_type=joint claim on a measured spin-only file
        claim = joinpath(tmp, "claim.extxyz")
        write(claim, replace(txt, "config_type=spin-only" => "config_type=joint"))
        @test_throws ArgumentError read_extxyz(claim)
        # (d) positions that differ from the reference crystal
        other = Crystal(Lattice(3.0 .* [1.0 0 0; 0 1 0; 0 0 1]),
                        [0.0 0.5; 0.0 0.5; 0.0 0.4], [1, 2], ["Fe", "Ge"])
        err = try; read_extxyz(f; reference = other); nothing; catch e; e; end
        @test err isa ArgumentError && occursin("reference", err.msg)
    end

    @testset "the reference-geometry gate bands round-off, not displacement" begin
        # This is an INTERCHANGE format, so the file's writer is routinely a different
        # build, and the reference side is a freshly recomputed `vectors * frac` whose
        # last bits depend on the arithmetic path. The distinction the gate exists to
        # make -- spin-only vs displaced (joint) -- lives at >~ 1e-3 A, orders above
        # either. Both scales are asserted so the band cannot drift into either one.
        data = [mkdat(i) for i = 1:2]
        f = joinpath(tmp, "band.extxyz")
        write_extxyz(f, data, xt)
        frac = Matrix(xt.frac_positions)
        frac[3, 2] += eps(1.5) / 3.0            # one ulp of a 1.5 A coordinate
        near = Crystal(xt.lattice, frac, xt.species, xt.species_labels)
        # a genuine last-bit difference, not a no-op the gate never sees
        @test cartesian_positions(near) != cartesian_positions(xt)
        @test length(read_extxyz(f; reference = near)) == 2
        frac[3, 2] = xt.frac_positions[3, 2] + 1e-3 / 3.0     # 1e-3 A displacement
        far = Crystal(xt.lattice, frac, xt.species, xt.species_labels)
        @test_throws ArgumentError read_extxyz(f; reference = far)
    end

    @testset "loud checks: claims never override measurements" begin
        data = [mkdat(i) for i = 1:2]
        f = joinpath(tmp, "claims.extxyz")
        write_extxyz(f, data, xt)
        txt = read(f, String)
        bad = joinpath(tmp, "claims_bad.extxyz")
        # units_field: the "T" header mislabel is refused
        write(bad, replace(txt, "units_field=eV/muB" => "units_field=T"))
        @test_throws ArgumentError read_extxyz(bad)
        # mconstr columns without a constraint_mode key
        write(bad, replace(txt, " constraint_mode=4" => ""))
        @test_throws ArgumentError read_extxyz(bad)
        # truncated file
        write(bad, join(split(txt, "\n")[1:3], "\n"))
        @test_throws ArgumentError read_extxyz(bad)
        # a bare token in the info line, an unterminated quote, a bad Properties
        write(bad, replace(txt, "pbc=\"T T T\"" => "pbc=\"T T T\" stray"))
        @test_throws ArgumentError read_extxyz(bad)
        write(bad, replace(txt, "pbc=\"T T T\"" => "pbc=\"T T T"))
        @test_throws ArgumentError read_extxyz(bad)
        write(bad, replace(txt, "species:S:1:pos:R:3" => "pos:R:3:species:S:1"))
        @test_throws ArgumentError read_extxyz(bad)
        # provenance keys this package has no field for are accepted and ignored
        write(bad, replace(txt, "config_type=spin-only" =>
                                "config_type=spin-only setup_id=s1 soc=false"))
        @test length(read_extxyz(bad)) == 2
        # a missing bcon column reads as a zero field (no constraint, zero torques)
        nob = joinpath(tmp, "nobcon.extxyz")
        rx = r"(?m)^((?:Fe|Ge)(?: \S+){6})(?: \S+){3}((?: \S+){6})$"
        write(nob, replace(txt, ":bcon:R:3" => "", rx => s"\1\2"))
        nb = read_extxyz(nob)
        @test all(iszero, nb[1].field) && all(iszero, nb[1].torques)
        @test nb[1].moments_bare == data[1].moments_bare
        # writer refuses mixed channel presence (one file = one observation set)
        nom = SpinDatum(-1.0, mkdirs(), zeros(3, nat))
        @test_throws ArgumentError write_extxyz(joinpath(tmp, "mix.extxyz"),
                                                [data[1], nom], xt)
        # ... and a mixed constraint_mode
        @test_throws ArgumentError write_extxyz(joinpath(tmp, "mix2.extxyz"),
                                                [data[1], mkdat(3; mode = 1)], xt)
        # ... and an atom-count mismatch with the crystal
        @test_throws ArgumentError write_extxyz(joinpath(tmp, "mix3.extxyz"),
                                                [SpinDatum(-1.0, randn(rng, 3, 3),
                                                           zeros(3, 3))], xt)
    end

    @testset "axis gates: generation-time and load-time, loud" begin
        zhat = repeat([0.0, 0.0, 1.0], 1, nat)
        mags = [1.0, 1.0]
        M = 0.9 .* zhat
        B0 = zeros(3, nat)
        mk(dirs; kw...) = SpinDatum(-1.0, mags' .* dirs, B0; kw...)
        # consistent mode-1 datum passes
        d_ok = mk(zhat; moments_bare = M, constraint_axes = zhat, constraint_mode = 1)
        @test SCEFitting.check_moment_gates([d_ok]) === nothing
        # a NEGATIVE readout with the converged direction on the negative side is
        # consistent (the sign gate tests agreement, not positivity)
        d_neg = mk(-zhat; moments_bare = -M, constraint_axes = zhat, constraint_mode = 1)
        @test SCEFitting.check_moment_gates([d_neg]) === nothing
        # sign flip: y > 0 but the converged direction points the other way
        d_flip = mk(-zhat; moments_bare = M, constraint_axes = zhat, constraint_mode = 1)
        err = try; SCEFitting.check_moment_gates([d_flip]); nothing; catch e; e; end
        @test err isa ArgumentError && occursin("sign-consistency", err.msg)
        # rows below the gate floor carry no sign information and are not gated
        d_small = mk(-zhat; moments_bare = 1e-4 .* zhat, constraint_axes = zhat,
                     constraint_mode = 1)
        @test SCEFitting.check_moment_gates([d_small]) === nothing
        # sign_gate_min must be positive
        @test_throws ArgumentError SCEFitting.check_moment_gates([d_ok];
                                                                sign_gate_min = 0.0)
        # nothing to gate (no mode, or mode 1 without bare moments): silent pass
        @test SCEFitting.check_moment_gates([mk(zhat)]) === nothing
        @test SCEFitting.check_moment_gates([mk(-zhat; constraint_axes = zhat,
                                                 constraint_mode = 1)]) === nothing
        # axis-angle p99 (mode 4): a 10° stale axis is loud at the 5° default
        th = deg2rad(10.0)
        rot = [cos(th) 0.0 sin(th); 0.0 1.0 0.0; -sin(th) 0.0 cos(th)]
        d_ang = mk(rot * zhat; moments_bare = M, constraint_axes = zhat,
                   constraint_mode = 4)
        err = try; SCEFitting.check_moment_gates([d_ang]); nothing; catch e; e; end
        @test err isa ArgumentError && occursin("axis-angle", err.msg)
        @test SCEFitting.check_moment_gates([d_ang]; axis_angle_p99_max = 15.0) ===
              nothing
        # mode 4 without bare moments: magmoms stand in for |y|, the angle still gates
        d_ang4 = mk(rot * zhat; constraint_axes = zhat, constraint_mode = 4)
        @test_throws ArgumentError SCEFitting.check_moment_gates([d_ang4])
        # the gate is a PERCENTILE: one collapse-row outlier among many clean rows
        # passes (the measured FeGe τ0.5 situation), a systematic offset does not
        many = [mk(zhat; moments_bare = M, constraint_axes = zhat, constraint_mode = 4)
                for _ = 1:100]
        @test SCEFitting.check_moment_gates(vcat(many, [d_ang])) === nothing
        @test_throws ArgumentError SCEFitting.check_moment_gates(
            vcat(many[1:10], [d_ang for _ = 1:10]))
        # the writer runs the gates: a violating set never becomes a file
        @test_throws ArgumentError write_extxyz(joinpath(tmp, "gate.extxyz"),
                                                [d_flip], xt)
        # ... and the reader re-runs them: a file corrupted after generation is
        # caught at load. Two corruptions, one gauge: per-atom columns are
        # species pos(3) mw(3) bcon(3) mint(3) mconstr(3).
        f = joinpath(tmp, "gated.extxyz")
        write_extxyz(f, [d_ok], xt)
        lines = split(read(f, String), "\n")
        # (a) flipping the WHOLE mconstr axis flips y with it — a mode-1 axis sign
        # is a gauge, and the gate must NOT fire on it (that is the physics: the
        # transverse-penalty constraint prescribes an axis, not an orientation)
        toka = split(lines[3])
        @test toka[16] == "1.0"                 # mconstr z of atom 1 is where we think
        toka[16] = "-1.0"
        gauge = joinpath(tmp, "gated_gauge.extxyz")
        write(gauge, join([lines[1], lines[2], join(toka, " "), lines[4]], "\n"))
        @test length(read_extxyz(gauge)) == 1
        # (b) flipping the BARE MOMENT alone breaks the sign consistency between
        # the converged direction and the readout — loud at load
        tokb = split(lines[3])
        @test tokb[13] == "0.9"                 # mint z of atom 1
        tokb[13] = "-0.9"
        bad = joinpath(tmp, "gated_bad.extxyz")
        write(bad, join([lines[1], lines[2], join(tokb, " "), lines[4]], "\n"))
        err = try; read_extxyz(bad); nothing; catch e; e; end
        @test err isa ArgumentError && occursin("sign-consistency", err.msg)
    end

    @testset "EMBSET pair: loud sibling checks, energies uncompared" begin
        e1 = joinpath(tmp, "EMBSET")
        e2 = joinpath(tmp, "EMBSET_mint")
        wemb(p, mz, e0; nconf = 3, bfy = 0.02) = open(p, "w") do io
            for c = 1:nconf
                println(io, e0 - c)
                for a = 1:2
                    println(io, "$a 0.0 0.0 $mz 0.01 $bfy 0.0")
                end
            end
        end
        wemb(e1, 1.2, -1.0)
        wemb(e2, 1.1, -2.0)                     # energies differ ON PURPOSE
        pd = read_embset_pair(e1, e2; constraint_mode = 4)
        @test length(pd) == 3
        @test pd[1].magmoms[1] == 1.2           # MW file: decomposition source
        @test pd[1].moments_bare[3, 1] == 1.1   # mint file: bare target
        @test pd[1].energy == -2.0              # the MW-file energy is used
        @test pd[1].constraint_mode == 4
        # the plain reader is unchanged by the refactor that made the pair possible
        single = read_embset(e1)
        @test [d.energy for d in single] == [d.energy for d in pd]
        @test all(single[c].directions == pd[c].directions for c = 1:3)
        @test all(single[c].torques == pd[c].torques for c = 1:3)
        # config-count mismatch (the EMBSET_mint_100 provenance-bug class)
        wemb(e2, 1.1, -2.0; nconf = 2)
        @test_throws ArgumentError read_embset_pair(e1, e2; constraint_mode = 4)
        # field-block mismatch: not siblings
        wemb(e2, 1.1, -2.0; bfy = 0.03)
        err = try; read_embset_pair(e1, e2; constraint_mode = 4); nothing; catch e; e; end
        @test err isa ArgumentError && occursin("field", err.msg)
        # per-config constraint axes attach, and the gates run over them
        wemb(e2, 1.1, -2.0)
        ax = repeat([0.0, 0.0, 1.0], 1, 2)
        pd1 = read_embset_pair(e1, e2; constraint_mode = 1, constraint_axes = ax)
        @test all(d.constraint_axes == ax for d in pd1)
        axv = [copy(ax) for _ = 1:3]
        pd2 = read_embset_pair(e1, e2; constraint_mode = 1, constraint_axes = axv)
        @test all(d.constraint_axes == ax for d in pd2)
        @test_throws ArgumentError read_embset_pair(e1, e2; constraint_mode = 1,
                                                    constraint_axes = axv[1:2])
        # mode 1 without axes dies at the datum ctor (loud, not a silent fallback)
        @test_throws ArgumentError read_embset_pair(e1, e2; constraint_mode = 1)
        # the gates run over the pair: a flipped bare moment against the axis is loud
        wemb(e2, -1.1, -2.0)
        @test_throws ArgumentError read_embset_pair(e1, e2; constraint_mode = 1,
                                                    constraint_axes = ax)
    end

    @testset "ExtxyzFile source → SCEDataset (E/T paths inert to the trio)" begin
        data = [mkdat(i) for i = 1:4]
        f = joinpath(tmp, "src.extxyz")
        write_extxyz(f, data, xt)
        b = SCEBasis(xt, BasisSpec(; nbody = 2, lmax = [1, 1], cutoff = 3.0))
        ds = SCEDataset(b, ExtxyzFile(f))
        @test length(ds.configs) == 4
        # identical to the direct-datum path
        ds2 = SCEDataset(b, data)
        @test maximum(abs, ds.X_E .- ds2.X_E) < 1e-14
        @test ds.y_E == ds2.y_E
        @test maximum(abs, ds.y_T .- ds2.y_T) < 1e-15        # mw re-derivation, ulps
    end

    @testset "saboteur closures: gate boundaries, frame checks, writer round-trips" begin
        zhat = [0.0, 0.0, 1.0]
        one(dirs; kw...) = SpinDatum(-1.0, reshape(dirs, 3, 1), zeros(3, 1); kw...)
        # mode 4 with an ANTIPARALLEL axis: the angle is 180° and must fire — there
        # is no gauge freedom in the direction-pinning mode (a sign(y) leak into
        # mode 4 would turn this into 0°)
        d_anti = one(zhat; moments_bare = reshape(0.9 .* zhat, 3, 1),
                     constraint_axes = reshape(-zhat, 3, 1), constraint_mode = 4)
        err = try; SCEFitting.check_moment_gates([d_anti]); nothing; catch e; e; end
        @test err isa ArgumentError && occursin("axis-angle", err.msg) &&
              occursin("180", err.msg)
        # ... whereas in mode 1 the same geometry is the gauge: y = ê_c·M < 0 and
        # ê_MW·ê_c < 0 agree, the signed axis is +ê_MW, angle 0
        d_g = one(zhat; moments_bare = reshape(0.9 .* zhat, 3, 1),
                  constraint_axes = reshape(-zhat, 3, 1), constraint_mode = 1)
        @test SCEFitting.check_moment_gates([d_g]) === nothing
        # the percentile is ⌊n/100⌋ tolerated rows, with ceil: 100 rows tolerate 1
        # stale row and refuse 2; 101 rows tolerate 1 and refuse 2 (floor/round would
        # tolerate 2 of 101)
        th = deg2rad(10.0)
        rot = [cos(th) 0.0 sin(th); 0.0 1.0 0.0; -sin(th) 0.0 cos(th)]
        clean = one(zhat; constraint_axes = reshape(zhat, 3, 1), constraint_mode = 4)
        stale = one(rot * zhat; constraint_axes = reshape(zhat, 3, 1),
                    constraint_mode = 4)
        rows(nc, ns) = vcat([clean for _ = 1:nc], [stale for _ = 1:ns])
        @test SCEFitting.check_moment_gates(rows(99, 1)) === nothing       # n = 100
        @test_throws ArgumentError SCEFitting.check_moment_gates(rows(98, 2))
        @test SCEFitting.check_moment_gates(rows(100, 1)) === nothing      # n = 101
        @test_throws ArgumentError SCEFitting.check_moment_gates(rows(99, 2))
        @test_throws ArgumentError SCEFitting.check_moment_gates(rows(98, 1))  # n = 99
        # a zero axis column is "no axis" at the gate too: mode 4 without bare
        # moments would otherwise read |y| = magmoms and a 90° angle on that atom
        d_z = SpinDatum(-1.0, hcat(zhat, zhat), zeros(3, 2);
                        constraint_axes = hcat(zhat, zeros(3)), constraint_mode = 4)
        @test SCEFitting.check_moment_gates([d_z]) === nothing
        # the knobs forward through every wrapper
        f_stale = joinpath(tmp, "stale.extxyz")
        xt1 = Crystal(Lattice(3.0 .* [1.0 0 0; 0 1 0; 0 0 1]), zeros(3, 1), [1], ["Fe"])
        @test_throws ArgumentError write_extxyz(f_stale, [stale, clean], xt1)
        write_extxyz(f_stale, [stale], xt1; axis_angle_p99_max = 15.0)
        @test_throws ArgumentError read_extxyz(f_stale)
        @test length(read_extxyz(f_stale; axis_angle_p99_max = 15.0)) == 1
        @test_throws ArgumentError read_configs(ExtxyzFile(f_stale))
        @test length(read_configs(ExtxyzFile(f_stale; axis_angle_p99_max = 15.0))) == 1
        other1 = Crystal(Lattice(3.0 .* [1.0 0 0; 0 1 0; 0 0 1]), fill(0.1, 3, 1), [1],
                         ["Fe"])
        @test_throws ArgumentError read_configs(ExtxyzFile(f_stale; reference = other1,
                                                           axis_angle_p99_max = 15.0))
        @test length(read_configs(ExtxyzFile(f_stale; reference = xt1,
                                             axis_angle_p99_max = 15.0))) == 1

        # frame-2 violations, each refused by name (one file = one structure family,
        # one observation set, one constraint scheme, one claim)
        data = [mkdat(i) for i = 1:2]
        f = joinpath(tmp, "frames.extxyz")
        write_extxyz(f, data, xt)
        lines = split(read(f, String), "\n")
        # frame 2 = lines 5..8 (count, info, Fe, Ge)
        function frame2(edit)
            L = copy(lines); edit(L)
            q = joinpath(tmp, "frame2.extxyz"); write(q, join(L, "\n")); q
        end
        @test_throws ArgumentError read_extxyz(frame2(L -> (L[7], L[8]) =
            (replace(L[7], r"^Fe" => "Ge"), replace(L[8], r"^Ge" => "Fe"))))  # species
        @test_throws ArgumentError read_extxyz(frame2(L -> L[6] =
            replace(L[6], "Lattice=\"3.0" => "Lattice=\"3.1")))                # lattice
        @test_throws ArgumentError read_extxyz(frame2(L -> L[6] =
            replace(L[6], "constraint_mode=4" => "constraint_mode=1")))         # mode
        @test_throws ArgumentError read_extxyz(frame2(L -> L[6] =
            replace(L[6], "config_type=spin-only" => "config_type=joint")))     # claim
        @test_throws ArgumentError read_extxyz(frame2(L -> begin                # mint gone
            L[6] = replace(L[6], ":mint:R:3" => "")
            rx = r"^((?:Fe|Ge)(?: \S+){9})(?: \S+){3}((?: \S+){3})$"
            L[7] = replace(L[7], rx => s"\1\2"); L[8] = replace(L[8], rx => s"\1\2")
        end))
        # info-line and Properties duplicates are refused, not arbitrated
        @test_throws ArgumentError read_extxyz(frame2(L -> L[2] = L[2] * " energy=0.0"))
        @test_throws ArgumentError read_extxyz(frame2(L -> begin
            L[2] = replace(L[2], ":mconstr:R:3" => ":mconstr:R:3:mw:R:3")
            L[3] = L[3] * " 0.0 0.0 0.0"; L[4] = L[4] * " 0.0 0.0 0.0"
        end))
        @test_throws ArgumentError read_extxyz(frame2(L -> begin
            L[2] = replace(L[2], ":mconstr:R:3" => ":mconstr:R:3:tag:S:1")
            L[3] = L[3] * " x"; L[4] = L[4] * " y"
        end))

        # a triclinic cell with asymmetric positions: the lattice string is
        # column-major (columns = lattice vectors, the ASE / SLCE convention) and
        # the reference check holds — a transposed writer or reader would fail here
        # (by whole angstroms, far outside the round-off band) where the cubic
        # fixtures cannot tell
        Atri = [3.0 0.3 0.2; 0.0 3.4 0.5; 0.0 0.0 3.9]
        xtri = Crystal(Lattice(Atri), [0.1 0.6; 0.2 0.7; 0.3 0.9], [1, 2], ["Fe", "Ge"])
        ftri = joinpath(tmp, "tri.extxyz")
        write_extxyz(ftri, data, xtri)
        ttxt = read(ftri, String)
        @test occursin("Lattice=\"3.0 0.0 0.0 0.3 3.4 0.0 0.2 0.5 3.9\"", ttxt)
        @test length(read_extxyz(ftri; reference = xtri)) == 2
        @test_throws ArgumentError read_extxyz(ftri; reference = xt)
        posline = split(split(ttxt, "\n")[3])
        @test parse.(Float64, posline[2:4]) == vec(Matrix(cartesian_positions(xtri))[:, 1])

        # a datum without the trio writes a shorter header and reads back
        plain = [SpinDatum(-1.0 - 0.1i, 1.2 .* mkdirs(), 0.01 .* randn(rng, 3, nat))
                 for i = 1:2]
        fp = joinpath(tmp, "plain.extxyz")
        write_extxyz(fp, plain, xt)
        @test occursin("Properties=species:S:1:pos:R:3:mw:R:3:bcon:R:3 ", read(fp, String))
        bp = read_extxyz(fp)
        @test all(d.moments_bare === nothing && d.constraint_mode === nothing for d in bp)
        @test bp[2].field == plain[2].field

        # free-text values: whitespace is quoted so the file reloads; a quote or
        # a line break inside a value is refused up front
        fq = joinpath(tmp, "quoted.extxyz")
        write_extxyz(fq, data, xt; source = "C:\\Program Files\\x", field_sign = "vasp +B",
                     comment = "two words")
        @test length(read_extxyz(fq)) == 2
        @test occursin("source=\"C:\\Program Files\\x\"", read(fq, String))
        @test_throws ArgumentError write_extxyz(fq, data, xt; comment = "say \"hi\"")
        @test_throws ArgumentError write_extxyz(fq, data, xt; source = "a\nb")

        # EMBSET pair: per-config axes are per config (a stale config 2 is named),
        # and the same moments in both files are not a pair
        e1 = joinpath(tmp, "EMBSET_b"); e2 = joinpath(tmp, "EMBSET_mint_b")
        wemb2(p, mz) = open(p, "w") do io
            for c = 1:3
                println(io, -1.0 - c)
                for a = 1:2
                    println(io, "$a 0.0 0.0 $mz 0.01 0.02 0.0")
                end
            end
        end
        wemb2(e1, 1.2); wemb2(e2, 1.1)
        axs = [hcat(zhat, zhat) for _ = 1:3]
        axs[2] = hcat(-zhat, zhat)               # antiparallel on config 2 (mode 4)
        err = try
            read_embset_pair(e1, e2; constraint_mode = 4, constraint_axes = axs)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError && occursin("axis-angle", err.msg)
        pdv = read_embset_pair(e1, e2; constraint_mode = 4,
                               constraint_axes = axs, axis_angle_p99_max = 180.0)
        @test pdv[2].constraint_axes == axs[2] && pdv[1].constraint_axes == axs[1]
        @test_throws ArgumentError read_embset_pair(e1, e1; constraint_mode = 4)
        cp(e1, joinpath(tmp, "EMBSET_copy"); force = true)
        @test_throws ArgumentError read_embset_pair(e1, joinpath(tmp, "EMBSET_copy");
                                                    constraint_mode = 4)
    end
end
