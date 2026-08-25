# Extended-XYZ (ASE dialect) training-data container — the canonical on-disk format
# for constrained-noncollinear spin training sets, shared with SLCE.jl (same dialect,
# so a spin-only file moves between the two packages unchanged). One frame per
# configuration: a `Lattice`/`Properties` info line plus one line per atom. The
# design decisions (recorded in the adiabatic-moment design record, D8):
#
#   * The structure (lattice + positions) is ALWAYS stored, even for spin-only data —
#     self-containment structurally removes the "which POSCAR pairs with which
#     EMBSET" provenance-bug class. The redundancy costs a few MB.
#   * spin-only vs joint is decided FROM THE DATA (positions agreeing across
#     frames → spin-only), never from a flag; a `config_type` claim is allowed but a
#     mismatch with the measured answer is a loud error. THIS package is pure spin:
#     a joint file (displaced frames, `forces` columns, or a `joint` claim) is
#     refused by name — SLCE.jl reads it — never silently flattened to its spins.
#   * Per-atom columns exist only when observed: `species:S:1:pos:R:3:mw:R:3`
#     [`:bcon:R:3`][`:mint:R:3`][`:mconstr:R:3`], 1:1 with `SpinDatum`'s channels
#     (`mw` = smoothed moment vectors `magmoms .* directions`; `bcon` = constraining
#     field; `mint` = bare moments; `mconstr` = constraint axes). `SpinDatum` always
#     carries a field, so the writer always emits `bcon`; a file without it reads as
#     a zero field (no constraint, zero torques).
#   * A generated file is self-contained INCLUDING the constraint axes: the axis
#     gates ([`check_moment_gates`](@ref)) run at generation and again at every load,
#     so archived axes are re-verified, never believed.
#   * DFT-code vocabulary stays in the generator (SCETools' `oszicar_to_extxyz`);
#     this reader/writer is code-neutral. Upstream's provenance keys (`setup_id`,
#     `soc`, `reference_id`) are accepted and ignored: `SpinDatum` has no
#     provenance field.

"""
    ExtxyzFile(path; reference = nothing, zero_moment_atol = 1e-10,
               sign_gate_min = 5e-3, axis_angle_p99_max = 5.0)

An [`AbstractDFTSource`](@ref) for an extended-XYZ training-set file, so it drops
straight into the pipeline:

```julia
dataset = SCEDataset(basis, ExtxyzFile("train.extxyz"))
```

`read_configs` on it is [`read_extxyz`](@ref)`(path; ...)` — see there for the format
and the keyword arguments.
"""
struct ExtxyzFile <: AbstractDFTSource
    path::String
    reference::Union{Crystal,Nothing}
    zero_moment_atol::Float64
    sign_gate_min::Float64
    axis_angle_p99_max::Float64
end

function ExtxyzFile(path::AbstractString; reference::Union{Crystal,Nothing} = nothing,
                    zero_moment_atol::Real = 1e-10, sign_gate_min::Real = 5e-3,
                    axis_angle_p99_max::Real = 5.0)::ExtxyzFile
    return ExtxyzFile(String(path), reference, Float64(zero_moment_atol),
                      Float64(sign_gate_min), Float64(axis_angle_p99_max))
end

read_configs(src::ExtxyzFile)::Vector{SpinDatum} =
    read_extxyz(src.path; reference = src.reference,
                zero_moment_atol = src.zero_moment_atol,
                sign_gate_min = src.sign_gate_min,
                axis_angle_p99_max = src.axis_angle_p99_max)

# ── info-line lexing ───────────────────────────────────────────────────────────────

# Split an extxyz info line into `key = value` pairs; values may be double-quoted
# (quotes stripped, spaces preserved). Insertion order is irrelevant to the reader.
function _xyz_info(line::AbstractString, path::AbstractString,
                   frame::Int)::Dict{String,String}
    out = Dict{String,String}()
    i = firstindex(line)
    n = lastindex(line)
    while i <= n
        c = line[i]
        if isspace(c)
            i = nextind(line, i)
            continue
        end
        j = i                      # key runs to '='
        while j <= n && line[j] != '='
            isspace(line[j]) &&
                throw(ArgumentError("extxyz $path frame $frame: bare token " *
                                    "\"$(line[i:prevind(line, j)])\" in the info " *
                                    "line (expected key=value pairs)"))
            j = nextind(line, j)
        end
        j <= n || throw(ArgumentError("extxyz $path frame $frame: key without " *
                                      "'=' at the end of the info line"))
        key = line[i:prevind(line, j)]
        isempty(key) && throw(ArgumentError("extxyz $path frame $frame: empty key " *
                                            "in the info line"))
        i = nextind(line, j)       # past '='
        if i <= n && line[i] == '"'
            i = nextind(line, i)
            k = i
            while k <= n && line[k] != '"'
                k = nextind(line, k)
            end
            k <= n || throw(ArgumentError("extxyz $path frame $frame: unterminated " *
                                          "quote in the info line (key \"$key\")"))
            val = line[i:prevind(line, k)]
            i = nextind(line, k)
        else
            k = i
            while k <= n && !isspace(line[k])
                k = nextind(line, k)
            end
            val = line[i:prevind(line, k)]
            i = k
        end
        # a repeated key would silently let the last value win (Dict semantics);
        # a file that says two things is refused, not arbitrated
        haskey(out, key) &&
            throw(ArgumentError("extxyz $path frame $frame: key \"$key\" appears " *
                                "twice in the info line"))
        out[key] = val
    end
    return out
end

# Parse `Properties=species:S:1:pos:R:3:…` into an ordered `(name, type, ncols)` list.
function _xyz_properties(spec::AbstractString, path::AbstractString,
                         frame::Int)::Vector{Tuple{String,Char,Int}}
    tok = split(spec, ':')
    length(tok) % 3 == 0 && !isempty(tok) ||
        throw(ArgumentError("extxyz $path frame $frame: malformed Properties " *
                            "\"$spec\" (need name:type:count triples)"))
    props = Tuple{String,Char,Int}[]
    for t = 1:3:length(tok)
        name = String(tok[t])
        ty = tok[t + 1]
        length(ty) == 1 && ty[1] in ('S', 'R', 'I', 'L') ||
            throw(ArgumentError("extxyz $path frame $frame: unsupported column " *
                                "type \"$ty\" for \"$name\""))
        nc = tryparse(Int, tok[t + 2])
        (nc === nothing || nc < 1) &&
            throw(ArgumentError("extxyz $path frame $frame: bad column count " *
                                "\"$(tok[t + 2])\" for \"$name\""))
        any(p -> p[1] == name, props) &&
            throw(ArgumentError("extxyz $path frame $frame: property \"$name\" " *
                                "appears twice in Properties (the second block " *
                                "would silently overwrite the first)"))
        # the only string column this reader knows is the leading species; a second
        # one would be read INTO species (last wins) — refused instead
        ty[1] == 'S' && !isempty(props) &&
            throw(ArgumentError("extxyz $path frame $frame: string column " *
                                "\"$name\" after species is not supported"))
        push!(props, (name, ty[1], nc))
    end
    return props
end

_xyz_number(s::AbstractString, what::String, path::AbstractString)::Float64 = begin
    v = tryparse(Float64, s)
    (v === nothing || !isfinite(v)) &&
        throw(ArgumentError("extxyz $path: cannot parse $what from \"$s\"" *
                            (v === nothing ? "" : " (non-finite)")))
    v
end

# One parsed frame: species labels, per-column 3 × nat blocks, and the info dict.
struct _XyzFrame
    species::Vector{String}
    cols::Dict{String,Matrix{Float64}}
    info::Dict{String,String}
end

function _read_xyz_frames(path::AbstractString)::Vector{_XyzFrame}
    isfile(path) || throw(ArgumentError("no such extxyz file: $path"))
    lines = readlines(path)
    frames = _XyzFrame[]
    i = 1
    while i <= length(lines)
        if isempty(strip(lines[i]))
            i += 1
            continue
        end
        frame = length(frames) + 1
        nat = tryparse(Int, strip(lines[i]))
        (nat === nothing || nat < 1) &&
            throw(ArgumentError("extxyz $path frame $frame: expected an atom count " *
                                "line, got \"$(lines[i])\""))
        i + 1 <= length(lines) ||
            throw(ArgumentError("extxyz $path frame $frame: missing info line"))
        info = _xyz_info(lines[i + 1], path, frame)
        haskey(info, "Properties") ||
            throw(ArgumentError("extxyz $path frame $frame: no Properties key"))
        props = _xyz_properties(info["Properties"], path, frame)
        props[1][1] == "species" && props[1][2] == 'S' && props[1][3] == 1 ||
            throw(ArgumentError("extxyz $path frame $frame: the first property " *
                                "must be species:S:1"))
        ntok = sum(p[3] for p in props)
        i + 1 + nat <= length(lines) ||
            throw(ArgumentError("extxyz $path frame $frame: truncated — $nat atom " *
                                "lines declared, file ends early"))
        species = Vector{String}(undef, nat)
        cols = Dict{String,Matrix{Float64}}(p[1] => Matrix{Float64}(undef, p[3], nat)
                                            for p in props if p[2] != 'S')
        for a = 1:nat
            tok = split(lines[i + 1 + a])
            length(tok) == ntok ||
                throw(ArgumentError("extxyz $path frame $frame atom $a: expected " *
                                    "$ntok columns, got $(length(tok))"))
            off = 0
            for (name, ty, nc) in props
                if ty == 'S'
                    species[a] = String(tok[off + 1])
                else
                    m = cols[name]
                    for k = 1:nc
                        m[k, a] = _xyz_number(tok[off + k],
                                              "frame $frame atom $a column $name",
                                              path)
                    end
                end
                off += nc
            end
        end
        push!(frames, _XyzFrame(species, cols, info))
        i += 2 + nat
    end
    isempty(frames) && throw(ArgumentError("extxyz $path: no frames"))
    return frames
end

_xyz_lattice(info::Dict{String,String}, path::AbstractString, frame::Int) = begin
    haskey(info, "Lattice") ||
        throw(ArgumentError("extxyz $path frame $frame: no Lattice key — the " *
                            "structure is always stored (self-containment rule)"))
    v = split(info["Lattice"])
    length(v) == 9 ||
        throw(ArgumentError("extxyz $path frame $frame: Lattice needs 9 numbers"))
    A = Matrix{Float64}(undef, 3, 3)
    for c = 1:3, r = 1:3
        A[r, c] = _xyz_number(v[3 * (c - 1) + r], "frame $frame Lattice", path)
    end
    A                                       # columns = lattice vectors
end

# ── reader ─────────────────────────────────────────────────────────────────────────

# Absolute band (Å) for every geometry comparison a reader makes: frame against frame,
# and file against a reference `Crystal`. This is an INTERCHANGE format, so none of the
# three quantities involved is exact — the writer is routinely a different build, the
# reference cartesians are a recomputed `vectors * frac`, and the file carries only as
# many decimals as its writer printed. That last one sets the scale: ASE's extxyz writer
# defaults to `%16.8f`, a text grid of 1e-8 Å whose rounding error reaches 5e-9 Å, so a
# band at 1e-8 Å would leave a factor of two and a producer printing six decimals would
# reproduce exactly the false verdict the band exists to prevent. 1e-6 Å clears an
# eight-decimal grid by ~200x and still sits three orders below the ≳ 1e-3 Å displacement
# the spin-only / joint distinction is about — at 1e-6 Å a Φ ~ 10 eV/Å² force is
# ~1e-5 eV/Å, far under DFT noise. It is therefore also a FLOOR: a displacement smaller
# than this is not representable through a file, and reads as sitting at the reference.
const _REF_GEOM_ATOL = 1e-6


"""
    read_extxyz(path; reference = nothing, zero_moment_atol = 1e-10,
                sign_gate_min = 5e-3, axis_angle_p99_max = 5.0) -> Vector{SpinDatum}

Read an extended-XYZ training-set file (the format [`write_extxyz`](@ref) emits; ASE
dialect, shared with SLCE.jl) into [`SpinDatum`](@ref)s. Per-atom columns map 1:1 onto
the datum's channels: `mw` → `magmoms .* directions` (the smoothed moment
decomposition), `bcon` → `field` (absent ⇒ zero field, zero torques), `mint` →
`moments_bare`, `mconstr` → `constraint_axes`. Info keys consumed: `energy` (required,
eV), `constraint_mode` (`1`/`4`, required when `mconstr` columns are present),
`units_field` (must be `eV/muB` when present), `config_type` (an optional claim),
`pbc` (must be fully periodic). Upstream's provenance keys (`setup_id`, `soc`,
`reference_id`) are accepted and ignored — a `SpinDatum` carries no provenance.

**spin-only vs joint is measured, not read**, and this package holds spin-only data
only: every frame's positions must agree, a `forces` column must be absent, and a
`config_type=joint` claim is refused — all by name, pointing at SLCE.jl, which reads
joint files. A file is never silently flattened to its spins.

Every geometry comparison — frame against frame, and, with a `reference::Crystal`,
the stored lattice and positions against it — holds to an absolute **1e-6 Å** band,
not exactly. This is an interchange format: the writer is routinely a different
build printing a finite number of decimals, and the reference cartesians are a
recomputed `vectors * frac`, so an exact test would answer a question about the
physics with a question about arithmetic. The band sits three orders below the
≳ 1e-3 Å displacement the refusal is about, and is therefore also a floor — a
structure displaced by less than 1e-6 Å reads as sitting at the reference. Species
must match exactly. Without a reference the file is taken on its own terms — which for a
**single-frame** file means a displaced structure cannot be told from a reference
one (there is no second frame to differ from): pass `reference` whenever that
distinction matters.

Every load re-runs the axis-consistency gates ([`check_moment_gates`](@ref)) — the
archived constraint axes are re-verified against the converged moment directions,
never believed.
"""
function read_extxyz(path::AbstractString;
                     reference::Union{Crystal,Nothing} = nothing,
                     zero_moment_atol::Real = 1e-10,
                     sign_gate_min::Real = 5e-3,
                     axis_angle_p99_max::Real = 5.0)::Vector{SpinDatum}
    frames = _read_xyz_frames(path)
    nat = length(frames[1].species)
    joint_msg = " — displaced (joint spin–lattice) data are not representable in " *
                "this pure-spin package; SLCE.jl reads them. Nothing is flattened " *
                "to its spins silently"

    # cross-frame consistency: one file = one structure family
    A1 = _xyz_lattice(frames[1].info, path, 1)
    for (f, fr) in enumerate(frames)
        length(fr.species) == nat ||
            throw(ArgumentError("extxyz $path frame $f: $(length(fr.species)) " *
                                "atoms, frame 1 has $nat"))
        fr.species == frames[1].species ||
            throw(ArgumentError("extxyz $path frame $f: species differ from frame 1"))
        _xyz_lattice(fr.info, path, f) == A1 ||
            throw(ArgumentError("extxyz $path frame $f: Lattice differs from " *
                                "frame 1 (varying cells are not supported here)"))
        haskey(fr.cols, "pos") && size(fr.cols["pos"], 1) == 3 ||
            throw(ArgumentError("extxyz $path frame $f: no pos:R:3 columns"))
        haskey(fr.cols, "mw") ||
            throw(ArgumentError("extxyz $path frame $f: no mw:R:3 columns (the " *
                                "smoothed moment channel is required)"))
        haskey(fr.cols, "forces") &&
            throw(ArgumentError("extxyz $path frame $f: forces columns" * joint_msg))
        for key in ("mw", "bcon", "mint", "mconstr")
            haskey(fr.cols, key) == haskey(frames[1].cols, key) ||
                throw(ArgumentError("extxyz $path frame $f: column \"$key\" " *
                                    "presence differs from frame 1 (one file = " *
                                    "one observation set)"))
            haskey(fr.cols, key) && size(fr.cols[key], 1) != 3 &&
                throw(ArgumentError("extxyz $path frame $f: column \"$key\" must " *
                                    "be :R:3"))
        end
        pbc = get(fr.info, "pbc", "T T T")
        pbc == "T T T" ||
            throw(ArgumentError("extxyz $path frame $f: pbc = \"$pbc\" — only " *
                                "fully periodic cells are supported"))
        uf = get(fr.info, "units_field", "eV/muB")
        uf == "eV/muB" ||
            throw(ArgumentError("extxyz $path frame $f: units_field = \"$uf\" — " *
                                "the constraining field must be in eV/muB (the " *
                                "generator converts; a \"T\" here is the header " *
                                "mislabel this key exists to correct)"))
    end

    # constraint_mode: uniform across frames, required iff mconstr columns exist
    modes = [get(fr.info, "constraint_mode", nothing) for fr in frames]
    allequal(modes) ||
        throw(ArgumentError("extxyz $path: constraint_mode differs across frames " *
                            "(one file = one constraint scheme)"))
    cmode = nothing
    m1 = modes[1]              # bound local: the !== nothing narrowing must be
    if m1 !== nothing          # inference-visible (JET: tryparse(Int, ::Nothing))
        cmode = tryparse(Int, m1)
        cmode === nothing &&
            throw(ArgumentError("extxyz $path: constraint_mode = \"$m1\" " *
                                "is not an integer"))
    end
    haskey(frames[1].cols, "mconstr") && cmode === nothing &&
        throw(ArgumentError("extxyz $path: mconstr columns without a " *
                            "constraint_mode info key — the axis rule is keyed by " *
                            "the mode, declare it (1 = transverse-penalty type, " *
                            "4 = direction-pinning type)"))

    # spin-only vs joint: measured from positions — and only spin-only is admissible
    ref_pos = frames[1].cols["pos"]
    # Banded like every other geometry comparison here (`_REF_GEOM_ATOL`). "One writer,
    # one file" does not make this one exact: a producer that recomputes or re-wraps
    # positions per frame prints a different last decimal for the same structure, and an
    # exact test then answers a question about the physics with a question about
    # formatting. Banding it also keeps the two verdicts consistent — a file must not be
    # spin-only with a `reference` and joint without one.
    frame_dev = maximum(fr -> maximum(abs, fr.cols["pos"] - ref_pos), frames)
    frame_dev <= _REF_GEOM_ATOL ||
        throw(ArgumentError("extxyz $path: positions differ across frames " *
                            "(max |Δ| = $frame_dev Å > $_REF_GEOM_ATOL)" * joint_msg))
    if reference !== nothing
        size(reference.frac_positions, 2) == nat ||
            throw(ArgumentError("extxyz $path: $nat atoms per frame, reference " *
                                "crystal has $(size(reference.frac_positions, 2))"))
        reflab = [reference.species_labels[s] for s in reference.species]
        reflab == frames[1].species ||
            throw(ArgumentError("extxyz $path: species differ from the reference " *
                                "crystal ($(frames[1].species[1]) … vs " *
                                "$(reflab[1]) …)"))
        lat_dev = maximum(abs, Matrix(reference.lattice.vectors) - A1)
        lat_dev <= _REF_GEOM_ATOL ||
            throw(ArgumentError("extxyz $path: Lattice differs from the reference " *
                                "crystal's lattice (max |Δ| = $lat_dev Å > " *
                                "$_REF_GEOM_ATOL)"))
        # Banded, not exact: an exact comparison answered a question about the physics
        # ("these are joint spin–lattice data") with a question about arithmetic — see
        # `_REF_GEOM_ATOL`. The deviation is reported so a user just outside the band
        # can tell a different structure from a band that is too tight.
        refc = Matrix(cartesian_positions(reference))
        pos_dev = maximum(abs, ref_pos - refc)
        pos_dev <= _REF_GEOM_ATOL ||
            throw(ArgumentError("extxyz $path: positions differ from the reference " *
                                "crystal's (max |Δ| = $pos_dev Å > $_REF_GEOM_ATOL)" *
                                joint_msg))
    end

    # config_type: an optional claim, checked against the measurement (spin-only)
    claims = unique(get(fr.info, "config_type", nothing) for fr in frames)
    claim = length(claims) == 1 ? claims[1] :
            throw(ArgumentError("extxyz $path: config_type differs across frames"))
    if claim !== nothing
        claim in ("spin-only", "joint") ||
            throw(ArgumentError("extxyz $path: config_type = \"$claim\" (expected " *
                                "spin-only or joint)"))
        claim == "spin-only" ||
            throw(ArgumentError("extxyz $path: config_type claims \"joint\" but the " *
                                "positions say \"spin-only\" — the flag never " *
                                "overrides the measurement, and a joint set is not " *
                                "readable here in any case (SLCE.jl reads it)"))
    end

    data = Vector{SpinDatum}(undef, length(frames))
    for (f, fr) in enumerate(frames)
        haskey(fr.info, "energy") ||
            throw(ArgumentError("extxyz $path frame $f: no energy key"))
        energy = _xyz_number(fr.info["energy"], "frame $f energy", path)
        field = get(fr.cols, "bcon", nothing)
        data[f] = SpinDatum(energy, fr.cols["mw"],
                            field === nothing ? zeros(3, nat) : field;
                            zero_moment_atol = zero_moment_atol,
                            moments_bare = get(fr.cols, "mint", nothing),
                            constraint_axes = get(fr.cols, "mconstr", nothing),
                            constraint_mode = cmode)
    end
    check_moment_gates(data; sign_gate_min = sign_gate_min,
                       axis_angle_p99_max = axis_angle_p99_max,
                       label = "read_extxyz($path)")
    return data
end

# ── writer ─────────────────────────────────────────────────────────────────────────

"""
    write_extxyz(path, data, crystal; field_sign = nothing, source = nothing,
                 comment = nothing, sign_gate_min = 5e-3,
                 axis_angle_p99_max = 5.0) -> Nothing

Write [`SpinDatum`](@ref)s as an extended-XYZ training-set file (the format
[`read_extxyz`](@ref) — and SLCE.jl — reads back; numbers are printed
shortest-round-trip, so every stored value survives the file bit-exactly — the one
non-bitwise step in a datum round-trip is re-deriving `directions`/`magmoms` from the
written moment vectors, exact up to the unit normalization). The structure comes from
`crystal`, written verbatim into every frame — the structure is **always** stored,
self-containment being the point of the format — with `config_type=spin-only`.

Channel presence must be uniform across `data` (one file = one observation set) and
the `constraint_mode` must be uniform. The constraining field is always written
(`SpinDatum` always carries one). The axis-consistency gates
([`check_moment_gates`](@ref)) run before anything is written — a file that would
fail its own load gate is never produced. `field_sign` / `source` are provenance
strings recorded verbatim in the info line (the generator documents there which
sign convention it normalized from, and what it read).
"""
function write_extxyz(path::AbstractString, data::AbstractVector{SpinDatum},
                      crystal::Crystal;
                      field_sign::Union{Nothing,AbstractString} = nothing,
                      source::Union{Nothing,AbstractString} = nothing,
                      comment::Union{Nothing,AbstractString} = nothing,
                      sign_gate_min::Real = 5e-3,
                      axis_angle_p99_max::Real = 5.0)::Nothing
    isempty(data) && throw(ArgumentError("write_extxyz: no data"))
    nat = size(crystal.frac_positions, 2)
    for (c, d) in enumerate(data)
        size(d.directions, 2) == nat ||
            throw(ArgumentError("write_extxyz: config $c has " *
                                "$(size(d.directions, 2)) atoms, crystal has $nat"))
    end
    for (name, get_ch) in (("mint", d -> d.moments_bare),
                           ("mconstr", d -> d.constraint_axes))
        present = get_ch(data[1]) !== nothing
        all((get_ch(d) !== nothing) == present for d in data) ||
            throw(ArgumentError("write_extxyz: channel \"$name\" is present on " *
                                "some configs and absent on others — one file = " *
                                "one observation set"))
    end
    allequal(d.constraint_mode for d in data) ||
        throw(ArgumentError("write_extxyz: constraint_mode differs across configs " *
                            "(one file = one constraint scheme)"))
    check_moment_gates(data; sign_gate_min = sign_gate_min,
                       axis_angle_p99_max = axis_angle_p99_max,
                       label = "write_extxyz($path)")

    # Free-text info values: the reader's lexer has no escape sequence, so a
    # double quote or a line break inside a value can never round-trip — refused;
    # a value with whitespace is written quoted (the reader accepts quoted values),
    # so the file this writes is always one its own reader loads.
    for (name, v) in (("field_sign", field_sign), ("source", source),
                      ("comment", comment))
        v === nothing && continue
        (occursin('"', v) || occursin('\n', v) || occursin('\r', v)) &&
            throw(ArgumentError("write_extxyz: $name contains a double quote or a " *
                                "line break, which the extxyz info line cannot carry"))
    end
    _q(v) = any(isspace, v) ? "\"" * v * "\"" : v

    cmode = data[1].constraint_mode
    has_mint = data[1].moments_bare !== nothing
    has_mconstr = data[1].constraint_axes !== nothing

    A = Matrix(crystal.lattice.vectors)
    latstr = join(string.(vec(A)), " ")         # columns = lattice vectors
    refc = Matrix(cartesian_positions(crystal))
    labels = [crystal.species_labels[s] for s in crystal.species]
    props = "species:S:1:pos:R:3:mw:R:3:bcon:R:3" * (has_mint ? ":mint:R:3" : "") *
            (has_mconstr ? ":mconstr:R:3" : "")

    open(path, "w") do io
        for d in data
            println(io, nat)
            print(io, "Lattice=\"", latstr, "\" Properties=", props,
                  " energy=", string(d.energy), " pbc=\"T T T\" config_type=spin-only")
            cmode === nothing || print(io, " constraint_mode=", cmode)
            print(io, " units_field=eV/muB")
            field_sign === nothing || print(io, " field_sign=", _q(field_sign))
            source === nothing || print(io, " source=", _q(source))
            comment === nothing || print(io, " comment=\"", comment, "\"")
            println(io)
            for a = 1:nat
                print(io, labels[a])
                for k = 1:3
                    print(io, " ", string(refc[k, a]))
                end
                for k = 1:3
                    print(io, " ", string(d.magmoms[a] * d.directions[k, a]))
                end
                for k = 1:3
                    print(io, " ", string(d.field[k, a]))
                end
                for (flag, m) in ((has_mint, d.moments_bare),
                                  (has_mconstr, d.constraint_axes))
                    flag || continue
                    for k = 1:3
                        print(io, " ", string(m[k, a]))
                    end
                end
                println(io)
            end
        end
    end
    return nothing
end
