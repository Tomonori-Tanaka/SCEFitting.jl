"""
TOML input files.

A human-authored `input.toml` collects an SCE run's *setup* parameters — the
crystal, the interaction (cluster) spec, and the symmetry settings — in one
readable file, so a basis can be built with `SCEBasis("input.toml")` instead of
constructing `Crystal` / `BasisSpec` in Julia. Training data and the choice of
estimator are intentionally kept out of this file (load data and fit in Julia),
mirroring the basis/data separation.

Schema:

    [structure]
    lattice        = [[ax,ay,az],[bx,by,bz],[cx,cy,cz]]  # each entry = one lattice vector (Å)
    positions      = [[fx,fy,fz], ...]                    # fractional, one per atom
    species        = [1,1,2, ...]                         # per-atom species index (1-based)
    species_labels = ["Fe","Pt"]
    pbc            = [true,true,true]                      # optional, default all true

    [interaction]
    nbody    = 2
    cutoff   = 3.0                # Å; `inf` = the full Wigner–Seitz cell; or a table, below
    lmax     = [2,2]              # per species (index order), or a label table, below
    lsum     = 4                  # optional Σl cap per body order (scalar or table)
    isotropy = false              # optional, default false
    images   = "minimum_image"    # optional: "minimum_image" (default) or "all_images"
    tie_tol  = 1e-8               # optional: relative same-distance band (see `SCEBasis`)

    # Label-keyed alternatives ("*" = fallback; pair keys are unordered, resolved by
    # specificity: concrete > "A-*" > "*-*"; body orders outside nbody are errors):
    #
    #     [interaction.lmax]
    #     "*" = 3
    #     B   = 0
    #
    #     [interaction.lsum]        # keys = body orders (bare-integer TOML keys)
    #     1 = 0
    #     2 = 4
    #
    #     [interaction.cutoff]      # body-keyed: scalar per order ...
    #     2 = 8.0
    #     [interaction.cutoff.3]    # ... or a species-pair table per order
    #     "Fe-*" = 6.0
    #     "*-*"  = 8.0
    #
    # A species-pair table directly under [interaction.cutoff] (keys like "Fe-Fe")
    # applies to every body order.

    [symmetry]                                            # optional section
    backend = "spglib"                                    # "none" (default) or "spglib"
    tol     = 1.0e-5                                       # optional, default 1e-5

    [moment]                      # optional section: the pointed site-moment basis
    nbody       = 3           # optional, default 3 (1 to 4). Star members grow as
                              #   C(z, N-1)*N! -- measured 36x from 3 to 4 bodies on a
                              #   54-atom cell at a 3NN star radius
    lmax_mark   = 2           # optional, default 2: cap on the marked site's own ê factor
    lmax_env    = [2]         # per species (index order), or a label table (below)
    sampled     = ["Fe"]      # REQUIRED: species the downstream consumer samples —
                              #   labels, ["*"] (= every species), or per-species booleans
    marked      = ["Fe"]      # optional (default every species): whose moments are expanded
    cutoff_pair = 4.1         # REQUIRED: mark–environment bond radius (Å) of 2-body
                              #   clusters — scalar (`inf` = whole WS cell) or pair table
    cutoff_star = 4.1         # optional, default = cutoff_pair: the N-1 mark-environment
                              #   spokes of a star, nbody >= 3 (env-env edges free).
                              #   Per star order via a body-keyed table (below)
    lsum        = 4           # optional, default uncapped: total spin rank per label
    isotropy    = true        # optional, default true (L_S = 0 only) — NOTE the default
                              #   differs from [interaction].isotropy (false)

    # Label-keyed alternatives (same rules as [interaction]). `cutoff_pair` takes no
    # body-order table — it is the 2-body radius. `cutoff_star` does: keys must cover
    # exactly 3:nbody, and each value is a scalar or a species-pair table.
    #
    #     [moment.lmax_env]
    #     "*" = 2
    #     Rh  = 0
    #
    #     [moment.cutoff_pair]
    #     "Fe-Fe" = 4.1
    #     "*-*"   = 3.0
    #
    #     [moment.cutoff_star]      # per star order: scalar per order ...
    #     3 = 4.1
    #     [moment.cutoff_star.4]    # ... or a species-pair table for one order
    #     "Fe-Fe" = 2.5
    #     "*-*"   = 0.0

The `[moment]` section is read into a `MomentSpec` (every value is handed to the
`MomentSpec` keyword constructor, which does all the validation); unknown keys are
refused (unlike `[interaction]`, which ignores keys it does not read), and so is the
upstream spelling `soc` (use `isotropy`; the polarities are opposite). The moment
basis is always minimum-image and shares `[interaction].tie_tol`. See
`MomentBasis(path)`.
"""

_input_require(d, key, ctx) =
    haskey(d, key) ? d[key] :
    throw(ArgumentError("[$ctx]: required key \"$key\" is missing"))

function _crystal_from_input(d)::Crystal
    lat = _input_require(d, "lattice", "structure")
    length(lat) == 3 ||
        throw(ArgumentError("[structure].lattice must list 3 lattice vectors (got $(length(lat)))"))
    A = MMatrix{3,3,Float64}(undef)
    @inbounds for k = 1:3
        v = lat[k]
        length(v) == 3 ||
            throw(ArgumentError("[structure].lattice vector $k must have 3 components"))
        for i = 1:3
            A[i, k] = Float64(v[i])   # entry k = k-th lattice vector = column k of the matrix
        end
    end
    pbc = if haskey(d, "pbc")
        p = d["pbc"]
        length(p) == 3 || throw(ArgumentError("[structure].pbc must have 3 entries"))
        (Bool(p[1]), Bool(p[2]), Bool(p[3]))
    else
        (true, true, true)
    end
    lattice = Lattice(SMatrix(A); pbc = pbc)

    pos = _input_require(d, "positions", "structure")
    nat = length(pos)
    nat > 0 || throw(ArgumentError("[structure].positions is empty"))
    fr = Matrix{Float64}(undef, 3, nat)
    @inbounds for a = 1:nat
        p = pos[a]
        length(p) == 3 ||
            throw(ArgumentError("[structure].positions[$a] must have 3 components"))
        for i = 1:3
            fr[i, a] = Float64(p[i])
        end
    end

    species = Int[Int(s) for s in _input_require(d, "species", "structure")]
    length(species) == nat ||
        throw(ArgumentError("[structure].species has $(length(species)) entries for $nat atoms"))
    labels = String[String(s) for s in _input_require(d, "species_labels", "structure")]
    return Crystal(lattice, fr, species, labels)   # Crystal validates species range etc.
end

# TOML sub-tables arrive as Dict{String,Any}; body-order keys are digit strings
# ("2"), species/pair keys anything else. Convert to the BasisSpec sugar forms.
_is_bodykey(k::AbstractString) = !isnothing(match(r"^\d+$", k))

function _lmax_from_input(x)
    x isa AbstractDict && return [String(k) => Int(v) for (k, v) in x]
    return Int[Int(v) for v in x]
end

function _lsum_from_input(x)
    x isa Real && return Int(x)
    x isa AbstractDict || throw(ArgumentError(
        "[interaction].lsum must be an integer or a body-order table"))
    return [(_is_bodykey(k) ? parse(Int, k) :
             throw(ArgumentError("[interaction].lsum: key $(repr(k)) is not a " *
                                 "body order"))) => Int(v) for (k, v) in x]
end

_pairtable_from_input(x::AbstractDict, ctx::String) =
    [String(k) => (v isa Real ? Float64(v) :
                   throw(ArgumentError("$ctx: $(repr(k)) must be a number"))) for (k, v) in x]

function _cutoff_from_input(x)
    x isa Real && return Float64(x)
    x isa AbstractDict ||
        throw(ArgumentError("[interaction].cutoff must be a number or a table"))
    ks = collect(keys(x))
    if all(_is_bodykey, ks)          # body-keyed: scalar or pair table per order
        return [parse(Int, k) => (v isa Real ? Float64(v) :
                                  _pairtable_from_input(v, "[interaction].cutoff.$k"))
                for (k, v) in x]
    elseif !any(_is_bodykey, ks)     # one species-pair table for every order
        return _pairtable_from_input(x, "[interaction].cutoff")
    end
    throw(ArgumentError("[interaction].cutoff mixes body-order keys with pair keys"))
end

function _interaction_from_input(d, labels::Vector{String})::BasisSpec
    haskey(d, "pair_cutoff") &&
        throw(ArgumentError("[interaction].pair_cutoff was replaced by `cutoff` " *
                            "(a scalar is equivalent; see the input-schema docstring " *
                            "for per-body / per-pair tables)"))
    nbody = Int(_input_require(d, "nbody", "interaction"))
    lmax = _lmax_from_input(_input_require(d, "lmax", "interaction"))
    cutoff = _cutoff_from_input(_input_require(d, "cutoff", "interaction"))
    lsum = haskey(d, "lsum") ? _lsum_from_input(d["lsum"]) : nothing
    isotropy = haskey(d, "isotropy") ? Bool(d["isotropy"]) : false
    return BasisSpec(labels; nbody = nbody, lmax = lmax, cutoff = cutoff, lsum = lsum,
                     isotropy = isotropy)
end

# ── [moment] ───────────────────────────────────────────────────────────────────────
# The reader only collects and converts; every range / length / consistency check is
# `MomentSpec`'s (one validation locus). Errors raised here are the purely syntactic
# ones a keyword constructor cannot see: unknown keys, wrong TOML value kinds, labels.

const _MOMENT_KEYS = ("nbody", "lmax_mark", "lmax_env", "sampled", "marked",
                      "cutoff_pair", "cutoff_star", "lsum", "isotropy")

# A species list: labels (with "*" = every species) or per-species booleans. Converts
# to the `Vector{Bool}` the spec takes; the length check of the boolean form is the
# spec's. TOML hands mixed arrays over as `Vector{Any}`, so the element kind is tested
# explicitly — `[1, 0]` must not pass as booleans.
function _species_list_from_input(x, labels::Vector{String}, what::String)::Vector{Bool}
    (x isa AbstractVector && !isempty(x)) ||
        throw(ArgumentError("$what must be a non-empty array of species labels or of " *
                            "booleans"))
    nkd = length(labels)
    if all(v -> v isa Bool, x)
        return Bool[v for v in x]
    elseif all(v -> v isa AbstractString, x)
        out = falses(nkd)
        seen = Set{String}()
        for v in x
            k = String(v)
            k in seen && throw(ArgumentError("$what: duplicate entry $(repr(k))"))
            push!(seen, k)
            i = _species_key_index(k, labels, what)
            i === nothing ? (out .= true) : (out[i] = true)
        end
        return out
    end
    throw(ArgumentError("$what must be an array of species labels (strings, \"*\" = " *
                        "every species) or an array of booleans, not $(repr(x))"))
end

# TOML integers arrive as `Int`; `Bool <: Integer`, so `isa Integer` would let `true`
# through as 1 — the kind test is on the concrete type.
_is_toml_int(v) = v isa Int

function _moment_lmax_env_from_input(x, labels::Vector{String})::Vector{Int}
    what = "[moment].lmax_env"
    (x isa AbstractVector || x isa AbstractDict) ||
        throw(ArgumentError("$what must be a per-species integer array or a label " *
                            "table (a bare scalar is not accepted, as for " *
                            "[interaction].lmax)"))
    all(_is_toml_int, x isa AbstractDict ? values(x) : x) ||
        throw(ArgumentError("$what entries must be integers"))
    return _resolve_species_table(_lmax_from_input(x), length(labels), labels, what)
end

function _moment_cutoff_from_input(x, labels::Vector{String}, key::String)
    what = "[moment].$key"
    x isa Bool && throw(ArgumentError("$what must be a number; got $(repr(x))"))
    x isa Real && return Float64(x)
    x isa AbstractDict ||
        throw(ArgumentError("$what must be a number or a species-pair table"))
    any(_is_bodykey, keys(x)) &&
        throw(ArgumentError("$what: body-order keys are not accepted at this level. " *
                            "`cutoff_pair` is one radius for the 2-body clusters; a " *
                            "per-order star table belongs directly under " *
                            "`[moment.cutoff_star]`, and each of ITS entries is a " *
                            "scalar or a species-pair table"))
    return _resolve_pair_table(_pairtable_from_input(x, what), length(labels), labels, what)
end

# `cutoff_star` is the one pointed cutoff that IS per body order (a 4-body probe is
# only affordable on a shell narrower than the 3-body one). A body-keyed table is
# converted to the body-keyed PAIR form `MomentSpec` takes, so the "keys must cover
# exactly 3:nbody" rule stays where every other range check is — in the constructor —
# rather than being enforced twice with two messages. Syntax is this function's only
# business: body keys must not be mixed with species-pair keys.
function _moment_star_cutoff_from_input(x, labels::Vector{String})
    what = "[moment].cutoff_star"
    (x isa AbstractDict && any(_is_bodykey, keys(x))) ||
        return _moment_cutoff_from_input(x, labels, "cutoff_star")
    all(_is_bodykey, keys(x)) ||
        throw(ArgumentError("$what mixes body-order keys with species-pair keys"))
    ks = sort([parse(Int, String(k)) for k in keys(x)])
    return [N => _moment_cutoff_from_input(x[string(N)], labels, "cutoff_star.$N")
            for N in ks]
end

function _moment_scalar(d, key::String, ::Type{T}, kind::String, default) where {T}
    haskey(d, key) || return default
    v = d[key]
    v isa T || throw(ArgumentError("[moment].$key must be $kind; got $(repr(v))"))
    return v
end

function _moment_lsum_from_input(d)::Union{Nothing,Int}
    haskey(d, "lsum") || return nothing
    d["lsum"] isa AbstractDict &&
        throw(ArgumentError("[moment].lsum is one total spin rank per label — no " *
                            "body-order table here (unlike [interaction].lsum)"))
    return _moment_scalar(d, "lsum", Int, "an integer", nothing)
end

function _moment_from_input(d, labels::Vector{String})::MomentSpec
    d isa AbstractDict || throw(ArgumentError("[moment] must be a table"))
    haskey(d, "soc") &&
        throw(ArgumentError("[moment].soc is the upstream spelling — use `isotropy` " *
                            "(opposite polarity: `isotropy = true` keeps L_S = 0 only, " *
                            "which upstream calls `soc = false`)"))
    for k in keys(d)
        k in _MOMENT_KEYS ||
            throw(ArgumentError("[moment]: unknown key $(repr(k)) (allowed: " *
                                join(_MOMENT_KEYS, ", ") * ")"))
    end
    lmax_env = _moment_lmax_env_from_input(_input_require(d, "lmax_env", "moment"), labels)
    sampled = _species_list_from_input(_input_require(d, "sampled", "moment"), labels,
                                       "[moment].sampled")
    marked = haskey(d, "marked") ?
        _species_list_from_input(d["marked"], labels, "[moment].marked") : nothing
    cutoff_pair = _moment_cutoff_from_input(_input_require(d, "cutoff_pair", "moment"),
                                            labels, "cutoff_pair")
    nbody = _moment_scalar(d, "nbody", Int, "an integer", 3)
    cutoff_star = haskey(d, "cutoff_star") ?
        _moment_star_cutoff_from_input(d["cutoff_star"], labels) : nothing
    lmax_mark = _moment_scalar(d, "lmax_mark", Int, "an integer", 2)
    lsum = _moment_lsum_from_input(d)
    isotropy = _moment_scalar(d, "isotropy", Bool, "a boolean", true)
    try
        return MomentSpec(; lmax_env, sampled, lmax_mark, marked, nbody, cutoff_pair,
                          cutoff_star, lsum, isotropy)
    catch err
        err isa ArgumentError || rethrow()
        throw(ArgumentError("[moment]: " * err.msg))
    end
end

function _backend_from_name(name)::AbstractSymmetryBackend
    n = lowercase(String(name))
    n == "none" && return NoSymmetry()
    (n == "spglib" || n == "spg") && return SpglibBackend()
    throw(ArgumentError("[symmetry].backend = $(repr(name)) is not recognized " *
                        "(use \"none\" or \"spglib\")"))
end

function _image_selection_from_name(name)::AbstractImageSelection
    n = lowercase(String(name))
    n == "minimum_image" && return MinimumImage()
    n == "all_images" && return AllImages()
    throw(ArgumentError("[interaction].images = $(repr(name)) is not recognized " *
                        "(use \"minimum_image\" or \"all_images\")"))
end

"""
    read_setup(path) -> (; crystal, spec, backend, tol, images, tie_tol, moment)

Parse a human-authored TOML input file (schema in the file-level docstring of
`src/io/input.jl`) into the in-memory `crystal::Crystal`, `spec::BasisSpec` (from the
file's `[interaction]` section), symmetry `backend::AbstractSymmetryBackend`,
`tol::Float64`, the periodic-image selection `images::AbstractImageSelection`, the
same-distance band `tie_tol::Float64` (`[interaction].tie_tol`, defaulting to the
`SCEBasis` default), and `moment::Union{Nothing,MomentSpec}` — the pointed
site-moment truncation from the optional `[moment]` section (`nothing` when the
section is absent). Training data and the estimator are **not** part of the file
(see [`SCEDataset`](@ref) / [`fit`](@ref)). See also `SCEBasis(path)` and
`MomentBasis(path)`.
"""
function read_setup(path::AbstractString)::@NamedTuple{crystal::Crystal,
                                                       spec::BasisSpec,
                                                       backend::AbstractSymmetryBackend,
                                                       tol::Float64,
                                                       images::AbstractImageSelection,
                                                       tie_tol::Float64,
                                                       moment::Union{Nothing,MomentSpec}}
    doc = TOML.parsefile(path)
    haskey(doc, "structure") ||
        throw(ArgumentError("input file is missing the [structure] section"))
    haskey(doc, "interaction") ||
        throw(ArgumentError("input file is missing the [interaction] section"))
    crystal = _crystal_from_input(doc["structure"])
    spec = _interaction_from_input(doc["interaction"], crystal.species_labels)
    length(spec.lmax) == length(crystal.species_labels) ||
        throw(ArgumentError("[interaction].lmax has $(length(spec.lmax)) entries for " *
                            "$(length(crystal.species_labels)) species"))
    images = haskey(doc["interaction"], "images") ?
        _image_selection_from_name(doc["interaction"]["images"]) : MinimumImage()
    # `tie_tol` changes the emitted basis just like `images` does, so a setup that
    # needed a widened band must be reproducible from its own file — it rides in
    # `[interaction]` next to `images` (validated by the `SCEBasis` constructor).
    tie_tol = haskey(doc["interaction"], "tie_tol") ?
        Float64(doc["interaction"]["tie_tol"]) : _SAME_DIST_RTOL
    sym = get(doc, "symmetry", Dict{String,Any}())
    backend = haskey(sym, "backend") ? _backend_from_name(sym["backend"]) : NoSymmetry()
    tol = haskey(sym, "tol") ? Float64(sym["tol"]) : 1e-5
    moment = haskey(doc, "moment") ?
        _moment_from_input(doc["moment"], crystal.species_labels) : nothing
    return (; crystal, spec, backend, tol, images, tie_tol, moment)
end

"""
    SCEBasis(path::AbstractString; backend = nothing, tol = nothing, images = nothing,
             tie_tol = nothing) -> SCEBasis

Build an [`SCEBasis`](@ref) directly from a TOML input file ([`read_setup`](@ref)).
The file's `[symmetry]` backend/tol and `[interaction]` `images`/`tie_tol` are used
unless overridden by the keyword arguments (e.g. `backend = SpglibBackend()` forces
Spglib regardless of the file). Using the Spglib backend requires `using Spglib`.
"""
function SCEBasis(path::AbstractString;
                  backend::Union{Nothing,AbstractSymmetryBackend} = nothing,
                  tol::Union{Nothing,Real} = nothing,
                  images::Union{Nothing,AbstractImageSelection} = nothing,
                  tie_tol::Union{Nothing,Real} = nothing)::SCEBasis
    inp = read_setup(path)
    be = backend === nothing ? inp.backend : backend
    tl = tol === nothing ? inp.tol : Float64(tol)
    im = images === nothing ? inp.images : images
    tt = tie_tol === nothing ? inp.tie_tol : Float64(tie_tol)
    return SCEBasis(inp.crystal, inp.spec; backend = be, tol = tl, images = im,
                    tie_tol = tt)
end

"""
    MomentBasis(path::AbstractString; backend = nothing, tol = nothing, tie_tol = nothing)
        -> MomentBasis

Build a pointed [`MomentBasis`](@ref) directly from a TOML input file
([`read_setup`](@ref)): the crystal from `[structure]`, the truncation from the
`[moment]` section (required here — its absence is an `ArgumentError` naming the
section), symmetry from `[symmetry]`, and the same-distance band from
`[interaction].tie_tol` (the file's `[interaction]` section is still required, as for
every setup file). The symmetry settings (`backend`, `tol`) and `tie_tol` are each
overridable by the keyword of the same name; the crystal and the truncation always
come from the file. The moment basis is always minimum-image — an
`[interaction].images = "all_images"` setting is ignored here, with a warning. Using
the Spglib backend requires `using Spglib`.
"""
function MomentBasis(path::AbstractString;
                     backend::Union{Nothing,AbstractSymmetryBackend} = nothing,
                     tol::Union{Nothing,Real} = nothing,
                     tie_tol::Union{Nothing,Real} = nothing)::MomentBasis
    inp = read_setup(path)
    inp.moment === nothing &&
        throw(ArgumentError("input file is missing the [moment] section — " *
                            "MomentBasis(path) needs the pointed truncation there " *
                            "(see the input-schema docstring)"))
    inp.images isa AllImages &&
        @warn "[interaction].images = \"all_images\" does not apply to the moment " *
              "basis, which is always minimum-image; ignored" path = path
    be = backend === nothing ? inp.backend : backend
    tl = tol === nothing ? inp.tol : Float64(tol)
    tt = tie_tol === nothing ? inp.tie_tol : Float64(tie_tol)
    return MomentBasis(inp.crystal, inp.moment; backend = be, tol = tl, tie_tol = tt)
end
