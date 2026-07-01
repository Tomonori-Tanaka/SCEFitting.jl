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
    nbody       = 2
    pair_cutoff = 3.0                                     # Å, or `inf` for the full Wigner–Seitz cell
    lmax        = [2,2]                                    # per species
    isotropy    = false                                   # optional, default false
    images      = "minimum_image"                          # optional: "minimum_image" (default) or "all_images"

    [symmetry]                                            # optional section
    backend = "spglib"                                    # "none" (default) or "spglib"
    tol     = 1.0e-5                                       # optional, default 1e-5
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

function _interaction_from_input(d)::BasisSpec
    nbody = Int(_input_require(d, "nbody", "interaction"))
    pair_cutoff = Float64(_input_require(d, "pair_cutoff", "interaction"))
    lmax = Int[Int(x) for x in _input_require(d, "lmax", "interaction")]
    isotropy = haskey(d, "isotropy") ? Bool(d["isotropy"]) : false
    return BasisSpec(; nbody = nbody, pair_cutoff = pair_cutoff, lmax = lmax, isotropy = isotropy)
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
    read_setup(path) -> (; crystal, spec, backend, tol, images)

Parse a human-authored TOML input file (schema in the file-level docstring of
`src/io/input.jl`) into the in-memory `crystal::Crystal`, `spec::BasisSpec` (from the
file's `[interaction]` section), symmetry `backend::AbstractSymmetryBackend`,
`tol::Float64`, and the periodic-image selection `images::AbstractImageSelection`.
Training data and the estimator are **not** part of the file (see
[`SCEDataset`](@ref) / [`fit`](@ref)). See also `SCEBasis(path)`.
"""
function read_setup(path::AbstractString)::@NamedTuple{crystal::Crystal,
                                                       spec::BasisSpec,
                                                       backend::AbstractSymmetryBackend,
                                                       tol::Float64,
                                                       images::AbstractImageSelection}
    doc = TOML.parsefile(path)
    haskey(doc, "structure") ||
        throw(ArgumentError("input file is missing the [structure] section"))
    haskey(doc, "interaction") ||
        throw(ArgumentError("input file is missing the [interaction] section"))
    crystal = _crystal_from_input(doc["structure"])
    spec = _interaction_from_input(doc["interaction"])
    length(spec.lmax) == length(crystal.species_labels) ||
        throw(ArgumentError("[interaction].lmax has $(length(spec.lmax)) entries for " *
                            "$(length(crystal.species_labels)) species"))
    images = haskey(doc["interaction"], "images") ?
        _image_selection_from_name(doc["interaction"]["images"]) : MinimumImage()
    sym = get(doc, "symmetry", Dict{String,Any}())
    backend = haskey(sym, "backend") ? _backend_from_name(sym["backend"]) : NoSymmetry()
    tol = haskey(sym, "tol") ? Float64(sym["tol"]) : 1e-5
    return (; crystal, spec, backend, tol, images)
end

"""
    SCEBasis(path::AbstractString; backend = nothing, tol = nothing, images = nothing)
        -> SCEBasis

Build an [`SCEBasis`](@ref) directly from a TOML input file ([`read_setup`](@ref)).
The file's `[symmetry]` backend/tol and `[interaction].images` are used unless
overridden by the keyword arguments (e.g. `backend = SpglibBackend()` forces Spglib
regardless of the file). Using the Spglib backend requires `using Spglib`.
"""
function SCEBasis(path::AbstractString;
                  backend::Union{Nothing,AbstractSymmetryBackend} = nothing,
                  tol::Union{Nothing,Real} = nothing,
                  images::Union{Nothing,AbstractImageSelection} = nothing)::SCEBasis
    inp = read_setup(path)
    be = backend === nothing ? inp.backend : backend
    tl = tol === nothing ? inp.tol : Float64(tol)
    im = images === nothing ? inp.images : images
    return SCEBasis(inp.crystal, inp.spec; backend = be, tol = tl, images = im)
end
