"""
Tabular access to fitted SCE coefficients.

`coeftable(fit_or_model)` returns an [`SCECoefficients`](@ref): a labeled, ordered
view of the SALC coefficients `Jϕ`, one row per design-matrix column, keyed by the
structural fields of each [`SALCKey`](@ref). It implements the **Tables.jl**
interface, so it drops straight into `DataFrame`, `CSV.write`, `Arrow.write`, … — the
library supplies the data source (it owns the mapping from internal storage to
meaningful rows); the caller brings the table / IO package. The reference energy `j0`
is the model intercept (not a per-row quantity); read it with [`intercept`](@ref).
"""

const _COEF_NAMES = (:body, :orbit_id, :decors, :L_S, :Lf, :block, :J, :alias_group, :split)
const _COEF_COLTYPES = (Int, Int, String, Int, Int, Int, Float64, Int, Symbol)
const _CoefRow = NamedTuple{_COEF_NAMES,
                            Tuple{Int,Int,String,Int,Int,Int,Float64,Int,Symbol}}

# The sorted decoration label as a flat string. A pure-spin decor renders as its
# bare `l` (so a v4-shaped key reads "1,1,2" exactly as before); a displacement
# factor renders as `u(k:l)`, a combined decor as `l+u(k:l)`. The colon inside
# the token is deliberate: the sites are comma-joined, so a token may not
# contain a comma or the column stops being splittable (`"2+u(0,1),u(1,0)"`
# would read as four sites).
function _decor_string(decors::AbstractVector{SiteDecor})::String
    token(d) = is_pure_spin(d) ? string(d.spin_l) :
               has_spin(d) ? "$(d.spin_l)+u($(d.disp_k):$(d.disp_l))" :
               "u($(d.disp_k):$(d.disp_l))"
    return join((token(d) for d in decors), ",")
end

"""
    SCECoefficients

A Tables.jl-compatible table of fitted SCE coefficients: parallel `keys` /
`jphi` (column `J`) / `alias_group` (index into [`alias_groups`](@ref), `0` for none)
/ `split` (`:free`, `:convention`, `:legacy` — see [`SCEPredictor`](@ref)), plus the
intercept `j0`. Build it with [`coeftable`](@ref); iterate it for `NamedTuple` rows,
or hand it to any Tables.jl sink.
"""
struct SCECoefficients
    keys::Vector{SALCKey}
    jphi::Vector{Float64}
    j0::Float64
    alias_group::Vector{Int}
    split::Vector{Symbol}

    function SCECoefficients(keys::Vector{SALCKey}, jphi::Vector{Float64}, j0::Real,
                             alias_group::Vector{Int}, split::Vector{Symbol})
        length(keys) == length(jphi) ||
            throw(ArgumentError("got $(length(keys)) keys for $(length(jphi)) coefficients"))
        length(alias_group) == length(jphi) && length(split) == length(jphi) ||
            throw(ArgumentError("alias_group / split must have one entry per coefficient"))
        return new(keys, jphi, Float64(j0), alias_group, split)
    end
end

# Without provenance columns: no alias group, every coefficient `:free`.
SCECoefficients(keys::Vector{SALCKey}, jphi::Vector{Float64}, j0::Real) =
    SCECoefficients(keys, jphi, j0, zeros(Int, length(jphi)), fill(:free, length(jphi)))

@inline _coef_row(c::SCECoefficients, i::Int)::_CoefRow =
    (body = c.keys[i].body, orbit_id = c.keys[i].orbit_id,
     decors = _decor_string(c.keys[i].decors), L_S = c.keys[i].L_S,
     Lf = c.keys[i].Lf, block = c.keys[i].block, J = c.jphi[i],
     alias_group = c.alias_group[i], split = c.split[i])

"""
    coeftable(f::SCEFit) -> SCECoefficients
    coeftable(m::SCEPredictor) -> SCECoefficients

A Tables.jl-compatible table of the fitted coefficients — one row per SALC with
columns `body`, `orbit_id`, `decors` (the sorted decoration label as a string:
a pure-spin key reads like `"1,1,2"`, displacement factors as `u(k:l)`), `L_S`,
`Lf`, `block`, `J` (the coefficient `Jϕ`), `alias_group` (the index into
[`alias_groups`](@ref) of the cross-orbit alias group the SALC belongs to, `0` for
none) and `split` (`:free`, or `:convention` when `J` was read back from a tied
column with the equal per-bond weight and is therefore a convention rather than a
measurement; `:legacy` for a model loaded from a pre-v7 file). For an `SCEFit`
the table is the SALC-space expansion of `coef(f)` (see [`SCEPredictor`](@ref)`(f)`).
The intercept `j0` is available via [`intercept`](@ref).
Example: `using DataFrames; DataFrame(coeftable(f))`.
"""
coeftable(f::SCEFit)::SCECoefficients = coeftable(SCEPredictor(f))
function coeftable(m::SCEPredictor)::SCECoefficients
    # A predictor may carry a key/coefficient subset of its basis (a display or
    # hand-built slice); the alias index is per basis SALC, so it applies only to a
    # full-length coefficient vector.
    ties = _column_ties(m.basis)
    ag = length(m.jphi) == length(ties.col_of) ? _alias_index(ties) : zeros(Int, length(m.jphi))
    return SCECoefficients(copy(m.keys), copy(m.jphi), m.j0, ag, copy(m.split))
end
# `copy` duplicates the `jphi` and key vectors so the table is independent of later
# refits; the `SALCKey`s themselves are shared (they are treated as immutable — used
# as `Dict` keys — so their `ls` is never mutated).

# --- collection / accessor interface ---
Base.length(c::SCECoefficients) = length(c.keys)
Base.eltype(::Type{SCECoefficients}) = _CoefRow
Base.getindex(c::SCECoefficients, i::Integer) = _coef_row(c, Int(i))
Base.firstindex(::SCECoefficients) = 1
Base.lastindex(c::SCECoefficients) = length(c)
function Base.iterate(c::SCECoefficients, i::Int = 1)
    i > length(c) && return nothing
    return (_coef_row(c, i), i + 1)
end

coef(c::SCECoefficients)::Vector{Float64} = c.jphi
intercept(c::SCECoefficients)::Float64 = c.j0

# --- Tables.jl row-access source ---
Tables.istable(::Type{SCECoefficients}) = true
Tables.rowaccess(::Type{SCECoefficients}) = true
Tables.rows(c::SCECoefficients) = c
Tables.schema(::SCECoefficients) = Tables.Schema(_COEF_NAMES, _COEF_COLTYPES)

# --- display ---
_fmtj(x::Float64)::String = string(round(x; sigdigits = 6))

Base.show(io::IO, c::SCECoefficients) =
    print(io, "SCECoefficients(", length(c), " terms, j0=", c.j0, ")")

function Base.show(io::IO, ::MIME"text/plain", c::SCECoefficients)
    n = length(c)
    print(io, "SCECoefficients: ", n, " term", n == 1 ? "" : "s", ", j0 = ", c.j0)
    n == 0 && return
    nshow = min(n, 20)
    cols = (["body"; [string(c.keys[i].body) for i = 1:nshow]],
            ["orbit_id"; [string(c.keys[i].orbit_id) for i = 1:nshow]],
            ["decors"; [_decor_string(c.keys[i].decors) for i = 1:nshow]],
            ["L_S"; [string(c.keys[i].L_S) for i = 1:nshow]],
            ["Lf"; [string(c.keys[i].Lf) for i = 1:nshow]],
            ["block"; [string(c.keys[i].block) for i = 1:nshow]],
            ["J"; [_fmtj(c.jphi[i]) for i = 1:nshow]],
            ["alias"; [c.alias_group[i] == 0 ? "" : string(c.alias_group[i]) for i = 1:nshow]],
            ["split"; [c.split[i] === :free ? "" : string(c.split[i]) for i = 1:nshow]])
    w = map(col -> maximum(length, col), cols)
    for r = 1:(nshow + 1)
        println(io)
        for k = 1:length(cols)
            print(io, "  ", lpad(cols[k][r], w[k]))
        end
    end
    n > nshow && print(io, "\n  ⋮ (", n - nshow, " more)")
    return
end
