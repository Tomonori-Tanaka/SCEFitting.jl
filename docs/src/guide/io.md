# Persistence and I/O

```@meta
CurrentModule = SCEFitting
```

Building a SALC basis is the expensive step, and a fitted model is worth keeping. This
page covers three I/O seams: saving and reloading models, building a basis from a
human-authored `input.toml`, tabulating coefficients, and reading DFT training data
through the code-agnostic source interface.

## Saving and reloading a model

`SCEFitting.save` / `SCEFitting.load` serialize a basis or a model to a
self-contained, human-readable **TOML** document. (They are intentionally *not* exported —
the names clash with FileIO / JLD2 / CSV — so qualify them.)

```julia
SCEFitting.save("model.toml", SCEPredictor(sce_fit))      # or save("basis.toml", basis)
model = SCEFitting.load(SCEPredictor, "model.toml")
predict_energy(model, configs)
```

The document stores the crystal, the space-group operations, the basis spec, the full
SALC basis, and (for a model) `j0` and the per-`SALCKey` coefficients. On reload the basis
is rebuilt *verbatim* (no re-projection), and coefficients are re-paired to it **by key**,
not by position — so a reloaded model predicts identically even if the basis were rebuilt
in a different order. TOML is chosen over JSON deliberately: stdlib `TOML` round-trips
`Float64` exactly and expresses the deep nested SALC document, so input and dump share one
zero-dependency format.

## A human-authored `input.toml`

Instead of constructing `Crystal` / `BasisSpec` in Julia, you can describe the *setup*
(crystal + interaction + optional symmetry) in a TOML file and build the basis from it.
Training data and the estimator stay in Julia.

```toml
# input.toml
[structure]
lattice   = [[8.0, 0.0, 0.0], [0.0, 8.0, 0.0], [0.0, 0.0, 10.0]]   # each entry = one cell vector
positions = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.25], [0.0, 0.0, 0.5], [0.0, 0.0, 0.75]]
species   = [1, 1, 1, 1]
species_labels = ["Fe"]

[interaction]
nbody    = 2
cutoff   = 2.6              # or `inf` for the whole Wigner–Seitz cell
lmax     = [1]
isotropy = true

[symmetry]
backend = "spglib"          # or "none"
tol     = 1.0e-5            # symprec: a Cartesian distance in Å, in (0, 0.1]
```

```julia
basis = SCEBasis("input.toml")     # reads [symmetry] backend/tol from the file
```

The `[interaction]` section also takes the label-keyed / per-body forms (see
[The interaction specification](basis.md#The-interaction-specification)) — species
tables with a `"*"` fallback, a per-body-order `lsum`, and cutoffs per body order
and species pair:

```toml
[interaction]
nbody    = 3
isotropy = true

[interaction.lmax]
"*" = 3
B   = 0

[interaction.lsum]          # keys = body orders; omitted orders are uncapped
1 = 0
2 = 4
3 = 4

[interaction.cutoff]        # body-keyed scalars ...
2 = 8.0

[interaction.cutoff.3]      # ... or a species-pair table per order
"Fe-*" = 6.0
"*-*"  = 8.0
```

[`read_setup`](@ref) returns the parsed setup (including the image selection) if you want
to inspect it before building.

### The pointed site-moment basis

The truncation of the adiabatic site-moment channel ([Site moments](moment.md)) can
live in the same file as an optional `[moment]` section, so a fit script builds both
bases from one setup:

```toml
[moment]
nbody       = 3            # optional, default 3 (1 to 4; members grow as C(z,N-1)*N!)
lmax_mark   = 2            # optional, default 2: the marked site's own ê factor
lmax_env    = [2]          # per species (index order), or a label table
sampled     = ["Fe"]       # REQUIRED: species the consumer samples (labels, ["*"], or booleans)
marked      = ["Fe"]       # optional, default every species: whose moments are expanded
cutoff_pair = 4.1          # REQUIRED: mark–environment bond radius (Å); `inf` = whole WS cell
cutoff_star = 4.1          # optional, default = cutoff_pair: a star's N-1 mark spokes
                           #   (or a body-keyed table, one entry per star order)
lsum        = 4            # optional, default uncapped: Σl cap per label
                           #   (or a PARTIAL body-keyed table; see below)
isotropy    = true         # optional, default true (L_S = 0 only)
```

```julia
basis        = SCEBasis("input.toml")
moment_basis = MomentBasis("input.toml")   # same [symmetry] and [interaction].tie_tol
```

Every value is handed to the [`MomentSpec`](@ref) keyword constructor, which owns the
validation (`sampled` is required there too, and every species with `lmax_env > 0`
must be sampled). `lmax_env` takes the label-table form (`[moment.lmax_env]` with a
`"*"` fallback) and `cutoff_pair` / `cutoff_star` the species-pair-table form
(`[moment.cutoff_pair]` with `"Fe-*"` / `"*-*"` keys), exactly as in `[interaction]`.
`cutoff_pair` takes no body-order table — it is the 2-body radius — but `cutoff_star`
does, one entry per star order, and the keys must cover exactly `3:nbody` so that no
order is left silently at the default:

```toml
[moment.cutoff_star]        # per star order: a scalar per order ...
3 = 4.1
[moment.cutoff_star.4]      # ... or a species-pair table for one order
"Fe-Fe" = 2.5
"*-*"   = 0.0
```

`lsum` takes a body-keyed table too, but keyed by the body order itself (from 1, not
from 3) and **partial**: an order the table does not name stays uncapped.

```toml
[moment.lsum]               # 3-body up to Σl = 6, 4-body up to 4; 1- and 2-body uncapped
3 = 6
4 = 4
```

The two tables differ on purpose. A missing radius would silently leave a star order at
`cutoff_pair` — a number with no relation to what that order needs — so `cutoff_star`
must cover `3:nbody` exactly; a missing cap just means "no cap", which is both safe and
the documented default. Sub-tables must come after the bare `[moment]` keys, as always
in TOML.

Unknown keys are refused, and so is the upstream spelling `soc` (use `isotropy`; the
polarities are opposite). Two things to keep in mind:

- `[moment].isotropy` defaults to **`true`** — the `MomentSpec` default — whereas
  `[interaction].isotropy` defaults to `false`. Write both explicitly when they matter.
- The moment basis is always minimum-image (`[interaction].images` does not apply to
  it) and shares `[interaction].tie_tol`, so a widened same-distance band applies to
  both channels. `MomentBasis(path)` still needs the `[interaction]` section, like
  every setup file.

## Tabular coefficients

[`coeftable`](@ref) returns an [`SCECoefficients`](@ref) — a **Tables.jl** source with one
row per SALC. The columns are read straight off each [`SALCKey`](@ref) (the stable
design-matrix-column identity), plus the fitted coefficient:

| Column | Meaning |
|---|---|
| `body` | body order ``N`` of the cluster (2 = pair, 3 = triplet, …) |
| `orbit_id` | index of the cluster symmetry orbit at that body order |
| `decors` | the per-site decoration multiset, as a comma-joined string. A pure-spin site renders as its bare ``l``, so this package's keys read `"1,1,2"` exactly as the old `ls` column did; a displacement factor would render as `u(k:l)` (colon inside the token, so the column always splits on commas back into its sites) |
| `L_S` | total coupled **spin** rank of the label. Every key this package builds is pure spin, so `L_S == Lf` throughout; the column exists because the key layout is shared with the spin–lattice expansion, where the two differ |
| `Lf` | final coupled angular momentum ``L_f`` of the invariant |
| `block` | disambiguates independent ``l``-orderings / coupling paths sharing the same `(body, orbit_id, decors, L_S, Lf)` |
| `J` | the fitted coefficient ``j_\varphi`` for this SALC (DFT energy unit, e.g. eV) |
| `alias_group` | index into `alias_groups(basis)` of the cross-orbit alias group the SALC belongs to (`0` for none) — distinct orbits the training cell cannot tell apart, tied into one design column |
| `split` | `:free`, or `:convention` when `J` was read back from a tied column with the equal per-bond weight (a convention, not a measurement); `:legacy` for a model loaded from a pre-v7 file |

The rows are in design-matrix column order, the same order as [`coef`](@ref) and
`basis.salc_basis.keys`. It drops straight into any table or IO package:

```julia
using DataFrames
coef_df = DataFrame(coeftable(sce_fit))   # or CSV.write("J.csv", coeftable(sce_fit))
intercept(sce_fit)                        # the reference energy j0 (not a row)
```

The boundary is deliberate: the library owns the mapping from internal storage
(`SALCKey` + coefficient position) to meaningful, labeled rows; the caller brings the
table / IO / plotting package.

## Reading DFT data

DFT-code I/O is isolated at the *training-data boundary*. The core owns only the
**abstract seam**: a source implements [`read_configs`](@ref)`(src) -> Vector{SpinDatum}`,
and the SCE pipeline only ever sees the code-agnostic [`SpinDatum`](@ref) /
[`SCEDataset`](@ref) — once you have the data, the originating code is irrelevant.

```julia
basis   = SCEBasis(crystal, interaction)
src     = some_source                                    # any AbstractDFTSource
dataset = SCEDataset(basis, src)                         # read_configs(src) under the hood
fit(SCEFit, dataset, OLS(); torque_weight = 0.5)
```

The torque target from a constrained calculation is
``\boldsymbol\tau_a = \boldsymbol m_a \times \boldsymbol B_a`` (the physical /
Landau–Lifshitz torque), the *same* physical quantity, sign, and layout as the model's
[`predict_torque`](@ref) ``= -\hat{\boldsymbol e}_a \times \partial E/\partial\hat{\boldsymbol e}_a``,
so the co-fit is consistent.

The **concrete DFT-code adapters** live in the companion
[SCETools.jl](https://github.com/Tomonori-Tanaka/SCETools.jl) package, not in the core. For
VASP:

```julia
using SCEFitting, SCETools
using SCETools.VASP: read_poscar, Oszicar

crystal = read_poscar("POSCAR")                          # → Crystal
basis   = SCEBasis(crystal, interaction)
src     = Oszicar(["run1/OSZICAR", "run2/OSZICAR"])      # an AbstractDFTSource (constrained NCL)
dataset = SCEDataset(basis, src)
fit(SCEFit, dataset, OLS(); torque_weight = 0.5)
```

Adding another code is one more sibling adapter in SCETools — the core and its exports do
not change. SCETools.VASP also writes the *inverse* direction (sampled configurations →
constrained-noncollinear INCAR / input sets) for generating new training data.

### Legacy Magesty training sets (EMBSET)

The one concrete format the core ships is Magesty's **EMBSET** training set — it is
DFT-code-agnostic (energies plus per-atom moment and constraining-field vectors, exactly
what [`SpinDatum`](@ref) stores), so an existing Magesty data set drops straight in:

```julia
dataset = SCEDataset(basis, EmbsetFile("EMBSET"))        # or read_embset("EMBSET")
```

See [`read_embset`](@ref) for the format and its validation rules (the reader is
cross-checked against Magesty's own in the oracle suite).

Next: [Sunny export](sunny.md).
