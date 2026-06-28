# Persistence and I/O

```@meta
CurrentModule = MagestyRebuild
```

Building a SALC basis is the expensive step, and a fitted model is worth keeping. This
page covers three I/O seams: saving and reloading models, building a basis from a
human-authored `input.toml`, tabulating coefficients, and reading DFT training data
through the code-agnostic source interface.

## Saving and reloading a model

`MagestyRebuild.save` / `MagestyRebuild.load` serialize a basis or a model to a
self-contained, human-readable **TOML** document. (They are intentionally *not* exported —
the names clash with FileIO / JLD2 / CSV — so qualify them.)

```julia
MagestyRebuild.save("model.toml", SCEModel(f))      # or save("basis.toml", basis)
model = MagestyRebuild.load(SCEModel, "model.toml")
predict_energy(model, configs)
```

The document stores the crystal, the space-group operations, the interaction, the full
SALC basis, and (for a model) `j0` and the per-`SALCKey` coefficients. On reload the basis
is rebuilt *verbatim* (no re-projection), and coefficients are re-paired to it **by key**,
not by position — so a reloaded model predicts identically even if the basis were rebuilt
in a different order. TOML is chosen over JSON deliberately: stdlib `TOML` round-trips
`Float64` exactly and expresses the deep nested SALC document, so input and dump share one
zero-dependency format.

## A human-authored `input.toml`

Instead of constructing `Crystal` / `Interaction` in Julia, you can describe the *setup*
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
nbody       = 2
pair_cutoff = 2.6           # or `inf` for the whole Wigner–Seitz cell
lmax        = [1]
isotropy    = true

[symmetry]
backend = "spglib"          # or "none"
tol     = 1.0e-5
```

```julia
basis = SCEBasis("input.toml")     # reads [symmetry] backend/tol from the file
```

[`read_input`](@ref) returns the parsed setup (including the image selection) if you want
to inspect it before building.

## Tabular coefficients

[`coeftable`](@ref) returns an [`SCECoefficients`](@ref) — a **Tables.jl** source with one
row per SALC (`body`, `orbit_id`, `ls`, `Lf`, `block`, `J`). It drops straight into any
table or IO package:

```julia
using DataFrames
df = DataFrame(coeftable(f))        # or CSV.write("J.csv", coeftable(f))
intercept(f)                        # the reference energy j0 (not a row)
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
using MagestyRebuild, SCETools
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

Next: [Sunny export](sunny.md).
