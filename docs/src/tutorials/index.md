# Tutorials

Narrated, end-to-end runs of the full workflow. Each is self-contained and executable, and
each is also available as a script under [`examples/`](https://github.com/Tomonori-Tanaka/Magesty_rebuild.jl)
in the repository.

```@contents
Pages = ["heisenberg_chain.md", "kagome_threebody.md"]
Depth = 1
```

- **[Heisenberg chain](heisenberg_chain.md)** — the canonical first run: build a 2-body
  isotropic basis, fit synthetic Heisenberg energies, recover the coupling `J`, then add a
  torque co-fit. Establishes the `basis → dataset → fit → predict` rhythm.
- **[Kagome three-body](kagome_threebody.md)** — beyond pairs: a 3-body interaction on the
  kagome lattice, where symmetry-equivalent sites force several `l`-orderings into one
  multi-term SALC. Recovers a synthetic 3-body model from energies and torques.

For shorter, feature-by-feature snippets see the [Guide](../guide/basis.md); for the
fastest possible first fit see [Getting started](../getting_started.md).
