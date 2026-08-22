# Theory

This section explains *what* the package computes and *why the rebuild is shaped the way
it is*. It has four chapters, readable independently:

```@contents
Pages = ["sce.md", "resolvability.md", "moment.md", "architecture.md"]
Depth = 2
```

- **[The spin-cluster expansion](sce.md)** — the formalism: the energy as a linear model in
  symmetry-adapted invariants of real spherical harmonics, the design matrix, the fit, and
  the torque as the analytic derivative of the same surface.
- **[Periodic resolvability](resolvability.md)** — why a finite supercell can only resolve
  the minimum-image (Wigner–Seitz-cell) interactions, the face/edge/corner boundary ties,
  and the compact-cluster criterion that extends this to three- and four-body clusters.
- **[Predicting site moments](moment.md)** — the pointed (site-marked) expansion of the
  adiabatic moment ``m_i(\boldsymbol e)``: covariance instead of invariance, the mark as a
  decoration, why the same projector works (Frobenius reciprocity), the site sum rule
  that ties it to the energy expansion, the signed-projection regression and its gate.
- **[Architecture](architecture.md)** — the design choices that distinguish this rebuild
  from Magesty.jl: pluggable seams, the core/extension split, canonical `SALCKey`
  addressing, the combined-space projection for `N ≥ 3`, and the numerical-oracle
  validation methodology.

For the prose-heavy rationale behind individual refinements, the repository's
`docs/design-notes.md` is the long-form companion to this section.
