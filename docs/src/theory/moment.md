# Predicting site moments: the pointed expansion

```@meta
CurrentModule = SCEFitting
```

[The spin-cluster expansion](sce.md) fits the *energy* of a spin configuration. This chapter
explains how the same machinery predicts the *moment magnitude of every site* as a function
of the configuration — the "pointed" (site-marked) expansion behind [`MomentBasis`](@ref),
[`MomentDataset`](@ref) and [`predict_moment`](@ref). The [guide](../guide/moment.md) shows
how to run it; this page explains what it computes and why it is built the way it is. It is
written to be read top to bottom, starting from a three-atom example that already contains
every idea.

## The problem: the adiabatic map

A constrained noncollinear DFT calculation fixes the spin *directions*
``\boldsymbol e = \{\hat{\boldsymbol e}_a\}`` and lets the electrons relax. Two things come
out: the total energy ``E(\boldsymbol e)``, and the magnetic moment of every site — a vector
``\boldsymbol M_a`` whose *magnitude* the electrons chose for that particular direction
pattern. The correspondence

```math
\boldsymbol e \;\longmapsto\; m_a(\boldsymbol e), \qquad a = 1, \dots, n
```

is the **adiabatic map**: the moment a site carries once the electrons have relaxed around a
given configuration. It is a function of the configuration in exactly the same sense as
``E(\boldsymbol e)`` is. The cluster expansion was invented to write ``E(\boldsymbol e)``
with a few symmetry-determined coefficients; the claim of this chapter is that
``m_a(\boldsymbol e)`` can be expanded with the *same* machine, and that its expansion
carries information the energy expansion cannot.

Why one wants it: a Monte-Carlo consumer needs ``m_a`` for the Zeeman energy
``-\mu_B\, m_a\, \hat{\boldsymbol e}_a \cdot \boldsymbol B`` and for ``M(T)``; a dynamics
consumer divides by ``m_a`` in the precession rate and in the thermal-noise amplitude; and
the moment of an *induced* site (Rh in FeRh, Ge in FeGe) is not a variable of its own at all
— it is entirely a function of its neighbours' directions, and a model that treats it as a
constant is simply wrong.

## Invariant versus covariant

The energy is a scalar of the whole configuration: moving the configuration with a
space-group operation ``g`` leaves it unchanged,
``E(g \cdot \boldsymbol e) = E(\boldsymbol e)``. This is **invariance**, and the energy
basis functions are built to have it.

A site moment *names a site*. Moving the configuration with ``g`` moves the site too:

```math
m_{g(a)}(g \cdot \boldsymbol e) \;=\; m_a(\boldsymbol e) \qquad (\forall\, g \in G).
```

This is **covariance**. In words: *symmetry-equivalent sites, placed in symmetry-equivalent
environments, carry the same moment.* That sentence is the entire symmetry content of the
pointed expansion, and the difference between it and invariance is why the group used for
symmetrization will shrink from the cluster's stabilizer ``H`` to a subgroup ``H_i``.

## Why an evaluation axis

Constrained DFT does not return a number ``m_a``; it returns a vector ``\boldsymbol M_a``,
and the *constraint* decided along which axis that moment was measured — the prescribed
direction for a direction-pinning constraint, the penalty axis for a transverse-penalty one.
The quantity expanded is therefore

```math
m_a \;=\; m\bigl(\hat{\boldsymbol e}_a;\ \boldsymbol e_{\text{env}}\bigr),
```

"the magnitude measured along the axis ``\hat{\boldsymbol e}_a``, as a function of the
environment", and the observable is the **signed projection**
``y = \hat{\boldsymbol e}_a \cdot \boldsymbol M_a``. In direction-pinning data
``\hat{\boldsymbol e}_a`` coincides with the converged direction, so one may ask whether the
distinction is needed. It is: (i) for a transverse penalty the axis is the constraint axis
and differs from the converged direction; (ii) at run time the question is "given a
direction, what is the magnitude?", so the axis is best kept as an independent variable;
(iii) in the basis functions the axis is the marked site's own factor and ``\boldsymbol e``
the environment's, and the [sum rule](@ref "The site sum rule") below is precisely the
identity obtained by setting the two equal.

## The smallest example: a three-site star

Before the general construction, one star cluster shows every concept. Take a central atom
``i`` and its two neighbours ``j, k``, so the cluster is ``S = (i, j, k)``, and let the
crystal have a mirror ``\sigma`` that swaps ``j \leftrightarrow k`` and fixes ``i``: the
cluster's stabilizer is ``H = \{1, \sigma\}``.

**(a) Energy side — symmetrize with ``H``.** Put ``l = 1`` on every site: the product of
three ``l = 1`` harmonics is a cubic in the spin components, odd under time reversal, so it
is dropped. The next candidate, ``\boldsymbol l = (0, 1, 1)`` (``l = 0`` on the centre),
would be ``\boldsymbol e_j \cdot \boldsymbol e_k`` — but the **energy basis forbids
``l = 0`` on any site**: ``\bar Z_{00} \equiv 1``, so site ``i`` vanishes from the function
and the column becomes indistinguishable from the pair cluster ``(j, k)``. The energy-side
star therefore starts at ``\boldsymbol l = (2, 1, 1)``, whose isotropic component is

```math
(\boldsymbol e_i \cdot \boldsymbol e_j)(\boldsymbol e_i \cdot \boldsymbol e_k)
  - \tfrac13\, \boldsymbol e_j \cdot \boldsymbol e_k .
```

**(b) Pointed side — mark site ``i``, symmetrize with ``H_i``.** A *mark* on site ``i``
declares "this function describes the moment of site ``i``". It is a label, not a value:
whether ``\boldsymbol e_i`` appears in the function is a separate question. Marking does two
things:

- Only operations that *fix the mark* may be used for symmetrization. With the mark on
  ``i``, ``\sigma`` fixes it and ``H_i = H``. With the mark on ``j``, ``\sigma`` would carry
  it to ``k``, so ``H_j = \{1\}`` — and "the function marked at ``j``" and "the function
  marked at ``k``" are one channel, related by ``\sigma``.
- **``l = 0`` becomes legal on the marked site.** The mark remembers which site the function
  is about, so the site does not vanish from the expansion when its own harmonic is the
  constant.

Thus with the mark on ``i`` and ``(l_{\mathrm{mark}}, l_j, l_k) = (0, 1, 1)``,

```math
\Psi^{(i)}(\boldsymbol e, \hat{\boldsymbol e}_i) \;=\; \boldsymbol e_j \cdot \boldsymbol e_k
```

is a valid pointed basis function. It does **not** depend on the evaluation axis: it says
"the moment on ``i`` is large when its neighbours are parallel and small when they are
antiparallel" — the **induced moment**, the local Stoner condition, and the leading term for
Rh in FeRh and for the ``|\boldsymbol h|^2`` response in FeGe. With ``(2, 1, 1)`` instead,

```math
\Psi^{(i)}(\boldsymbol e, \hat{\boldsymbol e}_i) \;\propto\;
  (\hat{\boldsymbol e}_i \cdot \boldsymbol e_j)(\hat{\boldsymbol e}_i \cdot \boldsymbol e_k)
  - \tfrac13\, \boldsymbol e_j \cdot \boldsymbol e_k ,
```

which *does* depend on the axis (larger when the axis aligns with the neighbours).

**(c) Counting channels.** The mark can sit on ``i``, ``j`` or ``k``, but ``j`` and ``k``
are exchanged by ``\sigma``, so the *orbits of the mark position* are ``\{i\}`` and
``\{j, k\}``: two pointed channels, not three. Their lengths are ``[H : H_i] = 1`` and
``[H : H_j] = 2`` and add up to ``N = 3``. In general **one ``N``-body energy orbit splits
into at most ``N`` pointed channels.**

**(d) The sum rule.** Build the ``(2, 1, 1)`` function for all three mark positions, set
each evaluation axis back to the configuration's own direction
(``\hat{\boldsymbol e}_i = \boldsymbol e_i``, …) and add them: the energy-type ``(2, 1, 1)``
function
``\sum_{\text{centre}}[(\boldsymbol e_c \cdot \boldsymbol e_j)(\boldsymbol e_c \cdot \boldsymbol e_k) - \frac13 \boldsymbol e_j \cdot \boldsymbol e_k]``
is recovered *exactly* (up to an ordering multiplicity). Adding the ``(0, 1, 1)`` functions
instead gives a ``\boldsymbol e_j \cdot \boldsymbol e_k``-type function of *fewer bodies*.
That degeneracy is why the energy side forbids ``l = 0``; on the pointed side the mark keeps
the decomposition apart, so each term is independent. **This is the precise sense in which
the moment expansion contains information the energy expansion cannot provide.**

**(e) One row of the regression.** One constrained configuration ``\boldsymbol e^{(c)}``
yields one design-matrix row *per marked atom*,

```math
y_{c,i} = \hat{\boldsymbol e}^{(c)}_i \cdot \boldsymbol M^{(c)}_i, \qquad
X_{(c,i),\varphi} = \Psi^{(i)}_\varphi\bigl(\boldsymbol e^{(c)}, \hat{\boldsymbol e}^{(c)}_i\bigr),
```

whereas the energy side has one row per configuration. Symmetry-equivalent sites (``j`` and
``k``) share the same columns, so covariance is satisfied by the *structure* of the
coefficients, not by a constraint added afterwards.

In summary: the mark is a label that (1) shrinks the symmetrization group ``H \to H_i`` and
(2) allows ``l = 0`` on the marked site; ``l_{\mathrm{mark}} = 0`` channels are axis-free
induced moments, ``l_{\mathrm{mark}} \ge 1`` channels are axis-dependent responses; there is
one channel per orbit of mark positions (``\le N``); and summing over marks returns an
energy-type function — of the same body order when ``l_{\mathrm{mark}} \ge 1``, and of lower
body order (the "new information") when ``l_{\mathrm{mark}} = 0``.

## The pointed expansion in general

### Marked clusters and the subgroup ``H_i``

Distinguish one site ``i`` of a cluster ``S`` as the mark. The ``G``-orbit of the marked
cluster ``(i, S)`` is a *pointed orbit*, with stabilizer

```math
H_i \;=\; \{\, h \in H : \pi_h(i) = i \,\} \;\subseteq\; H,
```

where ``\pi_h`` is the site permutation ``h`` induces on the cluster. The mark positions
fall into orbits ``O_1, \dots, O_r`` under ``H``, of lengths ``[H : H_{i_k}]`` summing to
``N``; the number of pointed channels is the number of orbits ``r \le N``. Concretely, for a
pair ``(i, j)`` the mark on ``i`` and the mark on ``j`` give "the effect of ``j`` on the
moment of ``i``" and "the effect of ``i`` on the moment of ``j``" **separate coefficients**
— for an Fe–Ge pair, the Fe-marked and Ge-marked columns are different channels. On B20 FeGe
the Fe site symmetry is ``C_3``, and the single nearest-neighbour pair orbit of the energy
basis splits into two pointed orbits (measured: 12 members times 2 channels, covariance
residual ``1.3 \times 10^{-14}``).

### The mark as a decoration

The mark must label a site without contributing a value, but the projection engine only
knows products of per-site factors. The solution is to add one more per-site factor whose
value is always 1 or 0. Let ``\boldsymbol u`` be an *indicator field* that is
``\hat{\boldsymbol x}`` on the marked atom and ``\boldsymbol 0`` elsewhere, and use as the
mark factor the ``(k, l) = (1, 0)`` Racah solid harmonic of that field,

```math
\text{mark factor} \;=\; |\boldsymbol u_b|^2 R_{00}(\boldsymbol u_b) \;=\; |\boldsymbol u_b|^2
\;=\; \begin{cases} 1 & b = i \\ 0 & b \ne i, \end{cases}
```

direction-independent because ``R_{00} \equiv 1``. In the code this is the displacement
decor `SiteDecor(disp = (1, 0))` of the decor SALC engine, which is why the pointed basis
needs no projection code of its own. The factor plays two roles:

1. **As a label.** The marked site carries a *different* decoration from the others, so the
   assignment set the projector works on records which site is marked — the mark breaks the
   permutation symmetry automatically, exactly as ``\sigma`` became unusable in the example.
2. **As a selector at evaluation.** Every term carries exactly one mark factor, so a member
   not marked at the row's atom is *exactly* zero. Building the row ``(c, a)`` means raising
   the indicator field on atom ``a`` and nothing else.

### Why the same projector is correct (Frobenius reciprocity)

The engine's Reynolds projector averages over the full stabilizer ``H``, which moves the
mark around; the requirement stated above is symmetrization under ``H_i``. The two agree,
and that is what licenses reusing the energy-side engine unchanged.

The intuition: inside the space the projector acts on, collect the assignments whose mark
sits on site ``i`` into a subspace ``V_i``. An element of ``H`` moves the mark, so it maps
``V_i`` to ``V_{\pi_h(i)}`` — ``H`` permutes the blocks ``\{V_i\}`` over the orbit ``O`` of
mark positions. The ``V_i`` component of an ``H``-averaged vector is therefore nothing but
an ``H_i``-average, and the remaining components are copies of it carried around the orbit
by ``H``. An ``H``-invariant vector is completely determined by its ``V_i`` component, and
there are as many of them as there are ``H_i``-invariant vectors in ``V_i``.

Formally, when the decoration multiset contains exactly one mark and
``V_O = \bigoplus_{i \in O} V_i`` is one permutation orbit, ``H`` permutes the blocks
transitively, the stabilizer of ``V_i`` is exactly ``H_i``, hence
``V_O \cong \operatorname{Ind}_{H_i}^{H} V_i`` and Frobenius reciprocity gives

```math
\operatorname{Inv}_H(V_O) \;\cong\; \operatorname{Inv}_{H_i}(V_i) \qquad (i \in O),
```

*including dimensions*, the isomorphism being restriction to the ``V_i`` component. When the
mark positions fall into several orbits,
``\operatorname{Inv}_H(V) = \bigoplus_O \operatorname{Inv}_{H_{i_O}}(V_{i_O})`` (the engine
projects one permutation orbit at a time, so this is its normal mode of operation). The
consequence worth holding on to: this is *not* "build the energy basis, then split it". One
``H``-invariant vector contains the terms of *all* mark positions simultaneously, and the
row's atom selects one at evaluation — so the physical requirement that all
symmetry-equivalent sites share one coefficient is guaranteed by construction, not imposed.

### The basis function and the expansion

Collecting the terms whose mark sits on site ``i``,

```math
\Psi^{(i)}_\varphi(\boldsymbol e, \hat{\boldsymbol e})
  \;=\; \sum_{\substack{m \in \mathcal O \\ \text{mark} = i}}\ \sum_t\ \sum_{\boldsymbol\mu}
  F^{(\varphi)}_{m,t}[\boldsymbol\mu]\;
  \underbrace{\bar Z_{l_{\mathrm{mark}} \mu_0}(\hat{\boldsymbol e})}_{\text{marked site's axis factor}}
  \prod_{j \in \text{env}(t)} \bar Z_{l_j \mu_j}\bigl(\boldsymbol e_{a_j}\bigr),
```

where ``\bar Z_{lm}`` are the RMS-normalized real harmonics (``\bar Z_{00} \equiv 1``, so
for ``l_{\mathrm{mark}} = 0`` the axis factor disappears, as in the ``(0, 1, 1)`` example)
and ``F`` is the folded coefficient tensor of the SALC. The moment expansion is

```math
m_i(\hat{\boldsymbol e}; \boldsymbol e) \;=\; \sum_\varphi c_\varphi\, \Psi^{(i)}_\varphi(\boldsymbol e, \hat{\boldsymbol e}),
```

the same shape as the energy expansion. The intercept is **not** a separate parameter: the
one-body, ``l_{\mathrm{mark}} = 0`` SALC is identically 1 and lives inside the basis, so a
per-orbit intercept ``\mu^{(0)}`` appears once per Wyckoff orbit automatically. This is why
the moment fit neither centers its columns nor adds a global intercept (a regularized
estimator shrinks ``\mu^{(0)}`` like any other column — choose it deliberately).

### Selection rules

1. **Time reversal.** The target ``y = \hat{\boldsymbol e} \cdot \boldsymbol M`` is even
   under ``\boldsymbol e \to -\boldsymbol e``, so the total spin rank — **including the
   mark's own ``l_{\mathrm{mark}}``** — must be even. Odd mark ranks pair with odd
   environment content.
2. **Isotropy.** Without spin–orbit coupling only the total spin rank ``L_S = 0`` survives
   (`isotropy = true`, the analogue of the energy basis' ``L_f = 0`` screen).
3. **``l_{\mathrm{mark}} \ge 0`` on the mark, ``l_j \ge 1`` in the environment.** The energy
   side's ban on ``l = 0`` is lifted *on the marked site only*, because the mark factor
   keeps the site from vanishing. Physically these are the axis-independent,
   environment-determined channels; on FeGe the ``(0, 1, 1)`` star is the dominant term.
4. **Environment species must be sampled.** A species the downstream consumer does not
   sample may not appear as an *environment* factor (the basis could not be evaluated at run
   time), so `lmax_env` must be 0 for it. A *mark* may sit on any species — the induced
   moment of a species that never appears in the energy model is exactly what the channel is
   for (`marked`, default: every species).
5. **Mark-aware cutoffs.** A three-body star is cut on its two *mark–environment* bonds
   only; the environment–environment edge is free. Applying the energy side's all-edge
   compact-cluster rule at the same radius keeps 3 of the 15 nearest-neighbour star pairs of
   FeGe and costs 20–32 % in ``\sigma`` (measured).

## Relation to the energy expansion

### The site sum rule

Set each marked site's evaluation axis to its own configuration direction,
``\hat{\boldsymbol e}_i = \boldsymbol e_i``. Then, summed over all marked atoms of the
reference cell,

```math
\sum_i \Psi^{(i)}_\varphi(\boldsymbol e, \boldsymbol e_i) \;=\; \widehat\Phi_\varphi(\boldsymbol e),
```

where the right-hand side is an ``H``-invariant, spin-only cluster function with the same
cluster orbit and the same multiset of spin ranks. *Proof sketch:* raise the indicator field
on every atom at once; every mark factor becomes 1, every term has exactly one mark, so the
evaluation equals the sum over mark positions; a uniform-amplitude field is ``G``-invariant
and the mark factor depends on ``|\boldsymbol u|`` only, so what remains is a
``G``-invariant function of the spins alone. The identity can be checked independently of
the implementation (design rows summed over marked atoms against a uniform-field SALC
evaluation), which makes it a physical-identity gate.

The two cases of the rule are the general form of the example's step (d):

- ``l_{\mathrm{mark}} \ge 1``: ``\widehat\Phi_\varphi`` lies in the span of the energy SALCs
  with ``\boldsymbol l = (l_{\mathrm{mark}}, l_j, \dots)``. The pointed basis is the energy
  basis decomposed by "which site is myself".
- ``l_{\mathrm{mark}} = 0``: the marked site drops out and ``\widehat\Phi_\varphi`` is an
  energy function of *lower body order*. The energy side forbade ``l = 0`` precisely to
  avoid this degeneracy; on the pointed side the mark retains the position, so the
  decomposed terms are independent and only become degenerate once ``\sum_i`` is taken.

### Coefficients are related, but not derivable

``E`` and ``m`` are two readouts of one adiabatic map, so the energy coefficients
``\{J_\alpha\}`` and the moment coefficients ``\{c_\varphi\}`` are not unrelated. A Landau
model with explicit magnitudes,

```math
\mathcal E(\{m_i\}, \boldsymbol e) = \sum_i \bigl(A m_i^2 + B m_i^4\bigr)
  - \sum_{ij} J_{ij}\, m_i m_j\, (\boldsymbol e_i \cdot \boldsymbol e_j),
```

defines ``m^*_i(\boldsymbol e)`` by stationarity, and the envelope theorem gives
``\partial E / \partial(\boldsymbol e_i \cdot \boldsymbol e_j) = -J_{ij}\, m^*_i(\boldsymbol e)\, m^*_j(\boldsymbol e)``
(ordered pairs). Left side from the energy fit, right side from the moment fit: a
**consistency diagnostic** between the two expansions. It is not an exact constraint (the
real energy is not of Landau form, and the higher orders are not small), and conversely
``\{c_\varphi\}`` cannot be derived from ``\{J_\alpha\}``. The two are fitted separately —
from the same constrained calculations, reading a second output. No additional DFT is
required.

### Where the ``l_{\mathrm{mark}} = 0`` channel comes from

Project the Landau picture onto the pointed channels with the local field
``\boldsymbol h_i = \sum_{j \in \mathrm{nn}(i)} \boldsymbol e_j``. The linear (Neumann)
response ``\sum_j \hat{\boldsymbol e}_i \cdot \boldsymbol e_j`` is the pointed pair
``(1, 1)``; the quadratic response splits exactly as

```math
(\hat{\boldsymbol e}_i \cdot \boldsymbol h_i)^2
  = \underbrace{\tfrac13 |\boldsymbol h_i|^2}_{(0,1,1)\ \text{channel}}
  + \underbrace{\sum_{j,k}\Bigl[(\hat{\boldsymbol e}_i \cdot \boldsymbol e_j)(\hat{\boldsymbol e}_i \cdot \boldsymbol e_k)
      - \tfrac13\, \boldsymbol e_j \cdot \boldsymbol e_k\Bigr]}_{(2,1,1)\ \text{channel}},
```

(using
``\langle (\hat{\boldsymbol e} \cdot \boldsymbol a)(\hat{\boldsymbol e} \cdot \boldsymbol b)\rangle_{\hat{\boldsymbol e}} = \frac13 \boldsymbol a \cdot \boldsymbol b``)
— the two functions of the example. On FeGe the ``|\boldsymbol h|^2`` form (the
environment-only local Stoner condition) matters more than the
``(\hat{\boldsymbol e} \cdot \boldsymbol h)^2`` form (``\sigma`` 0.035 vs 0.049 ``\mu_B``).

### Correspondence table

| | Energy ``E(\boldsymbol e)`` | Moment ``m_i(\hat{\boldsymbol e}; \boldsymbol e)`` |
|:--|:--|:--|
| Basis function | ``\Phi_\alpha(\boldsymbol e)`` | ``\Psi^{(i)}_\varphi(\boldsymbol e, \hat{\boldsymbol e})`` |
| Required symmetry | invariance | covariance |
| Symmetrization group | stabilizer ``H`` | mark-fixing subgroup ``H_i`` |
| Carrier | (ordering, path, ``M``) | the same, plus the mark position |
| Projector | ``P = \frac{1}{\lvert H \rvert}\sum_h U(h)`` | identical (Frobenius reciprocity) |
| Site factors | ``\bar Z_{l_i \mu_i}(\boldsymbol e_{a_i})``, ``l_i \ge 1`` | the same, plus the mark factor ``\lvert\boldsymbol u\rvert^2``, ``l_{\mathrm{mark}} \ge 0`` |
| Normalization | ``(4\pi)^{N/2}`` | ``(4\pi)^{n_{\mathrm{spin}}/2}`` (``n_{\mathrm{spin}}`` = number of spin factors) |
| Time reversal | ``\sum_i l_i`` even | ``l_{\mathrm{mark}} + \sum_j l_j`` even |
| No spin–orbit | ``L_f = 0`` | ``L_S = 0`` |
| Design-matrix row | configuration ``c`` | (configuration ``c``, marked atom ``i``) |
| Target | ``E^{(c)}`` | ``y_{c,i} = \hat{\boldsymbol e}_i \cdot \boldsymbol M_i`` |
| Intercept | ``j_0``, a separate parameter | the one-body ``l_{\mathrm{mark}} = 0`` column, per Wyckoff orbit |
| Cutoff | all edges compact | mark–environment bonds only (3-body) |
| Sharing | — | symmetry-equivalent sites share one coefficient |

### Explicit low-order forms (no spin–orbit)

With the addition theorem
``\sum_m \bar Z_{lm}(\boldsymbol a)\bar Z_{lm}(\boldsymbol b) = (2l+1) P_l(\boldsymbol a \cdot \boldsymbol b)``:

- **One body.** ``l_{\mathrm{mark}} = 0`` gives the constant 1 — the per-orbit intercept
  ``\mu^{(0)}``. (A single-site ``l_{\mathrm{mark}} = 2`` term has ``L_S = 2`` and is
  screened out without spin–orbit coupling.)
- **Two body.** ``L_S = 0`` forces ``l_{\mathrm{mark}} = l_j = l``, so
  ``\Psi^{(i)} \propto \sum_{j \in s(i)} P_l(\hat{\boldsymbol e}_i \cdot \boldsymbol e_j)``
  over a shell ``s(i)``: ``l = 1`` is the linear response, ``l = 2`` the ``P_2`` term. Same
  functional form as the energy's ``\sum_{ij} P_l(\boldsymbol e_i \cdot \boldsymbol e_j)``;
  the only differences are that ``\boldsymbol e_i`` has become the independent axis and that
  there is no sum over ``i`` (the sum rule).
- **Star ``(0, 1, 1)``.**
  ``\Psi^{(i)} \propto \sum_{(j,k)} \boldsymbol e_j \cdot \boldsymbol e_k``,
  axis-independent: the induced-moment / local-Stoner channel. For Rh in FeRh the
  free-moment magnitude
  ``|\boldsymbol V_{\mathrm{Rh}}|^2 = c^2\bigl(8 + 2\sum_{j<k} \boldsymbol e_j \cdot \boldsymbol e_k\bigr)``
  has exactly this form.
- **Star ``(2, 1, 1)``.**
  ``\Psi^{(i)} \propto \sum_{(j,k)}\bigl[(\hat{\boldsymbol e}_i \cdot \boldsymbol e_j)(\hat{\boldsymbol e}_i \cdot \boldsymbol e_k) - \frac13\, \boldsymbol e_j \cdot \boldsymbol e_k\bigr]``.

## The regression

### Target and design matrix

From one constrained configuration ``c``, one row per marked atom ``i``:

```math
y_{c,i} = \hat{\boldsymbol e}^{(c)}_i \cdot \boldsymbol M^{(c)}_{\mathrm{bare},i}, \qquad
X_{(c,i),\varphi} = \Psi^{(i)}_\varphi\bigl(\boldsymbol e^{(c)}, \hat{\boldsymbol e}^{(c)}_i\bigr), \qquad
\boldsymbol y \approx X \boldsymbol c .
```

``\boldsymbol M_{\mathrm{bare}}`` is the **bare integrated moment** (`moments_bare`; VASP's
`M_int`), distinct from the smoothed moment MW that supplies the *directions*: the ratio
``|\mathrm{MW}| / |\boldsymbol M_{\mathrm{bare}}|`` is configuration-dependent (FeGe:
``0.691 \pm 0.017``), so the two must not be mixed. The implementation evaluates ``\Psi`` by
substituting the **marked column only** of the configuration matrix by
``\hat{\boldsymbol e}_i``; members not marked at that atom vanish before they can read the
substituted column, so the substitution is exact, not an approximation.

### The evaluation-axis rule

One principle: the axis is **"the axis along which the constraint evaluated the adiabatic
map"**, and only its *source* depends on the constraint mode — the direction of the
converged moment for direction pinning (mode 4), the constraint axis for a transverse
penalty (mode 1). The axis is never *inferred* from the configuration: in mode-1 rows with
``\|\boldsymbol M\| \to 0`` the MW direction is meaningless, and using it shifts the
coefficients by a measured factor of 2.1.

### Why a signed projection

- ``|\boldsymbol M|`` is **non-analytic at zero** on an induced site (it behaves like
  ``c\,|\sum_j \boldsymbol e_j|``), which no polynomial basis can represent; ``m^2`` loses
  the sign and ``\ln m`` diverges.
- A mode-1 penalty ``E_p = \lambda \sum |\boldsymbol M_\perp|^2`` punishes the *transverse*
  component only, so the projection can be negative and genuinely is (FeRh Rh: 3848 of 7744
  rows have ``y < 0``). The signed projection is the only analytic readout that keeps the
  information of constrained data.
- **Equivalence to vector regression.**
  ``\|\boldsymbol M_i - m\,\hat{\boldsymbol e}_i\|^2 = (y_i - m)^2 + \text{const}``, so the
  projected least squares has the same normal equations as a vector-reconstruction loss. A
  naive vector regression of ``\boldsymbol M`` itself, however, sees observations collapsed
  onto ``(\boldsymbol V \cdot \hat{\boldsymbol e})\hat{\boldsymbol e}`` and, with isotropic
  sampling ``\mathbb E[\hat{\boldsymbol e} \hat{\boldsymbol e}^\top] = I/3``, attenuates
  every coefficient by ``1/3``.

### The row gate

Rows on which the adiabatic map is not well defined — a moment far off its axis — are
dropped. In a cancellation-free form,

```math
\boldsymbol M_\perp = \boldsymbol M - y\,\hat{\boldsymbol e}, \qquad
g = \frac{\|\boldsymbol M_\perp\|^2}{|\boldsymbol M|} = |\boldsymbol M|\sin^2\theta \;\le\; \varepsilon\ [\mu_B].
```

``g`` is the transverse residual weighted by ``|\boldsymbol M|`` and is *proportional* to
the penalty field projected on ``\hat{\boldsymbol M}``
(``-\boldsymbol B \cdot \hat{\boldsymbol M} / 2\lambda`` with
``\boldsymbol B = -2\lambda \boldsymbol M_\perp``) — proportional, not equal, because the
penalty acts on MW while ``g`` is computed from the bare moment. The tolerance
``\varepsilon`` (`gate_eps`) has no default: it is a statement about the data and is yours
to make. The per-orbit survival is reported, and a fit whose rows of some orbit fall below a
floor is *refused* rather than allowed to claim that orbit's columns.

### The induced moment vector

Writing
``m_i(\hat{\boldsymbol e}; \boldsymbol e) = a_0(\boldsymbol e) + \hat{\boldsymbol e} \cdot \boldsymbol V_i(\boldsymbol e) + O(l \ge 2)``,
the ``l = 1`` sector ``\boldsymbol V_i`` is the induced moment vector. It is read off by an
antisymmetrized six-point evaluation,

```math
V_\alpha = \tfrac12\bigl[m(+\hat{\boldsymbol\alpha}) - m(-\hat{\boldsymbol\alpha})\bigr],
\qquad \alpha = x, y, z,
```

(with [`predict_moment`](@ref) and its `axes` keyword). All even-``l_{\mathrm{mark}}``
columns cancel exactly, so this is exact for ``l_{\mathrm{mark}} \le 2``; an
``l_{\mathrm{mark}} = 3`` column would leak an ``l = 3`` component. A three-point evaluation
along ``\hat{\boldsymbol x}, \hat{\boldsymbol y}, \hat{\boldsymbol z}`` is *not* enough:
with the default ``l_{\mathrm{mark}} = 2``,
``m(\hat{\boldsymbol\alpha}) = a_0 + V_\alpha + Q(\hat{\boldsymbol\alpha})`` and the even
orders do not separate. "Direction ``\parallel \boldsymbol V``, magnitude
``|\boldsymbol V|``" is exact only for isotropic linear response; in general
``|\boldsymbol M_{\mathrm{free}}| = \max_{\hat{\boldsymbol e}} m_i(\hat{\boldsymbol e})``.

### Periodic resolvability

On cell-periodic data some pointed columns **vanish identically** and some are
**structurally dependent** (Wigner–Seitz boundary ties, images folding onto the same
reference atom) — the pointed form of [Periodic resolvability](resolvability.md). Because
each row's identity (marked atom, the mark's ``(l, \mu)`` on the independent variable
``\hat{\boldsymbol e}``, the environment factors on ``\boldsymbol e``) can be expanded
symbolically, the verdict is reached **without looking at any data**: vanishing columns are
frozen to exact zero, dependent directions are disclosed as min-norm representatives, and
the decision is a relative cut on the symbolic expansion, never "this value happens to be
small". A basis in which two environment factors land on the same reference-cell atom is
refused as unclassifiable, loudly ([`moment_resolvability`](@ref)).

## Scope and caveats

1. **Reference geometry only.** The expansion is ``m_i(\boldsymbol e)`` at the clamped-ion
   reference structure; this package carries no displacements, and mixing geometries would
   silently mix two structures.
2. **A true collapse to zero is not representable — and passes the gate.** Where the
   adiabatic map becomes singular (``m \to 0`` on a site that is not an induced site; FeGe
   at high constraint temperature shows a few such rows), the signed projection cannot
   represent it, and ``g = |\boldsymbol M|\sin^2\theta \to 0`` as ``|\boldsymbol M| \to 0``
   regardless of the angle, so such rows pass. The gate removes *large* moments off their
   axis, nothing else. A minimum-moment gate is deliberately absent (it would delete nearly
   every row of an induced species); the limitation is reported, not hidden (on FeGe, 0.24 %
   of the rows carry 14 % of the squared error).
3. **Extrapolation.** An induced site trained on an antiferromagnetic host covers
   ``|\hat{\boldsymbol e} \cdot \boldsymbol h_1| \lesssim 3.6`` for eight neighbours; the
   ferromagnetic value 8 is an extrapolation. [`moment_local_field`](@ref) /
   [`moment_coverage`](@ref) monitor the local-field coordinates
   ``(\|\boldsymbol h_1\|, \hat{\boldsymbol e} \cdot \hat{\boldsymbol h}_1)`` so an
   extrapolation is named before it is trusted.
4. **Neighbour enumeration is minimum-image.** An environment slot landing on the marked
   atom's own periodic image would let the marked-column substitution rewrite an environment
   factor; the neighbour list excludes ``i = j``, so this cannot happen by construction, and
   the two-factors-on-one-atom case is refused by the resolvability gate.
5. **Gauge (mode 1).** The penalty is even in ``\hat{\boldsymbol e}``, so
   ``\boldsymbol M(-\hat{\boldsymbol e}) = \boldsymbol M(\hat{\boldsymbol e})`` and the
   physical map
   ``m(\hat{\boldsymbol e}) = \hat{\boldsymbol e} \cdot \boldsymbol M(\hat{\boldsymbol e})``
   is odd: a global flip of the axes is a gauge, and only the even-``l_{\mathrm{mark}}``
   columns depend on that convention. The readers' sign-consistency gate is invariant under
   the flip, so it does not fix the gauge — its job is to detect stale or substituted axes.
   Gauge mixing is read off the per-orbit anti-alignment count instead (FeRh: Fe 0 of 7744 =
   single gauge; Rh 3833 of 7744 = the axis is an independent input and the sign is
   physics).

## Normalization bookkeeping

- ``Z_{00} = (4\pi)^{-1/2}`` and ``\bar Z_{00} = \sqrt{4\pi}\, Z_{00} = 1``.
- Racah:
  ``R_{00}(\boldsymbol u) = \sqrt{4\pi}\,|\boldsymbol u|^0 Z_{00}(\hat{\boldsymbol u}) = 1``
  exactly, so the mark factor is ``|\boldsymbol u|^{2k}`` itself.
- Energy: every one of the ``N`` sites is a spin factor, prefactor ``(4\pi)^{N/2}``.
- Pointed: prefactor ``(4\pi)^{n_{\mathrm{spin}}/2}`` with ``n_{\mathrm{spin}}`` the number
  of *spin* factors — an ``l_{\mathrm{mark}} = 0`` mark carries none and is not counted; an
  ``l_{\mathrm{mark}} \ge 1`` mark carries one.
- The two conventions (``Z_{00} = (4\pi)^{-1/2}`` and ``R_{00} = 1``) are the same "RMS = 1
  under the uniform measure" rule in two guises; a ``\sqrt{4\pi}`` is lost only when they
  are mixed, and the implementation pins ``R_{00} \equiv 1`` with an assertion.

## Symbols to code

| This page | Implementation |
|:--|:--|
| ``\bar Z_{lm}`` | `Harmonics.Zlm` (RMS normalization) |
| ``R_{lm}`` | `SolidHarmonics` (Racah convention) |
| mark factor ``\lvert\boldsymbol u\rvert^2`` | `SiteDecor(disp = (1, 0))` |
| ``\Psi^{(i)}_\varphi`` | the SALCs of a [`MomentBasis`](@ref), evaluated with the mark selector |
| truncation (``l_{\mathrm{mark}}``, ``l_j`` caps, cutoffs, `sampled`, `marked`) | [`MomentSpec`](@ref) |
| ``X``, ``y``, the axis rule, the gate | [`MomentDataset`](@ref) (`gate_eps`, `coverage_floor`) |
| ``\{c_\varphi\}`` | `fit(MomentFit, …)`, [`MomentModel`](@ref) |
| ``m_i(\hat{\boldsymbol e}; \boldsymbol e)`` | [`predict_moment`](@ref) (`axes` = the evaluation axes) |
| vanishing / dependent columns | [`moment_resolvability`](@ref) |
| coverage, residual profile, simple-feature floor | [`moment_local_field`](@ref), [`moment_band_profile`](@ref), [`moment_simple_floor`](@ref) |

## References

1. R. Drautz and M. Fähnle, *Phys. Rev. B* **69**, 104404 (2004) — the spin-cluster
   expansion of the energy.
2. M. Uhl and J. Kübler, *Phys. Rev. Lett.* **77**, 334 (1996); A. V. Ruban, S. Khmelevskyi,
   P. Mohn, and B. Johansson, *Phys. Rev. B* **75**, 054402 (2007) — configuration-dependent
   (longitudinal) moments in itinerant magnets, the physics the adiabatic map encodes.
3. J.-P. Serre, *Linear Representations of Finite Groups* (Springer, 1977), §7 — Frobenius
   reciprocity.
