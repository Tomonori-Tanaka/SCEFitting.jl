# Design: cross-orbit alias groups

Status: landed (2026-09-11)

## Summary

**設計列 = 結んだ列**。エイリアス群の SALC は設計行列で 1 列に重み付き加算され、
推定器・罰則計量・群ラベル・選択(`select_fit` / `select_support`)・`refit` の support・
CV・GCV/`effective_dof`・`residuals_torque` はすべて列空間(`n_columns(basis)`)で動く
(折った `X_E`/`X_T` を読むだけなので自動的に列空間)。SALC 空間(`n_salcs(basis)`)へは `SCEPredictor(f)` で
展開する:`jphi_salc[j] = w_j · jphi_col[col_of[j]]`。無エイリアス基底では
`col_of = 1:p`、`w = 1`、折り畳みをスキップするので既存挙動はバイト同一。

結び方は「列加算・各軌道の係数 = 結んだ係数」(重み `w_j = ±1` が通常)。
代表列 + `J/k` は列比が 1 のときしか関数を再現しない(3 重タイで点群が 2 像だけ融合すると
orbit A の列は B の 2 倍;`2J_A + J_B = S` が決まり、加算タイなら `J_A = J_B = S/3` で
ボンドあたり等分が自動、代表列方式は `1.5S` を返す)。一般の重みは
**メンバーあたりのテンソルノルム**で定義する:`ν_j = raw_j / √n_members(j)`、
`w_j = sign_j · ν_ref / ν_j`、`sign_j = sign(c_j)`、`c_j = ⟨v_ref, v_j⟩/‖v_ref‖²` は
最小二乗の比例係数(輸送コピーなら `|w_j| = 1` を厳密に取り、`|w ∓ 1| < 1e-10` で ±1 に
スナップ)。`raw_j` は `_function_vector` が返す**非集約**ノルム `√(Σ_members Σ_terms ‖T‖²)`
なので `raw_j/√n_members` は正確にメンバーあたりのテンソルノルム(`n_members` は
`_canonicalize_members` 後の数)。スナップされない重みは「同一表現の輸送コピー」という前提が
崩れた証拠であり、有効な重みではない(テストがピンする)。これが「ボンドあたり等分」の定義。

**却下した代替**:(i) SLCE 流の凍結 — 決定済みの和まで捨て、タイル MC で結合が欠落。
(ii) 基底で併合(引き継ぎ案)— 軌道が消えて開示不能、スパイラル拡張で作り直し。
(iii) pinv / ridge に任せる — 同一列は等分されるが、`salc_groups` が別群なので
`select_support` が片方だけ落とせば 1:0、多重度が違えば罰則計量が配分を決める。
ハードタイの意味は数値ではなく「選択ラベルの統合・再現性・警告の沈黙」。
(iv) `_assemble_problem` でだけ結ぶ(X·T)— 推定器の列ごとデータ(metric、column_groups)を
全呼び出し点で縮約する必要があり、`select_fit`/`cross_validate` の最も危険な経路に触る。

## Module layout

| Target | Change |
|---|---|
| `src/basis/aliases.jl`(新) | `AliasGroup`、`_ColumnTies`、`alias_groups(basis)`、`column_ties(salc_basis, crystal; alias_rtol)`、`_ALIAS_RTOL = 1e-6`、`@info` 整形 |
| `src/sce/model.jl` | `SCEBasis` に `ties::ColumnTies` フィールド(内部コンストラクタで計算、キーワード `alias_rtol`);外部コンストラクタに `alias_rtol` キーワード;`n_columns(basis)`;`SCEPredictor` に `split::Vector{Symbol}`;`SCEDataset` docstring |
| `src/fitting/design.jl` | `_design_energy` / `_design_torque` の末尾で `_fold_columns(X, ties)` |
| `src/fitting/fit.jl` | `SCEPredictor(f)` の展開;docstring(`jphi` は列空間) |
| `src/fitting/diagnostics.jl` | `coef` docstring(fit = 列空間、predictor = SALC 空間) |
| `src/fitting/selection.jl` | `salc_groups`(列ごと、群併合)、`group_costs`(`labels[col_of[j]]`)、`penalty_metric`(列ごとに Σ w_j Φ_j を評価)、`select_fit` の長さ検査を `n_columns` |
| `src/sce/coeftable.jl` | `alias_group::Int`(0 = なし)と `split::Symbol` 列;`coeftable(f)` は `SCEPredictor(f)` 経由 |
| `src/sce/introspect.jl`, `src/sce/bilinear.jl`, `src/interop/sunny.jl` | `_warn_convention_split(model, what)` を先頭で呼ぶ(maxlog = 1) |
| `src/io/persist.jl` | `PERSIST_SCHEMA_VERSION = 7`、`couplings[*].split`、読み込み時 v < 7 は `:legacy`(群内)/`:free` |
| `src/SCEFitting.jl` | `include("basis/aliases.jl")`(salcbasis の後、model の前);`public alias_groups, AliasGroup, n_columns` |
| `docs/src/theory/resolvability.md` | 新節 "When symmetry does not fuse the tie" + "How this package handles it"(P4₂/mnm 10 組の実例、規約の意味、q·ΔR = π 方向は規約依存) |
| `docs/src/api.md`, `SPEC.md`, `CHANGELOG.md`, `CLAUDE.md` | API 追加、SALC basis 節・fitting 節、Unreleased、coupled site + 発散台帳行 |
| `test/unit/test_aliases.jl`(新) | 下記 |

## API

```julia
struct AliasGroup
    salcs::Vector{Int}                 # SALC indices (design-order), ≥ 2, from ≥ 2 orbits
    orbit_ids::Vector{Int}
    body::Int
    atoms::Vector{Int}                 # representative member atoms of salcs[1]
    delta_shifts::Vector{SVector{3,Int}}   # per orbit: lattice-vector offset of its
                                       # representative member relative to orbit 1
                                       # (per-site difference collapsed to one vector for pairs)
    distance::Float64                  # pair separation (Å); NaN for N ≥ 3
    weights::Vector{Float64}           # w_j per SALC in `salcs` (±1 for transported copies)
    kind::Symbol                       # :proportional (tied) | :span_collapsed | :unequal_norm (reported only)
end

alias_groups(basis::SCEBasis) -> Vector{AliasGroup}     # public, unexported
n_columns(basis::SCEBasis) -> Int                       # public, unexported
SCEBasis(crystal, spec; ..., alias_rtol::Union{Nothing,Real} = 1e-6)
```

`_ColumnTies`(内部):`col_of::Vector{Int}`(長さ p)、`weight::Vector{Float64}`
(長さ p)、`columns::Vector{Vector{Int}}`(長さ p_t)、`groups::Vector{AliasGroup}`、
`trivial::Bool`。

`coeftable` 行:`(body, orbit_id, decors, L_S, Lf, block, J, alias_group, split)`。
`split ∈ (:free, :convention, :legacy)`。

`split` の決まり方:
- `SCEPredictor(basis, j0, jphi)` / 4 引数位置形(手設定モデル)= 全 `:free`(群内で異なる値を
  置くことは合法。合成真値の生成に使う)。
- `SCEPredictor(::SCEFit)` = 結んだ群の SALC は `:convention`、他は `:free`。`refit` で結んだ列が
  support から落ちた場合も同じ(係数 0 の `:convention`)。
- 永続化: v7 は `couplings[*].split` を保存・検証(未知の値・欠落は refuse)。v < 7 は
  `_split_status(basis.ties)` で群内 `:legacy`、他 `:free`。`_basis_doc` は ties を書かない。
  `_basis_from_doc` → 4 引数内部コンストラクタが既定 `alias_rtol` で再計算する。
  4 引数内部コンストラクタは `model.jl` / `persist.jl` のほか test 3 箇所から呼ばれる
  (キーワード既定値で互換)。

## Types and conventions

- 数値規約は不変。エイリアスの無い基底では `X_E`/`X_T`/`jphi` がバイト同一。
- `n_salcs(basis)` の意味を「SALC 基底関数の数」に限定し、「設計列数」は `n_columns`。
- `SCEFit.jphi` は列空間、`SCEPredictor.jphi` は SALC 空間。`coef(f)` は列空間の
  ベクトル(docstring で明記)。`coeftable(f)` は SALC 空間(展開済み)。
- 検出閾値 `_ALIAS_RTOL = 1e-6`:実測(Nd₂Fe₁₄B l02)で真のエイリアス 10 組のうち 4 組の
  比例残差が 1.49e-8(ゲージ固定の丸め)、`_AGG_DEP_RTOL = 1e-8` では取りこぼす。
  MnTe の非エイリアス生存者は O(1e-2)。閾値帯 `1.5e-8 ≪ 1e-6 ≪ 1e-2`。
- SLCE.jl との発散:SLCE は `unresolvable_columns` で共有軌道を丸ごと凍結(和も捨てる)、
  本パッケージは結んで等分(規約として開示)。台帳行を追加。

## Impact on coupled sites

- [x] `Harmonics` / `AngularMomentum` / SALC projection / normalization tests: 触らない。
- [x] `SALCKey` order / design columns / TOML persistence / `coeftable`: 設計列が
      SALC より少なくなり得る(`n_columns ≤ n_salcs`);TOML v7 に `split`;`coeftable` に 2 列。
- [x] Energy kernel ↔ torque kernel: 両設計行列を同じ `_fold_columns` で折る(`w` 共通)。
- [ ] Decor engine ↔ pointed moment basis ↔ `momentfit.jl` doors: `MomentBasis` は別経路、
      `_function_vector` が装飾 SALC を拒否するので、装飾 SALC を含むチャネルは検出をスキップ。
- [ ] Readers' `zero_moment_atol` ↔ dataset doors: 無関係。
- [x] Upstream divergence ledger (SLCE.jl): 行追加。
- [ ] `.claude/agents/` references: モジュール名・Makefile 不変。
- [x] `SPEC.md` / `docs/src/api.md` updates.

## Test strategy

`test/unit/test_aliases.jl`(すべて `NoSymmetry()` = 群が像を関係づけない状況を作る):

1. **P1 3 原子(a/2 タイ)**:a = 3 立方、原子 1 (0,0,0)、2 (0.5,0.2,0.1)、3 (0.25,0.3,0.42)、
   `lmax = 1`、pair cutoff 2.6。オラクル:手計算で原子 1–2 の ±x 像は等距離、2 軌道の
   SALC 値は `evaluate_salc = 2√3·(e₁·e₂)`(v0 done-line と同じ規格化)で一致。
   `alias_groups` = 1 群(`kind == :proportional`、`weights == [1, 1]`)、
   `n_columns == n_salcs − 1`。合成データ `E = J₊ φ₊ + J₋ φ₋`(`SCEPredictor(basis, 0, [J₊, J₋, …])`
   で生成)の OLS が `jphi_salc = (J₊+J₋)/2` を両 SALC に返す。`predict_energy`/
   `predict_torque` は生成モデルと一致(訓練セル上の完全再現)。
2. **正方 P セル Δz = c/2**:同上を c 方向で(`delta_shifts` が `(0,0,±1)` を指す)。
3. **合成回収(配分はセルが決める)**:1 の結晶の 2×1×1 超胞では ±x 像が別原子になり
   エイリアス 0 群;真値 J₊ ≠ J₋ を置いて乱数配置のエネルギー+トルクを生成し、OLS が
   両方を別々に回収(rtol 1e-8)。同じ真値で 1×1×1 周期の配置を生成(超胞エネルギー / 2)
   すると 1×1×1 基底の OLS は各軌道に (J₊+J₋)/2。
4. **変異**:`alias_rtol = nothing` で群 0、`n_columns == n_salcs`、OLS がランク警告
   (`@test_logs (:warn, r"rank deficient")`)。
5. **閾値帯ピン**(change detector と明記):真の群の最大比例残差 < `_ALIAS_RTOL`。
6. **persist**:v7 往復で `split` が保存・復元;v6 文書(`split` 無し)を読むと群内は `:legacy`。
7. **列空間の整合**:`salc_groups(basis)` の長さ = `n_columns`、エイリアス群の 2 軌道が同一
   ラベル;`penalty_metric(basis)` の長さ = `n_columns`、結んだ列の値 = 折った列の分散
   (`Var[Φ₊+Φ₋] = 4 Var[Φ₊]` を同一列で確認);`GroupAdaptiveRidge(basis; …)` +
   `select_fit` が走る。
8. **回帰**:既存スイート(bcc Fe、FeRh、MnTe 系フィクスチャ)で `alias_groups` が空、
   `make test-pin` 不変。
9. **部分融合タイ(Cm、手組み空間群)**:a=b=4, c=5、原子 (0,0,0), (0.5,0.5,0.25), (0.3,0.3,0.5)、
   対角鏡映 x↔y。1–2 ペアの 4 隅像のうち鏡映が 2 つを融合 → 軌道メンバー数 [1,2,1]。
   オラクル = ボンド重み付き平均 Σ nᵢJᵢ / Σ nᵢ(真値 [0.9,0.2,0.1] で 0.35;軌道平均 0.4 と異なる)。
   「ボンドあたり」と「軌道あたり」を分ける唯一のゲート(数値レビュー major 1)。
   `alias_rtol` の上限 `_ALIAS_RTOL_MAX = 1e-2`(`tie_tol` と同型、major 3)。スナップしない重みは
   `:unequal_norm` として結ばない(major 2)。

実系(Nd₂Fe₁₄B l02 = 10 組、結んだ X_T 満階数 169/169;l06 = 30;l044_c4.0 = 20)は CI 外。
`~/jijs/scefit/2-14-1/nd2fe14b/1x1x1/scefit/` のスクリプトで確認し、結果を CHANGELOG に記録。

## Risks and open items

- 異方チャネル(`Lf > 0`)では列が比例しない(張る空間の一致)。`:span_collapsed` の報告のみで
  結ばない;OLS 警告が残る。等分の定義自体が要議論(SLCE の凍結が正しい可能性)。
- `salc_groups` の併合で MC 選択の粒度が粗くなる(2 軌道の全チャネルが 1 群)。
  `group_costs` docstring の「粗い分割可」の範囲内。
- `penalty_metric` の結んだ列は分散を直接計算するので λ の意味が旧フィットと変わる
  (エイリアスのある基底のみ)。
- `dof(f)` は列数 + 1(結んだ後)。
- `alias_rtol` は永続化しない(`tie_tol`/`images` と同じ流儀)。ロード時は既定値で再計算。
  既定値以外で構築したモデルはロード後に群・`split` が変わり得る(記録)。
- `_fold_columns` は `_design_energy`/`_design_torque` の末尾(ホットパス)に入るが、無エイリアス
  基底では同じ行列オブジェクトを返す恒等パスなので測定対象の性能は不変。エイリアスがある基底での
  折り畳みは O(n·p) の列加算で設計行列の評価より桁で小さい。BENCH_LOG 記録は不要と判断。
- `SCEFit.jphi`/`SCEPredictor` のフィールド変更は「`SALCKey`/`SALC` field surface ↔ ALL test
  environments」の coupled site に当たる: `test/sunny`、`test/glmnet`、`test/oracle`、
  `test/parity`、`examples/*.jl`、`bench/*.jl`、`test/pin/payload.jl` は全て無エイリアス
  フィクスチャ(bench_nd2fe14b の既定 cutoff 4.0 Å は a/2, c/2 未満)なので有効なまま。
