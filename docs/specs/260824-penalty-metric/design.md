# Design: 罰則計量の正当化と λ 選択（エネルギー / モーメント両チャネル）

Status: draft (2026-08-24) — spec-reviewer 第 3 回の blocker 5 / major 10 を反映済み

## Summary

新しい推定器型は作らない。**列ごとの罰則スケール `m_j`（計量）を導入し、それ 1 つで
2 つの仕事を片付ける**。

1. **スケール不変にする**（両チャネル）— `m_j` を基底固有の参照ノルムにする。
2. **μ₀ 切片を罰則から外す**（モーメント側）— `m_j = 0` にする。`m_j = 0` は
   「この列は罰則を受けない」の意味。
3. **λ 選択**（モーメント側）— `cross_validate` / `effective_dof` / `gcv`。

②を群重み `v_g = 0` でやる案は**採らない**。計量で済むなら `Ridge` / `AdaptiveRidge` /
`GroupAdaptiveRidge` の**3 種すべてが同じ機構**で切片を外せるが、群重みは
`GroupAdaptiveRidge` にしか無いので `fit(MomentFit, ds, Ridge(λ))` が黙って μ₀ を縮める
穴が残る。`group_weights` の `≥ 0` 緩和は**不要**になったのでスコープから外す。

### 設計の要 1: 設計行列は触らない

計量を「設計列を標準化する」で実装してはいけない。設計を書き換えると `OLS` の解も
（数学的には同一でも）浮動小数点で動き、`fit` ↔ `refit` のバイト同一性
（`_assemble_problem` を共有している理由そのもの）が崩れる。罰則側だけに入れる。

### 設計の要 2: 計量は重みの**分母**に入れる（第 2 回レビュー B1 / B2）

罰則を `λ Σ_j m_j w_j β_j²` と書くだけでは**適応系がスケール不変にならない**。
列を `c_j` 倍（`β_j → β_j/c_j`、`m_j → c_j² m_j`）したとき不変であるためには
`w'_j = w_j` が要る:

| 推定器 | `w_j` | 再スケール後 | 不変 |
|---|---|---|---|
| `Ridge` | `1` | `1` | ✔ |
| `AdaptiveRidge`（現行式） | `1/(β_j² + ε)` | `1/(β_j²/c_j² + ε)` | ✘ |
| `GroupAdaptiveRidge`（現行式） | `v_g/(‖β_g‖² + p_g ε)` | `v_g/(Σβ_k²/c_k² + p_g ε)` | ✘ |

群内が一様な `c` でも `‖β_g‖²/c²` になるので破れる。⇒ **計量を分母に入れる**:

```
AdaptiveRidge:        w_j = 1 / (m_j·β_j² + ε)
GroupAdaptiveRidge:   w_j = v_g / (Σ_{k∈g} m_k·β_k² + p_g·ε)
罰則対角              D_j = m_j · w_j          (これが `_gar_weights!` の返り値)
```

`Σ m β²` が再スケール不変量なので `w'_j = w_j` が厳密に成立し、`D'_j = c_j² D_j` で
罰則項全体が不変になる。`ε` / `p_g·ε` を再スケールしなくてよいのはこのため。

**ただし分母だけでは足りない。** 適応リッジは非凸代理の不動点反復なので、目的関数の
不変性は収束先の不変性を意味しない。**3 つとも可換にする**:

1. **cold start も罰則対角で組む** — `AdaptiveRidge` の iteration 0（`estimators.jl:446`）は
   `Symmetric(XtX + λ·I(p))` なので `λ·Diagonal(m)` に、`_solve_gar` の cold start
   （`:504-508`、`w[j] = group_weights[column_groups[j]]`）は `D_j = m_j·v_g` に。
2. **停止則を計量座標で書く** — 現行 `rel = max_j|Δβ_j| / max(max|β_new|, eps)`
   （`:452-453`, `:509-513`）は `β_j → β_j/c_j` で不変でない。`√m_j β_j` が再スケール
   不変量なので

   ```
   rel = max_{j: D_j>0} √m_j·|Δβ_j| / max( max_{j: D_j>0} √m_j·|β^new_j|, eps )
   ```

   罰則列に限ることで μ₀ 支配の問題（第 1 回 M2）も同時に解ける。
3. `λ == 0` の OLS ルーティングは計量と無関係なのでそのまま。

**同じ 1 行が group-L0 の意味論も守る。** `docs/design-notes.md` §13（L437-441）は
不動点の性質を明示している:

> `λ·v_g·‖β_g‖²/(‖β_g‖² + p_g ε) → λ·v_g`: the fixed multiplier `v_g` **is** the
> group-L0 weight

計量を分母の**外**に掛けると収束時の群あたり罰則は `λ·v_g·⟨m⟩_g` になり **`v_g` が
group-L0 重みでなくなる** — `select_fit` の Pareto も θ 掃引も §13 の導出もその上に
乗っているので、「`cost_weights` を設計どおりに効かせる」という本 spec の目的と逆行する。
分母に入れれば `λ v_g (Σmβ²)/(Σmβ² + p_g ε) → λ v_g` で保存される。

### 罰則スケール `m_j`

**推定器が実際に見る列**に合わせる。`_assemble_problem`（`fit.jl:9-45`）が組む設計は
`w = 0` で中心化エネルギー列、`w > 0` で中心化エネルギーと**白色化トルク**の積み上げ:

```
エネルギー: m_j(w) = (1−w)·Var_ref[Φ_j]
                   + w·(1/(3·n_atoms))·E_ref[ Σ_a ‖(∂Φ_j/∂e_a) × e_a‖² ]
モーメント: m_j    = E_ref[Φ_j²]                    (中心化しない)
```

トルク項の `1/(3·n_atoms)` は必須。`_assemble_problem`（`fit.jl:16-20`）は energy 行を
`√((1−w)/n_E)`、torque 行を `√(w/n_T)` で白色化し、`n_T = n_E·3·n_atoms`
（`design.jl:29-31` の `block = 3*nat`）なので、推定器が見る列ノルム二乗の期待値は上式に
なる。**原子平均と 3 成分平均を落とすと 3×3×3 bcc Fe で 162 倍ずれる**。計量の一様な
定数倍は無害ではない（`w_j = 1/(m_jβ_j²+ε)` の `ε` 校正と実効 λ が動く）。
`√((1−w)/n_E)` / `√(w/n_T)` の `1/n` は**行平均を取ることで既に織り込まれている —
二重に掛けない**。`Φ_j` は `(4π)^(n_spin/2)` 込みの**設計列関数**を指す。

`Var` と `E[·²]` の使い分けは中心化の有無に対応する。**エネルギー側は常に分散**で定義し、
二次モーメントで代用しない — 一致するのは「メンバの原子がすべて相異なる」特殊ケース
（全サイト `l ≥ 1`（`salcbasis.jl:628`）から `E_ref[Φ_j] = 0`）にすぎない。
`AllImages` の自己像メンバは両因子が同じ球に乗るので `E_ref[Φ_j] ≠ 0` になりうる
（ただしその基底は `_refuse_self_image_basis` が `SCEDataset` の門で止めるので
罰則付きフィットには届かない — これは定義の**反例**であって根拠ではない）。

`m_j = 0` は「罰則を受けない」の意味を持たせる。`penalty_metric(mb)` は
μ₀ 列（切片）と `ds.vanishing` に載る恒等ゼロ列に対して厳密に `0` を返す。

**なぜ訓練データの列ノルムではないのか。**

1. **λ をセル間で比較したい。** 対の spec の 4 体実験は「3×3×3 の 108 列 vs +4 体」と
   「3×3×3 vs 4×4×4」の比較で成り立つ。`Var_ref[Φ_j] ∝ n_cells` かつ中心化後の `y` の
   分散も `∝ n_cells` なので RSS と罰則が同じ次数で伸び、比較が意味を持つ。
2. **CV が素直に書ける。** データ由来の標準化は fold ごとに再計算しないとホールドアウトが
   漏れる。基底固有の計量は **fold 不変**なので、その配慮が要らない。
3. **事前分布は関数に置くべきで、データの励起され方に置くべきではない。** 転移可能な
   物理モデルが目的なので、関数空間の量で罰する方が筋が通る。代償は Q5。

**計算**: 固定 seed の一様乱数配置 `K` 本から**列ごとの二乗和をチャンク累積**する。
`X_T` を実体化しない — `K·3·n_atoms × p` は 3×3×3 bcc Fe（54 原子・108 列）で
約 260 MB、4×4×4 で 4 倍。しかもこれを走らせるのは**構築子**（`Ridge(basis; …)`）。
`K → ∞` で解析値に収束し、その解析値がテストのオラクルになる。

**RNG はバージョン安定でなければならない。** 計量は**すべての罰則付き係数の入力**なので、
ストリームが動けば記録済み fit も回帰ピンも静かに動く。この package は乱数依存を明確に
避ける設計（`_grouped_folds`（`selection.jl:360-374`）は「no RNG dependency」を売りにし、
`fingerprint` は load 時に再計算する）なので、ここだけ例外にはできない。⇒ サンプラを
**明示**し（`MersenneTwister(seed)` + `randn` 正規化、または決定論的な golden-spiral 構成）、
「Julia のバージョン跨ぎで同一ストリームであることは保証されない」旨を docstring に書き、
**代表列 3 本の計量値を §14 のピンに含める**。

### 偶然の事前分布と意図的な事前分布を分ける

エネルギー側には既に**意図して設計された**群重み `v_g = √p_g·(c_g/c̄)^θ`
（`cost_weights`、式は `selection.jl:146`）があり、`select_fit` はこの θ を掃いて
(cost, error) の Pareto 前線を描く。そこに列ノルムの偶然の効果が重畳しているので、
**掃いている θ が何を意味するのかが現状では言えない**。

しかも `c_g`（`group_costs`、`selection.jl:96-113` — 群の distinct entry の union に
わたる member-site 数の和）と偶然の列ノルムは **軌道サイズ因子を共有する**ので、
「重畳」どころか**部分的に同一量**である。一方 `√p_g` はチャネル数（対称性が決める
不変量の本数）でメンバ多重度とは別物 — むしろ高対称軌道では反相関する。
⇒ 計量が偶然分（軌道サイズ × 畳み込み多重度）を消し、`v_g` が意図する階層を担う。
M6 の受け入れ確認「計量前後で `select_fit` の群の順序が変わる」は**期待された結果**。

### 代替案と却下理由

- **(a) `penalty_free::Vector{Int}` を推定器に足す。** 計量が既に列ごとの量なので二重。却下。
- **(b) 群重み `v_g = 0` で切片を外す。** `GroupAdaptiveRidge` にしか無く、
  `Ridge` / `AdaptiveRidge` 経路に穴が残る。却下（§Summary）。
- **(c) `fit` が推定器を差し替える。** action at a distance。却下。
- **(d) モーメント側も軌道ごとに中心化して μ₀ を解析的に外す（Frisch–Waugh–Lovell）。**
  **技術的には可能**である。μ₀ 列は行台が互いに素な**軌道インジケータ**なので互いに直交し、
  「掃き出す」は軌道ブロックごとの中心化に等しい。エネルギー側 `j0` の正確な一般化。
  却下の理由は「軌道別だから外せない」ではなく:
  1. 中心化は**すべての推定器が見る設計**を変えるので `OLS` が浮動小数点で動く。
  2. `ds.vanishing` で凍結された軌道では中心化してはならず、凍結集合と連動する分岐が増える。
  3. 掃き出しは「全 μ₀ を外す」しか表現できない。計量 0 は同じ推定量を、設計を触らずに、
     列ごとに得る。
  ⇒ 同じ推定量を得る 2 経路のうち機構の少ない方を採る。この理由を docstring に残す。

## Module layout

| Target | Change |
|---|---|
| `src/fitting/estimators.jl` | `Ridge` / `AdaptiveRidge` / `GroupAdaptiveRidge` に `metric` と `metric_provenance`（下記 API）。**計量が入る独立サイトは 6 つ**（下表）。`Ridge` の解を `Symmetric(X'X + λ·Diagonal(m)) \ X'y` に。**収束判定を計量座標かつ罰則列に限る**（§設計の要 2）。**PD ガードは 3 経路すべて**（`Ridge` / `AdaptiveRidge` / `_solve_gar`）: `free = findall(iszero, D)` が非空なら `XtX[free,free]` の Cholesky を先に取り、失敗は**列名指し**（`_intercept_columns` の軌道名）で拒否。`Symmetric \` は特異でも投げずに数値ゴミを返す（`estimators.jl:421-428` が `lambda == 0` について「‖β‖ ~ 1e16 vs OLS 0.48」と実測記録）。**長さ検査は `solve_coefficients` の入口**で `length(metric) == size(X, 2)`（ctor は設計を知らないので不可能。文言は `column_groups`（`:519-523`）に合わせ「built on a different basis」を含める）。docstring。 |
| `src/fitting/selection.jl` — `select_fit` | 推定器フィールドを展開する箇所は **5 つ**: `:547` パス本体、`:597` GCV 重み、`:629` `:cv` fold 内、**`:645` 最終 cold 再解の推定器再構築**（`GroupAdaptiveRidge(lams[sel], est.column_groups, …)` の**内側 ctor への位置引数展開**で `path.fit` を作る）、`:662` 選択点 GCV。**`:645` を落とすと λ パスは計量つき、返る `path.fit` だけ計量なし**になり、`path.score[selected] == gcv(path.fit)`（`:653-670` が明示的に守る性質）と `n_alive`/`cost` の再導出が全部壊れる。⇒ `metric` は内側 ctor に**既定なしの位置引数**で入れ、漏れが `MethodError` になるようにする（decor エンジンの `isotropy` と同じ「静かに落ちるより落ちろ」方針）。 |
| `src/fitting/selection.jl` — その他 | `penalty_metric(basis::SCEBasis; torque_weight, nconfig, seed)`。`Ridge(basis; …)` / `GroupAdaptiveRidge(basis; …)` が既定で載せる。`cost_weights` docstring に役割分担の節。`_edof` の罰則なし列対応（**`XtX` 経路含む**）。`cross_validate(::MomentDataset, …)` / `effective_dof(::MomentFit)` / `gcv(::MomentFit)` / `_gcv_neff(::MomentFit)` / `MomentCVResult`。 |
| `src/fitting/momentfit.jl` | `penalty_metric(mb::MomentBasis; free_intercepts = true, nconfig, seed)`。`_intercept_columns(mb)::Vector{Int}`（判別式 `key.body == 1 && length(key.decors) == 1 && key.decors[1].spin_l == 0`；`is_marked(d)` を併記して意図を明示）。`GroupAdaptiveRidge(mb; …)` が計量を既定で載せる。**`_reduce_to_active` に `metric` の縮約**を追加。 |
| `src/fitting/momentfit.jl` / `fit.jl` — 列構造推定器の縮約 | `metric` は `Ridge` / `AdaptiveRidge` を**列構造つき推定器**にする。CLAUDE.md が名指しする失敗モード（「a new column-structured estimator needs its own `_reduce_to_active` or it hits `DimensionMismatch`」）に該当するので、`_reduce_to_active(::Ridge, …)` / `(::AdaptiveRidge, …)` を追加。`refit`（`fit.jl:180` の `solve_coefficients(estimator, view(X, :, support), y)`）と `select_support` も support への縮約が要る（または名指しで拒否）。 |
| `docs/design-notes.md` §13（L406-…） | 重み写像に計量が入るので導出を更新（不動点 `→ λ v_g` の再導出、`ε` の意味、GCV の df）。**`CLAUDE.md:497` の「design-notes §13」参照は正しい**（`docs/design-notes.md:406` を指す）ので**触らない**。 |
| `src/SCEFitting.jl` | `export penalty_metric, MomentCVResult`。 |
| `test/unit/test_{selection,fit,momentfit}.jl` | §Test strategy。 |
| `docs/src/guide/{fitting,moment}.md` | 罰則計量の節。**実行される `@example` の λ 値（`fitting.md:110, 149, 189, 228, 246`）は意味が変わるので再校正**（strict ビルドは実行するので数値を印字している箇所は落ちうる）。 |
| `docs/src/api.md` / `SPEC.md` / `CHANGELOG.md`（**BREAKING**） / `CLAUDE.md`（coupled-site L490-512） / `docs/specs/README.md` / `bench/` | 更新。 |
| upstream `SLCE.jl` | 純スピン両チャネル（M5）。joint は Q4。 |

## API

```julia
Ridge(; lambda, metric = nothing, metric_provenance = nothing)
AdaptiveRidge(; lambda, epsilon = 1e-8, max_iter = 50, tol = 1e-6, metric = nothing, …)
GroupAdaptiveRidge(column_groups, group_weights; lambda, metric = nothing, …)
#   metric: 有限かつ ≥ 0 の列ごとの罰則スケール。0 = その列は罰則を受けない。
#           全列 0 は拒否。長さ検査は solve_coefficients の入口。
#   metric_provenance::Union{Nothing,@NamedTuple{torque_weight::Float64, nconfig::Int,
#                                                seed::Int, fingerprint::UInt64}}
#     `show` に出し、`fit` / `select_fit` / `cross_validate` が torque_weight と基底
#     fingerprint の不一致を**拒否**する（`SCEPredictor` の fingerprint と同じ
#     「信じずに再検査」）。m_j(w) が w 依存になった以上、これが無いと
#     「エネルギー専用計量 × トルク優勢フィット」が黙って通り、どのゲートにも掛からない
#     （スケール不変性は m ∝ c² なら何でも成り立つので検出できない）。

penalty_metric(basis::SCEBasis; torque_weight::Real = 0.0,
               nconfig::Integer = 2000, seed::Integer = 1)::Vector{Float64}
penalty_metric(mb::MomentBasis; free_intercepts::Bool = true,
               nconfig::Integer = 2000, seed::Integer = 1)::Vector{Float64}
#   参照アンサンブルは **mode-4 恒等代入（axes = e）** で組む。計量は基底固有量なので
#   訓練データの mode に依存させない（CLAUDE.md「評価軸は constraint_mode で決まる」）。
#   `ds.vanishing` 相当の恒等ゼロ列は `moment_resolvability(mb)` を呼んで**構造的に**
#   0 を返す。MC 推定が偶然 0 に落ちた列は 0 と扱わず、床を入れず、loud に拒否する
#   （「厳密ゼロは構造的にのみ許される」）。

Ridge(basis::SCEBasis; lambda, torque_weight = 0.0, metric_nconfig, metric_seed)
GroupAdaptiveRidge(basis::SCEBasis; lambda, theta = 1.0, torque_weight = 0.0, …)
GroupAdaptiveRidge(mb::MomentBasis; lambda, …)     # 計量が μ₀ を外す

cross_validate(ds::MomentDataset, estimator; nfolds = 5, seed = 1)::MomentCVResult
effective_dof(f::MomentFit)::Float64
gcv(f::MomentFit)::Float64
```

`metric` の**列順契約**: `SALCKey` ソート順 = 設計列順（CLAUDE.md「Design-matrix columns
are identified by `SALCKey`」）。長さ違いは `column_groups` と同文言で拒否する。

`free_intercepts = true` が既定。μ₀ を縮めることに物理的意味は無い。`false` は
意図的な比較のためだけに残す。

### 診断が再構成する設計（**必須**）

`fit(MomentFit, …)` は `ds.X[ds.keep, active]` を**縮約した推定器**（`momentfit.jl:413`）で
解き（`:415`）、格納するのは**未縮約**の推定器と凍結ゼロ込みの全長 `f.coeffs`
（`:421-424`）。診断が素直に `_penalty_diagonal(f.estimator, f.coeffs)` と
`_edof(ds.X, …)` を呼ぶと長さ検査は通ったうえで (i) 解いていない凍結列に自由度を課金し、
(ii) 凍結列込みの群サイズ `p_g` で**ソルバとは違う罰則対角**を使う。これは
`_refuse_refit_diagnostic` が「上流で df 40.0 対 正直な ≈ 5」と実測した失敗そのもの
（`selection.jl:245-249`）。

```julia
active = trues(size(ds.X, 2)); active[ds.vanishing] .= false
X = ds.X[ds.keep, active];  y = ds.y[ds.keep]
beta = f.coeffs[active];    est = _reduce_to_active(f.estimator, active)
```

`_gcv_neff(::MomentFit) = count(ds.keep)`。大域切片は無く、μ₀ は設計内の列として
`rank(X_F)` 項が自由度を課金する（エネルギー側の `+1` と同じ意味の項）。

### `MomentCVResult`

`CVResult`（`selection.jl:882-893`）は `torque_weight` / `rmse_energy` / `rmse_torque` を
持ち `show` もそれを印字する。モーメント CV には対応物が無く、fold ごとの自然な量は
`y = ê·M` の RMSE [μ_B] 一本。**流用せず別型**: `nfolds`, `seed`, `n_holdout`,
`score`（fold ごとのベクトル）, `pooled_score`（**λ 選択に使うスカラ** — `CVResult` の
`pooled_*` と同じく「whole-dataset numbers, not means of the per-fold columns」）,
`score_defined`（disclosure, Q1）, `rmse_moment`。

### 計量が入る 6 サイト（`_gar_weights!` だけではない）

requirements の初版は「`_gar_weights!` が罰則対角の唯一の定義」と書いていたが、
**これは事実に反する**。CLAUDE.md L494-497 自身が逆を述べている:

> `_gar_weights!` (the single definition of `wⱼ = v_g/(‖β_g‖² + p_g·ε)`;
> `_penalty_diagonal` **has one method per linear estimator**, and `AdaptiveRidge`'s
> `1/(β² + ε)` **must stay in sync with its solve loop**)

| # | サイト | 変更 |
|---|---|---|
| 1 | `solve_coefficients(::Ridge)` `estimators.jl:432` | `Diagonal(m)` |
| 2 | `solve_coefficients(::AdaptiveRidge)` iteration 0 `:446` | `λ·Diagonal(m)` |
| 3 | 同 重み更新 `:449` `@. w = 1/(beta^2 + ε)` | `1/(m·beta^2 + ε)`、返り値は `D = m·w` |
| 4 | `_solve_gar` cold start `:504-508` | `D_j = m_j·v_g` |
| 5 | `_gar_weights!` `:467-480` | 分母に `m`、返り値 `D = m·w` |
| 6 | `_penalty_diagonal` の 3 メソッド `selection.jl:174-193` | `(::Ridge) → (λ, m)`、`(::AdaptiveRidge) → (λ, m_j/(m_jβ_j²+ε))`、`(::GroupAdaptiveRidge)` は `_gar_weights!` 経由 |

正しい不変条件は「**`_gar_weights!` は群形の唯一の定義**。`Ridge` と `AdaptiveRidge` は
ソルバ側と `_penalty_diagonal` 側の**対になった 2 サイト**を持ち、計量は両方に入る」。
6 のいずれかを落とすと `gcv` / `effective_dof` / `select_fit` の GCV が**ソルバと違う
罰則対角**を使う — `_refuse_refit_diagnostic`（`selection.jl:245-249`）が名指しで潰したのと
同じ失敗クラス。

## Types and conventions

- **物理規約の変更なし。** 変わるのは推定量の定義（罰則対角）だけ。
- **新しい不変条件**: 罰則対角 `D_j = m_j w_j` は `_gar_weights!` が唯一の定義で、
  ソルバの反復・`_penalty_diagonal`・`select_fit` の 4 呼び出し点がすべてそこを通る。
  罰則なし列（`D_j = 0`）がある場合の有効自由度は `df = rank(X_F) + Σ_i s_i/(s_i + λ)`。
- **上流差分**: SLCE.jl は両チャネルに同じ欠陥。M7 / Q4 の段取り。

### `_edof` の拡張

`X = [X_F X_P]`、`D = diag(0_F, W)`、`A = X'X + λD` と置くと

```
tr(H) = tr(A⁻¹X'X) = p − λ·tr((A⁻¹)_PP W)
(A⁻¹)_PP = (X_P' M X_P + λW)⁻¹,   M = I − X_F(X_F'X_F)⁻¹X_F'
G = W^{-1/2}(X_P' M X_P)W^{-1/2}
⇒ tr(H) = p_F + Σ_i s_i/(s_i + λ)        (s_i = eig(G))
```

零固有値は寄与 0 なので非零 `s_i` に限ってよく、`p_F = rank(X_F)`。λ → ∞ で
`df → rank(X_F)`。

**前提**: (1) `A ≻ 0` ⟺ **`X_F` が列フルランク**（`v'Av = ‖Xv‖² + λΣ_P w_j v_j²` が 0 に
なるのは `v_P = 0` かつ `X_F v_F = 0` のときだけ。罰則ブロックは `λW ≻ 0` なので常に安全）
⇒ ガードは `A` の分解成否ではなく **`X_F` のランク検査**。(2) `λ > 0`。
(3) 密ハット参照が定義されるのは `p_F ≤ n` のとき。

**実装方針**: `free = findall(iszero, D)` の分岐を **`XtX` キーワードの有無より先**に
置く。`_edof` は常に `X` を受け取る（`XtX` はキーワード）ので双対形も使える。
`Q` = `X_F` の thin QR、`B = X̃_P − Q(Q'X̃_P)` と置くと `G` の非零固有値 =
`B'B`（`p_P×p_P`）= `BB'`（`n×n`）の非零固有値。4 分岐:

| `free` | 形状 | 使う形 |
|---|---|---|
| 空 | `p ≤ n` | 現行の `Symmetric(X̃'X̃)`（`XtX` があればそれを使う） |
| 空 | `n < p` | 現行の双対 `Symmetric(X̃X̃')` |
| 非空 | `p_P ≤ n` | `B'B`（Gram 形。`XtX` から `G = W^{-1/2}(X_P'X_P − X_P'X_F(X_F'X_F)⁻¹X_F'X_P)W^{-1/2}` としても同値） |
| 非空 | `n < p_P` | `BB'`（双対）。`p_F ≤ n` は前提 3 |

`XtX` キーワードが節約するのは **`p ≤ n` の Gram 形成だけ**。`M` を `n×n` で作らない。
`free` が空なら現行の式をそのまま通す（エネルギー側のビット等号の最短保証）。
なお `penalty_metric(::SCEBasis)` は 0 を返さない設計なので、エネルギー側で `free` が
非空になるのは**ユーザ供給の計量**を `select_fit` に渡した場合に限る。

**代替（却下）**: `D = 0` を `1e-12` に置換。条件数が `1e12` 悪化。`Ridge` の
`lambda == 0` 分岐（`estimators.jl:422-428`）と同じ方針で拒否。

### `_solve_gar` の収束判定

現行の停止則 `rel = max|Δβ| / max(maximum(abs, beta_new), eps)` は分母が**大域**。
μ₀ ≈ 2.2 μ_B が罰則なしで残ると `tol = 1e-6` が実質 2e-6 の**絶対**閾値になり、
~1e-3 の罰則係数が未収束のまま止まる。`_penalty_diagonal` は返された `beta` から `D` を
再計算するので、早期停止は「解かれていない平滑化行列」を診断に渡す。⇒ **停止則を
罰則列（`D_j > 0`）に限る**。

## Impact on coupled sites

- [x] **GCV ↔ `_assemble_problem` ↔ `islinear` ↔ GAR weight map**（CLAUDE.md L490-512）:
      本 spec の中心。`_gar_weights!` が唯一の定義（計量は分母と返り値の両方）、
      `_penalty_diagonal` はそこを通す、`_edof` が拡張、`_solve_gar` の停止則が変わる、
      **`select_fit` の 4 呼び出し点**が計量を運ぶ。密ハット行列テストが両チャネルで
      一致することがゲート。
- [x] **`docs/design-notes.md` §13**: 重み写像・不動点・`ε` の校正の導出が動く。
- [x] **`fit` ↔ `refit` は `_assemble_problem` を共有**: 設計を触らないのでバイト同一性は
      不変。ただし `refit` / `select_support` は列構造推定器の support 縮約が要る。
- [x] Decor エンジン ↔ pointed 基底 ↔ `momentfit.jl` の門: `_intercept_columns` の判別式は
      `_mark_decor(0) = SiteDecor(disp = (1,0))`（`momentbasis.jl:170-171`）と
      `SiteDecor` 内側構築子の `spin_l == 0 ⟹ has_disp`（`decor.jl:97-99`）から誤検出不能。
      `penalty_metric(mb)` は `_design_moment` を通すので評価規約に自動追従。
- [x] **エネルギーカーネル ↔ トルクカーネル**: 共フィットの計量が両ブロックに跨がる。
      スケール不変性自体は保たれる（`Φ_j → c_jΦ_j` なら勾配列も同じ `c_j` 倍）が、
      **尺度は `w` で変わる**ので `m_j(w)` を `w` の関数として定義する（§罰則スケール）。
- [ ] `SALCKey` 順序 / 設計列 / TOML 永続化: 影響なし（係数の値は変わるが列の同一性・
      順序・スキーマは不変）。`metric` の列順契約は API 節。
- [x] 上流差分 ledger: M5 / Q4。
- [ ] `.claude/agents/`: 変更なし。
- [x] `SPEC.md` / `docs/src/api.md`: `SPEC.md` の fitting 節、`api.md` に 2 行。

## Test strategy

1. **罰則なし列の厳密性（解析オラクル）** — 直交設計で `m_j = 0` の列の係数は**任意の λ・
   任意の重み・任意の計量**で `β_free = (x_F'y)/‖x_F‖²`（OLS 値）に厳密一致
   （`rtol = 1e-14`）。`A` がブロック対角になることの帰結。罰則列は適応系なので閉形式に
   ならない（不動点が三次式）ため、収束解が定義方程式を満たすこと（罰則列ごとの相対残差
   < tol）と λ 単調で `|β_pen|` が減ることを見る。
2. **スケール不変（解析オラクル、両チャネル、3 推定器すべて）** — 列を乱数正対角 `C` で
   再スケールし計量を `C²` 倍すると罰則付き**予測** `X·β̂` が一致（`rtol = 1e-10`）。
   群内**不均一**な `C` を必ず含める（B1 の破れ方はそこに出る）。
   **現行実装と「計量を分母の外に置いた実装」の両方で落ちることを先に確認**。
3. **group-L0 不動点の保存（解析オラクル）** — 計量ありで λ→大の極限で群あたり罰則寄与が
   `λ·v_g` に収束（`docs/design-notes.md` §13 の性質）。分母の外に置くと `λ·v_g·⟨m⟩_g` に
   なることを対照として示す。
4. **参照計量の解析検算（解析オラクル）** — 2 体 `l=1` の SALC について
   `E[Z_{1m}(e_i)Z_{1m'}(e_j)] = δ_{ij}δ_{mm'}/(4π)` から `E[Φ²]` を手で閉じた形で出し、
   `penalty_metric` が `|m̂/m_exact − 1| < 5/√K`（実測 σ を併記）で一致。
   併せて `MinimumImage` 基底で `E_ref[Φ_j] ≈ 0`、**`AllImages` の自己像列のうち
   `L_S = L_f = 0` の等方チャネル**で `E_ref[Φ_j] ≠ 0`（`L_f ≠ 0` ではテンソルのトレースが
   消えて空振りする）を確認 — エネルギー側の計量を分散で定義する理由の直接ゲート。
5. **`_edof` の密ハット参照（独立実装オラクル）** — `tr(X(X'X+λD)⁻¹X')` と一致。
   `D = 0` の列 0 / 1 / 複数本、`p ≤ n` と `n < p`（自由列は `p_F ≤ n`）、
   **`XtX` キャッシュ経路と非キャッシュ経路の両方**。λ → ∞ で `df → rank(X_F)`。
6. **`select_fit` の返り値が計量つき（独立実装オラクル）** — パス上の点は warm start
   （`selection.jl:546-551` の `beta0 = prev`）で `fit` は cold、`:648-650` が「warm and
   cold agree only within the IRLS tol」と明記しているので**ビット一致は主張しない**。
   ゲートは: `path.fit.jphi` が、**推定器を独立に組み直した**
   `fit(SCEFit, ds, GroupAdaptiveRidge(basis; lambda = path.lambdas[path.selected],
   theta, metric))` とビット一致。これが `:645` の計量落ち（B2）の直接検出器になる。
7. **μ₀ の λ 依存（解析オラクル）** — 実基底では μ₀ は λ に完全不感ではない
   （`β_F = (X_F'X_F)⁻¹X_F'(y − X_P β_P(λ))`）。不感なのは `X_F ⊥ X_P` の直交設計
   （ゲート 1）に限る。ゲートは極限: λ → ∞ で `Ridge(λ)`（計量なし）の μ₀ → 0、
   計量つきの μ₀ → `X_F` だけへの OLS 解に収束。有限 λ では `|μ₀(λ) − μ₀(∞)|` が単調減少。
   **3 推定器すべてで測る**（B の却下理由の直接ゲート）。
8. **CV の分割性質（実装非依存オラクル）** — (a) 全行を分割、(b) 各 config の行が単一 fold、
   (c) **各 fold の学習側で `ds.orbit_rep[train & ds.keep]` が全マーク軌道を覆う**。
   fold 数の扱いはエネルギー側（`selection.jl:960-965`）を踏襲。
9. **CV のリークガード** — `PrecomputedPilot` / それを積んだ `AdaptiveLasso` を拒否。
10. **診断の凍結列整合（独立実装オラクル）** — `ds.vanishing` が非空な基底で
    `effective_dof(f)` が縮約設計の密ハットと一致し、**全設計で組んだ値とは一致しない**。
11. **縮約の網羅** — `metric` つき `Ridge` / `AdaptiveRidge` / `GroupAdaptiveRidge` を
    `fit(MomentFit, …)`（凍結列あり）と `refit` / `select_support`（support 縮約）に
    通して `DimensionMismatch` が出ないこと、または名指しで拒否されること。
12. **拒否** — 全列 0、`islinear` false での `gcv`、非 PD（`X_F` がランク落ち）、
    `metric` 長さ違い。
13. **`OLS` の不変（実装非依存オラクル）** — 捕獲比較ではなく、正規方程式
    `X'(y − Xβ̂) ≈ 0`（`rtol` 明示）＋「列を `C` 倍しても `Xβ̂` が不変」で押さえる。
14. **エネルギー側診断のビット等号（回帰ピン、変更検出器と明示）** — `metric = nothing`
    かつ正の群重みでの `gcv(::SCEFit)` / `effective_dof(::SCEFit)` /
    `GroupAdaptiveRidge` の解が `==`。捕獲コミットと日付、再取得条件を併記。

受け入れ確認（テスト外）: tasklist M6。

## Risks and open items

- **Q1 CV スコアをゲート行だけで測るか。** 推奨 **`ds.keep` で学習・評価**（ゲートは
  物理的な定義域の問題で汎化の問題ではない）。`MomentCVResult.score_defined` を併記。
- **Q2 `refit` / `select_fit` の支持規則は据え置き。** `|jϕ_j|·‖X[:,j]‖ > threshold`
  （`fit.jl:169`, `selection.jl:566-572`）は列再スケールに対して**すでに厳密に不変**
  （`β_j → β_j/c_j`, `‖x_j‖ → c_j‖x_j‖`）なので、そもそも計量の対象外。
  「別の問い」ではなく「もともとスケール不変」が据え置きの根拠。docs にそう書く。
- **Q3 参照アンサンブルの既定 `K`。** `nconfig = 2000` は当て推量。ゲート 4 の収束測定で
  「相対誤差 < 1e-3 に必要な `K`」を実測して決める。`K` / `seed` / `torque_weight` は
  推定器の provenance として `show` に出す。
- **Q4 SLCE joint 経路（先送りの具体的理由）。** DISP 因子は `|u|^{2k} R_{lm}(u)`
  （SLCE `decor.jl:28`）で、`R_{lm}` は単位球正規化を持たない**立体調和**、`u` について
  次数 `2k+l` の斉次式。したがって
  (i) 参照分布に**長さスケール σ** が要る、(ii) **σ の選択自体が非調和次数にわたる階層
  事前分布**になる（本 spec が偶然分から切り離そうとしている当のもの）、
  (iii) `m_j` が `Å^{2(2k+l)}` の次元を持ちスカラー λ が単一の意味を失う。
  ⇒ SCEFitting（純スピン両チャネル）を先に完結させ、joint は別 spec。
- **Q5 参照分布と訓練分布の乖離。** 参照は一様乱数スピン、実データ（bcc Fe / FeRh / klm）は
  ほぼ共線的な低温配置。参照で弱く励起されるが訓練で強く効く列は `λ·m_j` が小さく
  **実質無罰則**になる。「事前分布は関数に置く」という立場の直接の代償。
  ⇒ 診断を出す: `m_j` と訓練列分散 `Var_train[Φ_j]` の比の max/min を
  `show` / `coeftable` の隣に出し、乖離が見えるようにする。
- **Q6a `p_g` に計量ゼロ列を数えるか。** `w_j = v_g/(Σ_{k∈g} m_kβ_k² + p_g·ε)` の分子側の
  和には `m_k = 0` の列が寄与しないのに `p_g` には残るので、§13 の「群サイズ非依存の
  per-coefficient floor」の議論が微妙にずれる。`p_g` を**罰則列数**
  `count(m_k > 0, k ∈ g)` にする案がある。不動点 `→ λv_g` はどちらでも保存されるので
  純粋に `ε` 校正の問題。
- **Q6 `ε` の校正。** 分母が `Σ m_k β_k² + p_g ε` になると `ε` は「係数の大きさの床」から
  **「計量重み付きの大きさの床」**に変わる。既定 `1e-8` の意味、`estimators.jl:257-260` の
  docstring、§13 の「群サイズ非依存の per-coefficient floor」議論、`AdaptiveRidge` との
  退化一致（テストでピン留め）が全部動く。据え置くか再校正するかを決める。
- **Q7 `m_j` の量的な大きさ。** 「メンバ多重度で `N!` 倍」という単純な見積もりは**正しくない**:
  `_canonical_basis`（`salcbasis.jl:246-268`）はアンカー空間で単位ノルムの係数ベクトルを
  作り、`_canonicalize_members` で `N!` の ordering が畳まれるとき同じ単項式キーに
  コヒーレントに加算されるのは `ls` 多重集合の**安定化群 S** の分だけ（異なる assignment は
  別キーに散る）。⇒ `E[Φ²] ~ 軌道サイズ × S`、`S ∈ [1, N!]` で、汚染は**体数だけでなく
  同一体数内のラベル間でも変わる**（3 体で `ls=(1,1,1)` は `S=6`、`(1,1,2)` は `S=2`、
  `(1,2,3)` は `S=1`）。動機の主張は強まるが具体数値は**ゲート 4 の実測で埋める**。
- **記録済み fit の再取得。** bcc Fe l02…l044 / `m_*`、FeRh、klm の罰則付き結果は λ の
  意味が変わる。`AdaptiveLasso(pilot = Ridge(...))`（`fitting.md:149`）も pilot 経由で
  変わるので CHANGELOG の BREAKING で名指しする。
- **数値結果への影響**: `OLS` は両チャネルでゼロ。罰則付き推定は**すべて変わる**。
  `_assemble_problem` の `w = 0` 枝を直したとき（`fit.jl:30-41`、SLCE `572cbe0`）と
  同格の BREAKING として扱う。
