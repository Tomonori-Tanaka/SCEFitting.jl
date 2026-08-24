# Requirements: 罰則計量の正当化と λ 選択（エネルギー / モーメント両チャネル）

Status: draft (2026-08-24) — spec-reviewer 第 3 回反映済み。後方互換は制約にしない、
エネルギー側にも波及させる方針（ユーザー指示 2026-08-24）

## Goal

罰則付き推定を**基底の正規化規約に依存しない形**に直し、λ を選ぶ層を揃える。

1. **罰則計量**（両チャネル）— λ‖β‖² は列の再スケールで不変でないため、現状は
   「どの規約で基底を書いたか」で罰則の強さが変わる。基底固有の計量で測る形に直す。
2. **μ₀ 切片を罰則から外す**（モーメント側のみ）— 切片を 0 に縮めることに物理的意味が無い。
   エネルギー側は `j0` が設計の外なので既に正しい。
3. **λ 選択層**（モーメント側のみ）— `cross_validate` / `gcv` / `effective_dof` の
   `MomentFit` 版。エネルギー側は既にある。

既存の呼び出しの数値が変わること、記録済みの fit を取り直す必要が生じることは許容する
（ユーザー指示）。変えてよい代わりに **CHANGELOG に BREAKING として立てる**。

## Background

### ① μ₀ 切片が罰則を受ける（モーメント側）

エネルギー側は設計を列中心化し `j0 = ybar − dot(xbar, jphi)` と解析的に外すので、罰則付き
設計に切片列が存在しない（`src/fitting/fit.jl:92`）。モーメント側は `l = 0` の 1 体
`[MARK]` 列がそのまま μ₀ で設計行列の中にいる。docstring も「the μ₀ columns are penalized
like every other column — v1 is OLS-first」と明記（`momentfit.jl:391-394`）。μ₀ 列は
**軌道インジケータ列**（該当軌道の行では同一の定数、他のマーク原子の行では項リストが空で
厳密に `0.0`）で、係数は裸のモーメント長（bcc Fe で ≈ 2.2 μ_B）。`Ridge` はこれを 0 に
向けて縮めるので **λ を上げるほど予測モーメントが系統的に短くなる**。μ_B = 0 は特別な点では
ないので、この縮小に意味は無い。

### ② 罰則が列スケールに依存する（**両チャネル**）

`Ridge` は `Symmetric(X'X + λI) \ X'y`（`src/fitting/estimators.jl:432`）で標準化しない。

**機構はメンバ多重度**であって調和関数のスケールではない。後者はむしろこの問題を消すための
因子で、`Zₗₘ` の per-site `(4π)^(−1/2)`（`Harmonics.jl:7`）が一様乱数スピンで
`E[Z²] = 1/(4π)` を与え、`n_spin` 個のスピン因子分の `(1/4π)^{n_spin}` を
`(4π)^{n_spin/2}` がちょうど打ち消して**単項式 1 個の RMS を 1** にする
（`n_spin` = スピン因子を持つ decor 数。エネルギー側は `= N`、pointed で rank-0 マークなら
`= N−1`；`salc.jl:279, 470-471`）。結合テンソルも射影器の直交固有ベクトルで正規化済み。

残るのは SALC が軌道のメンバ和であることで、`‖x_j‖² ~ 軌道サイズ × 畳み込み多重度`。
**「`N!` 倍」ではない**: `_canonicalize_members` で `N!` の ordering が畳まれるとき同じ
単項式キーにコヒーレントに加算されるのは `ls` 多重集合の**安定化群サイズ `S`** の分だけで、
異なる assignment は別キーに散る（`S ∈ [1, N!]`；3 体なら `ls=(1,1,1)` で `S=6`、
`(1,1,2)` で `S=2`、`(1,2,3)` で `S=1`）。⇒ 汚染は**体数だけでなく同一体数内のラベル間でも
変わる**。主張はむしろ強まるが、具体的な倍率は断定せず design §Test strategy 4 の実測で
埋める。純スピンのエネルギー基底も同じ全順序規約なので、**この汚染はエネルギー側にも
同じだけある**。

**向きが問題**: 直交近似で ridge の相対縮小率は `‖x_j‖²/(‖x_j‖² + λ)` なので、
**ノルムが大きい列ほど縮まない**。現状は「体数が上がるほど罰則が軽い」という**暗黙の
事前分布**が入っており、物理的に望ましい向き（高次の補正ほど強く抑える）の逆。

**エネルギー側では、この偶然の事前分布が意図的な事前分布を汚染している。**
`cost_weights(basis; theta)`（`selection.jl:132`）は
`v_g = √p_g · (c_g/c̄)^theta` という**意図して設計された**群重み（Yuan–Lin の群サイズ因子 ×
MC 収縮コスト）を与え、`select_fit` はこれを掃いて (cost, error) の Pareto 前線を描く。
列ノルムの偶然の効果がそこに重畳しているので、**掃いている θ が何を意味するのかが
現状では言えない**。計量で偶然分を消して初めて `cost_weights` が設計どおりに効く。

さらにモーメント側では、この汚染が対の spec
[`260824-pointed-nbody-general`](../260824-pointed-nbody-general/) の 4 体実験に直撃する。
標準化せずに 4 体列を足すと 4 体列は実効 λ が数倍小さい状態で入るので、**「4 体が効いた」
という偽陽性に系統的に傾く**。CV が 0.0217 μB を割っても物理か有利な事前分布か区別できない。

⇒ 二段に分ける: (i) 計量で**偶然の**事前分布を消し、(ii) 意図的な階層は**群重みで明示する**。
今は (i) と (ii) が混ざっている。

### ③ λ を選ぶ層がモーメント側に無い

`gcv` / `effective_dof` / `select_fit` / `cross_validate` はすべて `::SCEFit` /
`::SCEDataset` 型付け（`selection.jl` — `gcv` 347 / `effective_dof` 287 / `select_fit` 502 /
`cross_validate` 944）。`MomentFit` 版が無く、λ は手で決めるしかない。

### 前例

λ の意味を変える BREAKING は**この package で既に一度やっている**。`_assemble_problem` の
`w = 0` 枝を `√(1/n_E)` スケールに直した変更（`fit.jl:30-41` のコメント）は
「BREAKING: λ for an energy-only penalized fit is `n_E` times smaller than before」と
明記され、SLCE.jl `572cbe0` から逆移植された。本 spec は同種・同格の変更で、扱いも同じにする。

### エネルギー側の計量は分散で定義する

計量は**推定器が実際に見る列**に合わせる。エネルギー側は `_assemble_problem` が列中心化
するので `m_j = Var_ref[Φ_j]`、モーメント側は中心化しないので `m_j = E_ref[Φ_j²]`。

理由は素直に「`_assemble_problem` が列中心化するから」。二次モーメントと一致するのは
特殊ケースにすぎない: 純スピンのエネルギー SALC は全サイトが `l ≥ 1`（`salcbasis.jl:628`）
なので、**メンバの原子がすべて相異なる限り** `E_ref[Φ_j] = 0`。`MinimumImage` はこの条件を
満たす。反例として `AllImages` の自己像メンバ `(a,0)-(a,R)` は両因子が同じ球に乗り
`E[Z_{lm}(e_a)Z_{l'm'}(e_a)] = δδ/(4π) ≠ 0`（ただしその基底は `_refuse_self_image_basis` が
`SCEDataset` の門で止めるので罰則付きフィットには届かない — 一致が定義ではないことの
**反例**であって根拠ではない）。

## Scope

Includes:

- **罰則計量**（両チャネル）: `Ridge` / `AdaptiveRidge` / `GroupAdaptiveRidge` に列計量
  `metric` を持たせ、罰則対角に乗せる（design §罰則計量）。設計行列は**書き換えない**。
- 基底由来の計量の構築: `penalty_metric(basis::SCEBasis; …)` /
  `penalty_metric(mb::MomentBasis; …)`、および基底を取る推定器コンストラクタ。
- モーメント側: μ₀ 列を**構造的に**罰則から外す（`penalty_metric(mb)` が `0` を返す）。
  `group_weights` の `≥ 0` 緩和は**不要**（計量 1 つで 3 推定器すべてが片付く）。
- **`select_fit` の 4 呼び出し点**（`selection.jl:547, 597, 629, 662`）が推定器フィールドを
  展開して `_solve_gar` / `_gar_weights!` を呼ぶので、計量を明示的に通す。
- **列構造推定器の縮約**: `_reduce_to_active(::Ridge, …)` / `(::AdaptiveRidge, …)`、
  `refit` / `select_support` の support 縮約（または名指し拒否）。
- `docs/design-notes.md` §13（`:406-`）の導出更新（重み写像に計量、不動点 `→ λ v_g` の
  再導出、`ε` の意味）。
- `_edof` の罰則なし列対応（**`XtX` キャッシュ経路を含む** — `select_fit` が通る）。
- `_solve_gar` の収束判定を罰則列に限る。
- `cross_validate(::MomentDataset, …)` / `effective_dof(::MomentFit)` / `gcv(::MomentFit)` /
  `MomentCVResult`。
- テスト、docs、`SPEC.md`、`CHANGELOG.md`（**BREAKING**）、`CLAUDE.md` coupled-site
  L490-512 の更新。**`CLAUDE.md:497` の「design-notes §13」参照は正しい**
  （`docs/design-notes.md:406` を指す）ので触らない。
- **上流 SLCE.jl への移植**（両チャネル + joint 経路。design Q4 で段取り）。
- 記録済み fit の再取得（package 外、`~/jijs`）: bcc Fe l02…l044 / `m_*`、FeRh、klm。

Excludes:

- `ElasticNet` / `Lasso` / `AdaptiveLasso` の計量と罰則なし列。GLMNet の
  `standardize` / `penalty_factor` に乗る話で機構が別（本 spec は線形推定器 3 種に閉じる）。
- `refit` / `select_fit` / `select_support` の**支持規則**の変更。
  `|jϕ_j|·‖X[:,j]‖ > threshold` は「この fit の予測にどれだけ寄与しているか」を測る規則で、
  「どれだけ強く罰するか」とは別の問い。据え置き、2 つのスケールが共存することを文書化
  （design Q2）。
- `select_fit` / `select_support` のモーメント版（`group_costs` に対応物が無い）。
- 4 体列の実データ評価（対の spec）。

## Invariants

物理・数値の正しさに関わるものだけ。「既存の呼び出しと同じ値」は不変条件にしない。

- **`OLS` の結果は両チャネルで不変**（ビット一致）。最小二乗は列スケールに不変なので、
  計量を入れても変わってはいけない。**後方互換ではなく実装の正しさの検査** —
  変わったら設計行列を触ってしまっている。
- **罰則はスケール不変**: 設計列を任意の正の対角 `C` で再スケールし計量を `C²` 倍すると、
  罰則付き**予測** `X·β̂` が不変。②の直接の定義。
- **罰則なし列は直交設計で λ に厳密不感**（一般設計では不感でない — `β_F` は
  `β_P(λ)` を通じて λ に依存する）。
- **罰則対角 `D_j = m_j w_j` の定義は推定器ごとに対で存在する**。`_gar_weights!` は
  **群形の**唯一の定義だが、`Ridge` と `AdaptiveRidge` はソルバ側と `_penalty_diagonal`
  側の**対になった 2 サイト**を持つ（CLAUDE.md L494-497 が明記）。計量は
  **6 サイトすべて**に入る（design §計量が入る 6 サイト）。`select_fit` の
  **5 呼び出し点**（`selection.jl:547, 597, 629, 645, 662`）も同じ対角を運ぶ。
- **group-L0 の不動点が保存される**: 収束時の群あたり罰則寄与が `→ λ·v_g`
  （`docs/design-notes.md` §13 L437-441）。計量を重みの**分母**に入れることでこれが成立する
  （外に掛けると `λ·v_g·⟨m⟩_g` になり `v_g` が group-L0 重みでなくなる）。
- **`PrecomputedPilot` はモーメント CV でも拒否**（エネルギー側 `selection.jl:953` と同理由）。
- **モーメント CV の fold 単位は配置**（`ds.row_config`）。
- **`_assemble_problem` の中心化・白色化は不変**（`fit` ↔ `refit` のバイト同一性）。
- 物理規約（spin layout、`(4π)^(N/2)`、torque 符号 `τ = −e × ∂E/∂e`、mode rule、
  `ds.keep` の意味、`ds.vanishing` の厳密ゼロ凍結、`SALCKey` の順序と列アドレス、
  TOML 永続化）は一切触らない。

## Completion criteria

- [ ] `make test-all`（4 threads）/ `make test-pin` 緑。
- [ ] **スケール不変ゲート**（解析オラクル、両チャネル）: 列を乱数正対角で再スケールし
      計量を `C²` 倍しても罰則付き予測が一致（`rtol = 1e-10`）。**現行実装で落ちることを
      先に確認**（欠陥の存在証明）。
- [ ] **参照計量の解析検算**（解析オラクル）: 2 体 `l=1` の SALC で一様乱数スピン下の
      `E[Φ²]` を手で閉じた形（`E[Z_{1m}(e_i)Z_{1m'}(e_j)] = δ_{ij}δ_{mm'}/(4π)`）から出し、
      `|m̂/m_exact − 1| < 5/√K`（実測 σ を併記）で一致。
- [ ] **中心化が期待値で恒等**（解析オラクル、エネルギー側）: 参照アンサンブルで
      `E_ref[Φ_j] = 0`（全サイト `l ≥ 1` の帰結）を数値的に確認。
- [ ] **`cost_weights` の分離**: 計量導入前後で、同一 θ・同一 λ の `select_fit` が選ぶ
      群の順序が変わることを記録（偶然分が除かれたことの実測。テストではなく受け入れ確認）。
- [ ] **罰則なし列の厳密性**（解析オラクル）: 直交設計で罰則なし列の係数が任意の λ・
      任意の重み・任意の計量で OLS 値に厳密一致（`rtol = 1e-14`）。
- [ ] **μ₀ の λ 依存**（解析オラクル）: 実基底では μ₀ は λ に完全不感ではない
      （`β_F = (X_F'X_F)⁻¹X_F'(y − X_P β_P(λ))`）。ゲートは極限: λ → ∞ で `Ridge` の
      μ₀ → 0、本 spec の推定器の μ₀ → `X_F` だけへの OLS 解に収束し、有限 λ では
      `|μ₀(λ) − μ₀(∞)|` が単調に減る。
- [ ] **`_edof` の密ハット参照**（独立実装オラクル）: `tr(X(X'X+λD)⁻¹X')` と一致。
      罰則なし列 0 / 1 / 複数本、`p ≤ n` と `n < p`（自由列は `p_F ≤ n`）、
      **`XtX` キャッシュ経路と非キャッシュ経路の両方**。λ → ∞ で `df → rank(X_F)`。
- [ ] **`_solve_gar` の収束**: 収束解が定義方程式を**罰則列ごとの相対残差** `≤ tol` で
      満たす（停止則の分母が大域だと μ₀ に支配されて罰則列が未収束で止まる）。
- [ ] **`OLS` の不変**（実装非依存オラクル）: 捕獲比較ではなく、正規方程式
      `X'(y − Xβ̂) ≈ 0`（`rtol` 明示）＋「列を `C` 倍しても `Xβ̂` が不変」。
- [ ] **`fit` ≡ `select_fit`**（独立実装オラクル）: 計量つき `GroupAdaptiveRidge` で、
      `select_fit` のパス上の λ 点と `fit` の同 λ が係数ビット一致。
- [ ] **group-L0 不動点の保存**（解析オラクル）: λ→大で群あたり罰則寄与 → `λ·v_g`。
      計量を分母の外に置いた実装では `λ·v_g·⟨m⟩_g` になることを対照で示す。
- [ ] **縮約の網羅**: `metric` つき 3 推定器を `fit(MomentFit, …)`（凍結列あり）と
      `refit` / `select_support`（support 縮約）に通して `DimensionMismatch` が出ない、
      または名指しで拒否される。
- [ ] **SLCE 移植**: 純スピン両チャネルへ同型移植し SLCE 側 `make test-all` / docs strict 緑。
      joint 側が未移植である期間を ledger 行に明記。
- [ ] **診断の凍結列整合**（独立実装オラクル）: `ds.vanishing` が非空な基底で
      `effective_dof(f::MomentFit)` が `fit` と同じ三つ組（`active` 列 / `keep` 行 /
      縮約推定器）の密ハットと一致し、**全設計で組んだ値とは一致しない**。
- [ ] **CV の分割性質**（実装非依存オラクル）: fold 集合が全行を分割、各 config の行が
      単一 fold、**各 fold の学習側で全マーク軌道が覆われる**（`coverage_floor` は
      データセット全体の門で fold の門ではない）。fold 数の扱いはエネルギー側
      （`selection.jl:960-965`）を踏襲。
- [ ] **拒否**: `PrecomputedPilot`（CV）、`islinear` false（`gcv`）、全群 0、
      非 PD（ある軌道の行が解に 1 本も残らない）、`_reduce_to_active` が罰則群を残さない場合。
- [ ] `make docs` strict 緑。`SPEC.md` / `docs/src/api.md`（`MomentCVResult` の 1 行）/
      `CHANGELOG.md`（**BREAKING**）/ `CLAUDE.md` / `docs/specs/README.md` 更新。
- [ ] `bench/BENCH_LOG.md`: 計量構築（参照アンサンブルの設計行列 1 回分）と λ 掃引の実測。
- [ ] Tier 2 レビューパネル（4 軸）、`numerical-reviewer` の指摘は全件適用。
- [ ] 受け入れ確認（テスト外）: bcc Fe `m_full_l2` で `cross_validate` の λ 曲線・選択 λ・
      μ₀ を記録（**対の spec の 4 体実験のベースライン**）。エネルギー側は l044 系列で
      λ の意味の変化を記録。

## References

- `src/fitting/fit.jl`（`_assemble_problem` L9-45、解析 `j0` L92、`w=0` 枝の BREAKING
  前例 L30-41）。
- `src/fitting/estimators.jl`（`Ridge` L420-433、`_gar_weights!` L467-480、
  `_solve_gar` L487-518、`GroupAdaptiveRidge` 内側 ctor L273-311、`islinear` L343-347、
  `lambda == 0` の近ゼロ禁止 L422-428）。
- `src/fitting/selection.jl`（`cost_weights` L116-149、`GroupAdaptiveRidge(basis; …)` L158、
  `_penalty_diagonal` L174-196、`_rank_df` L202-208、`_edof` L214-235、`_gcv_neff` L241、
  `_refuse_refit_diagnostic` L250、`_grouped_folds` L366、`select_fit` L502、
  `cross_validate` L944、fold 数規則 L960-965）。
- `src/fitting/momentfit.jl`（`fit(MomentFit)` L368-425、`salc_groups` L614-670、
  `GroupAdaptiveRidge(mb; …)` L671-687、`_reduce_to_active` L688-718）。
- `src/basis/Harmonics.jl:7`（per-site `(4π)^(−1/2)`）、`src/basis/salc.jl:186, 279`
  （`(4π)^(N/2)`）、`src/basis/salcbasis.jl:628`（全サイト `l ≥ 1`）。
- CLAUDE.md coupled-site L490-512。
- 対の spec: [`260824-pointed-nbody-general`](../260824-pointed-nbody-general/)。
- 実測の動機: `~/jijs/scefit/for_package/fe_bcc/3x3x3cubic/scefit/INVESTIGATION-moment-residuals.md` §4d。
