# Tasklist: 罰則計量の正当化と λ 選択（エネルギー / モーメント両チャネル）

Status: landed (2026-08-25) — M1–M6 完了（実装コミット: 下記）。M6 の測定は
`~/jijs/scefit/for_package/fe_bcc/3x3x3cubic/scefit/M6-penalty-metric.md`

This file holds coarse-grained, commit-sized milestones. Day-to-day tracking
goes through `TaskCreate` in-session.

## Milestones

### M1 — `_edof` の罰則なし列対応（`fix(selection): effective dof with unpenalized columns`）

エネルギー側のビット等号を保ったまま `_edof` だけ先に直す。罰則なし列を作る経路が
まだ無いので、ゲートは人工設計で書く。

- [x] `free = findall(iszero, D)` の分岐を **`XtX` の有無より先**に。分割形は Gram だけで
      書く（`select_fit` のキャッシュ経路が `Inf` を作る穴を塞ぐ）。`M` は `n×n` で作らず
      `X_F` の thin QR で射影。前提（`X_F` 列フルランク / `λ > 0`）を導出コメントに明記。
- [x] テスト: design §5（密ハット参照、**4 分岐すべて**、λ→∞ で `df → rank(X_F)`）、
      §14（エネルギー側バイト等号ピン）。
- [x] `make test-all` 緑。

### M2 — 罰則スケール（`feat(estimators): a basis-intrinsic penalty scale`）

**計量は重み写像の分母に入れる**（外に掛けると適応系がスケール不変にならず、
group-L0 の不動点 `→ λ v_g` も壊れる）。

- [x] `estimators.jl`: 3 推定器に `metric`（`≥ 0`、全列 0 拒否、長さ検査）。
      **計量が入る 6 サイト全部**（design の表: `Ridge` の解 / `AdaptiveRidge` の
      iteration 0 と重み更新 / `_solve_gar` cold start / `_gar_weights!` /
      `_penalty_diagonal` の 3 メソッド）。返り値は罰則対角 `D_j = m_j w_j`。
      **PD ガードは 3 経路すべて**（`free` 非空なら `XtX[free,free]` の Cholesky、
      失敗は列名指し）、**収束判定は計量座標 `√m_j β_j` かつ罰則列に限る**、
      長さ検査は `solve_coefficients` の入口。`metric_provenance` と不一致拒否。
- [x] `selection.jl`: `penalty_metric(::SCEBasis; torque_weight, nconfig, seed)`、
      基底を取る構築子、`cost_weights` docstring の役割分担、
      **`select_fit` の 5 呼び出し点**（`:547, 597, 629, **645**, 662`）に計量を通す。
      内側 ctor の `metric` は**既定なしの位置引数**にして漏れを `MethodError` にする。
- [x] `momentfit.jl`: `penalty_metric(::MomentBasis; free_intercepts, …)`、
      `_intercept_columns(mb)`、`_reduce_to_active` の `metric` 縮約。
- [x] **列構造推定器の縮約**: `_reduce_to_active(::Ridge, …)` / `(::AdaptiveRidge, …)`、
      `refit` / `select_support` の support 縮約（または名指し拒否）。
- [x] `docs/design-notes.md` §13 の導出更新（重み写像に計量、不動点の再導出、`ε` の意味）。
- [x] テスト: design §1（罰則なし列の厳密性）、§2（スケール不変。**現行実装と
      「分母の外」実装の両方で落ちることを先に確認**、群内不均一な `C` を含む）、
      §3（group-L0 不動点）、§4（参照計量の解析検算）、§6（`select_fit` の返り値が
      計量つき）、
      §7（μ₀ の λ→∞ 極限、3 推定器）、§11（縮約の網羅）、§12（拒否）、§13（OLS 不変）。
- [x] `make test-all` / `make test-pin` 緑。

### M3 — λ 選択（`feat(selection): cross_validate / gcv / effective_dof for MomentFit`）

- [x] `selection.jl`: `MomentCVResult`、`cross_validate(::MomentDataset, …)`、
      `effective_dof` / `gcv`（**`fit` と同じ三つ組を再構成**）、`_gcv_neff(::MomentFit)`。
      `_grouped_folds` 再利用、fold 数の扱いはエネルギー側踏襲。
- [x] テスト: design §8（分割性質 + **fold 学習側の軌道被覆**）、§9（リークガード）、
      §10（凍結列整合）。
- [x] `make test-all` 緑。

### M4 — docs / 記録 / レビュー

- [x] `docs/src/guide/fitting.md` / `moment.md`（罰則計量、`cost_weights` との役割分担、
      支持規則との 2 スケール共存（Q2）、GCV が config 相関を無視する旨）。
- [x] `SPEC.md` fitting 節、`docs/src/api.md` に 2 行、
      `CHANGELOG.md`（**BREAKING** — `fit.jl:30-41` の前例と同格、記録済み fit の
      再取得が要る旨を名指し）、`CLAUDE.md` coupled-site L490-512 の更新 +
      L497 の「design-notes §13」参照先修正、`docs/specs/README.md`。
- [x] `bench/bench_solver.jl` に計量構築（**ピークメモリも**）と λ 掃引、
      `BENCH_LOG.md` にエントリ。
- [x] `test/glmnet/` / `test/sunny/` / `examples/*.jl` の推定器構築箇所を grep
      （field surface ↔ ALL test environments）。
- [x] `make docs` strict 緑。
- [x] **Tier 2 レビューパネル（4 軸）実施 (2026-08-24): blocker 4 / major 12 / minor 24。**
      `2c720fb` で全件適用（見送りは性能 m6/m8 の 2 件、理由はコミット本文）、
      `56f4dc5` で `nconfig` 既定を実測コストから 2048 に。

### M5 — SLCE.jl へ移植（Q4 の段取り）

- [x] SLCE.jl の**純スピン両チャネル**に同型移植（SLCE `57fc4fc`）: `fitting/metric.jl`
      新設、計量の 6 サイト、収束判定、pointed の μ₀ 免除、`_effective_dof_free`、
      `_reduce_to_active` の ridge 系、`with_lambda`、provenance の門 5 箇所。
      **群重み `≥ 0` の緩和は移植していない**（本 spec で落とした通り）。
- [x] **上流固有の追加**: ASR / freeze 再パラメータ化の下では罰則が `Z'·Diagonal(D)·Z`
      に圧縮されるので、計量は**基底列**で添字づけたまま（`column_groups` と同じ規約）。
      計量なしの一様罰則は γ 空間で厳密に `λ·I` を通す分岐を残し、無重み fit の
      ビット等号を保った。`_free_directions` は `null(diag(√m)·Z)` を SVD で取る
      （純スピンの freeze なら `Z` は選択行列なので添字版に厳密に一致）。
- [x] **移植しなかったもの（意図的）**: M3 のモーメント側 λ 選択 API
      （`cross_validate(::MomentDataset,…)` / `MomentCVResult` / `gcv` /
      `effective_dof(::MomentFit)` / `_cv_fold_count`）。M5 の本文が挙げていないため。
      ledger に行を書いた。その結果 SLCE 側の `_effective_dof_free` は内部からしか
      到達できない = **本 spec の M1 が上流で着地した時とまったく同じ状態**なので、
      ゲートも同じく人工設計＋密ハット参照で書いた。
- [ ] joint（変位）側は**参照アンサンブルの定義を決めてから**別途。
      `|u|^{2k} R_{lm}(u)` は一様乱数スピンでは定義できないので、変位の参照分布を
      spec に起こす（本 spec のスコープ外）。**上流では門で拒否**:
      `penalty_metric(::SLCEBasis)` は変位付き基底を名指しで拒み、
      `force_weight > 0` の fit も計量付きでは拒む。
- [x] SLCE 側 `TEST_MODE=all` 緑 53375（+124）、docs strict 緑。

### M6 — 受け入れ確認（package 外、コミットしない）

報告書: `~/jijs/scefit/for_package/fe_bcc/3x3x3cubic/scefit/M6-penalty-metric.md`
（生ログ `M6-report*.txt` を同じ場所に置いた）。

- [x] bcc Fe `m_full_l2`（108 列 / 5400 行）で `cross_validate` の λ 曲線・選択 λ・μ₀。
      **計量の一番はっきりした効き目がここに出た**: λ→∞（1e8）で計量ありは
      μ₀ = 2.168 を保ち他を 1.6e-06 まで潰して `effective_dof → 1.05`
      （= 罰則なし列 1 本、`_edof_free` の λ→∞ 極限の実データ確認）、計量なしは
      同じ λ で **μ₀ = 0.00006**（参照モーメント 2.09 μB をまるごと失う、CV 1.38 μB）。
      選択 λ は計量あり 1e+1（CV 0.021711、OLS 0.021712 をわずかに下回る）、
      計量なし 1e+0（CV 0.021730、OLS より悪い）。
- [x] **λ の意味の変化を 7 基底で定量化**（本命）。`λ_j = ‖X_j‖²/m_j`（計量あり）/
      `‖X_j‖²`（なし）の広がり = 「1 つの λ が列によって何倍違う意味を持つか」:
      計量なし 16×〜**7022×** → 計量あり 3.1×〜56.4×、**圧縮率 5.2×〜124.6×**。
      計量なしで最も離れた 2 列はどの基底でも「3 体 (1,1,2) 対 2 体 (1,1)」で、
      罰則が高体数を優遇していたという主張の直接確認。
- [ ] ~~計量前後で `select_fit` が選ぶ群の順序が変わる~~ → **測れなかった（否定的結果）**。
      記録済みのどの基底でも**群がひとつも死なない**: bcc Fe l044_c4.1（21 列/18 群、
      16300 行）は λ = 1e9…1e-2 の 23 点で `n_alive = 18/18`、FeRh l044_c5.0
      （82 列/62 群、77770 行）は 25 点で計量あり 2 群が一瞬落ちるだけ・計量なし 0 群。
      最小 λ でまとめて「死ぬ」のは相対しきい値（その λ の最大 scaled magnitude の
      1e-6）の副作用で、群スパース性ではない。理由は圧倒的な優決定。
      **含意: コスト考慮 Pareto 選択は実データ上でまだ一度も発火していない。**
      再挑戦には群が死ぬ設計（4 体星型 or 行数に迫る列数）が要る — 対の spec の
      4 体実験（`cutoff_star` の体数別化が前提）と同じ入口。
- [x] Q5 乖離診断: `m_j / ‖X_j‖²_train` は bcc Fe l044 で 0.103–1.178（11.5×）、
      pointed `m_full_l2` で 0.028–0.770（27.1×）、FeRh l044_c5.0 で 0.427–24.08
      （56.4×）。訓練集合（MFA m 掃引 / AFM–FM 2 集合）が一様ランダムスピンから
      遠いためで、spec が Q5 で予告した通り。計量ありで残る λ_j の広がりは
      **これがそのまま出ている**。
- [x] 記録済み罰則付き fit の再取得 → **対象ゼロ**。`~/jijs/scefit` 配下の `*.jl`
      全数走査で `OLS()` 113 箇所、罰則付き推定量・`select_fit` / `select_support` /
      `refit` は **0 箇所**。bcc Fe も FeRh も klm もすべて OLS で、λ はどこにも
      入っていない。CHANGELOG の BREAKING は「これから罰則を使う仕事への警告」で
      あって、再取得すべき成果物は存在しない。

## Exit checklist

Run through every item once implementation lands. ~~Strike through~~ items
that do not apply.

- [x] `make test-all` passes (4 threads) — 38795; SLCE 側 `TEST_MODE=all` 53375。
- [x] `make test-pin` passes — 104、-t 4 / -t 1 両方。ピン再取得なし。
- [x] `make docs` builds (strict) — 両パッケージ。
- [x] If results changed: regression or validation test added, oracle
      independent of the implementation（解析閉形式 / 密ハット参照 / 独立恒等式）。
- [x] If public API changed: `SPEC.md` and `docs/src/api.md` updated（両パッケージ）。
- [x] If a hot path was touched: before / after recorded in
      `bench/BENCH_LOG.md`（SCEFitting 側 `bench_solver.jl`。SLCE 側は計量構築の
      ベンチを新設していない — 同じカーネルで、λ パスへの上乗せは列あたり 1 乗算）。
- [x] Tier 2 review panel run (numerical / maintainability / performance /
      API axes) and findings resolved（`2c720fb`, `56f4dc5`）。
- [x] ~~If module names or Makefile targets changed: `.claude/agents/` swept.~~
- [x] If this diverges from SLCE.jl: divergence ledger row in `CLAUDE.md`
      （移植後に 2 行へ書き換え: モーメント λ 選択 API と綴りの差分。joint 側は
      両パッケージとも未対応で、上流は門で拒否する）。
- [x] `CHANGELOG.md` `[Unreleased]` updated（**BREAKING**、両パッケージ）.
- [x] `Status:` line in this file and the table in `docs/specs/README.md`
      updated in sync.
- [x] Implementation commit hash appended below.

## Implementation commits

- SCEFitting.jl: `4d7ed7c` (M1) → `c6ccdde` (M2) → `e24b45e` (M3) → `6680d9e` (M4 docs)
  → `2c720fb` (Tier 2 パネル全件適用) → `56f4dc5` (`nconfig` 既定 2048)。
- SLCE.jl: `57fc4fc` (M5)。
