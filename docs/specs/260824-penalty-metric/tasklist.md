# Tasklist: 罰則計量の正当化と λ 選択（エネルギー / モーメント両チャネル）

Status: in progress (2026-08-24) — M1–M4 着地（実装コミット: 下記）

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

- [ ] SLCE.jl の**純スピン両チャネル**に同型移植（`_edof` / 計量の 6 サイト /
      収束判定 / 計量 0 による pointed の μ₀ 免除）。**群重み `≥ 0` の緩和は本 spec で
      落としたので移植しない**。
- [ ] joint（変位）側は**参照アンサンブルの定義を決めてから**別途。
      `|u|^{2k} R_{lm}(u)` は一様乱数スピンでは定義できないので、変位の参照分布を
      spec に起こす（本 spec のスコープ外、ledger に段取りを残す）。
- [ ] SLCE 側の `make test-all` / docs strict 緑。

### M6 — 受け入れ確認（package 外、コミットしない）

- [ ] bcc Fe `m_full_l2`（108 列）で `cross_validate` の λ 曲線・選択 λ・μ₀ を記録。
      **対の spec の 4 体実験のベースライン**。
- [ ] エネルギー側 l044 系列で λ の意味の変化を記録し、**計量前後で `select_fit` が
      選ぶ群の順序が変わる**ことを確認（`c_g` は軌道サイズ因子を通じて偶然分と部分的に
      共線なので、これは**期待された結果**であって偶然ではない）。
- [ ] 参照分布と訓練分布の乖離診断（Q5）: `m_j / Var_train[Φ_j]` の max/min を記録。
- [ ] 記録済み罰則付き fit の再取得: bcc Fe l02…l044 / `m_*`、FeRh、klm。

## Exit checklist

Run through every item once implementation lands. ~~Strike through~~ items
that do not apply.

- [ ] `make test-all` passes (4 threads).
- [ ] `make test-pin` passes, or pins recaptured with the reason in
      `test/pin/PIN.md`.
- [ ] `make docs` builds (strict).
- [ ] If results changed: regression or validation test added, oracle
      independent of the implementation.
- [ ] If public API changed: `SPEC.md` and `docs/src/api.md` updated.
- [ ] If a hot path was touched: before / after recorded in
      `bench/BENCH_LOG.md`.
- [ ] Tier 2 review panel run (numerical / maintainability / performance /
      API axes) and findings resolved.
- [ ] ~~If module names or Makefile targets changed: `.claude/agents/` swept.~~
- [ ] If this diverges from SLCE.jl: divergence ledger row in `CLAUDE.md`
      （M7 の段取りと、joint 側が未移植である期間を明記）。
- [ ] `CHANGELOG.md` `[Unreleased]` updated（**BREAKING**）.
- [ ] `Status:` line in this file and the table in `docs/specs/README.md`
      updated in sync.
- [ ] Implementation commit hash appended below.
