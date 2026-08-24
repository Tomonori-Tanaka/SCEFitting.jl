# Tasklist: pointed moment basis の体数一般化（門は 4）

Status: in progress (2026-08-24) — M0 実測完了、requirements に書き戻し済み

This file holds coarse-grained, commit-sized milestones. Day-to-day tracking
goes through `TaskCreate` in-session.

**前提**: M1–M5 は単独で進められる。[`260824-penalty-metric`](../260824-penalty-metric/)
が要るのは **M6（実験）だけ**（4 体が効いたのか過学習したのかを区別する測定器）。
ただし pointed pin の L2 は `OLS()` のみで捕獲すること — 対の spec が触るのは
罰則付き推定量なので、そうしておけば先後どちらでも pin は動かない。

## Milestones

### M0 — measurement first（コミットしない、scratchpad）

門だけ外した使い捨てビルドで測り、結果を requirements / design に書き戻す。
**M0 が終わるまで実装に入らない**。

- [x] N=4 のラベル集合と各 `Σl` を列挙。mark `l=0` + env `(1,1,1)` の不在を確認
      （`Σl = 2⌈(N−1)/2⌉` の一般則の直接確認）。
- [x] `(1,1,2)` の `L_S = 0` 路が全サイト配置で一意であることの裏取り（design §3）。
- [x] **gate 3 のフィクスチャ探し**: `D = 1`（安定化群が自明）になる結晶と `MomentSpec` を
      特定する。`D` は `salcbasis.jl:428/569` が決めるので、候補ごとに
      `(assignment × path × Mf)` 次元を実測する。`D = 1` が取れないなら定数に `1/√D` を
      入れて手計算するか、代替オラクル（design §3 末尾）に落とす判断をする。
- [x] 27 原子 sc 超胞 / bcc Fe 3×3×3 で **(a) メンバ生成 / (b) SALC 射影 /
      (c) `moment_resolvability`** の時間を**別々に**実測。**列数の倍率**（requirements の
      「約 13 倍」は M0 前の目安）と軌道数も測る。
- [x] FeGe B20 は N=3 で既に `UnclassifiableBasis` なので、N=4 でも拒否されることを確認。
- [x] 結果を requirements「設計に効く物理事実」「列数の増加」と design Q4 に反映。

- [x] **Exit**: requirements / design への書き戻しが済み、gate 3 のフィクスチャが決まっている。

### M1 — pin 層の pointed フィクスチャ（`test(pin): pointed moment fixtures`）

一般化の**前に**捕獲する。現状 `test/pin/` は純スピンのエネルギーのみで、pointed 変更に
対して空虚に緑になる。

- [x] `capture.jl` に**フィクスチャ絞り込み**を足す（`:8` は `PIN_FIXTURES` を全再生成
      するので、そのままだと既存 5 本の値と `[meta]` を書き換えてしまう）。
      既存 5 本を再取得する場合は `PIN.md` 規則 5 / CLAUDE.md に従い**ユーザーの明示指示**
      のもとで行う。
- [x] `fixtures.jl` に pointed 項、`payload.jl` に pointed 分岐（L0 キーは `spin_ls` では
      なく**マーク印つき decor 多重集合**、L2 の 4 スロットはモーメント版で埋める）、
      必要なら `runtests.jl` の分岐、`pins/*.toml`、`PIN.md` の項。
      **`PIN_SCHEMA` は上げない**。L2 は **`OLS()` のみ**。
- [x] `make test-pin` 緑（既存 5 本が不変であることを含む）。

- [x] **Exit**: `make test-pin` 緑、既存 5 本のピン値と `[meta]` が不変。

### M2 — 一般化（`refactor(basis): general body order for the pointed enumeration`）

門は 3 のまま据え置いて中身だけ一般化。**この時点で数値は 1 bit も動かない**。

- [x] `_subsets` / `_permutations`、`_PERMS3` 削除。
- [x] `_moment_labels` 一般化、`_pointed_star_candidates(…, N)` 一般化、
      **`MomentBasis` ctor の星構築ループ**（L340-348）、`admit` の `body >= 3`。
- [x] `salc.jl:482`（`_eval_term_mixed` の直前）に decorated の `D`（rank-0 マークなら
      `N`、rank ≥ 1 なら `N+1`）を追記。**`:293` の「`D` == body order」は純スピン専用の
      カーネルについて正しいので触らない**。
- [x] `make test-all` / `make test-pin` 緑、ピン再取得なし（M1 のフィクスチャで担保）。

- [x] **Exit**: `make test-all` / `make test-pin` 緑、数値が 1 bit も動いていない。

### M3 — 門を 4 に開ける（`feat(basis): allow nbody = 4 in the pointed moment basis`）

- [ ] 門 `1 ≤ nbody ≤ 4` + 理由入りメッセージ。`test_input.jl` を `nbody = 5` に。
- [ ] N=4 ゲート 2–7, 9（独立ブルートフォース / 絶対正規化オラクル / 回転共変 /
      TR bitwise / resolvability / 4 体が生えている / 門）。
- [ ] `make test-all` 緑。

- [ ] **Exit**: `make test-all` 緑、N=4 ゲートが全て通る。

### M4 — SLCE.jl へ同時移植（`refactor(basis): …` + `feat(basis): …`）

- [ ] SLCE.jl に M2 / M3 と同一の変更。
- [ ] `test/parity/` の列 parity（同期の変更検出器）を両側で緑に。
- [ ] SLCE 側の `make test-all` / docs strict 緑。

- [ ] **Exit**: 両パッケージで `make test-all` 緑、parity ゲート緑。

### M5 — docs / bench / レビュー

- [ ] `docs/src/guide/moment.md`（`Σl = 2⌈(N−1)/2⌉` の一般則、星の辺規則の非対称と
      WS タイの注意、Q6 の 1 行、`:56`）、`docs/src/guide/io.md:99, 105`、
      `src/io/input.jl:63`、`src/basis/momentbasis.jl:30`、**`docs/src/api.md:329`**、
      **`docs/src/theory/moment.md:289, 381, 385-408`**、`SPEC.md`。
- [ ] `CLAUDE.md:542` の ledger 行に Cartesian projector の所在
      （`test/unit/test_mixedsalc.jl`）を追記。
- [ ] `bench/bench_moment.jl` + Makefile `bench-moment` を新設し、N=3 不変（回帰）と
      N=4 の **(a)(b)(c) 分離**実測 + TTFX を `BENCH_LOG.md` に。新設しない判断なら
      「M0 の一回計測」であることとフィクスチャ・ハーネスを BENCH_LOG に明記する。
      Makefile の `.PHONY`（L11-12）と `.claude/agents/` の Makefile ターゲット参照も掃く。
- [ ] `CHANGELOG.md` `[Unreleased]`、`docs/specs/README.md` の行。
- [ ] **Tier 2 レビューパネル（4 軸）実施、`numerical-reviewer` の指摘を全件適用。**

- [ ] **Exit**: `make docs` strict 緑、Tier 2 の指摘が全件解決。

### M6 — 実験（package 外、コミットしない）

- [ ] 同じ 3×3×3 bcc Fe データで `m_full_l2` + 4 体（l≤2, star 4.1）の CV を測り、
      **0.0217 μB を割るか**を `INVESTIGATION-moment-residuals.md` §4d に追記。

- [ ] **Exit**: §4d に測定結果（CV が 0.0217 を割るか）が追記されている。
## Exit checklist

Run through every item once implementation lands. ~~Strike through~~ items
that do not apply.

- [ ] `make test-all` passes (4 threads).
- [ ] `make test-pin` passes（**M1 の pointed フィクスチャを含めて**）, or pins
      recaptured with the reason in `test/pin/PIN.md`.
- [ ] `make docs` builds (strict).
- [ ] If results changed: regression or validation test added, oracle
      independent of the implementation.
- [ ] If public API changed: `SPEC.md` and `docs/src/api.md` updated.
- [ ] If a hot path was touched: before / after recorded in
      `bench/BENCH_LOG.md`（メンバ生成 / 射影 / resolvability を分けて）.
- [ ] Tier 2 review panel run (numerical / maintainability / performance /
      API axes) and findings resolved.
- [ ] If module names or Makefile targets changed: `.claude/agents/` swept.
- [ ] If this diverges from SLCE.jl: divergence ledger row in `CLAUDE.md`.
      （本 spec は同時移植なので新規行なし — ずれたまま片方だけ入れない）
- [ ] `CHANGELOG.md` `[Unreleased]` updated.
- [ ] `Status:` line in this file and the table in `docs/specs/README.md`
      updated in sync.
- [ ] Implementation commit hash appended below.
