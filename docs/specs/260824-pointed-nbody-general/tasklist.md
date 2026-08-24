# Tasklist: pointed moment basis の体数一般化（門は 4）

Status: in progress (2026-08-25) — M0–M4 着地（実装コミット: 下記）、残るは M5 の bench と M6

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

- [x] 門 `1 ≤ nbody ≤ 4` + 理由入りメッセージ。`test_input.jl` を `nbody = 5` に。
- [x] N=4 ゲート 2–7, 9（独立ブルートフォース / 絶対正規化オラクル / 回転共変 /
      TR bitwise / resolvability / 4 体が生えている / 門）。
- [x] `make test-all` 緑。

- [x] **Exit**: `make test-all` 緑、N=4 ゲートが全て通る。

### M4 — SLCE.jl へ同時移植（`refactor(basis): …` + `feat(basis): …`）

- [x] SLCE.jl に M2 / M3 と同一の変更（SLCE `5684572`）。ついでに 57a22de のパネル修正
      （自己像オラクル・多重度カウント・符号ゲージ・無言打ち切り警告・threading）も
      同時に入れた — 分けて入れる意味がない。
- [x] `test/parity/` に `nbody = 4` ケースを新設（FeGe 2×2×2 / `cutoff_star = 2.6` /
      43 列うち 4 体 8）。**worst relative column deviation = 0.00e+00**。
      スコープ注記は「受入数値は N=3、体数の門は専用ケースが見る」に書き換え。
- [x] SLCE 側の `make test-all` 相当（`TEST_MODE=all`）緑 53375、docs strict 緑。

- [x] **Exit**: 両パッケージで全スイート緑、parity ゲート緑。

### M5 — docs / bench / レビュー（`feat(basis): per-star-order cutoff_star` + docs）

- [x] `docs/src/guide/moment.md`（`Σl = 2⌈(N−1)/2⌉` の一般則、星の辺規則の非対称と
      WS タイの注意、Q6 の 1 行、`:56`）、`docs/src/guide/io.md:99, 105`、
      `src/io/input.jl:63`、`src/basis/momentbasis.jl:30`、**`docs/src/api.md:329`**、
      **`docs/src/theory/moment.md:289, 381, 385-408`**、`SPEC.md`。
- [x] `CLAUDE.md` の ledger 行に Cartesian projector の所在
      （`test/unit/test_mixedsalc.jl`、`_Q5`）を追記。
- [x] **`cutoff_star` の体数別化**（M6 の前提。スコープ移動の理由は requirements /
      design Q6 に記録）。`MomentSpec.cutoff_star::Vector{Matrix{Float64}}`、
      アクセサ `_star_cutoff` / `_star_cutoff_envelope`、TOML の body-keyed テーブル
      （鍵は `3:nbody` を過不足なく）、SLCE.jl へ同時移植（`6bd3faa`）。
      スカラー / 行列のブロードキャストで**既存の綴りはビット等号**（pin 104 緑・再取得なし）。
      ゲートは「per-order カットの意味」= 合成性を両方向で主張（両パッケージ）。
- [x] `bench/bench_moment.jl` + Makefile `bench-moment` を新設し、N=3（回帰基準）と
      N=4 の **(a) メンバ生成 /(b) 軌道縮約 /(c) 全ビルド /(d) resolvability /
      (e) 設計行列・計量** 分離実測 + (f) TTFX（子プロセス）を `BENCH_LOG.md` に。
      `.PHONY` と `.claude/agents/profiler.md`、`bench/README.md` も掃いた。
      **判明したこと**: 壁は SALC 射影で、(a)+(b) はどちらの次数でもビルドの 2 % 未満。
      4 体は列 12 本のために射影 +3.7 s。`moment_resolvability` は既定 rtol の結果を
      基底にキャッシュするので、ベンチは `rtol = 1e-10` を明示して未キャッシュで測る
      （さもないとキャッシュヒットを測ってしまう）。
- [x] `CHANGELOG.md` `[Unreleased]`、`docs/specs/README.md` の行。
- [ ] **Tier 2 レビューパネル（4 軸）実施、`numerical-reviewer` の指摘を全件適用。**

- [ ] **Exit**: `make docs` strict 緑、Tier 2 の指摘が全件解決。

### M6 — 実験（package 外、コミットしない）

**⚠ M0 の実測で、spec が書いた形のままでは走らないことが判明した（2026-08-25）。**
`m_full_l2` は `cutoff_star = 4.1`（3NN）で、そのまま `nbody = 4` にすると
bcc Fe 3×3×3 で**星メンバ 2,774,736 本 / P1 軌道 115,614**（N=3 の 36 倍 / 9.0 倍）。
さらに `lsum` 無指定・`lmax_mark = lmax_env = 2` だと N=4 のラベルだけで 6 本あり、
列数は数万に達する — 行数は 5,400（100 config × 54 原子）しかないので**統計的にも無意味**。

⇒ 4 体は**最近接殻だけ**（2.451 Å, C(8,3) = 56 星/原子）で入れるのが唯一意味のある
第一プローブだが、`MomentSpec` の `cutoff_star` は**全 N で 1 枚**なので
「3 体は 4.1 Å、4 体は 2.5 Å」が現状**表現できない**。

- [x] **前提**: `cutoff_star` を体数ごと（エネルギー側の `cutoff[N-1][a,b]` と同じ形）に
      する。これが無いと M6 は走らない。→ M5 で着地。
- [ ] `m_full_l2` に 4 体（1NN 星のみ, `lsum = 4`）を足した variant で CV を測り、
      **0.0217 μB を割るか**を `INVESTIGATION-moment-residuals.md` §4d に追記。
      **列数と行数の比を必ず併記する**（過学習と改善の区別に要る）。

- [ ] **Exit**: §4d に測定結果（CV が 0.0217 を割るか、列/行比つき）が追記されている。
## Exit checklist

Run through every item once implementation lands. ~~Strike through~~ items
that do not apply.

- [x] `make test-all` passes (4 threads) — 38795; SLCE 側 `TEST_MODE=all` 53375。
- [x] `make test-pin` passes（**M1 の pointed フィクスチャを含めて**）— 104、-t 4 / -t 1
      両方。ピン再取得なし。
- [x] `make docs` builds (strict) — 両パッケージ。
- [x] If results changed: regression or validation test added, oracle
      independent of the implementation.
- [x] If public API changed: `SPEC.md` and `docs/src/api.md` updated.
- [x] If a hot path was touched: before / after recorded in
      `bench/BENCH_LOG.md`（メンバ生成 / 軌道縮約 / 全ビルド / resolvability /
      設計行列・計量 / TTFX を分けて — `bench/bench_moment.jl`, 2026-08-25）。
- [x] Tier 2 review panel run (numerical / maintainability / performance /
      API axes) and findings resolved（`57a22de`）。
- [x] ~~If module names or Makefile targets changed: `.claude/agents/` swept.~~
- [x] If this diverges from SLCE.jl: divergence ledger row in `CLAUDE.md`.
      （同時移植で決着 — ledger の該当行は 2026-08-25 の closed 節に移した）
- [x] `CHANGELOG.md` `[Unreleased]` updated（両パッケージ）.
- [x] `Status:` line in this file and the table in `docs/specs/README.md`
      updated in sync.
- [x] Implementation commit hash appended below.

## Implementation commits

- SCEFitting.jl: `b808181` (M1 pin) → `f297df9` (M2) → `679b6fa` (M3) →
  `57a22de` (Tier 2 パネル全件適用) → parity `nbody = 4` ケース（本コミット）。
- SLCE.jl: `5684572`（M2 + M3 + 57a22de のパネル修正を一括）。
