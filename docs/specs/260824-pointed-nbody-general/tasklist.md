# Tasklist: pointed moment basis の体数一般化（門は 4）

Status: landed (2026-08-25) — M0–M6 完了（実装コミット: 下記）。M6 の測定は
`~/jijs/scefit/for_package/fe_bcc/3x3x3cubic/scefit/INVESTIGATION-moment-residuals.md` §4f

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
      4 体は列 12 本のために射影 +2.2 s。`moment_resolvability` は既定 rtol の結果を
      基底にキャッシュするので、ベンチは `rtol = 1e-10` を明示して未キャッシュで測る
      （さもないとキャッシュヒットを測ってしまう）。
      **レビューパネル後**: 3 つのスレッドループを `:greedy` に（既定は連続チャンク
      分割で、体数昇順に並ぶ高コスト項目が全部最後のチャンクに落ちる）。ビット等号、
      4 体ビルド 5.73 s → 4.36 s。TTFX の初回計測は無効（`setenv` が子の環境を
      置換して 1 スレッドになり、計時窓がパッケージロードの外側だった）→ 破棄して
      `addenv` + `-t` + 親側計時で再取得（9.5 / 13.0 s）。
- [x] `CHANGELOG.md` `[Unreleased]`、`docs/specs/README.md` の行。
- [x] **Tier 2 レビューパネル（4 軸）実施 (2026-08-25)**: blocker 0 / major 15 /
      minor 24。numerical の全件（major 1 = docstring の精度数値が BENCH_LOG と入れ替え、
      minor 6 件）を適用。maintainability の major 7 件（SLCE 側 ledger 不在、
      ledger 行 2 件が事実誤り、census 単一パス未移植、`setenv` 環境破壊、
      guide の例と数値の不一致）、performance の major 4 件のうち 3 件、api の
      major 3 件を適用。
      **見送り 1 件（理由付き）**: performance major 4 = SLCE の `Z'DZ` を
      IRLS 反復ごとに再構築する件のバッファ hoist。ビット等号ではあるが、(i) 到達
      するのは ASR/joint 基底 × 罰則付き推定量の `fit` のみで `select_fit` は通らず、
      (ii) 対の spec の M6 でこのプロジェクトの記録済み fit は全部 OLS と判明しており、
      (iii) この経路を測るベンチが無いので before/after を記録できない（CLAUDE.md の
      bench 規則）。同レビュアー自身も代数的縮約のほうは「ビット等号を失うので採るな」
      としている。
- [x] **Exit**: `make docs` strict 緑（両パッケージ）、Tier 2 の指摘が全件解決
      （上記 1 件は理由付き見送り）。

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
- [x] `m_full_l2` に 4 体を足した variant で CV を測り、`INVESTIGATION-moment-residuals.md`
      **§4f** に追記（生ログ `M6-nbody4-report.txt`）。ハーネスは記録済み実行そのまま
      （同じ分割 `randperm(MersenneTwister(1), 100)`）で、**ベースライン行は
      `m_full_l2/fit_summary.txt` の 0.020907 / 0.021689 を厳密に再現**。

      | 基底 | 列 | うち 4 体 | 列/行 | fit | CV |
      |---|---|---|---|---|---|
      | A `m_full_l2` (3 体, star 4.1) | 108 | 0 | 2.0 % | 0.020907 | **0.021689** |
      | B A + 4 体 1NN (2.5 Å) | 164 | 56 | 3.0 % | 0.020202 | **0.021258** |
      | C A + 4 体 2NN (2.9 Å) | 448 | 340 | 8.3 % | 0.019110 | **0.021781** |
      | F 3 体 star 4.1 lsum 4（D/E の対照） | 91 | 0 | 1.7 % | 0.021010 | 0.021692 |
      | D F + 4 体 1NN, lsum 4 | 103 | 12 | 1.9 % | 0.020614 | **0.021308** |
      | E F + 4 体 2NN, lsum 4 | 163 | 72 | 3.0 % | 0.020222 | **0.021212** |

      **答え: 割る。ただし 2 %**（最良 E で −2.2 %、フル 3 体を保った B で −2.0 %）。
      そして **C で過学習として飽和が見える**: in-sample は 0.020907 → 0.019110 と
      下がり続けるのに CV は A より悪化（列/行 8.3 %）。§4e が 1NN 固定対照で見た
      「12 → 47 列で伸びない」の先が撮れた。副産物として **F ≡ A**（Σl=6 の 3 体
      17 列は CV に何も足していない）を再確認し、同じ列数予算なら角運動量より体数の
      ほうがましだが幅は 2 % のまま、も出た。⇒ §4d の候補 (i)(ii) は棄却、残るは
      (iii) データ量。
      **`lsum` の体数別化があれば D/E はもっと素直に測れる**（現状 `MomentSpec.lsum`
      は 1 枚で全体数に効くので、4 体を lsum 4 に抑えると 3 体の Σl=6 も道連れになる）。
      エネルギー側の `[interaction].lsum` は体数別テーブルを取るので形はある。本 spec
      のスコープ外 — F 行が寄与ゼロを示しているので結論は動かない。

- [x] **Exit**: §4f に測定結果（CV が 0.0217 を割るか、列/行比つき）が追記されている。
## Exit checklist

Run through every item once implementation lands. ~~Strike through~~ items
that do not apply.

- [x] `make test-all` passes (4 threads) — **38869**; SLCE 側 `TEST_MODE=all` **53439**
      （SLCE は Aqua がテスト環境に無く 1 error、`git stash` したベースラインでも同じ
      = 本作業とは無関係の環境欠落。Pkg 操作は不可なので未対処）。
- [x] `make test-pin` passes（**M1 の pointed フィクスチャを含めて**）— 104、-t 4 / -t 1
      両方。ピン再取得なし（`:greedy` も per-order `cutoff_star` もビット等号）。
- [x] `make test-parity` — 全ケース worst relative column deviation `0.00e+00`。
- [x] `make docs` builds (strict) — 両パッケージ。
- [x] If results changed: regression or validation test added, oracle
      independent of the implementation.
- [x] If public API changed: `SPEC.md` and `docs/src/api.md` updated.
- [x] If a hot path was touched: before / after recorded in
      `bench/BENCH_LOG.md`（メンバ生成 / 軌道縮約 / 全ビルド / resolvability /
      設計行列・計量 / TTFX を分けて — `bench/bench_moment.jl`, 2026-08-25）。
- [x] Tier 2 review panel run (numerical / maintainability / performance /
      API axes) and findings resolved（`57a22de` = 一般化分、および 2026-08-25 の
      第 2 回 = 移植 + per-order `cutoff_star` 分。見送り 1 件の理由は M5 に記載）。
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
- SLCE.jl: `5684572`（M2 + M3 + 57a22de のパネル修正を一括）→ `6bd3faa`
  （per-order `cutoff_star`）→ 第 2 回パネル適用（本コミット）。
