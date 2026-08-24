# Design: pointed moment basis の体数一般化（門は 4）

Status: draft (2026-08-24) — spec-reviewer 第 3 回の blocker 3 / major 7 を反映済み

## Summary

3 つの直書きを一般 N に書き換え、門を 4 に置く。**コードは任意 N、門は 4** にするのは、
後で 5 に上げるのが定数 1 つで済み、かつ検証が届いていない領域を API で約束しないため。

一般化の中身:

1. **ラベル**: マーク 1 個（rank `0…lmax_mark`）+ 環境 `N−1` 個（rank `1…lmax_env` の
   非減少多重集合）、`Σl` 偶かつ `≤ lsum`。現行 3 分岐はこの一般式の N=1,2,3 展開。
2. **候補**: マーク原子の近傍から `N−1` 個の**組合せ**、`(atom, shift)` の重複を弾き、
   並進シグネチャで畳んでから **全 `N!` 再アンカー順序**に展開。現行の `x<y` +
   `_PERMS3` はこの N=3 展開。
3. **ctor の星構築ループ**: `if spec.nbody >= 3` の単発分岐を `for N = 3:spec.nbody` に
   （これを忘れると `nbody = 4` が黙って `nbody = 3` になる）。

その下（結合路・テンソル・SALC 射影・評価カーネル・分解可能性・群化）は既に任意 N。

### 星型の辺規則を N を上げても変えない理由

エネルギー側は `N ≥ 3` で **`C(N,2)` 本すべて**が同時に最小像であることを要求する
（compact-cluster 規準、CLAUDE.md L192-199）。pointed は **マーク–env の `N−1` 本だけ**を
縛り、env–env 辺は自由（`momentbasis.jl:207-212` の M2-5）。これを N を上げても維持する:

- 各 env サイトの位置は「マークのセル + そのスポークの最小像」で**一意に決まる**。
  スポークが 2 本でも 3 本でも同じで、`(mark; env 原子多重集合)` に対して星が一意なので、
  二つの軌道が同じ単項式を担うことがない。スポーク本数に依らない議論。
- エネルギー側の反例（CLAUDE.md L194-196「i–j と i–k を最小化する像が j–k を長い像に
  追いやる」）は pointed には**当たらない** — 星は j–k を一度も問わない。エネルギー側に
  compact 規準が要るのは、クラスターに区別された中心が無いから。

**ただしこの一意性は各スポークの最小像が一意な場合の話。** WS 境界のタイでは 1 スポークが
複数の最小像を持ち、タイ由来のメンバ多重度は `(タイ多重度)^(N−1)` で増える（N=3 の二乗から
一段上がる）。エネルギー側と違い pointed には `_reduce_orbit_salcs` が無い
（decorated SALC は拒否される — ledger 行）ので、この縮退は `moment_resolvability` と
fit 時の列凍結だけが受け止める。⇒ ゲート 2 のセルにタイのあるもの（`faces` / `fcc`）を
必ず含め、さらに小さなタイセルで N=4 が `UnclassifiableBasis` を投げることを明示的に押さえる。

これは上流差分ではなく **pointed と energy の役割差**なので、ledger 行ではなく
`_pointed_star_candidates` のコメントと `docs/src/guide/moment.md` に据える。

## Module layout

| Target | Change |
|---|---|
| `src/basis/momentbasis.jl`（`MomentSpec`） | 門 `1 ≤ nbody ≤ 4`（メッセージが上限の**理由**を名指し）。docstring の `nbody` / `cutoff_star` 行、`Σl = 2⌈(N−1)/2⌉` の一般則。 |
| `src/basis/momentbasis.jl`（`_moment_labels`） | `N == 1/2/3` 分岐 → 一般式（mark rank × env 多重集合、偶 `Σl`、`≤ lsum`）。 |
| `src/basis/momentbasis.jl`（`_pointed_star_candidates`） | 署名に `N` を追加。`_subsets(n, N−1)` + `allunique` の `(atom, shift)` 検査 + 全順列展開。`_PERMS3` 削除。**全順列は `clusters/enumerate.jl:19-29` の `_ordered_subsets(n, n)` がそのまま**（「distinct index の順序付き k タプル全部」なので `k = n` が全順列）— 新設するのは `_subsets` だけ。 |
| `src/basis/momentbasis.jl`（`MomentBasis` ctor, L340-348） | 星構築を `if spec.nbody >= 3` から `for N = 3:spec.nbody` ループに。`_pointed_star_candidates(crystal, nl3, spec, N)` / `_orbits_from_members(crystal, sg, stars, N)` / `push!(orbits, (N, k, O))` の直書き `3` を `N` に。**1–2 体側の `min(spec.nbody, 2)` は据え置き**（候補源が `build_clusters` で多重度規約が別。統合しない）。`nl3` / `dmin2_star` は `cutoff_star` 1 枚なので N をまたいで共有。 |
| `src/basis/momentbasis.jl`（`admit`） | `body == 3` → `body >= 3`（ループは既に「マーク以外の全サイト」を回している）。 |
| `src/io/input.jl` | スキーマ docstring の `nbody` コメント。リーダーのコードは不変。 |
| `src/basis/salc.jl:482` | pointed が通るのは `_eval_term_mixed` の側（`_eval_term` は純スピン専用で decorated を拒否するので `:293` の「`D` == body order」は**正しい。触らない**）。`_eval_term_mixed` に「decorated では `D` = スロット数 = 環境 `N−1` + マークの因子数（rank 0 なら 1、rank ≥ 1 なら 2）= `N` または `N+1`」を追記。 |
| `test/unit/test_momentbasis.jl` | N=4 ゲート群（§Test strategy）。 |
| `test/unit/test_ws_nbody.jl` | `_ptstar_brute` を `_wsnb_brute` の隣に追加。 |
| `test/unit/test_input.jl` | 拒否テストを `nbody = 5` に。 |
| `test/pin/{fixtures.jl,payload.jl,runtests.jl,capture.jl,pins/*.toml}` + `PIN.md` | **pointed フィクスチャ新設**。`runtests.jl` はペイロード鍵を直書きしている（`:73` の `("n_salcs","keys","members","terms")`、`:104-114` の L0′/L1/L2 スロット、`:57` の `PIN_SCHEMA` 一致検査）ので、pointed は L2 の 4 スロットをモーメント版で埋める（`sum(abs2, _design_moment(...))`、`MomentFit` の `r2`/`coef`、ホールドアウト `predict_moment`）か `runtests.jl` を分岐させる。**`PIN_SCHEMA` は上げない**（既存 5 本が全部赤になる）。`capture.jl:8` は `PIN_FIXTURES` を全再生成するので**フィクスチャ絞り込み**を足す。L0 キーは `spin_ls`（純 DISP のマークを落とす）ではなくマーク印つき decor 多重集合。L2 は `OLS()` のみ。 |
| `test/parity/runtests.jl` | SLCE 同期チェック（列ビット一致）。**CI 無し・外部データ依存**。 |
| `docs/src/guide/moment.md:56` / `io.md:99, 105` / `src/io/input.jl:63` / `src/basis/momentbasis.jl:30` | 「3-body star」の記述を一般化。 |
| `docs/src/api.md:329` | **編集必要**（第 2 回の「編集不要」は誤り）。`@docs` ブロックではなく**手書きの散文**「3-body stars are cut on the two mark–environment bonds」がある。 |
| `docs/src/theory/moment.md` | `:289`（「A three-body star is cut on its two mark–environment bonds」）、`:381`（比較表の行）、§"Explicit low-order forms" `:385-408`（1体/2体/star `(0,1,1)`/star `(2,1,1)` で止まっている）に N=4 の項を足す。**`:375` が既に `(4π)^{n_spin/2}` を、`:407` が gate 3 の不変量を published しているので、そこを引く**。 |
| `SPEC.md` / `CHANGELOG.md` / `docs/specs/README.md` | 更新。 |
| `bench/bench_moment.jl` + Makefile `bench-moment`（`.PHONY` L11-12 も） | モーメント側のベンチが存在しない（`bench/` は `bench_{clusters,design_matrix,end_to_end,nd2fe14b,salcbasis,solver}.jl` のみ）。新設して (a) メンバ生成 / (b) 射影 / (c) `moment_resolvability` を分けて測る。新設しない判断なら BENCH_LOG に「M0 の一回計測、フィクスチャと計測ハーネスは …」と明記する（Q7）。 |
| upstream `SLCE.jl` | 同一変更 + 列 parity。 |

## API

```julia
MomentSpec(; …, nbody = 3, …)   # 1 ≤ nbody ≤ 4（既定は 3 のまま）
```

既定は **3 のまま**。`nbody` は打ち切りの**宣言**であって、既定で 13 倍の列と
`Σl = 4` の内容を黙って持ち込むべきものではない（後方互換の話ではなく、既定値の意味の話）。

拒否メッセージ（案）:

```
nbody must be in 1:4; got 5. The enumeration and the SALC projection are written for
general N, but only N <= 4 is covered by the test oracles — raising the cap without
extending them would promise an unverified region.
```

private ヘルパ:

```julia
_subsets(n::Int, k::Int)::Vector{Vector{Int}}          # 1:n の昇順 k-部分集合
_permutations(n::Int)::Vector{Vector{Int}}             # 1:n の全順列（辞書順、決定的）
_moment_labels(spec::MomentSpec, N::Int)               # 署名不変、中身を一般化
_pointed_star_candidates(crystal, nl, spec, N::Int)    # ← N を追加
```

## Types and conventions

- **物理規約の変更はゼロ**。TR（偶 `Σl`）、`isotropy` の極性、`MinimumImage` 固定、
  `SALCKey` の順序 — すべて不変。
- **スケールは `(4π)^(n_spin/2)`**、`n_spin` = スピン因子を持つ decor 数
  （`salc.jl:470-471`）。**rank-0 マークの pointed ラベルでは `n_spin = N−1`** で、
  `(4π)^(N/2)` では**ない**。オラクルの定数を導く実装者はここで `√(4π)` ずれる。
- **記述は `L_S` に統一**。decor エンジンの `isotropy` が screen するのは `L_S`
  （`salcbasis.jl:396-406`）。pointed でも `Lf ≡ L_S` になるのは、マークの DISP 因子が
  `disp = (1,0)` = `|u|²R₀₀` で角運動量ランク 0（`decor.jl:28-30`）、射影器でも
  `slots[j].factor.l == 0 && continue` で飛ばされる（`salcbasis.jl:442`）ため。
- **新しい不変条件**: 候補源は「1 物理インスタンスあたり `N!` 再アンカー順序」を出す。
- **上流差分 ledger**: SLCE.jl と**同時**に同じ変更を入れるので新規行は立てない。
  ずれたまま片方だけ入れることは parity ゲートで禁止する。

## Impact on coupled sites

- [x] `Harmonics` / `AngularMomentum` / SALC 射影 / 正規化テスト: **射影層を N=4 で初めて
      踏む**。`coupling_paths` / `coeff_tensor_complex` / `complex_to_real_tensor` は一般再帰
      だが未検証 → 絶対正規化オラクル・回転共変・TR bitwise で押さえる。
- [x] `SALCKey` 順序 / 設計列 / TOML 永続化 / `coeftable`: N ≤ 3 バイト等号（pin 層）で不変を
      担保。N=4 の新キーは既存の並びの後ろに付く（`body` が第 1 キー）。
      **`salc_groups(::MomentBasis)`（`momentfit.jl:646-669`）は N 一般なので
      `GroupAdaptiveRidge(mb)` が新しい群を自動で拾う** — 対の spec と直結。
- [ ] エネルギーカーネル ↔ トルクカーネル: 影響なし（pointed は torque を持たない）。
- [x] Decor エンジン ↔ pointed 基底 ↔ `momentfit.jl` の門: (1) 乗数規約、(2) `admit` の
      置換軌道不変性、(3) マーク 1 個 — 3 つとも N 一般で維持。`moment_resolvability` の
      env 検査は `allunique` ベースで既に N 一般だが、**タイ多重度が `(タイ)^(N−1)` で
      増える**ので拒否の到達経路が広がる（§星型の辺規則）。
- [ ] Readers' `zero_moment_atol` ↔ dataset doors: 影響なし。
- [x] 上流差分 ledger: 新規行なし、ただし SLCE 同時変更が必須。
- [ ] `sce/introspect.jl` `multipole_terms` ↔ SCEMonteCarlo の `TiledHamiltonian` 取り込み:
      N/A（pointed 基底は永続化されず、`multipole_terms` は decorated 基底を拒否する）。
- [ ] `make test-downstream`: N/A（同上）。
- [x] `.claude/agents/`: Makefile に `bench-moment` を足すので**掃く**（agent 定義が
      Makefile ターゲットを名指ししている）。
- [x] `SPEC.md` / `docs/src/api.md`: `api.md` は bare binding が拾うので編集不要、
      `SPEC.md` の `nbody` 記述を更新。

## Test strategy

1. **N ≤ 3 バイト等号（pin 層）** — `test/pin/` に pointed フィクスチャを新設し、
   `keys` と `X` を `PIN.md` の protocol（`PIN_DATE` / `PIN_COMMIT` / 再取得規則）に
   従って固定。変更検出器であることを明示。**ad-hoc な `==` を `test_momentbasis.jl` に
   置かない**（pin の層構造を二重化しない）。
2. **星候補 vs 独立ブルートフォース（独立実装オラクル）** — `_ptstar_brute(cr, spec, N)` は
   `_wsnb_brute`（`test_ws_nbody.jl:22-63`）と同じく**生の格子並進から最小像を再導出し
   `NeighborList` を触らない**。定義は anchored + 順序付き N タプルの形で書く:
   「`shifts[1] = 0`、マーク種のサイトから出る N−1 本のスポークが全て最小像かつ
   `cutoff_star[species pair]` 内」。こう書くと `N!` 多重度が**定義から落ちてくる**。
   **ただし `_wsnb_brute` の「原子が相異」規則は継承しない** — `_wsnb_brute:56` は原子の
   重複を弾くが `_pointed_star_candidates:254` は**厳密な `(atom, shift)` の重複だけ**を
   弾き、同一近傍の 2 つの像は意図的に残す。逐語移植すると、本 spec が必須としている
   タイセル（`faces` / `fcc`）でちょうど赤になる。
   **相対タイ帯を production と揃える**: `_pointed_star_candidates` は cutoff 判定と
   最小像判定の**両方**に `fac = 1.0 + nl.tol` を掛ける（`momentbasis.jl:230, 232`、
   `_SAME_DIST_RTOL = 1e-8` は `geometry/neighborlist.jl:105`）。ブルートフォースが同じ帯を
   使わないと、必須としているタイセルでちょうど割れる。
   **N = 3, 4 のみ**（N=2 は候補源が `build_clusters` で多重度規約が別 — 混ぜるのは
   spec 自身が警告している罠）。エネルギー側の `N ≤ n_atoms(cr)` 制限も**当たらない**
   （`faces`/`generic`/`hex` は 3 原子だが 4 体星は張れる）。セルは 4 種 × cutoff 3 通り。
3. **N=4 絶対正規化オラクル（`test_normalization.jl` O1→O3 の流儀。`test_momentbasis.jl`
   に置く）** — 定数を捕獲せず**導出する**。比の恒常性だけでは**一様な**順序落ち
   （全軌道で 24 中 12 しか出さない等）を検出できない（比は不変で、ピンした定数だけが動く）
   ので、Invariant 2 の本ゲートはこちら。mark `l=0` + env `(1,1,2)`（`L_S = 0` の路が
   全サイト配置で**一意**）について:
   1. **閉形式を直接書く。** ランク `(1,1,2)` の唯一の `L_S = 0` 不変量はデカルトで
      `e_j·Q(e_l)·e_k = (e_j·e_l)(e_k·e_l) − (1/3)(e_j·e_k)`、`Q(e) = e eᵀ − I/3`。
      閉形式は内積だけで書けるので**追加部品は要らない**。同じ不変量が
      `docs/src/theory/moment.md:407` に star `(2,1,1)` として既に published されている
      （`l=2` が env でなくマークに載る形）ので、新規に導出するのではなく**それを引く**。
      `test_mixedsalc.jl:129-133` の `_Q5` は同じ空間だが、使うなら
      `test/unit/testutils.jl`（`runtests.jl:35` で最初に include）へ移すこと — 現状
      `test_momentbasis.jl` から見えるのは include 順の偶然（`runtests.jl:49` < `:62`）。
      （`test/oracle/` に射影器は無い。`CLAUDE.md:542` の ledger 行が
      「the ~45-line Cartesian projector」と場所を書いていないのが第 2 回の誤りの原因なので、
      `(test/unit/test_mixedsalc.jl)` を追記する。）
   2. **期待定数 = `N!`（順序規約） × `1/√D`（Reynolds 射影の単位固有ベクトル） ×
      `κ`（step 1 の正規化定数）**。`D` = そのブロックの `(assignment × path × Mf)` 次元
      （`salcbasis.jl:428`）で、`assignments = unique([t[p] for p in perms])`
      （`:569`、`perms` は `_stabilizer` のサイト置換）が決める。
      `κ` は tesseral 定数（`Harmonics.jl:30-32`）× 単位 Frobenius の結合テンソル
      × `(4π)^(n_spin/2)`（rank-0 マークなら `n_spin = N−1 = 3`）。

      **`1/√D` を落とすと定数が間違う。既存の 2 つのピンがそれを示す**:
      - O1/O3（P1 ペア、`D = 1`）: `2√3 = 2! × √3`。射影因子は無い — `D = 1` だから。
      - FeGe pointed star `(0,1,1)`（`test_momentbasis.jl:185` の `6.0`）: nn Fe₃ 三角形の
        C₃ 安定化群が 3 つのマーク配置を 1 つの assignment 集合に束ねるので `D = 3`、
        固有ベクトルは `(1,1,1)/√3`。`3! × (1/√3) × √3 = 6.0` — **ピン値と厳密一致**。
        `1/√D` を落とすと `3! × √3 = 6√3 ≈ 10.39` になる。
      ⇒ 「既存 N=3 オラクルと O3 は同じ形」は**誤り**。O3 に射影因子が無いのは
      `D = 1` だからで、pointed star には有る。

      **フィクスチャは `D = 1` になるものを選ぶ**（安定化群が自明な星）。そうすれば
      定数は `N! × κ` ちょうどで手計算が閉じる。ゲートは先に「そのラベル・その軌道の
      SALC が 1 本だけ」と「`s.key.block == 1`」を `@test` して `D = 1` を主張する。
      参考: `(1,1,2)` の `κ = (4π)^{3/2}·(1/√5)·(3/4π)·√(15/8π) = 3√(3/2)`
      （`1/√5` は `(1,1,2)→0` テンソルの単位 Frobenius 正規化 `‖T‖² = Σ_r‖Q_r‖_F² = 5`、
      `√(15/8π)` は加法定理が決める `Z_2m ↔ Q(e)` の定数）。`D = 1` なら列は
      `24 · 3√(3/2) · Σ_stars [...]`。
   3. **列 = 期待定数 × 手書き幾何和**、和はマーク原子 `a` を中心とする
      **（星 × その SALC の `assignments` 集合）**にわたる。
      **`N!` の順序展開はランク割り当てを配らない** — `_connect_all`
      （`salcbasis.jl:119-144`）の `perm` は `(atom, shift)` で整列するのでどの順序でも
      原子→decor 写像は同じであり、`_canonicalize_members`（`salc.jl:155-165`）が
      同一 canonical slot list に足し込む（＝ちょうど `N!` 倍、それ以外は無い）。
      1 列の中で env ランク割り当てを足し合わせるのは**安定化群だけ**。
      ⇒ `D = 1` のフィクスチャなら和は**星のみ**（`l=2` が載る env サイトは軌道が固定）に
      なり、手書き参照が一意に書ける。「(星 × env ランク割り当て)」で和を取ると 3 倍になる。
      全マーク原子 × 乱数配置 2 本で確認。**フィクスチャ（結晶と `MomentSpec` 全体）を
      spec に名指しする** — `D` を決めるので load-bearing。
   参照は `s.members` を一切読まない（O3 と同じ規律）。
   **代替**（(1,1,2) が使えない場合）: mark `l=1` + env `(1,1,1)` の列が 4 ベクトル不変式の
   3 次元空間 `Σ(a·b)(c·d)` に載ることを乱数配置の最小二乗で確認。`L_S=0` の路数 3 =
   空間の次元 3 が独立チェックになる。**ただしこれはスケール不変なので一様な順序落ちを
   検出できない**（全列を半分にしても同じ空間に載る）— gate 3 の歯は消える。
   採るなら絶対定数も併せて要求する（`1/√D` 込み、`n_spin = 4` なので `(4π)²`）か、
   「構造チェックであって gate 3 の代替ではない」と明記して歯を別の場所から持ってくる。
4. **回転共変（解析オラクル）** — **任意の SO(3) 回転**を `e` と `axes` の**両方**に掛けて
   `L_S = 0` 列が不変（`rtol = 1e-12`）。スピンだけ回すとマーク軸が取り残される。
   これは既存 `test_momentbasis.jl:213-231` の**空間群共変**（`sg.map_sym` の原子ラベル
   付け替えを伴う）とは**別の、より強い**主張なので、同じ testset に足すとしても
   別ゲートとして書く。
5. **TR（数値オラクル）** — `_design_moment(mb, [-e], [-ax]) == X`（bitwise、既存
   `test_momentbasis.jl:233` に N=4 を追加）。奇 `Σl` 列が混入すれば符号が反転する。
   ラベルの `iseven(Σl)` 確認と mark `l=0` + env `(1,1,1)` の不在確認は**構造の副次
   チェック**として残すが、これは `_moment_labels` の実装の言い換えであってオラクルではない、
   とコメントに明記する。
6. **resolvability（独立オラクル）** — 既存「pointed resolvability gate」
   （`test_momentbasis.jl:257-276`、既存 2 本はいずれも `nbody = 2`）に N=4 を追加:
   記号ランク ≡ 乱数設計の SVD ランク、スペクトルギャップ付き。
   **refuse しないことを求めるのは 27 原子 sc 超胞（`:463` の fixture）と bcc Fe 3×3×3**。
   **FeGe B20 基本胞は N=3 で既に拒否される**（`:302` の `@test_throws`、機構は
   `momentbasis.jl:615-621` の 108 対 98）ので、N=4 でも `@test_throws` を期待値にする
   （スポークが 3 本になれば 2 つの env が同一参照原子に落ちる確率は上がるだけ）。
   `rank` / `vanishing` / `null_combinations` を requirements に書き戻す。
7. **4 体が生えている** — `nbody = 4` の列数 > `nbody = 3` の列数、かつ `key.body == 4` の
   列が 1 本以上（B1 の直接ゲート）。
8. **SLCE 列 parity（同期の変更検出器 — 逐語移植なので独立実装オラクルではない）** —
   `test/parity/`。既存の許容（`runtests.jl:85-103`: 列ノルム `rtol = 1e-10`、
   max abs `1e-12`）を踏襲する（「ビット一致」ではない）。CI ジョブ無し・外部データ依存で、
   データ不在時は loud skip。
9. **門** — `nbody = 5` を両パッケージで拒否、`nbody = 4` は通る。
10. **既存 3 体オラクルの続投** — FeGe の 6.0 / 2√3 が不変。

## Risks and open items

- **Q1 門を 4 にするか無制限か。** 推奨 **4**。検証していない領域を API で約束しない。
- **Q2 既定 `nbody` を 3 のままにするか。** 推奨 **3 のまま**（§API）。
- **Q3 N=4 オラクルのラベル。** `(1,1,2)` の路一意性は確認済みなので M0 は裏取りのみ。
  代替も §3 に明記済み。
- **Q4 コスト — 律速は 2 系統ある。**
  (a) **メンバ数** `C(z, N−1) × N!`（transport と評価に効く）。`z ≈ 12` で N=3 が 396、
  N=4 が 5280。
  (b) **射影**（実際の律速）: Reynolds 射影の次元は `D = 配置数 × 路数 × (2Lf+1)`
  （`salcbasis.jl:428`。`L_S` と一致するのは rank-0 マークの偶然による）、配置数（`assignments`）は
  **メンバ安定化群の軌道**が上限（`salcbasis.jl:569`。`_multiset_arrangements` が決めるのは
  ブロック数であって射影器の `D` ではない — 混同すると `O((N!)³)` に読める）、`P` の組み立てに `|stab|` 因子
  （`:436`）、`eigen(Symmetric(P))` が `O(D³)`、内側テンソルは `Π(2lᵢ+1)·(2Lf+1)`
  （rank-0 マークなら `(2l+1)^{N−1}·(2Lf+1)`）。
  さらに (c) `moment_resolvability` は毎 `MomentDataset` 構築で走り、QR/SVD 項が伸びる
  （倍率は**見積もらず M0 で実測**する）。
  ⇒ M0 で (a)(b)(c) を**別々に**計測し、bench にも分けて載せる。想定外に重い場合は
  N=4 の既定 cutoff を下げる指針を docs に書く（実装は変えない）。
- **Q5 `_eval_term_mixed` の `Val(D)` 特殊化**（pointed が通るのは `_eval_term` ではなく
  `_eval_term_mixed`、`salc.jl:483`）。`D` = スロット数で、**rank-0 マークなら `N`**、
  rank ≥ 1 マークなら `N+1`。N=4 でコンパイル時間が増える。実測して bench に記録。
  TTFX 対策は本 spec では**やらない**（ビット一致が最優先）。
- **Q6 `cutoff_star` の精密化。** ここには**二つ別のもの**がある。
  - **スポーク別**（1 本目は近く 2 本目は遠く）は**定義できない** — サイトが置換で
    入れ替わる（ラベルは多重集合）。置換不変で意味のある精密化は「スポーク長の**総和**」や
    「クラスター直径」の上限。**この結論は不変で、実装しない**。docs に 1 行残す。
  - **体数別**（3 体は 4.1 Å、4 体は 2.5 Å）は置換不変で定義でき、エネルギー側の
    `cutoff[N-1][a,b]` と同じ形。当初スコープ外にしたが、**M0 の実測で M6 がこれ無しには
    走らないことが判明した**（3NN 半径のまま 4 体にすると bcc Fe 3×3×3 で
    2,774,736 メンバ / 115,614 P1 軌道 対 行数 5,400）ので **2026-08-25 にスコープへ移した**。
    実装: `MomentSpec.cutoff_star::Vector{Matrix{Float64}}`（体数 `N` は添字 `N-2`）、
    アクセサ `_star_cutoff` 1 箇所、星の近傍リストは `_star_cutoff_envelope` で 1 本。
    スカラー / 行列は全星次数へブロードキャストするので**既存の綴りはビット等号**。
    ゲートは「per-order カットの**意味**」= 合成性（`[r3, r4]` の体 b 内容 ≡ 単一半径
    `r_b` の体 b 内容）を両方向で主張する。
- **数値結果への影響**: N ≤ 3 では**ゼロ**（pin 層で担保）。N = 4 は新しい列で、既定
  `nbody = 3` なので明示的に上げた人だけが見る。
- **本 spec 単独では実用にならない。** 着地順は **罰則（対の spec）→ 一般化**。
  罰則側は単独で価値があり、かつ 4 体実験の測定器でもある。
