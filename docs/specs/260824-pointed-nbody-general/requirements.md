# Requirements: pointed moment basis の体数一般化（門は 4）

Status: draft (2026-08-24) — spec-reviewer 第 3 回の blocker 3 / major 7 を反映済み。
後方互換は制約にしない方針（ユーザー指示 2026-08-24）

## Goal

pointed モーメント基底の body order を「1/2/3 の直書き」から **任意 N の一般実装**に
書き換え、`MomentSpec` の門を `nbody ≤ 4` に置く。コードは一般、約束は検証が届いた
範囲だけ、という形にする。

## Background

- 現状 `MomentSpec` は `1 <= nbody <= 3` で拒否する（`src/basis/momentbasis.jl:137`）。
  上流 SLCE.jl も同一行（`SLCE.jl/src/basis/momentbasis.jl:87`）。
- 打ち切りは **3 箇所の直書き**に閉じている（`momentbasis.jl` 内で `body == 3` /
  `>= 3` が残るのは `:340, 344, 345, 374` のみで、下の 3 項がそれを尽くす）:
  1. `_moment_labels(spec, N)` の `N == 1/2/3` 分岐（`:180-203`、分岐は `:183-199`）。
     N ≥ 4 は空を返し、呼び出し側の `isempty(labels) && continue`（`:360`）が**黙って
     飛ばす**。
  2. `_pointed_star_candidates` の `for x, y = (x+1)` 二重ループと `_PERMS3`
     （`:226-270`、`_PERMS3` は `:224`）。
  3. **`MomentBasis` コンストラクタの星構築が body = 3 決め打ち**（`:340-347`）:
     `if spec.nbody >= 3` の単発分岐で `_orbits_from_members(…, 3)` /
     `push!(orbits, (3, k, O))`。①②だけ直して門を開けると **4 体軌道が 1 つも生えず、
     `nbody = 4` が `nbody = 3` と同一の基底になる**。SLCE.jl も同形（`:280-288`）。
- その下の機構は **すでに任意 N**（レビューで全経路を確認済み）:
  `AngularMomentum.coupling_paths`（`:182-210`）、`coeff_tensor_complex`、
  `CoupledBasis{R}`、`_orbit_salcs_decors(…, N::Int, …)`（`salcbasis.jl:518`）、
  `_member_sig(m, ::Val{N})`（`orbits.jl:102`）、`_orbits_from_members`、
  `_eval_term_mixed`（`salc.jl:483`）、`SALCScratch.z`（可変長）、
  `moment_resolvability` の env 検査（`allunique` ベース、`momentbasis.jl:632`）、
  `salc_groups(::MomentBasis)`（`momentfit.jl:646-669`）、`_design_moment`。
  射影・軌道・評価・分解可能性・群化・設計行列のどこにも `Val{N}` の上限も固定長バッファも
  帽子も無い。
- **エネルギー側の列挙層は N=4 まで独立検証済み**（`test/unit/test_ws_nbody.jl:131-150`、
  fcc の 4 体候補 192 個をオラクル固定）。一方 **SALC 射影層より先は N ≥ 4 が未検証** —
  スイート中の `nbody = 4` は 2 箇所（`test_momentbasis.jl:80`, `test_input.jl:303`）とも
  拒否テストで、`nbody ≥ 4` で基底を建てたテストは 1 本もない。
- **動機は実測**: bcc Fe 3×3×3 の pointed フィットは CV 0.0217 μB で飽和しており、床の
  原因が (a) 未実装の 4 体か (b) データ量かは**未確定**
  （`~/jijs/scefit/for_package/fe_bcc/3x3x3cubic/scefit/INVESTIGATION-moment-residuals.md`
  §4d）。周期境界条件の下では射程は表現可能性の制約にならないので「セルが小さいから」は
  候補から外れる。本 spec は **(a) を反証可能にする実験装置**でもある。
- 列数の増加は **M0 で実測済み（2026-08-24、scratchpad の使い捨てビルド）**:

  | 系 / 設定 | N=3 | N=4 | 倍率 |
  |---|---|---|---|
  | sc 27 原子, `cutoff_star = 1.1`（nn 殻のみ） | メンバ 2430 / 星軌道 2 / **全 11 列** | メンバ 12960 / 星軌道 2 / **全 17 列** | メンバ ×5.3、列 ×1.5 |
  | bcc Fe 3×3×3 (54 原子), `cutoff_star = 4.1` | メンバ 76,788 / P1 軌道 12,798 | メンバ **2,774,736** / P1 軌道 **115,614** | メンバ ×36、軌道 **×9.0** |
  | FeGe B20 (`lmax_env=[1,1]`, `cutoff_pair=3.0`) | 20 列 | 140 列 | ×7 |

  「約 13 倍」という第 1–3 回の目安は**軌道数で ×9.0** が実測値。sc 27 のメンバ数
  12960 は `27 × C(6,3) × 4! = 27 × 20 × 24` に厳密に一致する（列挙規約の裏取り）。
  時間は sc 27 で全ビルド 1.5 s / `moment_resolvability` 0.05 s なので
  `make test-all` のゲートとして成立する。

## 設計に効く物理事実

pointed ラベルは **マーク 1 個 + 環境 N−1 個**で、環境スロットは `SPIN` 因子なので
`l ≥ 1` が要る（`decor.jl:43`）。したがって `Σl ≥ N−1`。時間反転スクリーンは
マークのランクを含めた `Σl` の偶奇（`momentbasis.jl:185/190/195`）なので、切り上げて

> **pointed N 体セクターは `Σl = 2⌈(N−1)/2⌉` から始まる。**

N=3 は Σl=2（mark `l=0` + env `(1,1)`）、**N=4 は Σl=4**。素朴な対応物
mark `l=0` + env `(1,1,1)` は Σl=3 の奇で消える（整合性チェック: 3 本の `l=1` を
`L_S=0` に結合した唯一の不変量はスカラー三重積で擬スカラー）。

**結合は `L_S`（`Lf` ではない）**。decor エンジンの `isotropy` が screen するのは `L_S`
（`salcbasis.jl:396-406`）。pointed でも一致するのは、マークの DISP 因子が
`disp = (1,0)` = `|u|²R₀₀` で**角運動量ランク 0**（`decor.jl:28-30`）、射影器でも
`slots[j].factor.l == 0 && continue` で飛ばされる（`salcbasis.jl:442`）ため。

**スケールは `(4π)^(n_spin/2)`**、`n_spin` = スピン因子を持つ decor 数
（`salc.jl:470-471`）。**rank-0 マークなら `n_spin = N−1`**、rank ≥ 1 マークなら `N`。
`(4π)^(N/2)` ではない — オラクルの定数を導く実装者はここで `√(4π)` ずれる。

結合路数（`coupling_paths` で確認済み）: mark `l=0` + env `(1,1,2)` → `L_S = 0` の路は
**全サイト配置で一意**。mark `l=1` + env `(1,1,1)` → **3 本**（4 ベクトルの多重線形不変量が
`(a·b)(c·d)` の 3 組分けで張られることと一致）。

**M0 実測（2026-08-24）**: `lmax_mark = 2, lmax_env = [2], lsum = 4` での一般ラベル
列挙は N=1 → 2 本、N=2 → 3 本、N=3 → 4 本（現行の直書き分岐と一致）、
**N=4 → 2 本**（mark `l=0` + env `(1,1,2)` と mark `l=1` + env `(1,1,1)`）。
`Σl = 2⌈(N−1)/2⌉` は N=1…6 で `[0,2,2,4,4,6]` を返し、mark `l=0` + env `(1,1,1)` は
Σl=3 の奇で実際に落ちる。

**gate 3 のフィクスチャは P1（`NoSymmetry`）で取る**。`D` は Reynolds 射影が働く
`(assignment × path × Mf)` 空間の次元なので、サイト安定化群が自明なら `D = 1` になり、
期待定数は `N!·κ` ちょうど。既存ピンの非対称（P1 ペア `2!·√3 = 2√3` に射影因子が無く、
FeGe star `(0,1,1)` は `3!·(1/√3)·√3 = 6.0` と `1/√3` を持つ）はまさにこの差である。

## Scope

Includes:

- `src/basis/momentbasis.jl`: 門を `1 ≤ nbody ≤ 4` に。`_moment_labels` 一般化。
  `_pointed_star_candidates(crystal, nl, spec, N)` 一般化。**ctor の星構築を
  `for N = 3:spec.nbody` ループに**。`admit` の `body == 3` → `body >= 3`。
  `_subsets` / `_permutations` ヘルパ。`_PERMS3` 削除。ファイル冒頭の設計記録
  （`:30`「the 3-body cutoff is MARK–ENVIRONMENT-BOND based」）の一般化。
- `src/io/input.jl:63`、`docs/src/guide/io.md:99, 105`、`docs/src/guide/moment.md:56`、
  **`docs/src/api.md:329`**（`@docs` ブロックではなく手書き散文）、
  **`docs/src/theory/moment.md:289, 381, 385-408`**（低次形式の列挙が star `(2,1,1)` で
  止まっている）: 「3-body star」の記述を一般化。なお `theory/moment.md:375` は既に
  `(4π)^{n_spin/2}` を、`:407` は gate 3 の不変量を published している。
- テスト: `test/unit/test_momentbasis.jl`（N=4 ゲート群）、`test/unit/test_ws_nbody.jl`
  （`_ptstar_brute`）、`test/unit/test_input.jl`（拒否を `nbody = 5` に）、
  **`test/pin/`（pointed フィクスチャ新設 — `fixtures.jl` / `payload.jl` / `runtests.jl` /
  `capture.jl` / `pins/*.toml` / `PIN.md`）**、`test/parity/`（SLCE 同期）。
- `bench/bench_moment.jl` と `bench-moment` ターゲット（または M0 の一回計測を
  BENCH_LOG に記録する旨の明記 — design Q7）。
- docs / `SPEC.md` / `CHANGELOG.md` / `docs/specs/README.md`。
- 上流 SLCE.jl への**同時**移植。

Excludes:

- `nbody ≥ 5`（門で拒否。門の定数 1 つで後から開く）。
- エネルギー側 `BasisSpec` / `candidate_clusters` の挙動変更。
- 1–2 体経路の候補源変更（`build_clusters` のまま。多重度規約が別なので統合しない）。
- `cutoff_star` の精密化（design Q6）。
- 罰則計量・λ 選択（[`260824-penalty-metric`](../260824-penalty-metric/)）。

## Invariants

- **N ≤ 3 の設計列はビット一致**（pin 層で担保）。N = 1, 2 も不変。
- **メンバ乗数規約 = 全 `N!` 再アンカー順序**（CLAUDE.md coupled-site (1)）。
  **一様な順序落ちは比の恒常性では検出できない**ので、絶対正規化オラクルで押さえる。
- **星型の辺規則**: マーク–env の N−1 本だけを「最小像かつ `cutoff_star` 内」に縛り、
  env–env 辺は自由。エネルギー側の compact-cluster 規準とは意図的に違う。
- **`admit` の判定は置換軌道不変**（coupled-site (2)）。
- **各 pointed ラベルはマークをちょうど 1 個**（coupled-site (3)）。
- **`(4π)^(n_spin/2)` スケール**（`N/2` ではない）、`SALCKey` の順序と列アドレス、
  TR スクリーン、`isotropy` の極性、`MinimumImage` 固定、既存 TOML モデルの再読込・再予測。

## Completion criteria

- [ ] `make test-all`（4 threads）緑。
- [ ] **`make test-pin` 緑（pointed フィクスチャを追加した上で）**。現状 `test/pin/` は
      純スピンのエネルギー fixture のみで pointed 変更に対して**空虚に緑**になる。
      `PIN_SCHEMA` は**上げない**（上げると既存 5 本が全部赤）。`capture.jl` は
      `PIN_FIXTURES` を全再生成するので、**フィクスチャ絞り込みを足す**か、既存 5 本の
      再取得をユーザーの明示指示のもとで行う（`PIN.md` 規則 5 / CLAUDE.md「Always
      confirm — Recapturing regression pins」）。pointed の L2 は **`OLS()` のみ**で捕獲する
      （対の spec が触るのは罰則付き推定量なので先後どちらでも動かない）。
      L0 キーは `spin_ls` ではなく**マーク印つきの decor 多重集合全体**を印字する
      （`spin_ls` は純 DISP のマークを落とすので L0 が単射でなくなる）。
- [ ] **4 体が本当に生えている**: `nbody = 4` の列数 > `nbody = 3`、かつ `key.body == 4` の
      列が 1 本以上（③の直接ゲート）。
- [ ] **N=4 星候補 vs 独立ブルートフォース**: `_ptstar_brute` とメンバ集合が一致。
      セルにはタイのあるもの（`faces` / `fcc`）を必ず含める。**N = 3, 4 のみ**
      （N=2 は候補源が違う）。
- [ ] **N=4 絶対正規化オラクル**: 定数を**捕獲せず導出**して列と一致（design §3）。
- [ ] **N=4 TR（数値オラクル）**: `_design_moment(mb, [-e], [-ax]) == X`（bitwise）。
- [ ] **N=4 回転共変**: 共通回転を **`e` と `axes` の両方**に掛けて `L_S = 0` 列が不変。
- [ ] **N=4 resolvability**: **27 原子 sc 超胞**（`test_momentbasis.jl:463-505`）で
      refuse せず、`rank` が独立な乱数設計の SVD ランクと一致。**N=4 の spec を明示する**:
      既存 fixture は `nbody = 2, cutoff_pair = 1.8` で `a = 1.0` なら 26 近傍（6+12+8）を
      含み、そのまま 4 体星にすると `C(26,3)×4!×27 ≈ 1.7M` メンバで `make test-all` の
      ゲートにならない。⇒ **`cutoff_star = 1.1`（nn 殻のみ、`C(6,3) = 20`）と `lsum = 4`**
      を使う。マークの 6 nn は `L = 3` で 6 個の相異なる参照胞原子なので
      `allunique(env_atoms)`（`momentbasis.jl:632`）が成り立ち、ゲートは分類できる。
      **FeGe B20 基本胞では N=3 と同様に `@test_throws UnclassifiableBasis`**
      （`test_momentbasis.jl:302`、機構は `momentbasis.jl:615-621` の 108 対 98）。
      **bcc Fe 3×3×3 は完了条件にしない** — `test/` から到達できない
      （`bcc_fe(n)` は `bench/fixtures.jl:29` の別環境、`test/pin/fixtures.jl:18` は
      2 原子の conventional 胞）。M0 / M6 のインライン測定に回す。
      `rank` / `vanishing` / `null_combinations` の実測を本ファイルに書き戻す。
- [ ] **SLCE 列 parity（同期の変更検出器）**: `test/parity/` の既存許容
      （列ノルム `rtol = 1e-10`、max abs `1e-12`）で一致。**CI 無し・外部データ依存**、
      データ不在時は loud skip。
- [ ] `nbody = 5` が両パッケージで拒否され、メッセージが上限と理由を名指しする。
- [ ] `make docs` strict 緑。`SPEC.md` / `CHANGELOG.md` / `docs/specs/README.md` 更新。
- [ ] `bench`: N=3 不変（回帰）と N=4 の **(a) メンバ生成 / (b) 射影 /
      (c) `moment_resolvability`** を分けた時間・列数・TTFX。
- [ ] **既存 3 体オラクルの続投**: FeGe の `6.0` と P1 ペアの `2√3` が不変
      （design §3 が示すとおり、この 2 つが正規化の算術そのもの）。
- [ ] Tier 2 レビューパネル（4 軸）実施、`numerical-reviewer` の指摘は全件適用。
- [ ] **実験の受け入れ**（テスト外、対の spec 着地後）: 同じ 3×3×3 bcc Fe データで
      `m_full_l2`（108 列）+ 4 体（l≤2, star 4.1）の CV を測り、**0.0217 μB を割るか**を記録。

## References

- `src/basis/momentbasis.jl`（`MomentSpec` L79-170、`_moment_labels` L180-203、
  `_pointed_star_candidates` L207-270、ctor L311-412、星構築 L340-347、
  `admit` の body 分岐 L374、`isempty(labels)` L360、resolvability の機構 L615-621、
  env 検査 L632、ファイル冒頭の設計記録 L30）。
- `src/basis/{AngularMomentum.jl:175-245, salcbasis.jl:396-406/428/436/442/518,
  decor.jl:28-30/43/161-166, salc.jl:470-471/483, Harmonics.jl}`、
  `src/clusters/{enumerate.jl:19-29, orbits.jl:45-105}`、`src/fitting/momentfit.jl:646-669`。
- `test/unit/test_ws_nbody.jl:22-63`（`_wsnb_brute`）、`:131-150`、
  `test/unit/test_momentbasis.jl:149-213, 213-231, 233, 257-276, 302, 305-315, 463`、
  `test/unit/test_mixedsalc.jl:89-148`（デカルト不変量、特に `_Q5` L127-135）、
  `test/unit/test_normalization.jl:170-193`（O1→O3 の流儀）、
  `test/pin/{PIN.md, runtests.jl, payload.jl, capture.jl}`、`test/parity/runtests.jl:85-103`。
- CLAUDE.md coupled-site L293-320、L192-199。
- 対の spec: [`260824-penalty-metric`](../260824-penalty-metric/)。
- 実測の動機: `~/jijs/scefit/for_package/fe_bcc/3x3x3cubic/scefit/INVESTIGATION-moment-residuals.md` §4d。
