# Requirements: cross-orbit alias groups (Wigner–Seitz boundary ties across orbits)

Status: landed (2026-09-11)

## Goal

同じ参照セル原子集合を、空間群が関係づけない複数の周期像で結ぶ**別々の軌道**
(cross-orbit alias)を基底構築時に構造的に検出し、フィットでは 1 つの設計列に結び、
読み出しでは各軌道に**規約(convention)として等分**した係数を返し、その事実を
`coeftable` / TOML / 読み出し警告で開示する。

## Background

- Nd₂Fe₁₄B 1×1×1(68 原子、P4₂/mnm 16 操作、等方 2 体 cutoff `inf`)で、別々のペア軌道
  10 組の設計列が完全に一致し、X_T の階数が l チャネルごとに 10 不足する
  (`~/jijs/scefit/2-14-1/nd2fe14b/1x1x1/scefit/SUMMARY.md`、`pair_aliases.jl`)。
  各組は同じ原子対集合を同じ距離で、Δz = ±c/2(9 組)または面内 Δx = Δy = ±a/2(1 組)の
  2 像で結ぶ。「+ 像」と「− 像」を入れ替える操作が群に無いため別軌道になっている。
- 現状の唯一の門は `OLS` の相対ランク警告(`estimators.jl` `_OLS_RANK_RTOL = 1e-10`)。
  罰則付き推定器はランク落ちを黙って通し、`salc_groups` が 2 軌道を別群にするので
  `select_support` / 罰則計量が配分を**黙って決める**(片方だけ落とせば 1:0)。
- 引き継ぎ `/tmp/HANDOFF-2026-09-11-cross-orbit-tie-aliases.md`(ユーザー決定 = 1 相互作用と
  して推定し等分で読み出す)。3 体議論の記録 `~/Packages/STUDY-2026-09-11-cross-orbit-ties.md`。
- SLCE.jl は逆方針(共有軌道を丸ごと凍結、`unresolvable_columns`)。本 spec は SCE 系のみ。
  発散台帳に行を足す。
- 将来のスピンスパイラル(一般化 Bloch 定理)訓練データは L_S = 0 のタイを分離できる
  (条件 q·ΔR ∉ 2πZ かつ d₊ + d₋ ≠ 0)。そのため**基底では併合も削除もしない**(明示 R の
  軌道を保持)。Phase 2(spiral datum、回転対応カーネル)は本 spec の範囲外。

## Scope

Includes:

- `alias_groups(basis)`:構造検出(データ非依存)。主判定 = 同一チャネル
  `(body, decors, L_S, Lf)` 内で同じ原子多重集合の集合を持つ別軌道;比例性(集約
  shift-blind 単項式ベクトル、`_function_vector`)は「等分が定義できる」証明書。
  非比例(張る空間だけ潰れる)は `:span_collapsed` として報告のみ、結ばない。
- 設計列 = 結んだ列。`_design_energy` / `_design_torque` は SALC 列を重み付き加算して
  `n_columns(basis)` 本にする。無エイリアス基底では**バイト同一**(折り畳みをスキップ)。
- `SCEFit.jphi` は列空間(長さ `n_columns`)、`SCEPredictor.jphi` は SALC 空間(長さ
  `n_salcs`)。`SCEPredictor(f)` が展開(各 SALC の係数 = 結んだ係数 × 重み)。
- `SCEPredictor.split::Vector{Symbol}`(SALC ごと `:free | :convention | :legacy`)、
  `coeftable` に `alias_group` / `split` 列、TOML スキーマ v7(`couplings[*].split`)。
- `salc_groups(basis)` は列ごとのラベル(エイリアス群をまたぐ軌道群を併合)、
  `group_costs` / `cost_weights` / `penalty_metric` は列空間で計算。
- `multipole_terms` / `bilinear_terms` / `to_sunny` は `:convention` / `:legacy` の非零係数
  があれば一度警告。
- `SCEBasis(...; alias_rtol = 1e-6)`、`alias_rtol = nothing` で検出無効(変異テスト用)。
- 基底構築時の `@info`(群の数、軌道対、原子対、距離、像の格子ベクトル差、倍化すべき方向)。
- docs: `theory/resolvability.md` に SLCE の「When symmetry does not fuse the tie」相当節
  +本パッケージの扱い(P4₂/mnm 実例)、`SPEC.md`、`build_salc_basis` docstring、
  `_reduce_orbit_salcs` コメントの「not mergeable」改訂、`CLAUDE.md` 発散台帳・coupled site。

Excludes:

- スピンスパイラル datum、回転対応評価カーネル、WS 半径外の明示 R 列挙(Phase 2/3)。
- SLCE.jl 側の方針変更。
- 異方チャネル(`Lf > 0`)で列が比例しない場合の結び(検出・報告のみ)。
- SCEMonteCarlo のコード変更(`MultipoleTerm` 契約不変。両軌道に係数が入るので無変更で流れる)。

## Invariants

- Spin layout stays `3 × n_atoms` with unit columns; `magmoms` separate.
- Real tesseral `Zₗₘ`; the `(4π)^(N/2)` design scale applied exactly once.
- `SALCKey` column addressing and ordering; existing TOML models reload and
  predict identically(v2–v6 は `split = :legacy`(群内)/ `:free` で読める)。
- Torque sign `τ = −e × ∂E/∂e`; co-fit whitening.
- Minimum-image / Wigner–Seitz resolvability(no alias folding): `> L/2` の相互作用は
  列挙も折り畳みもしない(不変)。結ぶのは `MinimumImage` が既に保持している WS **境界**タイの
  像のうち、空間群が関係づけないため別軌道に置かれたもの**だけ**(軌道横断)。軌道内の縮約
  `_reduce_orbit_salcs` とその MnTe ゲートは不変。→ CLAUDE.md「Always confirm」の
  resolvability rule に該当。ユーザーは 2026-09-11 の議論(STUDY-2026-09-11-cross-orbit-ties.md)で
  「1 相互作用として推定し等分で読み出す」を決定済み。
- Threaded ≡ serial bitwise gates keep passing.
- エイリアスの無い基底(bcc Fe、FeRh、MnTe、全 pin フィクスチャ)では設計行列・係数・
  読み出しがバイト同一。
- 訓練セル上の `predict_energy` / `predict_torque` は、結ぶ前の min-norm OLS 解と一致
  (同一関数の再配分だから)。

## Completion criteria

- [ ] `make test-all` passes(4 threads)。`make test-pin` unchanged。
- [ ] 新ゲート `test/unit/test_aliases.jl`(オラクルは実装非依存):P1 3 原子 a/2 タイ、
      正方 P セル Δz = c/2、合成回収(倍化セルで J₊ ≠ J₋ を別々に回収;1×1×1 に畳むと
      各軌道 (J₊+J₋)/2)、変異(`alias_rtol = nothing` で群が消え OLS 警告が戻る)、
      閾値帯ピン、persist v7 往復と v6 読み込み。
- [ ] `SPEC.md` / `docs/src/api.md` / `theory/resolvability.md` 更新、`make docs` 通過。
- [ ] Tier 2 レビュー実施。
- [ ] `CLAUDE.md` 発散台帳行(SLCE = 凍結、本パッケージ = 結んで等分)。
- [ ] 既存フィクスチャ(bcc Fe、FeRh、MnTe、pin)で `alias_groups` が空で数値バイト同一。
- [ ] 実系 Nd₂Fe₁₄B(CI 外): l02 = 10 組・結んだ X_T rank 169/169(旧 169/179)、l06 = 30、l044_c4.0 = 20 を
      CHANGELOG に記録。

## References

- 議論記録 `STUDY-2026-09-11-cross-orbit-ties.md`(`~/Packages/` 直下、引き継ぎ文書の内容を含む)
- 再現: Nd₂Fe₁₄B 1×1×1 の `scefit/SUMMARY.md` と `pair_aliases.jl`(データ置き場、リポジトリ外)
- SLCE.jl `docs/src/theory/resolvability.md` "When symmetry does not fuse the tie"、
  `src/basis/resolvability.jl` `unresolvable_columns`、CHANGELOG の face (b) 節
- memory `sce-cross-orbit-ties`, `slce-zero-salc-defect`, `scefitting-anisotropic-defects`
