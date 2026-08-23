# Requirements: `[moment]` section in the TOML setup file

Status: landed (2026-08-24) — Q1–Q4 resolved as recommended; commit hash in tasklist.md

## Goal

モーメント予測（pointed 展開）の打ち切り仕様 `MomentSpec` を、エネルギー基底と同じ
`input.toml` に `[moment]` セクションとして書けるようにし、
`MomentBasis("input.toml")` の 1 行で pointed 基底が建つようにする。
ユーザーの `fit.jl` から `MOMENT_*` 定数群と `MomentSpec(; ...)` の手書きを消すのが目的。

## Background

- 現状 `input.toml` は `[structure]` / `[interaction]` / `[symmetry]` のみ
  （`src/io/input.jl` のファイル docstring がスキーマ）。`SCEBasis(path)` は
  `read_setup(path)` 経由で建つが、`MomentSpec` は Julia 側で
  `MomentSpec(; lmax_env, sampled, lmax_mark, nbody, cutoff_pair, cutoff_star,
  lsum, isotropy)` と書くしかない。
- 実運用（bcc Fe 3×3×3 の `~/jijs/scefit/for_package/fe_bcc/3x3x3cubic/scefit/*/fit.jl`、
  FeRh の `ferh/fit/*/fit.jl`）では、エネルギー側の打ち切りは `input.toml` に、
  モーメント側は `fit.jl` 冒頭の `const MOMENT_LMAX_MARK = 2` 等 6 個の定数に分かれて
  いて、variant を作るたびに 2 ファイルを同期して編集している（`m_p5.0_s4.1` などは
  `input.toml` は同一で `fit.jl` の定数だけが違う）。「セットアップは 1 ファイル」
  という `input.toml` の設計意図（basis/data 分離）からも、打ち切り仕様は TOML 側に
  あるべきもの。
- `MomentSpec` の検証規則（`sampled` 必須・`lmax_env > 0` なら sampled・対称な
  cutoff 行列・`nbody ∈ 1:3`）は既に内側コンストラクタ相当の keyword コンストラクタに
  集約されている（外側 keyword コンストラクタ；位置引数の既定コンストラクタは無検証）
  ので、TOML リーダーは **値を集めて同じコンストラクタに渡すだけ** にできる
  （検証の二重実装はしない）。

## Scope

Includes:

- `src/io/input.jl`: 任意セクション `[moment]` の読み取り（`_moment_from_input`）。
  `read_setup` の返り値に `moment::Union{Nothing,MomentSpec}` を追加。
- `MomentBasis(path::AbstractString; backend, tol, tie_tol)` コンストラクタ
  （`SCEBasis(path)` と同じ上書き規則）。
- 種ラベル・種ペア表の sugar（`"*"` fallback、`"Fe-*"` 等）は `[interaction]` と
  同じ `src/sce/truncation.jl` のリゾルバを再利用。
- `[moment]` の未知キーは **エラー**（タイポ保護。新設セクションなので後方互換の
  制約がない）。
- テスト（`test/unit/test_input.jl` 拡張）、docs（`docs/src/guide/io.md`、
  `docs/src/guide/moment.md`；`api.md` は自動）、`SPEC.md`、`CLAUDE.md` ledger、`CHANGELOG.md`。
- ユーザー側 `fit.jl` の簡略化（`~/jijs/.../fe_bcc/3x3x3cubic/scefit/` は
  ユーザー依頼の書き込み対象）は本 spec の完了後に別途行う（パッケージ外）。

Excludes:

- `MomentDataset` の `gate_eps` / `coverage_floor` を TOML に入れること
  → design.md の open item Q2 で扱う（推奨は「入れる、ただし別サブテーブル」）。
- `MomentModel` の TOML 永続化（現状 `persist.jl` に無い）。別 spec。
- `[interaction]` リーダーの未知キー拒否（既存挙動の変更、別件）。
- `MomentSpec` 本体・`MomentBasis(crystal, spec)` の挙動変更。
- 上流 SLCE.jl への逆移植（SLCE は `soc` 綴りで同名セクションが無い。必要なら
  別途 ledger 行を起こす）。

## Invariants

- `MomentSpec` の検証規則は変えない：TOML から建てた `MomentSpec` は同じ値を
  keyword コンストラクタに渡したものと **フィールド単位で等しい**（`==`）。
- `MomentBasis(path)` は `MomentBasis(read_setup(path).crystal, spec; backend, tol,
  tie_tol)` と **SALC キー列・`marked_atoms`・`records` が同一**（同じコードパスを
  呼ぶだけで、新しい基底構築経路は作らない）。
- `[structure]` / `[interaction]` / `[symmetry]` の既存スキーマと既存ファイルの
  読み取り結果は不変（`[moment]` が無いファイルは `moment === nothing` で、
  `read_setup` の他フィールドはバイト同一）。
- `isotropy` の極性は `MomentSpec` のまま（`true` = `L_S = 0` のみ）。upstream の
  `soc` 綴りは受け付けない（`soc` キーは「`isotropy` を使え」と名指しで拒否）。
- モーメント基底は常に `MinimumImage`（`[interaction].images` は影響しない）。
  同距離帯 `tie_tol` は `[interaction].tie_tol` の値を **そのまま共有**する
  （エネルギー側と別の帯を持たせない：片側だけ殻が割れる事故を防ぐ）。
- Spin layout `3 × n_atoms`、`SALCKey`、`(4π)^(N/2)`、torque 符号、mode rule、
  gate は一切触らない（本 spec は I/O 層のみ）。

## Completion criteria

- [ ] `make test-all` が通る（4 threads）。`make test-pin` 不変。
- [ ] `test/unit/test_input.jl` に `[moment]` のゲートが入る：
  (a) 手書き keyword `MomentSpec` とのフィールド等価（オラクル = 手書き値）、
  (b) `MomentBasis(path)` ≡ `MomentBasis(crystal, spec)` のキー列一致、
  (c) 省略時既定値（`nbody 3` / `lmax_mark 2` / `cutoff_star = cutoff_pair` /
  `lsum` 無し / `isotropy true` / `marked` 全種）が docstring に書いた値と一致、
  (d) 拒否：必須キー欠落、未知キー、`soc` キー、未知ラベル、無順序ペアキーの重複
  （`"Fe-Rh"` と `"Rh-Fe"`；リゾルバは正規化して対称行列しか作らないので
  `MomentSpec` の対称性検査には到達しない）、`sampled` 不整合（`MomentSpec` 由来の
  メッセージがそのまま出る）、`[moment]` 無しでの `MomentBasis(path)`。
- [ ] `make docs` が strict で通り、`io.md` にスキーマ、`moment.md` に TOML 経路が載る
  （`api.md` は bare binding `MomentBasis` の `@docs` が新メソッドの docstring を自動で
  拾う — `SCEBasis(path)` と同じ — ので編集不要；明示シグネチャを足すと重複警告で
  strict ビルドが落ちる）。
- [ ] `SPEC.md` の I/O 節（L216 は既に `tie_tol` が抜けて古い — `tie_tol` と `moment` の
  両方を足す）と `CHANGELOG.md` `[Unreleased]`、`docs/specs/README.md` の行を更新。
- [ ] bcc Fe `l024_c4.1` の `input.toml` に `[moment]` を足して
  `MomentBasis("input.toml")` が既存 `fit.jl` と同じ 83 SALC を出す（受け入れ確認、
  テストではない）。

## References

- `src/io/input.jl`（現スキーマ）、`src/basis/momentbasis.jl` L63–152
  （`MomentSpec`）、`src/sce/truncation.jl`（`_resolve_species_table` /
  `_resolve_pair_table`）。
- `docs/src/guide/moment.md` §Spec and basis（`sampled` / `marked` の設計意図）。
- 実運用の例：`~/jijs/scefit/for_package/fe_bcc/3x3x3cubic/scefit/l024_c4.1/{input.toml,fit.jl}`。
