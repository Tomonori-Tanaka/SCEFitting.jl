# Tasklist: pointed `lsum` の体数別化

Status: landed (2026-08-25)

This file holds coarse-grained, commit-sized milestones. Day-to-day tracking
goes through `TaskCreate` in-session.

## Milestones

### M1 — `MomentSpec.lsum` を体数別に

- [x] フィールドを `Vector{Int}` に、キーワード型を広げて `_resolve_lsum(x, nbody)` へ委譲
- [x] `_moment_labels` の遮蔽を `spec.lsum[N]` に;`N` の範囲を確認
- [x] 空セクタ警告が `lsum[$body]` とその値を名指し
- [x] docstring(添字規約と `cutoff_star` との非対称の理由を含む)

### M2 — テスト

- [x] 手計算ラベル数オラクル(導出をコメントに)
- [x] 合成ゲート(体数別 ≡ 単一値、両方向・非空主張つき)
- [x] スカラーのビット等号ゲート
- [x] 検証ゲート(範囲外 / 重複 / 負値 / 位置指定ベクトル)

### M3 — TOML `[moment].lsum`

- [x] `_moment_lsum_from_input` を body-keyed 受理に
- [x] スキーマ docstring に表の例
- [x] `test_input.jl` の「Dict は拒否」テストを反転、受理ゲートを追加

### M4 — ドキュメント

- [x] `SPEC.md` の `MomentSpec` 記述
- [x] `docs/src/guide/moment.md`(`Σl` 床と体数別キャップの関係)
- [x] `docs/src/guide/io.md`(公開 TOML スキーマページ、部分指定可の理由つき)
- [x] `CLAUDE.md` の sugar 解決 coupled-site 行に `basis/momentbasis.jl`
- [x] `CHANGELOG.md` `[Unreleased]` に **breaking** 見出し;コミット本文に
      `BREAKING CHANGE:`(`MomentSpec.lsum` の読み出し側が壊れる)

### M5 — SLCE.jl 移植

- [x] 同一の struct / 解決 / 読み出し / 警告 / docstring 変更(TOML なし)
- [x] 同一のテスト群を上流に
- [x] 上流スイート緑

## Exit checklist

Run through every item once implementation lands. ~~Strike through~~ items
that do not apply.

- [x] `make test-all` passes (4 threads).
- [x] `make test-pin` passes, or pins recaptured with the reason in
      `test/pin/PIN.md`.
- [x] `make test-parity` passes (needs a sibling SLCE.jl; no CI job). This is the
      only real-data gate on the scalar spelling's bit-identity — the parity cases
      spell `lsum = 2` / `lsum = 4`.
- [x] `make docs` builds (strict).
- [x] If results changed: regression or validation test added, oracle
      independent of the implementation.
- [x] If public API changed: `SPEC.md` and `docs/src/api.md` updated.
- [x] ~~If a hot path was touched: before / after recorded in `bench/BENCH_LOG.md`.~~
      (遮蔽述語 1 行、列数は変わらない)
- [ ] Tier 2 review panel run (numerical / maintainability / performance /
      API axes) and findings resolved.
- [x] ~~If module names or Makefile targets changed: `.claude/agents/` swept.~~
- [x] ~~If this diverges from SLCE.jl: divergence ledger row in `CLAUDE.md`.~~
      (発散しない — 同時移植。既存の pointed 行に添字規約と parity 検出器を追記)
- [x] `CHANGELOG.md` `[Unreleased]` updated.
- [x] `Status:` line in this file and the table in `docs/specs/README.md`
      updated in sync.
- [ ] Implementation commit hash appended below.

## Commits

<!-- filled in as milestones land -->
