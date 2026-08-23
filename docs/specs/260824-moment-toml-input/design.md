# Design: `[moment]` section in the TOML setup file

Status: landed (2026-08-24) — Q1–Q4 resolved as recommended; commit hash in tasklist.md

## Summary

`input.toml` に任意セクション `[moment]` を足し、`read_setup` がそれを
`MomentSpec` に解決して返す。`MomentBasis(path)` は `SCEBasis(path)` の写し
（ファイルの `[symmetry]`・`[interaction].tie_tol` を既定値に、keyword で上書き可）。
リーダーは **値を集めて既存の `MomentSpec(; ...)` keyword コンストラクタに渡すだけ**
で、検証は一切複製しない。種ラベル表・種ペア表は `[interaction]` と同じリゾルバ
（`_resolve_species_table` / `_resolve_pair_table`）を再利用するので、書式も
エラー文も両セクションで揃う。

ユーザーの `fit.jl` は

```julia
basis = SCEBasis("input.toml")
mb    = MomentBasis("input.toml")          # ← 6 個の MOMENT_* 定数と MomentSpec(...) が消える
mds   = MomentDataset(mb, data; gate_eps = 2e-3)
```

になる。

### スキーマ（`src/io/input.jl` のファイル docstring に追記する本文）

```toml
[moment]                      # optional section: the pointed site-moment basis
nbody       = 3               # optional, default 3 (1, 2, or 3)
lmax_mark   = 2               # optional, default 2: cap on the marked site's own ê factor
lmax_env    = [2]             # per species (index order), or a label table (below)
sampled     = ["Fe"]          # REQUIRED: species the downstream consumer samples —
                              #   labels, "*" (= every species), or per-species booleans
marked      = ["Fe"]          # optional (default every species): whose moments are expanded
cutoff_pair = 4.1             # REQUIRED: mark–environment bond radius (Å) of 2-body
                              #   clusters — scalar (`inf` = whole WS cell) or a pair table
cutoff_star = 4.1             # optional, default = cutoff_pair: the two mark bonds of
                              #   3-body stars (environment–environment edge is free)
lsum        = 4               # optional, default uncapped: total spin rank per label
isotropy    = true            # optional, default true (L_S = 0 only) — NOTE the default
                              #   differs from [interaction].isotropy (false)

# Label-keyed alternatives (same rules as [interaction]):
#     [moment.lmax_env]
#     "*" = 2
#     Rh  = 0
#     [moment.cutoff_pair]          # pair keys, unordered, specificity-resolved
#     "Fe-Fe" = 4.1
#     "*-*"   = 3.0
#     [moment.cutoff_star]
#     "*-*"   = 4.1
```

規則：

- `[moment]` が無い → `read_setup(path).moment === nothing`。他フィールドは不変。
- **未知キーはエラー**（`[moment]: unknown key "lmax_enviroment" (allowed: ...)`）。
  `soc` は upstream 綴りとして名指しで拒否（「`isotropy` を使え。極性は逆」）。
- `sampled` / `marked`: 文字列配列（ラベル；配列内の `"*"` は全種、`sampled = "*"`
  の裸文字列は不可）または `Bool` 配列（種インデックス順）。リーダーは **変換だけ**
  （ラベル → `Vector{Bool}`）し、長さ・整合性の検査は `MomentSpec` に任せる。
  `eltype` が `String` でも `Bool` でもない配列（混在 `["Fe", true]`、整数 `[1, 0]` —
  TOML は混在配列を `Vector{Any}` で通すので明示的に弾く）、重複・未知ラベルはエラー。
  `marked` 省略 = 全種（`MomentSpec` の既定と同じ）。
- `lmax_env`: 整数配列（種インデックス順）またはラベル表（`"*"` fallback あり）
  → `_resolve_species_table(x, nkd, labels, "[moment].lmax_env")`。裸のスカラーは
  `[interaction].lmax` と同じく不可（リゾルバ自体は broadcast を受けるが、リーダーで
  弾いて両セクションの書式を揃える）。
- `cutoff_pair` / `cutoff_star`: 実数（`inf` 可）またはペア表 → 値を
  `_pairtable_from_input(x, ctx)` で `Float64` に揃えてから
  `_resolve_pair_table(x, nkd, labels, "[moment].cutoff_pair")`（文字列値で
  `_check_cutoff_value` が `MethodError` を出すのを防ぐ）。body-keyed 表（`2 = …`）は
  `_is_bodykey` で検出して **名指しで拒否**（pointed 基底の cutoff は body 別でなく
  役割別 pair / star なので、キー名がその役割；`_pair_key_parts` の "not of the form
  A-B" に落とさない）。
- `lsum`: 整数のみ（`[interaction].lsum` の body 表は受け付けない — `MomentSpec.lsum`
  はラベル全体の総ランク 1 個）。
- `lsum` / `nbody` / `lmax_mark`: `isa Integer`、`isotropy`: `isa Bool` を明示検査して
  `ArgumentError`（`Int(Dict)` / `Bool("x")` の `MethodError` に落とさない）。
- `read_setup` は従来どおり `[structure]` と `[interaction]` を必須とする：
  `MomentBasis(path)` はモーメント基底だけが欲しい場合でも `[interaction]` を要る
  （`tie_tol` の供給元）。docstring と io.md に明記。
- 集めた値は `MomentSpec(; lmax_env, sampled, lmax_mark, marked, nbody, cutoff_pair,
  cutoff_star, lsum, isotropy)` に渡す。`sampled` 不整合・非対称・範囲外は
  `MomentSpec` 自身が投げる（メッセージはそのまま、文脈 `[moment]` を前置して再送）。

## Module layout

| Target | Change |
|---|---|
| `src/SCEFitting.jl` | `include("io/input.jl")` を `basis/momentbasis.jl` の **後**に移す（L74–76 の I/O 節コメントも合わせて更新）（`read_setup` の返り値型 `Union{Nothing,MomentSpec}` と `MomentBasis(path)` の定義に型が要る）。`io/input.jl` は `Crystal` / `BasisSpec` / backend 型にしか依存しないので移動は安全。`dftsource.jl` / `embset.jl` / `extxyz.jl` は `read_setup` を参照しない（grep で確認済み）。 |
| `src/io/input.jl` | ファイル docstring にスキーマ追記。`_moment_from_input(d, labels)::MomentSpec`、`_species_list_from_input(x, labels, what)::Vector{Bool}` を追加。`read_setup` の NamedTuple に `moment` を追加。`MomentBasis(path; backend, tol, tie_tol)` を `SCEBasis(path)` の直後に定義。 |
| `src/sce/truncation.jl` | 変更なし（リゾルバ再利用のみ）。 |
| `src/basis/momentbasis.jl` | 変更なし。`MomentSpec` docstring の冒頭に "or from a TOML `[moment]` section, see `read_setup`" の 1 行だけ追加。 |
| `test/unit/test_input.jl` | `[moment]` の testset 群を追加（下記）。 |
| `docs/src/guide/io.md` | `[moment]` スキーマの小節。 |
| `docs/src/guide/moment.md` | §Spec and basis に TOML 経路を追記。 |
| `docs/src/api.md` | **編集不要**：bare binding `MomentBasis` の `@docs`（L333–338）が全メソッドの docstring を自動で描く（`SCEBasis(path)` が L110 の bare `SCEBasis` に載るのと同じ）。明示シグネチャを足すと重複警告 → `warnonly = false` で strict ビルド失敗。 |
| `SPEC.md` | I/O 節 L216 の `read_setup` シグネチャ（既に `tie_tol` が抜けて古い）に `tie_tol` と `moment` を足し、`[moment]` の 1 行。 |
| `CLAUDE.md` | 上流差分 ledger に 1 行；coupled-site 「BasisSpec sugar resolution ↔ canonical consumers」の bullet に `_moment_from_input` を読者として追記（「sugar 形を足したら TOML リーダーも更新」の対象に含める）。 |
| `docs/specs/README.md` | 本 spec の行。 |
| `CHANGELOG.md` | `[Unreleased]` → Added。 |

## API

```julia
# read_setup: one more field; callers destructuring by name are unaffected
function read_setup(path::AbstractString)::@NamedTuple{crystal::Crystal,
                                                       spec::BasisSpec,
                                                       backend::AbstractSymmetryBackend,
                                                       tol::Float64,
                                                       images::AbstractImageSelection,
                                                       tie_tol::Float64,
                                                       moment::Union{Nothing,MomentSpec}}

"""
    MomentBasis(path::AbstractString; backend = nothing, tol = nothing, tie_tol = nothing)

Build a pointed [`MomentBasis`](@ref) from a TOML input file: the crystal from
`[structure]`, the truncation from `[moment]` (required here — `ArgumentError` names
the missing section), symmetry from `[symmetry]`, and the same-distance band from
`[interaction].tie_tol`, each overridable by the keyword. The moment basis is always
minimum-image; `[interaction].images` does not apply to it.
"""
function MomentBasis(path::AbstractString;
                     backend::Union{Nothing,AbstractSymmetryBackend} = nothing,
                     tol::Union{Nothing,Real} = nothing,
                     tie_tol::Union{Nothing,Real} = nothing)::MomentBasis
```

内部ヘルパは `_`-prefix、export 変更なし（`MomentBasis` は既に export 済み、
`read_setup` も export 済み）。

## Types and conventions

- 物理・数値規約への影響ゼロ（I/O 層のみ）。`MomentSpec` の検証を TOML 側で
  複製しない＝**検証の所在は 1 箇所**（coupled-site を増やさない）。
- `isotropy` 既定値の非対称（`[interaction]` false / `[moment]` true）は
  `MomentSpec` の既定値を尊重する。理由：モーメント側の既定はパッケージ API で
  既に `true` と公開済みで、TOML だけ変えると同じ名前が 2 つの既定を持つ。
  docstring・guide の両方に明記する（open item Q3）。
- SLCE.jl との差分：upstream に `[moment]` TOML は無い（upstream は `soc` 綴り）。
  ledger 行は「SCEFitting 独自の I/O 拡張、逆移植するなら `soc` 極性に注意」。

## Impact on coupled sites

- [ ] `Harmonics` / `AngularMomentum` / SALC projection / normalization tests: 影響なし。
- [ ] `SALCKey` order / design columns / TOML persistence / `coeftable`: 影響なし
      （`model.toml` の永続化フォーマットは触らない）。
- [ ] Energy kernel ↔ torque kernel: 影響なし。
- [x] Decor engine ↔ pointed moment basis ↔ `momentfit.jl` doors: `MomentBasis(path)` は
      既存 `MomentBasis(crystal, spec)` を呼ぶだけ。ゲート (b) でキー列同一を確認。
- [ ] Readers' `zero_moment_atol` ↔ dataset doors: 影響なし。
- [x] Upstream divergence ledger (SLCE.jl): 1 行追加（上記）。
- [ ] `.claude/agents/` references: モジュール名・Makefile 変更なし → 不要。
- [x] `SPEC.md` / `docs/src/api.md` updates: 上記表。

## Test strategy

すべて `test/unit/test_input.jl` に追加。オラクルは **手書きの keyword 呼び出しと
docstring に書いた既定値**（TOML パーサと同じコードを通らない）。

1. **フィールド等価**（オラクル = 手書き `MomentSpec(; ...)`）: 2 種系
   `["Fe","Rh"]` で全キーを与えたファイルを読み、`read_setup(p).moment` の 9 フィールドが
   `MomentSpec(; lmax_env=[2,0], sampled=[true,false], lmax_mark=2, marked=[true,true],
   nbody=3, cutoff_pair=[4.1 3.0; 3.0 3.0], cutoff_star=…, lsum=4, isotropy=true)` と
   `==`。ラベル表形（`"*"` fallback、`"Fe-*"`）と添字形の両方で同じ結果。
2. **既定値**（オラクル = docstring の値）: `lmax_env`・`sampled`・`cutoff_pair` だけの
   最小ファイルで `nbody == 3`、`lmax_mark == 2`、`cutoff_star == cutoff_pair`、
   `lsum == typemax(Int)`、`isotropy == true`、`marked == fill(true, nkd)`。
3. **基底同一**（オラクル = 既存 Julia 経路）: `_INPUT_FULL` 相当の 2 原子鎖 +
   `[moment]` で `MomentBasis(p; backend = NoSymmetry())` と
   `MomentBasis(read_setup(p).crystal, read_setup(p).moment)` の
   `salc_basis.keys`・`marked_atoms`・`n_salcs` が一致。keyword 上書き（`tie_tol`）が
   `mb.tie_tol` に反映される。
4. **拒否**: `sampled` 欠落 / `cutoff_pair` 欠落 / `lmax_env` 欠落、未知キー
   （`lmax_enviroment`）、`soc = false`（メッセージに `isotropy` を含む）、未知ラベル
   `"Co"`、`sampled = ["Fe", true]`（混在）、`sampled = [1, 0]`（整数配列）、`sampled = [true]` を 2 種系で（長さ — `MomentSpec` の
   メッセージ）、`[moment.cutoff_pair]` に `"Fe-Rh"` と `"Rh-Fe"`（重複ペアキー）、
   `lmax_env = 2`（裸スカラー）、
   `lmax_env = [2, 2]` + `sampled = ["Fe"]`（`MomentSpec` の M3-1 拒否がそのまま出る）、
   `[moment]` 無しで `MomentBasis(p)`（メッセージに `[moment]`）、
   `lsum` に body 表、`cutoff_pair` に body 表。
5. **既存ファイル不変**: `_INPUT_FULL` / `_INPUT_MINIMAL` で `moment === nothing` かつ
   他フィールドは従来のテストがそのまま通る。

受け入れ確認（テスト外）: bcc Fe `l024_c4.1/input.toml` に `[moment]` を足して
`MomentBasis("input.toml")` が 83 SALC、`Im-3m` で建つこと。

## Risks and open items

- **Q1 `sampled` を必須のままにするか。** 単一種系ではいつも `sampled = ["Fe"]` と
  書くことになる。既定「`lmax_env > 0` の種を sampled とみなす」は `MomentSpec` が
  わざわざ拒否している組合せを黙って直すことになるので **必須のまま**を推奨
  （M3-1 決定 A の趣旨：消費側が何をサンプルするかは宣言事項）。
- **Q2 `gate_eps`（`MomentDataset`）を TOML に入れるか。** 入れない案を推奨：
  `input.toml` は「セットアップ＝結晶+打ち切り+対称性」で、`gate_eps` は
  **データの品質に対する宣言**（docstring: "the tolerance is a statement about the
  data, the caller must state it"）。入れる場合は `[moment.dataset] gate_eps = 2e-3`
  を `read_setup(...).moment_dataset::NamedTuple` で返し
  `MomentDataset(mb, data; setup.moment_dataset...)` と splat する形が最小
  （`MomentDataset` 自体は触らない）。ユーザー判断。
- **Q3 `isotropy` 既定値の非対称。** 上記のとおり `MomentSpec` 既定 `true` を尊重。
  代替 (b)「`[moment].isotropy` を必須にする」は安全だが冗長、(c)「省略時は
  `[interaction].isotropy` に従う」は SOC エネルギー fit で意図せず異方 pointed 基底
  （列数が数倍）になるので却下。
- `include` 順の移動は `io/input.jl` の依存が `Crystal`/`BasisSpec`/backend 型と
  `_SAME_DIST_RTOL` だけであること（grep 済み：`read_setup` と input.jl のヘルパは
  src 内で他から参照されない；`momentfit.jl` 末尾が `MomentBasis` 用 `save` を io 節の
  後に定義している前例あり）に依る；`make test-all` で担保。
- **Q4 レビュー tier。** I/O 層のみなので Tier 2 パネルを api-reviewer +
  maintainability-reviewer の 2 軸（numerical / performance は対象外）に縮約する案を
  推奨。ユーザー判断。
- 数値結果への影響なし（新しい構築経路を作らない）。
