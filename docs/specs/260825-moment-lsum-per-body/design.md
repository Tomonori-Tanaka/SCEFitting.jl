# Design: pointed `lsum` の体数別化

Status: landed (2026-08-25)

## Summary

`MomentSpec.lsum` を `Int` から `Vector{Int}`(長さ `nbody`、添字 = 体数 `N`)に
変える。キーワードの解決はエネルギー側の `_resolve_lsum(x, nbody)` を**そのまま
再利用**する — スカラーは `fill(x, nbody)`、`nothing` は `fill(LSUM_UNCAPPED, nbody)`、
body-keyed(`Dict` / `Pair` ベクトル)は指定された体数だけ上書きして残りは無制限。
範囲・重複・負値のエラー文言もエネルギー側と 1 文字同じになる。

読み出しは 1 箇所しかない。`_moment_labels(spec, N)` の

```julia
(iseven(t) && t <= spec.lsum) || continue   # →  t <= spec.lsum[N]
```

この関数はすでに `N` を引数に取っているので、シグネチャも呼び出し側も変わらない。

**却下した代替**: `lsum` 用に独自の解決器を書く案。`_resolve_lsum` は既に必要な
3 形式とすべての範囲チェックを持っており、独自版は文言が割れるだけで得が無い。
**却下した代替**: 位置指定ベクトル `[0, 4, 4]` の受理。`cutoff_star` は
`N − 2` オフセットのため位置指定に意味があったが、`lsum` は添字 = 体数なので
body-keyed が位置指定の上位互換であり、2 綴りを持つ理由が無い。`_resolve_lsum` の
catch-all がすでに body-keyed を名指しで案内して拒否する。

## Module layout

| Target | Change |
|---|---|
| `src/basis/momentbasis.jl` | 構造体フィールド `lsum::Vector{Int}`;キーワード型を広げ `_resolve_lsum` を呼ぶ;`_moment_labels` の遮蔽を `spec.lsum[N]` に;空セクタ警告が `lsum[$body]` とその値を名指し;docstring |
| `src/io/input.jl` | `_moment_lsum_from_input` を body-keyed 受理に(現在は `AbstractDict` を名指しで拒否);`[moment]` スキーマ docstring に表の例 |
| `SPEC.md` / `docs/src/guide/moment.md` / `docs/src/guide/io.md` | `lsum` の体数別綴りと `Σl` 床の関係;公開 TOML スキーマページ |
| `CLAUDE.md` | sugar 解決の coupled-site 行に `basis/momentbasis.jl` を追加(`_resolve_lsum` に 2 人目の consumer) |
| `test/unit/test_momentbasis.jl` | 体数別ゲート(手計算オラクル + 合成)・スカラー等号・検証ゲート |
| `test/unit/test_input.jl` | body-keyed TOML 受理;既存の「Dict は拒否」テストを反転 |
| SLCE.jl `src/basis/momentbasis.jl` + tests | 同一変更(TOML なし) |

## API

```julia
MomentSpec(; …, lsum::Union{Nothing,Integer,AbstractDict,AbstractVector} = nothing, …)

# 綴り(すべて nbody = 4 のとき)
lsum = nothing                     # 全体数 無制限
lsum = 4                           # 全体数 Σl ≤ 4        ← 既存と完全に同じ基底
lsum = [3 => 6, 4 => 4]            # 3 体は 6 まで、4 体は 4 まで、1・2 体は無制限
lsum = Dict(4 => 4)                # 4 体だけ 4、他は無制限
```

TOML(**別々のファイル**。同じ `[moment]` に裸の `lsum` とサブテーブル
`[moment.lsum]` を両方書くのは TOML の重複キー。サブテーブルは `[moment]` の裸キーの
**後**に置く — `cutoff_star` と同じ落とし穴):

```toml
# スカラー版
[moment]
lsum = 4
```

```toml
# 体数別版(キーは裸の整数 = 体数、[interaction].lsum と同じ綴り)
[moment]
cutoff_pair = 4.1
[moment.lsum]
3 = 6
4 = 4
```

`cutoff_star` の表は `3:nbody` を**過不足なく**覆う必要があるが、`lsum` の表は
**部分指定可**で未指定は無制限。半径が欠けると黙って `cutoff_pair`(その体数と
無関係な数)に落ちるのに対し、キャップが欠けるのは「制限なし」という安全で
文書化済みのデフォルトだから。

内部アクセスは **`_label_lsum(spec, N)`**(範囲チェック込み)。生の
`spec.lsum[N]` は残さない。理由は「規約を 1 箇所に閉じ込める」だけではない:
同じ struct に規約の違うベクトルが 2 本並び、`cutoff_star` 側は既に
`_star_cutoff(spec, N)` 経由で読めとコメントで警告している。`spec.lsum[N - 2]` と
書き間違えても**例外は出ず、別の体数のキャップが黙って効く**(N = 3 なら 1 体の
キャップ)。読み出しが 1 箇所しかない今なら実質コストゼロで、2 つの規約が隣り合う
アクセサとして並ぶ。

## Types and conventions

- `LSUM_UNCAPPED` = `typemax(Int)`(`src/basis/salcbasis.jl`)。無制限の哨兵は
  エネルギー側と共有し、pointed 側で別の哨兵を作らない。
- 添字規約: `lsum[N]` は体数 `N`。`BasisSpec.lsum` と同一、`cutoff_star[N−2]` とは
  意図的に異なる。この非対称は両方の docstring に理由付きで書く。
- 物理規約の変更なし。`Σl` の偶数遮蔽(時間反転)と床 `2⌈(N−1)/2⌉` は不変。
- 数値結果の変更なし(既存綴りはビット等号)。新しい表現力が増えるだけ。

## Impact on coupled sites

- [x] `Harmonics` / `AngularMomentum` / SALC projection / normalization tests: 影響なし
- [x] `SALCKey` order / design columns / TOML persistence / `coeftable`: 影響なし
      (pointed モデルは TOML 永続化を持たない — `io/persist.jl` は `BasisSpec` のみ)
- [x] Energy kernel ↔ torque kernel: 影響なし
- [ ] Decor engine ↔ pointed moment basis ↔ `momentfit.jl` doors: `_moment_labels` の
      遮蔽が唯一の読み出し点。`momentfit.jl` は `spec.lsum` を読まない(確認済み)
- [x] Readers' `zero_moment_atol` ↔ dataset doors: 影響なし
- [ ] **BasisSpec sugar resolution ↔ canonical consumers**: `_resolve_lsum` に 2 人目の
      consumer(`basis/momentbasis.jl`)が付く。受理形式やエラー文言を触ると
      `[interaction]` と `[moment]` の両チャネルが同時に動く — CLAUDE.md の当該行に
      `basis/momentbasis.jl` を名指しで足す
- [ ] Upstream divergence ledger (SLCE.jl): 同時移植するので新しい発散行は増えない。
      既存の `[moment]` TOML 行(こちらにしか無い)に body-keyed `lsum` を含める
- [x] `.claude/agents/`: 影響なし
- [ ] `SPEC.md` / `docs/src/api.md` updates: `MomentSpec` の記述を更新

## Test strategy

独立オラクルは **1 だけ**(手計算)。2 は実装の中核ルーチンを共有するので
オラクルではなく **property / invariant ゲート**(体数間の独立性)、3 は等号ゲート。

1. **手計算ラベル数(オラクル = 手で数えた閉形式)**。既存の手計算ラベル列挙
   testset(`test/unit/test_momentbasis.jl` の "(4) labels against a hand
   enumeration")を拡張する。`lmax_mark = 2, lmax_env = [2]`
   のとき、体数 `N` のラベルは「マーク階数 `lm ∈ 0:2` × 非減少環境多重集合
   `e ∈ 1:2^(N−1)`、`Σl = lm + Σe` が偶数かつ ≤ cap」。N = 2 は
   `(0,2) (1,1) (2,2)` の 3 本(Σl = 2, 2, 4)なので cap = 2 で 2 本、cap = 4 で 3 本。
   N = 4 の床は 4。これを `lsum = [2 => 2, 4 => 4]` のような混在 spec で確認する。
   導出はテスト内のコメントに書く。
2. **合成 property ゲート(オラクルではない — 実装を共有する自己整合ゲート)**。体数別 spec
   `[b1 => c1, b2 => c2]` の体数 `b` の列集合は、単一値 spec `lsum = c_b` の体数 `b`
   の列集合と**キー単位で一致**しなければならない。両方向・両体数で確認する
   (`cutoff_star` の体数別化と同じ形のゲート)。片方向だけだと「両方空」で
   受かるので、非空であることも同時に主張する。
3. **スカラーのビット等号**。`lsum = 4` と `lsum = [1 => 4, 2 => 4, 3 => 4]`
   (nbody = 3)が同一のキー列と同一の設計行列を返す。
4. **検証ゲート**。範囲外の体数キー・重複キー・負値・位置指定ベクトルが
   `ArgumentError` で落ちる(文言はエネルギー側と共有)。
5. **TOML ゲート**。`[moment.lsum]` の表が読める;`test_input.jl` の
   「Dict は拒否」テストを「Dict は受理され、体数別に効く」に反転。既存の
   `@test m.lsum == typemax(Int)` は `fill(typemax(Int), nbody)` に更新が必要
   (`Vector{Int}` 化で必ず落ちる)。`lsum = -1` の `@test_throws ArgumentError` は
   文言非依存なので通るが、文言は `"lsum must be ≥ 0"` → `"lsum: must be ≥ 0"` に
   変わる(エネルギー側と共有する意図的な変化)。
6. **既存の pin / parity**。`lsum = 2` / `lsum = 4` のスカラー綴りを使っている
   fixture(`test/pin/fixtures.jl`, `test/pin/pins/B2_FeRh_pointed.toml`,
   `test/parity/runtests.jl`)は**再取得なしで通る**こと自体が等号の証拠。

## Risks and open items

- **`_moment_labels` の `N` 範囲(確認済み・危険なし)**。本番の呼び出しは
  `for b in 1:spec.nbody` の 2 箇所のみ、テストも `N ≤ nbody`。それでも
  `_label_lsum` が `1 <= N <= spec.nbody` を明示チェックし、`BoundsError` ではなく
  「その体数は spec に無い」と言う。
- **未決なし**。API・添字規約・オラクルはすべて上で確定している。
- **後続作業(この spec の外)**: M6 変種 C(448 列・4 体 340 列)を罰則付き推定量で
  再測する。体数別 `lsum` が入ると「3 体無制限 + 4 体 `lsum = 4`」が初めて
  1 本の spec で書けるので、再測の設計はそれを使う。
