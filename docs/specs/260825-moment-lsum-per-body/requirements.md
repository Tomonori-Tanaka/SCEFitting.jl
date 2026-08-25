# Requirements: pointed `lsum` の体数別化

Status: landed (2026-08-25)

## Goal

`MomentSpec.lsum` を単一の `Int` から**体数別**(`Vector{Int}`、添字 = 体数 `N`)へ
拡張し、`[interaction].lsum` と同じ綴り(`lsum = 4` / `lsum = [1 => 0, 2 => 4]`)で
書けるようにする。スカラーは全体数へブロードキャストし、既存の綴りは**ビット等号**を
保つ。

## Background

**構造的な理由。** pointed の `N` 体セクタは `Σl = 2⌈(N−1)/2⌉` から始まる
(環境スロットは各々 `l ≥ 1`、時間反転で偶 `Σl` のみ)。したがって床は
N=1:0 / N=2:2 / N=3:2 / N=4:4。単一の `lsum` は高体数を**構造的に飢えさせる**:
`lsum = 4` のとき N=2・N=3 は偶レベルを 2 段(Σl = 2, 4)使えるのに N=4 は 1 段
(Σl = 4)しか使えない。実測ラベル数(`lmax_mark = 2, lmax_env = [2], lsum = 4`)は
N=1→2, N=2→3, N=3→4, **N=4→2**。

**実測上の理由。** 直前の spec の M6(4 体の効きの測定、
`for_package/fe_bcc/.../INVESTIGATION-moment-residuals.md` §4f)で
「3 体は無制限のまま 4 体だけ `lsum = 4`」が**表現できなかった**。`lsum = 4` を
置くと 3 体側の `Σl = 6` 列 17 本も道連れになり(93 → 76)、4 体の寄与と 3 体の
剥離を切り分けるために追加の対照 run(F)を回す羽目になった。測定設計の解像度が
API 側で潰れている。

**コスト。** 読み出しは `_moment_labels` の `t <= spec.lsum` **1 箇所**だけで、
その関数はすでに `N` を引数に取っている。解決器 `_resolve_lsum(x, nbody)` は
`src/sce/truncation.jl` にエネルギー側用として**既にあり**、`Nothing` /
`Integer` / body-keyed の 3 形式と範囲・重複・負値チェックを実装済み。

## Scope

Includes:

- `MomentSpec.lsum::Int` → `Vector{Int}`(長さ `nbody`、添字 = 体数)
- キーワード `lsum` が body-keyed table / Dict / Pair ベクトルを受理(`_resolve_lsum` 再利用)
- `_moment_labels` の遮蔽を `spec.lsum[N]` へ
- 空セクタ警告が責任を持つ半径・キャップを名指しする(`lsum[$body]` とその値)
- TOML `[moment].lsum` の body-keyed 受理(現在は Dict を名指しで拒否している)
- docstring / 入力スキーマ / `SPEC.md` / guide の更新
- SLCE.jl への同時移植(TOML を除く。上流に `[moment]` セクションは無い)

Excludes:

- `lmax_mark` / `lmax_env` の体数別化(必要になった実測がまだ無い)
- `cutoff_star` の再設計(直前の spec で着地済み)
- 位置指定ベクトル `lsum = [0, 4, 4]` の受理(下の Invariants を参照)

## Invariants

- **これは読み出し側の破壊的変更**。`MomentSpec` は export された公開型で、
  `spec.lsum` を `Int` として読んでいるコード(`spec.lsum == 4` 等)は壊れる。
  ビット等号なのは**構築される基底**であって、フィールドの型ではない。
  `CHANGELOG.md` は breaking 見出し、コミット本文に `BREAKING CHANGE:` を入れる。
- **スカラー・`nothing` の綴りはビット等号**: `lsum = 4` は全体数に 4、`nothing` は
  全体数に `LSUM_UNCAPPED`(= `typemax(Int)`)。遮蔽述語 `t <= cap` は不変なので
  既存の基底・列順・pin は 1 ビットも動かない。
- 添字の規約は **`lsum[N]` = 体数 `N` そのもの**(`BasisSpec.lsum` と同一、
  `salcbasis.jl` の `lsum_by_body[N]` と同一)。`cutoff_star[N − 2]` のオフセットとは
  意図的に**異なる**(`cutoff_star` は N ≥ 3 にしか存在しないが `lsum` は N = 1 から
  効く)。
- `[moment].lsum` は `[interaction].lsum` と**同じ綴り**を受ける。エネルギー側と
  異なる形式を発明しない。
- 位置指定ベクトルは受理せず、body-keyed を名指しで案内して**拒否**する
  (`_resolve_lsum` の catch-all がすでにそうしている)。`cutoff_star` は位置指定を
  受けるが、あちらは添字にオフセットがあり位置と体数が食い違うため別枠。
- `Σl` の床 `2⌈(N−1)/2⌉` と偶 `Σl` の時間反転遮蔽は不変。
- SALCKey の列順・設計行列・トルク符号・分解可能性ゲートは不変。
- 直列 ≡ 並列のビット等号ゲートは通り続ける。
- `(4π)^(N/2)` の設計スケールは 1 回だけ、実タセラル `Zₗₘ`、トルク符号
  `τ = −e × ∂E/∂e`、最小像の分解可能性 — いずれも不変。
- **表の部分指定可否は 2 つの表で異なってよい**。`[moment].cutoff_star` は
  `3:nbody` を過不足なく覆うことが必須(半径が欠けると黙って `cutoff_pair` に
  落ちる = 意味の無い数)。`[moment].lsum` は部分指定可(欠けた体数は無制限 =
  安全かつ文書化されたデフォルト)。この非対称は意図であり、両方の文書に理由を書く。

## Completion criteria

- [ ] `make test-all` / `make test-pin`(再取得なし)/ `make test-parity` / `make docs` が緑
- [ ] 体数別キャップの受入ゲートに**実装非依存のオラクル**がある
      (手計算のラベル数 + 合成ゲート: 体数別 spec の体数 `b` 内容 ≡ 単一値 spec の
      体数 `b` 内容)
- [ ] スカラー綴りのビット等号が明示ゲートで守られている
- [ ] TOML の body-keyed 受理ゲート(および現在の「Dict は拒否」テストの反転)
- [ ] `SPEC.md` / `docs/src/` / 入力スキーマ docstring 更新
- [ ] SLCE.jl 側に同一の変更が入り、上流のスイートが緑
- [ ] `CHANGELOG.md` `[Unreleased]` 更新

## References

- 直前の spec: `docs/specs/260824-pointed-nbody-general/`(`cutoff_star` の体数別化、
  同じ形の変更)
- M6 の実測: `~/jijs/scefit/for_package/fe_bcc/3x3x3cubic/scefit/INVESTIGATION-moment-residuals.md` §4f
- 再利用する解決器: `src/sce/truncation.jl` の `_resolve_lsum`
- エネルギー側の綴り: `docs/src/guide/basis.md`(`lsum = [2 => 4]`)
