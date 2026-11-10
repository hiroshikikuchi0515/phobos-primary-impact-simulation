# 図表コード監査 v3 — Icarus primary impacts on Phobos

2026-11-10受領の `phobos_figure_code_review_v3.tar.gz` に含まれたPython 10本を静的確認し、構文解析10/10に成功しました。**入力群を揃えた完全な実行再現性の検証ではありません。** 原稿は本文v25、Supplementary v26を基準にしています。

## 対応が明確なもの

- `plot_bottom_legends_v3.py`: 36条件のL/T比、D≥1 km比較、直径別観測L/T比、全球4象限の旧版図、Stickney25観測分率という**5系統の旧作図**を作成します。現行の本文Figure 1–7をそのまま再生成するものではありません。
- `build_available_results_v2.py`: Ernst v004の頂点と三角面から面積・Stickney25を計算し、観測の直径階級周辺分布と全球4象限モデル表を生成します。独立した形状OBJ、`sector_summary_v1.csv`、観測クレーター径別CSVが必要で、面積とN=1235のassertがあります。
- `plot_large_crater_global_v2.py`: 旧版の大型クレーターL/T比を作図します。
- `recompute_rotation_robustness_v30.py`: `rotation_gof_full_v25.csv`（144セル×19姿勢＝2736行）から回転頑健性の集計CSVと2×2図を出力します。最終Supplementary Figure S11に**対応する作図候補**ですが、原図とのピクセル・ラベル完全一致までは未検証です。
- `analyze_benchmark_v26.py`、`analyze_convergence_v26.py`、`plot_validation_v23.py`、`summarize_campaign_v25.py`: 旧検証バージョンの計算・作図です。最終10 seed図(S12)や改訂版S10を再生成できるとは主張しません。
- `build_pg1_figure_index_v2.py`、`index_existing_figures_v1.py`: 過去の画像を探索する補助ツールです。科学的結果を新たに生成しません。

## 本文Figure 1–7

| 図 | 現時点のコード確認 |
|---|---|
| Fig.1 | 模式図。今回の10本に生成元はありません |
| Fig.2 | 既登録 `workflow/as_used/.../make_v004_pg1_impact_maps_v1.py` は関連しますが、現行図との一致検証は未完了 |
| Fig.3 | 旧L/T図コードはあり。ただし現行4パネルのL/T・M/A図とは別版 |
| Fig.4 | 36条件×2カタログ×2マスクの最終p値ヒートマップ作図コードは今回未提供 |
| Fig.5 | 最終4象限モデル・観測図の作図コードは今回未提供 |
| Fig.6 | 直径別(G−3)/N図の最終作図コードは今回未提供 |
| Fig.7 | 旧Stickney図コードはあり。現行の36条件曲線版とは一致未確認 |

## Supplementary Figure S1–S12

- S1–S3: 旧版に関連図がありますが、最終改訂図そのものとは未照合です。
- S4–S6: 関連する本番マップコードは既登録。ただし横縦2:1や凡例を改訂した最終v26との一致は未確認です。
- S7/S8: 最終図の生成元は今回の10本にありません。
- S9: `plot_validation_v23.py` は**旧回転診断図**で、現在のS9ではありません。
- S10: `analyze_convergence_v26.py` は**旧注入半径図**で、現在のS10ではありません。
- S11: `recompute_rotation_robustness_v30.py` は該当する144セル図の候補です。CSV未同梱です。
- S12: `analyze_benchmark_v26.py` は初期ベンチマークであり、最終10 seed図の作図元ではありません。

## 公開前の注意点

1. `plot_bottom_legends_v3.py` と他の旧作図は元の `input/`, `data/`, `figures/` 階層を前提にしています。この追加ZIP内の置き場では**そのまま実行できません**。再現コードとしてではなく出典管理用のアーカイブです。
2. 外部形状モデルや観測カタログを公開するには配布条件を確認してください。新規に算出した集計値と第三者作成の観測値は区別します。
3. 最終図を再生成できると主張するには、図ごとの入力CSV、作図スクリプト、実行方法、生成結果の確認が必要です。
4. `rotation_gof_full_v25.csv` が手元にあれば `python3 analysis/validation/recompute_rotation_robustness_v30.py --input /path/to/rotation_gof_full_v25.csv --outdir output/rotation_v30` で出力可能です。

## 次の候補（新しい粒子計算は不要）

- 本文Fig.2とS4–S6の元になる15°格子CSV、マップ作図時の出力一式。
- 本文Fig.3–7とS7–S12の**最終版**作図スクリプト（ChatGPT側で生成したものも含む）。
- 座標回転全比較CSV `rotation_gof_full_v25.csv` とその観測カタログの帰属情報。
