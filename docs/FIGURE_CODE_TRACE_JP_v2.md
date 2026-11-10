# 図表生成コードの対応（監査用、現時点で判明している範囲）

| 対象 | 確認できた入力・処理 | 状況 |
|---|---|---|
| 36条件の元計算 | `workflow/as_used/.../tools/run_lwss_campaign_v1.py` → C++ v5 | **条件・seedを監査済**、外部OBJが必要 |
| L/T・M/Aの36条件集計 | `tools/summarize_lwss_campaign_v1.py` → 36行CSV | **集計CSVと整合**、26列イベントが必要 |
| 全球4象限・Stickney25 | `scripts/rebuild_exact_four_sector_from_raw_v2.py` ／ `rebuild_masked_model_only_v2.py` | **スクリプトを保存**、再実行は生イベント・観測座標が必要 |
| 衝突分布の円筒図 | `tools/make_v004_pg1_impact_maps_v1.py` | v∞=7.95298 km/s、a=9376/20000 km、15°ビンの図に対応。**最終掲載画像との一致は未検証** |
| それ以外のMain Fig. 1–7 / Supplementary Fig. S1–S12 | 追加の作図・統計処理が必要 | **今回の10ファイルでは生成元を特定し切れない** |

ラベルやファイル名だけで図の出所を断定せず、掲載したPNG/PDFと出力結果のチェックサム・元スクリプトを逐次照合してください。
