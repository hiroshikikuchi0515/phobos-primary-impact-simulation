# GitHubへ登録する手順（v2）

1. `Phobos_GitHub_additions_v2.zip` を**Windows側でダウンロードし**、「すべて展開」します。
2. GitHubの非公開リポジトリ `phobos-primary-impact-simulation` の Code → Add file → Upload files を開きます。
3. このZIP内の **repo_additions_v2フォルダの中身**をドラッグ＆ドロップしてアップロードします。ルートに `data/`、`docs/`、`scripts/`、`workflow/` が並ぶようにします。**ZIPそのものをアップロードしないでください。**
4. `Commit changes`（例：`Add audited 36-case provenance and analysis scripts`）で反映します。
5. GitHub上で `docs/FINAL36_CODE_AUDIT_JP_v2.md` を開き、内容を確認します。依然として **Private** を推奨します。

既存の `README.md`、`src/`、`CITATION.cff`、`.gitignore`、`data/derived/model_288_reference.csv` は**上書きしません**。

## WSLの元環境で検算したい場合（元データの再計算は行わない）

ZIPはWindowsで展開してGitHubへアップロードするためのもので、WSLのプロジェクト直下にコピーする必要はありません。

公開後、Gitでcheckoutしたリポジトリで、36条件のCSVを次のように検証できます（pandas使用）。

```bash
python3 scripts/audit_final36_aggregates_v2.py
```

この検算は36条件の**配布済み集計CSVの整合性**を確認するだけです。火星重力で軌道を新たに積分したり、全Figureを再計算したりはしません。

公開DOIを付けるのは、形状・観測カタログのライセンス、必要ファイル、README、リリースタグを点検した後で構いません。
