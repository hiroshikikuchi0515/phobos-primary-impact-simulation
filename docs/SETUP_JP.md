# GitHub登録手順（論文用の別リポジトリ）

1. GitHubの右上「+」から **New repository**、または https://github.com/new を開く。
2. Ownerが自分のアカウントであることを確認し、名前を `phobos-primary-impact-simulation` にする。
3. Descriptionに `Simulation and analysis code for primary impacts on Phobos (Icarus manuscript)` と入力する。
4. 公開する場合は **Public** を選ぶ。権利・内容を確認する前に公開するのが不安な場合は Privateで始め、後でPublicに変更可能（ただし公開前の履歴に注意）。
5. この雛形をアップロードする場合、GitHub上でREADME/ライセンスを先に自動作成する必要はない。
6. コードの権利や依存物の公開可否、README、実行条件、第三者データの扱いを点検する。
7. 版を確定したらGitHub Releaseを作成する。将来DOIが必要ならGitHub連携したZenodoでそのReleaseをアーカイブする。

## Gitコマンド例（リポジトリを空で作成した場合）

ZIPを展開してからそのディレクトリ内で実行する（リモートの名前は作成したGitHubリポジトリに合わせる）。

```bash
git init -b main
git add README.md CITATION.cff .gitignore src scripts analysis data docs
git commit -m 'Initialize Phobos primary-impact code archive'
git remote add origin git@github.com:hiroshikikuchi0515/phobos-primary-impact-simulation.git
git push -u origin main
```

SSH鍵を設定していない場合はHTTPSのURLを使う。**上記は公開可否の確認と、GitHubでの新規リポジトリ作成が済んでから実行すること。**
