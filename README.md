# RAE-ICCS-BO
## 概要
- SCALE-RMを用いた1800km*1800kmの気象シミュレーション並びに初期値介入実験を行える。
- .ncファイルはファイルサイズの都合上、自分で用意する or 著者からもらう。
- 著者はscale-5.5.3/scale-rm/test/tutorial/real/に配置しており、異なる階層に配置する場合はシンボリックリンクの張り替えが必要である。
- python fileが散乱しているが、最低限no_control_t=24.py (Free run用), gradio2.py (描画用), sim_1run.py, sim_ICCS_BO.pyだけあればよい。
## 下準備
### sim_1run.py
- init/にinit*.ncを用意
- run/にNo-control*, restart*, history*を用意
- run/cira.nc -> ../../../../../../data/rad/cira.nc
- config ファイルの確認
- Python3など実行に必要なパッケージのインストール
### sim_ICCS_BO.py
- 20 core以上のPCを想定しています。
- run/0000~run/0004にinit*.ncとhistory*.ncを用意する
- configファイルの確認

## 結果
### sim_1run.py
- 東北地方における6h累積降水量（目的関数値）に関する.pngと、結果の.ncが生成されます。
### sim_ICCS_BO.py
