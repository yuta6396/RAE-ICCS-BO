import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
import gradsio2 as gio  # 使用しているモジュール
from matplotlib.animation import FuncAnimation, PillowWriter

# 各種設定
outtype = 'fcst'                     # outtype in convert_letkfout.py
ncsuf = 'pe______.nc'
cmap = plt.cm.jet                   # カラーマップを指定
cmap.set_under('lightgray')         # 下限値の色
plt_extent = [125, 146, 26, 43]       # 描画範囲
norm = mcolors.Normalize(vmin=5, vmax=35)
TIME_INTERVAL = 900                # 例: 1時間を秒単位で表す（必要に応じて調整）
Objective_T_period = 24*int(3600/TIME_INTERVAL)              # 例: 6時間分のデータ
start_t = 10*int(3600/TIME_INTERVAL)   
base_dir = f"interval={TIME_INTERVAL}" # データのあるディレクトリ change!!
fny = 2
fnx = 2            # 格子数（例）

# --- 更新用の関数を定義 ---
def update(frame):
    """
    各フレーム（時刻 frame）ごとに降水量データを読み込み、描画内容を更新する関数
    history_name
    CTRL(base_dir = f"interval={TIME_INTERVAL}"): {base_dir}/no-control_interval={TIME_INTERVAL}_t={frame}
    TESTの例: ../test_result_1run/Min_Sum6hPREC_MPC_BO300_min_03-11-11-05  
                f"{base_dir}/interval={TIME_INTERVAL}_t={frame}"

    initial_history_name
    CTRL: {base_dir}/no-control_interval={TIME_INTERVAL}_t={start_t}
    """
    # 各時刻 frame に対応するディレクトリパス等を設定
    history_name = f"{base_dir}/no-control_interval={TIME_INTERVAL}_t={frame}"
    ncfilebase = f"{history_name}.{ncsuf}"
    
    try:
        # 例として 'PREC' 変数を読み込み、1hごとの降水量（単位変換済み）を取得
        var = gio.readvar2d('PREC', ncfilebase, outtype, fnx, fny)
    except ZeroDivisionError:
        print("ゼロ除算が発生しました。デフォルト値を使用します。")
        var = [np.zeros_like(slon)]
    except ValueError as e:
        print(f"値エラーが発生しました: {e}")
        var = [np.zeros_like(slon)]
    except Exception as e:
        print(f"予期しないエラーが発生しました: {e}")
        var = [np.zeros_like(slon)]
    
    # 単位変換 (例: kg/m/s から mm/h に変換する場合)
    this_prec = var[1] * 3600

    # タイトルを更新（例：時刻表示）
    time = frame - start_t
    hour = time//int(3600/TIME_INTERVAL)
    minite = time%int(3600/TIME_INTERVAL)
    if minite < int(3600/TIME_INTERVAL) - 1:
        ax.set_title(f"{hour}:{minite*15:02d}-{hour}:{(minite+1)*15}", fontsize=20, fontweight='bold', loc='center')
    else:
        ax.set_title(f"{hour}:{minite*15}-{hour+1}:{0:02d}", fontsize=20, fontweight='bold', loc='center')

    
    # 以前の描画を削除して新たな pcolormesh を描画
    global im
    im.remove()
    im = ax.pcolormesh(slon, slat, this_prec, cmap=cmap, norm=norm,
                       shading='auto', transform=ccrs.PlateCarree())
    return im

# --- 図と軸の設定 ---
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
ax.set_extent(plt_extent)
ax.coastlines(resolution='50m', linewidth=0.5)
ax.set_title("Total Precipitation", fontsize=20, fontweight='bold', loc='center')

# グリッド線の設定
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, linestyle='--', color='gray')
gl.xlocator = mticker.FixedLocator(range(125, 146))
gl.ylocator = mticker.FixedLocator(range(26, 43))
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 15}
gl.ylabel_style = {'size': 15}

# --- 初期フレームのデータ読み込み ---
initial_history_name = f"{base_dir}/no-control_interval={TIME_INTERVAL}_t={start_t}"
# initial_history_name = f"{base_dir}/interval={TIME_INTERVAL}_t=0"
ncfilebase_initial = f"{initial_history_name}.{ncsuf}"
slon = gio.readvar2d('lon', ncfilebase_initial, outtype, fnx, fny)
slat = gio.readvar2d('lat', ncfilebase_initial, outtype, fnx, fny)

# 初期フレームの precipitation データ（ゼロ行列）
prec_sum_matrix = np.zeros_like(slon)

# pcolormesh の初期描画
im = ax.pcolormesh(slon, slat, prec_sum_matrix, cmap=cmap, norm=norm,
                   shading='auto', transform=ccrs.PlateCarree())

# カラーバーの追加
cbar = fig.colorbar(im, ax=ax, orientation='vertical', pad=0.05, aspect=30)
cbar.set_label('mm', fontsize=18)
cbar.ax.tick_params(labelsize=15)

# --- アニメーションの作成 ---
anim = FuncAnimation(fig, update, frames=range(start_t, start_t+6*int(3600/TIME_INTERVAL)), blit=False)
# anim = FuncAnimation(fig, update, frames=range(0, 6*int(3600/TIME_INTERVAL)), blit=False)
# --- GIF として保存 ---
gif_filename = os.path.join(base_dir, "PREC-heatmap", f"LonLat_PREC_6h.gif")
# ディレクトリが存在しない場合は作成
os.makedirs(os.path.dirname(gif_filename), exist_ok=True)
writer = PillowWriter(fps=8)  # 例: 1秒間に2フレーム
anim.save(gif_filename, writer=writer, dpi=300)

# 動作確認用（必要に応じて）
# plt.show()
