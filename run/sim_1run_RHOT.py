import os
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import subprocess

import cartopy.crs as ccrs
import gradsio2 as gio
import numpy.ma as ma
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
import imageio  # GIF作成用
from datetime import datetime

# 時刻を計測するライブラリ
import time
import pytz
from datetime import datetime
from zoneinfo import ZoneInfo

import warnings

import random
import requests
matplotlib.use('Agg')

# 警告を抑制
warnings.filterwarnings('ignore', category=DeprecationWarning)

jst = pytz.timezone('Asia/Tokyo')# 日本時間のタイムゾーンを設定
current_time = datetime.now(jst).strftime("%m-%d-%H-%M")

"""
目的関数：0~t時間後までのある領域のPRECの総和
最適化変数：介入位置（x, y)とMOMX, MOMYそれぞれへの 介入量
"""

nofpe = 4
fny = 2
fnx = 2
X_size = 90
Y_size = 90

TIME_INTERVAL = 3600 #TIME_INTERVAL[sec]ごとに降水強度を出力できる
varname = 'PREC'

init_file = "../init/init_d01_20070714-180000.000.pe######.nc"  
org_file = "restart_t=10.pe######.nc"
history_file = "history_d01.pe######.nc" 

file_path = '/home/yuta/scale-5.5.3/scale-rm/test/tutorial/real/experiment_init_BO/run'

# 降水強度を最小化したい領域

Objective_X_low = 65
Objective_X_high = 75
Objective_Y_low = 60
Objective_Y_high = 70
Area_size=(Objective_X_high-Objective_X_low)*(Objective_Y_high-Objective_Y_low)
# 降水強度を最小化したい時間
Objective_T_period = 6 # この場合0~6までの7回加算


# 制御対象範囲
# Control_X_low = 45
# Control_X_high =90#AX =90
# Control_Y_low = 0
# Control_Y_high = 90 #AX =90
Control_Z_high = 4#MAX =36 36層存在 2の時429mが中心

# 介入領域の大きさ
Control_X_size = 5 #あんま変えられない　いくつか同時に変更する地点ありrandom_samples1 とか
Control_Y_size = 5
Control_Z_size = Control_Z_high+1
#　介入幅
Control_Var = 1.01

"""
"""

optimal_inputs=[62,50]

Control_purpose = "Min_Sum6hPREC" #Min_Max6hPREC, Min_Sum6hPREC
base_dir = f"../test_result_1run_RHOT/{Control_Var}_{optimal_inputs}_{current_time}"



def control_var_operation(var, pe, LB_Control_X, LB_Control_Y,  B_Control_Z):
    if 0 <= LB_Control_X <= 40 and 0 <= LB_Control_Y <= 40:
        if pe == 0:
            var[LB_Control_Y + 2 : LB_Control_Y + Control_Y_size + 2, LB_Control_X + 2 : LB_Control_X + Control_X_size + 2, B_Control_Z: B_Control_Z + Control_Z_size]*= Control_Var # (y, x, z) 3:5なら3, 4のみ
    elif  45 <= LB_Control_X <= 85 and 0 <= LB_Control_Y <= 40: # 右下
        if pe == 1:
            var[LB_Control_Y + 2 : LB_Control_Y + Control_Y_size + 2, LB_Control_X - 45 : LB_Control_X + Control_X_size - 45, B_Control_Z: B_Control_Z + Control_Z_size] *= Control_Var
    elif   0 <= LB_Control_X <= 40 and 45 <= LB_Control_Y <= 85: # 左上
        if pe == 2:
            var[LB_Control_Y - 45 : LB_Control_Y + Control_Y_size - 45, LB_Control_X + 2 : LB_Control_X + Control_X_size + 2, B_Control_Z: B_Control_Z + Control_Z_size] *= Control_Var
    elif   45 <= LB_Control_X <= 85 and 45 <= LB_Control_Y <= 85: # 右上
        if pe == 3:
            var[LB_Control_Y - 45 : LB_Control_Y + Control_Y_size - 45, LB_Control_X -45 : LB_Control_X + Control_X_size -45, B_Control_Z: B_Control_Z + Control_Z_size] *= Control_Var

    elif 0 <= LB_Control_X <= 40: # Y=41~44
        for Y_i in range(LB_Control_Y, 45): # pe == 0
            if pe == 0:
                var[Y_i + 2, LB_Control_X + 2 : LB_Control_X + Control_X_size + 2, B_Control_Z: B_Control_Z + Control_Z_size]*= Control_Var
        for Y_i in range(45, LB_Control_Y + Control_Y_size): # pe == 2
            if pe == 2:
                var[Y_i - 45, LB_Control_X + 2 : LB_Control_X + Control_X_size + 2, B_Control_Z: B_Control_Z + Control_Z_size]*= Control_Var

    elif 0 <= LB_Control_Y <= 40: # X=41~44
        for X_i in range(LB_Control_X, 45): # pe == 0
            if pe == 0:
                var[LB_Control_Y + 2 : LB_Control_Y + Control_Y_size + 2, X_i + 2, B_Control_Z: B_Control_Z + Control_Z_size]*= Control_Var
        for X_i in range(45, LB_Control_X + Control_X_size): # pe == 1
            if pe == 1:
                var[LB_Control_Y + 2 : LB_Control_Y + Control_Y_size + 2, X_i - 45, B_Control_Z: B_Control_Z + Control_Z_size]*= Control_Var

    elif 45 <= LB_Control_X <= 85: # Y=41~44
        for Y_i in range(LB_Control_Y, 45): # pe == 1
            if pe == 1:
                var[Y_i + 2, LB_Control_X - 45 : LB_Control_X + Control_X_size - 45, B_Control_Z: B_Control_Z + Control_Z_size]*= Control_Var
        for Y_i in range(45, LB_Control_Y + Control_Y_size): # pe == 3
            if pe == 3:
                var[Y_i - 45, LB_Control_X - 45 : LB_Control_X + Control_X_size - 45, B_Control_Z: B_Control_Z + Control_Z_size]*= Control_Var

    elif 45 <= LB_Control_Y <= 85: # X=41~44
        for X_i in range(LB_Control_X, 45): # pe == 2
            if pe == 2:
                var[LB_Control_Y - 45 : LB_Control_Y + Control_Y_size - 45, X_i + 2, B_Control_Z: B_Control_Z + Control_Z_size]*= Control_Var
        for X_i in range(45, LB_Control_X + Control_X_size): # pe == 3
            if pe == 3:
                var[LB_Control_Y - 45 : LB_Control_Y + Control_Y_size - 45, X_i - 45, B_Control_Z: B_Control_Z + Control_Z_size]*= Control_Var

    elif 41 <= LB_Control_X <= 44 and 41 <= LB_Control_Y <= 44:
        for X_i in range(LB_Control_X, 45): 
            for Y_i in range(LB_Control_Y, 45): # pe == 0
                if pe == 0:
                    var[Y_i + 2, X_i + 2, B_Control_Z: B_Control_Z + Control_Z_size]*= Control_Var
            for Y_i in range(45, LB_Control_Y + Control_Y_size): # pe == 2
                if pe == 2:
                    var[Y_i - 45, X_i + 2, B_Control_Z: B_Control_Z + Control_Z_sizeZ]*= Control_Var
        for X_i in range(45, LB_Control_X + Control_X_size):
            for Y_i in range(LB_Control_Y, 45): # pe == 0
                if pe == 1:
                    var[Y_i + 2, X_i - 45, B_Control_Z: B_Control_Z + Control_Z_size]*= Control_Var
            for Y_i in range(45, LB_Control_Y + Control_Y_size): # pe == 2
                if pe == 3:
                    var[Y_i - 45, X_i  - 45, B_Control_Z: B_Control_Z + Control_Z_size]*= Control_Var
    else:
        print("インデックスがリストの範囲外です。")
        sis.exit("Error")
    return var

def sim(LB_Control_X,LB_Control_Y, f):
    """
    得られた最適解を用いて目的関数の値を再度計算する。
    制御しない場合と制御した場合における、ある時刻のある領域の降水強度の値を返す。
    """

    # 目的の降水強度
    TEST_prec_matrix=np.zeros((Objective_T_period+1, Objective_X_high-Objective_X_low,Objective_Y_high-Objective_Y_low)) #ベイズ最適化の累積降水量  
    CTRL_prec_matrix=np.zeros((Objective_T_period+1, Objective_X_high-Objective_X_low,Objective_Y_high-Objective_Y_low)) # 制御前の累積降水量
    TEST_CTRL_prec_matrix = np.zeros((Objective_T_period+1, Objective_X_high-Objective_X_low,Objective_Y_high-Objective_Y_low)) #制御あり -制御なし　の各地点のある時刻の降水強度　負の値ほど、良い制御
    TEST_prec_sum=0
    CTRL_prec_sum=0

    for pe in range(nofpe):
        f.write(f"\n\n\n{pe=}\n")
        output_file = f"out.pe######.nc"
        # input file
        init = init_file.replace('######', str(pe).zfill(6))
        org = org_file.replace('######', str(pe).zfill(6))
        history = history_file.replace('######', str(pe).zfill(6))
        output = output_file.replace('######', str(pe).zfill(6))
        history_path = file_path+'/'+history
        if(os.path.isfile(history_path)):
            subprocess.run(["rm",history])
        subprocess.run(["cp", org, init]) #初期化

        with netCDF4.Dataset(init) as src, netCDF4.Dataset(output, "w") as dst:
            # global attributes のコピー
            dst.setncatts(src.__dict__)
            
            # dimensions のコピー
            for name, dimension in src.dimensions.items():
                dst.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))
            
            # 各変数のコピーと、最大・最小値およびインデックスの出力
            for name, variable in src.variables.items():
                if pe == 0:
                    f.write(f"Processing variable:{name}\n")
                # 出力用の変数を作成
                x = dst.createVariable(name, variable.datatype, variable.dimensions)
                dst[name].setncatts(src[name].__dict__)
                
                # 変数データの取得
                data = src[name][:]
                
                # 数値型でデータが存在する場合、最大値・最小値とその位置を計算して出力
                if np.issubdtype(data.dtype, np.number) and data.size > 0:
                    try:
                        max_val = np.nanmax(data)
                        min_val = np.nanmin(data)
                        # np.nanargmax/np.nanargmin はフラットなインデックスを返すので、元のshapeに戻す
                        max_index = np.unravel_index(np.nanargmax(data), data.shape)
                        min_index = np.unravel_index(np.nanargmin(data), data.shape)
                        f.write(f"Variable {name}: max = {max_val} at index {max_index}, min = {min_val} at index {min_index}\n")
                    except Exception as e:
                        f.write(f"Variable {name}: Error computing min/max: {e}\n")
                else:
                    f.write(f"Variable {name} is non-numeric or empty; skipping min/max check\n")
                
                # 変数名が特定の場合は操作を実施、それ以外はそのままコピー
                if name == "RHOT":
                    dst[name][:] = control_var_operation(data, pe, LB_Control_X, LB_Control_Y, 0)
                else:
                    dst[name][:] = data
        subprocess.run(["cp", output, init])

    for loop in range(0, Objective_T_period):
        TEST_prec_sum, CTRL_prec_sum, TEST_prec_matrix, CTRL_prec_matrix, TEST_CTRL_prec_matrix = state_update(loop, TEST_prec_matrix, CTRL_prec_matrix, TEST_CTRL_prec_matrix, TEST_prec_sum, CTRL_prec_sum)

    return TEST_prec_sum, CTRL_prec_sum, TEST_prec_matrix, CTRL_prec_matrix ,TEST_CTRL_prec_matrix
   
def state_update(loop, TEST_prec_matrix, CTRL_prec_matrix, TEST_CTRL_prec_matrix, TEST_prec_sum, CTRL_prec_sum):
    orgfile = f"no-control_t={loop+11}.pe######.nc"
    gpyoptfile=f"{base_dir}/t={loop}.pe######.nc" 

    for pe in range(nofpe):
        history = history_file.replace('######', str(pe).zfill(6))
        history_path = file_path+'/'+history
        if (os.path.isfile(history_path)):
            subprocess.run(["rm", history])
    subprocess.run(["mpirun", "-n", str(nofpe), "./scale-rm","run.d02.conf"])

    for pe in range(nofpe):
        gpyopt = gpyoptfile.replace('######', str(pe).zfill(6))
        history = history_file.replace('######', str(pe).zfill(6))
        subprocess.run(["cp", history,gpyopt])

    for pe in range(nofpe):  # history処理
        fiy, fix = np.unravel_index(pe, (fny, fnx))
        nc = netCDF4.Dataset(history_file.replace('######', str(pe).zfill(6)))
        onc = netCDF4.Dataset(orgfile.replace('######', str(pe).zfill(6)))
        nt = onc.dimensions['time'].size
        nx = nc.dimensions['x'].size
        ny = nc.dimensions['y'].size
        nz = nc.dimensions['z'].size
        gx1 = nx * fix
        gx2 = nx * (fix + 1)
        gy1 = ny * fiy
        gy2 = ny * (fiy + 1)
        if pe == 0:
            dat = np.zeros((nt, nz, fny*ny, fnx*nx))
            odat = np.zeros((nt, nz, fny*ny, fnx*nx))
        dat[:, 0, gy1:gy2, gx1:gx2] = nc[varname][:]
        odat[:, 0, gy1:gy2, gx1:gx2] = onc[varname][:]

    # データから目的関数の値を計算する
    # 目的関数に該当する領域以外もPRECは計算しない
    for j in range(Objective_X_low,Objective_X_high):
        for k in range(Objective_Y_low,Objective_Y_high):
            TEST_prec_matrix[loop, j-Objective_X_low, k-Objective_Y_low] += dat[1,0,k,j]*TIME_INTERVAL #(y,x,z)=(time,z,y,x)　j-Objective_X_low,k-Objective_Y_lowは[0,0]->[5,5]とか
            CTRL_prec_matrix[loop,j-Objective_X_low, k-Objective_Y_low] += odat[1,0,k,j]*TIME_INTERVAL
            TEST_CTRL_prec_matrix[loop, j-Objective_X_low, k-Objective_Y_low] = TEST_prec_matrix[loop, j-Objective_X_low, k-Objective_Y_low] - CTRL_prec_matrix[loop, j-Objective_X_low, k-Objective_Y_low]
            TEST_prec_sum+=TEST_prec_matrix[loop, j-Objective_X_low, k-Objective_Y_low]
            CTRL_prec_sum+=CTRL_prec_matrix[loop, j-Objective_X_low, k-Objective_Y_low]
    print(f"{TEST_prec_sum=}")
    print(f"{CTRL_prec_sum=}")
    return TEST_prec_sum, CTRL_prec_sum, TEST_prec_matrix, CTRL_prec_matrix ,TEST_CTRL_prec_matrix


def plot_PREC(TEST_prec_matrix, TEST_CTRL_prec_matrix, Objective_ratio):
    """
    最適解（制御入力）でシミュレーションした結果の可視化関数。
    ・目的時刻までの降水強度総和（シミュレーションと対照実験）のマップを描画  
    ・両者の差分マップ  
    ・グリッド上での累積降水量  
    ・各ループ（時刻）ごとの累積値（シミュレーション/対照）をまとめたマップ
    を作成する。
    """

    # ----- 各種設定 -----
    outtype = 'fcst'
    ncsuf = 'pe______.nc'
    cmap = plt.cm.jet
    cmap.set_under('lightgray')
    plt_extent  = [125, 146, 26, 43]      # 全体の描画範囲
    plt_extent_sub = [138, 144, 35, 41]     # 詳細表示用の範囲
    # plt_extent_sub = plt_extent
    norm      = mcolors.Normalize(vmin=15, vmax=165)
    norm_diff = mcolors.Normalize(vmin=-20, vmax=20)
    norm_1h = mcolors.Normalize(vmin=5, vmax=35)
    norm_diff_1h = mcolors.Normalize(vmin=-10, vmax=10)

    # ループごとの累積値を保存するリスト（各時刻ごとの状態を後で可視化）
    prec_list = []      # シミュレーション結果
    ctrl_prec_list = [] # 対照実験結果

    prec_sum = None
    ctrl_prec_sum = None

    # ----- 時刻ループ：各時刻ごとにデータを読み込み累積する -----
    for loop in range(Objective_T_period):
        # シミュレーション・対照それぞれのファイル名を生成
        history_name      = f"{base_dir}/t={loop}"
        CTRL_history_name = f"/home/yuta/scale-5.5.3/scale-rm/test/tutorial/real/experiment_init_BO/run/no-control_t={loop+11}"
        ncfile_sim  = f"{history_name}.{ncsuf}"
        ncfile_ctrl = f"{CTRL_history_name}.{ncsuf}"

        # シミュレーションの降水量データを読み込み
        try:
            var = gio.readvar2d('PREC', ncfile_sim, outtype, fnx, fny)
        except Exception as e:
            print(f"Error reading {ncfile_sim}: {e}")
            # エラー時は、適切な形状のゼロ配列を用いる（初回なら仮のサイズ fnx×fny）
            var = [np.zeros((fnx, fny)), np.zeros((fnx, fny))]
        # 対照実験の降水量データ
        try:
            ctrl_var = gio.readvar2d('PREC', ncfile_ctrl, outtype, fnx, fny)
        except Exception as e:
            print(f"Error reading {ncfile_ctrl}: {e}")
            ctrl_var = [np.zeros((fnx, fny)), np.zeros((fnx, fny))]

        # 単位変換 (例：kg/m/s → mm/h)
        var      = [v * TIME_INTERVAL for v in var]
        ctrl_var = [v * TIME_INTERVAL for v in ctrl_var]

        # 初回ループでは緯度経度情報も取得し、累積値を初期化
        if loop == 0:
            slon = gio.readvar2d('lon', ncfile_sim, outtype, fnx, fny)
            slat = gio.readvar2d('lat', ncfile_sim, outtype, fnx, fny)
            prec_sum      = var[1].copy()
            ctrl_prec_sum = ctrl_var[1].copy()
        else:
            prec_sum      += var[1]
            ctrl_prec_sum += ctrl_var[1]

        # 各時刻ごとの累積値をリストに保存（後でループごとに可視化）
        prec_list.append(var[1].copy())
        ctrl_prec_list.append(ctrl_var[1].copy())

    # 出力ディレクトリが存在しない場合は作成
    out_dir = os.path.join(base_dir, "PREC-heatmap")
    os.makedirs(out_dir, exist_ok=True)

    # ----- 1. 最終時刻の累積降水量マップ（TEST） -----
    fig1 = plt.figure(figsize=(12, 8))
    ax1 = fig1.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax1.set_extent(plt_extent_sub)
    ax1.coastlines(resolution='50m', linewidth=0.5)
    gl1 = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5,
                        linestyle='--', color='gray')
    gl1.xlocator = mticker.FixedLocator(range(125, 146))
    gl1.ylocator = mticker.FixedLocator(range(26, 43))
    gl1.top_labels = gl1.right_labels = False
    im1 = ax1.pcolormesh(slon, slat, prec_sum, cmap=cmap, norm=norm,
                         shading='auto', transform=ccrs.PlateCarree())
    cbar1 = fig1.colorbar(im1, ax=ax1, orientation='vertical', pad=0.05, aspect=30)
    cbar1.set_label('mm', fontsize=18)
    plt.savefig(os.path.join(out_dir, f"LonLat_{Objective_ratio}%.pdf"), dpi=1200)
    # plt.savefig(os.path.join(out_dir, f"LonLat_{Objective_ratio}%.pdf"), dpi=1200, bbox_inches = 'tight')

    # ----- 2. 累積降水量の差分マップ（TEST - CTRL） -----
    mask = np.ones(prec_sum.shape, dtype=bool)
    mask[60:69, 65:74] = False
    data = prec_sum - ctrl_prec_sum
    # masked_data = ma.array(data, mask=mask)
    masked_data = data
    cmap.set_bad(color='none')

    fig2 = plt.figure(figsize=(10, 8))
    ax2 = fig2.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax2.set_extent(plt_extent_sub)
    ax2.coastlines(resolution='50m', linewidth=0.5)
    gl2 = ax2.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5,
                        linestyle='--', color='gray')
    gl2.xlocator = mticker.FixedLocator(range(139, 143))
    gl2.ylocator = mticker.FixedLocator(range(36, 40))
    gl2.top_labels = gl2.right_labels = False
    im2 = ax2.pcolormesh(slon, slat, masked_data, cmap=cmap, norm=norm_diff,
                         shading='auto', transform=ccrs.PlateCarree())
    cbar2 = fig2.colorbar(im2, ax=ax2, orientation='vertical', pad=0.05, aspect=30)
    cbar2.set_label('mm', fontsize=18)
    plt.savefig(os.path.join(out_dir, f"Sub_LonLat_{Objective_ratio}%.pdf"), dpi=1200, bbox_inches = 'tight')
    # plt.savefig(os.path.join(out_dir, f"Sub_LonLat_{Objective_ratio}%.pdf"), dpi=1200, bbox_inches = 'tight')


    # ----- 5. ループごとの累積値の変化を可視化（シミュレーション/対照） -----
    # 2行 Objective_T_period列のサブプロットを作成
    fs3 = 12
    fig3, axs = plt.subplots(3, Objective_T_period, figsize=(4 * Objective_T_period, 12),
                              subplot_kw={'projection': ccrs.PlateCarree()})
    for i in range(Objective_T_period):
        ax_sim = axs[0, i]
        ax_sim.set_extent(plt_extent_sub)
        ax_sim.coastlines(resolution='50m', linewidth=0.5)
        im = ax_sim.pcolormesh(slon, slat, prec_list[i], cmap=cmap, norm=norm_1h,
                               shading='auto', transform=ccrs.PlateCarree())
        ax_sim.set_title(f"TEST t={i}-{i+1}h", fontsize=fs3)
    for i in range(Objective_T_period):
        ax_ctrl = axs[1, i]
        ax_ctrl.set_extent(plt_extent_sub)
        ax_ctrl.coastlines(resolution='50m', linewidth=0.5)
        im = ax_ctrl.pcolormesh(slon, slat, ctrl_prec_list[i], cmap=cmap, norm=norm_1h,
                                shading='auto', transform=ccrs.PlateCarree())
        ax_ctrl.set_title(f"CTRL t={i}-{i+1}h", fontsize=fs3)

    for i in range(Objective_T_period):
        ax_diff = axs[2, i]
        ax_diff.set_extent(plt_extent_sub)
        ax_diff.coastlines(resolution='50m', linewidth=0.5)
        im = ax_diff.pcolormesh(slon, slat, prec_list[i] - ctrl_prec_list[i], cmap=cmap, norm=norm_diff_1h,
                                shading='auto', transform=ccrs.PlateCarree())
        ax_diff.set_title(f"TEST - CTRL t={i}-{i+1}h", fontsize=fs3)    
    fig3.subplots_adjust(right=0.9)
    # 上段（TEST, CTRL）のカラーバー用 Axes（全体の上半分、例：上から55%位置から35%の高さ）
    cbar_ax_top = fig3.add_axes([0.92, 0.40, 0.02, 0.50])
    # 下段（TEST-CTRL）のカラーバー用 Axes（例：下から15%位置から35%の高さ）
    cbar_ax_bot = fig3.add_axes([0.92, 0.1, 0.02, 0.25])

    # ScalarMappable を作成（各カラーバーのため）
    sm_top = plt.cm.ScalarMappable(norm=norm_1h, cmap=cmap)
    sm_top.set_array([])  # 値は不要
    sm_bot = plt.cm.ScalarMappable(norm=norm_diff_1h, cmap=cmap)
    sm_bot.set_array([])

    # カラーバーを配置
    fig3.colorbar(sm_top, cax=cbar_ax_top, label='mm')
    fig3.colorbar(sm_bot, cax=cbar_ax_bot, label='mm')
    plt.savefig(os.path.join(out_dir, f"Loop_by_loop_{Objective_ratio}%.png"), dpi=1200)
    # plt.savefig(os.path.join(out_dir, f"Loop_by_loop_{Objective_ratio}%.pdf"), dpi=1200, bbox_inches = 'tight')
    # 不要なウィンドウを閉じる
    plt.close('all')

def make_directory(base_dir):
    os.makedirs(base_dir, exist_ok=False)
    # 階層構造を作成
    sub_dirs = ["PREC-heatmap", "summary"]
    for sub_dir in sub_dirs:
        path = os.path.join(base_dir, sub_dir)
        os.makedirs(path, exist_ok=True) 
    return


### 実行
def main():
    make_directory(base_dir)
    filename = f"config.txt"
    config_file_path = os.path.join(base_dir, filename) 
    # ファイルに書き込む
    with open(config_file_path,  "w") as f:
        # 降水強度を最小化したい時刻
        f.write("\n# 降水強度を最小化したい時間\n")
        f.write(f"{Objective_T_period=}\n")
        f.write(f"Control_Z_high = {Control_Z_high}  # MAX =35?\n")
        
        # 介入領域の大きさ
        f.write("\n# 介入領域の大きさ\n")
        f.write(f"Control_X_size = {Control_X_size}\n")
        f.write(f"Control_Y_size = {Control_Y_size}\n")
        f.write(f"Control_Z_size = {Control_Z_size}\n")
        f.write(f"{optimal_inputs=}")
        
    filename = f"variable_maxmin.txt"
    config_file_path = os.path.join(base_dir, filename) 
    with open(config_file_path,  "w") as f:
        # シミュレーションの実行
        for pe in range(nofpe):
            org = org_file.replace('######', str(pe).zfill(6))
            init = init_file.replace('######', str(pe).zfill(6))
            subprocess.run(["cp", org, init])

        TEST_prec_sum, CTRL_prec_sum, TEST_prec_matrix, CTRL_prec_matrix , TEST_CTRL_prec_matrix=sim(optimal_inputs[0],optimal_inputs[1],f)#  change
    # シミュレーション結果
    print(f"{TEST_prec_sum=}")
    print(f"{CTRL_prec_sum=}")
    Objective_ratio= TEST_prec_sum/CTRL_prec_sum
    print(f"%={Objective_ratio*100}%")
    plot_PREC( TEST_prec_matrix, TEST_CTRL_prec_matrix, Objective_ratio)

    print(base_dir)


def notify_slack(webhook_url, message, channel=None, username=None, icon_emoji=None):
    """
    Slackに通知を送信する関数。

    :param webhook_url: SlackのWebhook URL
    :param message: 送信するメッセージ
    :param channel: メッセージを送信するチャンネル（オプション）
    :param username: メッセージを送信するユーザー名（オプション）
    :param icon_emoji: メッセージに表示する絵文字（オプション）
    """
    payload = {
        "text": message
    }

    # オプションのパラメータを追加
    if channel:
        payload["channel"] = channel
    if username:
        payload["username"] = username
    if icon_emoji:
        payload["icon_emoji"] = icon_emoji

    try:
        response = requests.post(webhook_url, json=payload)
        response.raise_for_status()  # エラーがあれば例外を発生させる
        print("Slackへの通知が送信されました。")
    except requests.exceptions.RequestException as e:
        print(f"Slackへの通知に失敗しました: {e}")

def get_script_name():
    return os.path.basename(__file__)

if __name__ == "__main__":
    webhook_url = os.getenv("SLACK_WEBHOOK_URL")  # export SLACK_WEBHOOK_URL="OOOO"としておく
    try:
        main()
    except Exception as e:
        error_message = f"❌ {get_script_name()}でエラーが発生しました: {e}"
        notify_slack(webhook_url, error_message, channel="webhook")
        # エラー発生時はその後の処理も中断するため再度例外を送出する
        raise
    else:
        success_message = f"✅ {get_script_name()}の処理が完了しました。"
        notify_slack(webhook_url, success_message, channel="webhook")