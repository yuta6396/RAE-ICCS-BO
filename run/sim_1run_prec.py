import os
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import subprocess

import cartopy.crs as ccrs
import gradsio2 as gio
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
目的関数：t時間後のある領域のPRECの総和
最適化変数：介入位置（x, y, z)と介入量
"""

nofpe = 4
fny = 2
fnx = 2
X_size = 90
Y_size = 90

TIME_INTERVAL = 3600 #TIME_INTERVAL[sec]ごとに降水強度を出力できる
varname = 'PREC'

init_file = "../init/init_d01_20070714-180000.000.pe######.nc"  #使っていない
# org_file = "../init/init_d01_20070714-180000.000.pe######.org.nc"
org_file = "restart_t=10.pe######.nc"
history_file = "history_d01.pe######.nc"

orgfile = 'no-control_24h.pe######.nc'
file_path = '/home/yuta/scale-5.5.3/scale-rm/test/tutorial/real/experiment_init_BO/run'
gpyoptfile=f"BO_init_4var.pe######.nc"

# 降水強度を最小化したい領域

Objective_X_low = 65
Objective_X_high = 75
Objective_Y_low = 60
Objective_Y_high = 70
Area_size=(Objective_X_high-Objective_X_low)*(Objective_Y_high-Objective_Y_low)
# 降水強度を最小化したい時間
Objective_T_period = 6 # この場合0~6までの7回加算


# 制御対象範囲
Control_X_low = 45
Control_X_high =90#AX =90
Control_Y_low = 0
Control_Y_high = 90 #AX =90
Control_Z_high = 4#MAX =36 36層存在 2の時429mが中心

# 介入領域の大きさ
Control_X_size = 5 #あんま変えられない　いくつか同時に変更する地点ありrandom_samples1 とか
Control_Y_size = 5
Control_Z_size = Control_Z_high+1
#　介入幅
Control_Var_low = -20   #30~30でエラー確認済み(MOMY)
Control_Var_high = 20

Objective_T = 5

optimal_inputs=[71, 57, 0,0]





base_dir = f"../test_result_1run/prec_{optimal_inputs}_{current_time}"


def control_var_operation(var, pe, LB_Control_X, LB_Control_Y, Control_Y_size, Control_X_size, B_Control_Z, Control_Var):
    if 0 <= LB_Control_X <= 40 and 0 <= LB_Control_Y <= 40:
        if pe == 0:
            var[LB_Control_Y + 2 : LB_Control_Y + Control_Y_size + 2, LB_Control_X + 2 : LB_Control_X + Control_X_size + 2, B_Control_Z: B_Control_Z + Control_Z_size]+= Control_Var # (y, x, z) 3:5なら3, 4のみ
    elif  45 <= LB_Control_X <= 85 and 0 <= LB_Control_Y <= 40: # 右下
        if pe == 1:
            var[LB_Control_Y + 2 : LB_Control_Y + Control_Y_size + 2, LB_Control_X - 45 : LB_Control_X + Control_X_size - 45, B_Control_Z: B_Control_Z + Control_Z_size] += Control_Var
    elif   0 <= LB_Control_X <= 40 and 45 <= LB_Control_Y <= 85: # 左上
        if pe == 2:
            var[LB_Control_Y - 45 : LB_Control_Y + Control_Y_size - 45, LB_Control_X + 2 : LB_Control_X + Control_X_size + 2, B_Control_Z: B_Control_Z + Control_Z_size] += Control_Var
    elif   45 <= LB_Control_X <= 85 and 45 <= LB_Control_Y <= 85: # 右上
        if pe == 3:
            var[LB_Control_Y - 45 : LB_Control_Y + Control_Y_size - 45, LB_Control_X -45 : LB_Control_X + Control_X_size -45, B_Control_Z: B_Control_Z + Control_Z_size] += Control_Var

    elif 0 <= LB_Control_X <= 40: # Y=41~44
        for Y_i in range(LB_Control_Y, 45): # pe == 0
            if pe == 0:
                var[Y_i + 2, LB_Control_X + 2 : LB_Control_X + Control_X_size + 2, B_Control_Z: B_Control_Z + Control_Z_size]+= Control_Var
        for Y_i in range(45, LB_Control_Y + Control_Y_size): # pe == 2
            if pe == 2:
                var[Y_i - 45, LB_Control_X + 2 : LB_Control_X + Control_X_size + 2, B_Control_Z: B_Control_Z + Control_Z_size]+= Control_Var

    elif 0 <= LB_Control_Y <= 40: # X=41~44
        for X_i in range(LB_Control_X, 45): # pe == 0
            if pe == 0:
                var[LB_Control_Y + 2 : LB_Control_Y + Control_Y_size + 2, X_i + 2, B_Control_Z: B_Control_Z + Control_Z_size]+= Control_Var
        for X_i in range(45, LB_Control_X + Control_X_size): # pe == 1
            if pe == 1:
                var[LB_Control_Y + 2 : LB_Control_Y + Control_Y_size + 2, X_i - 45, B_Control_Z: B_Control_Z + Control_Z_size]+= Control_Var

    elif 45 <= LB_Control_X <= 85: # Y=41~44
        for Y_i in range(LB_Control_Y, 45): # pe == 1
            if pe == 1:
                var[Y_i + 2, LB_Control_X - 45 : LB_Control_X + Control_X_size - 45, B_Control_Z: B_Control_Z + Control_Z_size]+= Control_Var
        for Y_i in range(45, LB_Control_Y + Control_Y_size): # pe == 3
            if pe == 3:
                var[Y_i - 45, LB_Control_X - 45 : LB_Control_X + Control_X_size - 45, B_Control_Z: B_Control_Z + Control_Z_size]+= Control_Var

    elif 45 <= LB_Control_Y <= 85: # X=41~44
        for X_i in range(LB_Control_X, 45): # pe == 2
            if pe == 2:
                var[LB_Control_Y - 45 : LB_Control_Y + Control_Y_size - 45, X_i + 2, B_Control_Z: B_Control_Z + Control_Z_size]+= Control_Var
        for X_i in range(45, LB_Control_X + Control_X_size): # pe == 3
            if pe == 3:
                var[LB_Control_Y - 45 : LB_Control_Y + Control_Y_size - 45, X_i - 45, B_Control_Z: B_Control_Z + Control_Z_size]+= Control_Var

    elif 41 <= LB_Control_X <= 44 and 41 <= LB_Control_Y <= 44:
        for X_i in range(LB_Control_X, 45): 
            for Y_i in range(LB_Control_Y, 45): # pe == 0
                if pe == 0:
                    var[Y_i + 2, X_i + 2, B_Control_Z: B_Control_Z + Control_Z_size]+= Control_Var
            for Y_i in range(45, LB_Control_Y + Control_Y_size): # pe == 2
                if pe == 2:
                    var[Y_i - 45, X_i + 2, B_Control_Z: B_Control_Z + Control_Z_sizeZ]+= Control_Var
        for X_i in range(45, LB_Control_X + Control_X_size):
            for Y_i in range(LB_Control_Y, 45): # pe == 0
                if pe == 1:
                    var[Y_i + 2, X_i - 45, B_Control_Z: B_Control_Z + Control_Z_size]+= Control_Var
            for Y_i in range(45, LB_Control_Y + Control_Y_size): # pe == 2
                if pe == 3:
                    var[Y_i - 45, X_i  - 45, B_Control_Z: B_Control_Z + Control_Z_size]+= Control_Var
    else:
        print("インデックスがリストの範囲外です。")
        sis.exit("Error")
    return var

def sim(LB_Control_X,LB_Control_Y,Control_MOMX, Control_MOMY):
    """
    得られた最適解を用いて目的関数の値を再度計算する。
    制御しない場合と制御した場合における、ある時刻のある領域の降水強度の値を返す。
    """

    # 目的の降水強度
    TEST_prec_matrix=np.zeros((Objective_X_high-Objective_X_low,Objective_Y_high-Objective_Y_low)) #ベイズ最適化の累積降水量  
    CTRL_prec_matrix=np.zeros((Objective_X_high-Objective_X_low,Objective_Y_high-Objective_Y_low)) # 制御前の累積降水量
    TEST_CTRL_prec_matrix = np.zeros((Objective_X_high-Objective_X_low,Objective_Y_high-Objective_Y_low)) #制御あり -制御なし　の各地点のある時刻の降水強度　負の値ほど、良い制御
    TEST_prec_sum=0
    CTRL_prec_sum=0

    for pe in range(nofpe):
        output_file = f"out_d01.pe######.nc"
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
            # copy global attributes all at once via dictionary
            dst.setncatts(src.__dict__)
            # copy dimensions
            for name, dimension in src.dimensions.items():
                dst.createDimension(
                    name, (len(dimension) if not dimension.isunlimited() else None))
            # copy all file data except for the excluded
            for name, variable in src.variables.items():
                x = dst.createVariable(
                    name, variable.datatype, variable.dimensions)
                # copy variable attributes all at once via dictionary
                dst[name].setncatts(src[name].__dict__)
                if name == "MOMX":
                    var = src[name][:]
                    dst[name][:] = control_var_operation(var, pe, LB_Control_X, LB_Control_Y, Control_Y_size, Control_X_size, 0, Control_MOMX)
                elif name == "MOMY":
                    var = src[name][:]
                    dst[name][:] = control_var_operation(var, pe, LB_Control_X, LB_Control_Y, Control_Y_size, Control_X_size, 0, Control_MOMY)
                else:
                    dst[name][:] = src[name][:]
        subprocess.run(["cp", output, init])

    subprocess.run(["mpirun", "-n", str(nofpe), "./scale-rm","run.d01.conf"])

    # run　の結果から、目的関数の計算のために必要なデータを取ってくる。
    # time.sleep(0.3)
    for pe in range(nofpe):
        gpyopt = gpyoptfile.replace('######', str(pe).zfill(6))
        history = history_file.replace('######', str(pe).zfill(6))
        subprocess.run(["cp", history,gpyopt])
    for pe in range(nofpe):  # history処理
        fiy, fix = np.unravel_index(pe, (fny, fnx))
        nc = netCDF4.Dataset(history_file.replace('######', str(pe).zfill(6)))
        onc = netCDF4.Dataset(orgfile.replace('######', str(pe).zfill(6)))
        nt = nc.dimensions['time'].size
        nx = nc.dimensions['x'].size
        ny = nc.dimensions['y'].size
        nz = nc.dimensions['z'].size
        gx1 = nx * fix
        gx2 = nx * (fix + 1)
        gy1 = ny * fiy
        gy2 = ny * (fiy + 1)
        # print(f"onc[varname][:] = {onc[varname][:]}")
        if pe == 0:
            dat = np.zeros((nt, nz, fny*ny, fnx*nx))
            odat = np.zeros((nt, nz, fny*ny, fnx*nx))
        dat[:, 0, gy1:gy2, gx1:gx2] = nc[varname][:]
        odat[:, 0, gy1:gy2, gx1:gx2] = onc[varname][:]

    # データから目的関数の値を計算する
    # 目的関数に該当する領域以外もPRECは計算しない
    for j in range(Objective_X_low,Objective_X_high):
        for k in range(Objective_Y_low,Objective_Y_high):
            TEST_prec_matrix[j-Objective_X_low,k-Objective_Y_low] += dat[Objective_T,0,k,j]*TIME_INTERVAL #(y,x,z)=(time,z,y,x)　j-Objective_X_low,k-Objective_Y_lowは[0,0]->[5,5]とか
            CTRL_prec_matrix[j-Objective_X_low,k-Objective_Y_low] += odat[Objective_T,0,k,j]*TIME_INTERVAL
            TEST_CTRL_prec_matrix[j-Objective_X_low,k-Objective_Y_low] = TEST_prec_matrix[j-Objective_X_low,k-Objective_Y_low] - CTRL_prec_matrix[j-Objective_X_low,k-Objective_Y_low]
            TEST_prec_sum+=TEST_prec_matrix[j-Objective_X_low,k-Objective_Y_low]
            CTRL_prec_sum+=CTRL_prec_matrix[j-Objective_X_low,k-Objective_Y_low]
    return TEST_prec_sum, CTRL_prec_sum, TEST_prec_matrix, CTRL_prec_matrix ,TEST_CTRL_prec_matrix

def plot_PREC(TEST_prec_matrix, TEST_CTRL_prec_matrix):
    """
    得られた最適解（制御入力）でシミュレーションをした結果を、"可視化する"関数。
    目的時刻t=Objective_T の降水強度を可視化する
    """
    outtype = 'fcst'                    # outtype in convert_letkfout.py
    ncsuf = 'pe______.nc'

    # params for plot
    cmap = plt.cm.jet
    cmap.set_under('lightgray')
    plt_extent = [125,146,26,43]
    history_name = 'history_d01'
    norm = mcolors.Normalize(vmin=5, vmax=35)
    ncfilebasetmp = '{:s}.{:}'.format(history_name, ncsuf)

    # get lon lat
    slon = gio.readvar2d('lon',ncfilebasetmp,outtype,fnx,fny)
    slat = gio.readvar2d('lat',ncfilebasetmp,outtype,fnx,fny)
    try:
        var = gio.readvar2d('PREC',ncfilebasetmp,outtype,fnx,fny)
    except ZeroDivisionError:
        print("ゼロ除算が発生しました。デフォルト値を使用します。")
        var = 0  # または適切なデフォルト値
    except ValueError as e:
        print(f"値エラーが発生しました: {e}")
        var = 0  # または適切なデフォルト値
    except Exception as e:
        print(f"予期しないエラーが発生しました: {e}")
        var = 0  # または適切なデフォルト値

    # unit conversion (kg/m/s --> mm/d)
    var *=  60 * 60
    print(f"{len(var)=}")
    for t in range(len(var)):
        print(t)
        # plot
        fig = plt.figure(figsize=(12, 8))

        # scale
        ax = fig.add_subplot(1,1,1,projection=ccrs.PlateCarree())
        ax.set_extent(plt_extent)
        ax.coastlines(resolution='50m',linewidth=0.5)

        # grid line
        gl = ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=True,linewidth=0.5,linestyle='--',color='gray')
        gl.xlocator = mticker.FixedLocator(range(125, 146))
        gl.ylocator = mticker.FixedLocator(range(26, 43)) # 44
        gl.top_labels = False
        gl.right_labels = False
        gl.xlabel_style = {'size': 15}
        gl.ylabel_style = {'size': 15}

        im = ax.pcolormesh(slon,slat,var[t],cmap=cmap, norm=norm,shading='auto', transform=ccrs.PlateCarree())

        # pos = ax.get_position()
        # fig.subplots_adjust(right=0.85)
        # cax = fig.add_axes([pos.x1,pos.y0,0.02,pos.y1-pos.y0])
        # clb = fig.colorbar(im, orientation='vertical', cax=cax)
        cbar = fig.colorbar(im, ax=ax, orientation='vertical', pad=0.05, aspect=30)
        cbar.set_label('mm/h', fontsize=18)
        cbar.ax.tick_params(labelsize=15)
        # filename = os.path.join(base_dir, "PREC-heatmap", f"LonLat_t={t}.pdf")
        # plt.savefig(filename ,dpi=1200, bbox_inches = 'tight')
        filename = os.path.join(base_dir, "PREC-heatmap", f"LonLat_t={t}.png")
        plt.savefig(filename ,dpi=1200)
        plt.close()

    # ヒートマップの作成 緯度経度でなくグリッド
    ## 制御後の降水強度
    vmin = 5  # 最小値
    vmax = 35  # 最大値
    plt.figure(figsize=(10, 8))
    plt.grid(True)
    plt.imshow(TEST_prec_matrix.T, cmap='viridis', aspect='auto',origin='lower',vmin = vmin, vmax=vmax) #       カラーマップは好みに応じて変更可能
    plt.colorbar(label="precipitation (mm/h)")  # カラーバーのラベル
    plt.title(f"Time={Objective_T}h_X={Objective_X_low}-{Objective_X_high-1}_Y={Objective_Y_low}-{Objective_Y_high-1} ")
    plt.xlabel("X")
    plt.ylabel("Y")
    filename = os.path.join(base_dir, "PREC-heatmap", f"Grid.png")
    plt.savefig(filename ,dpi=300)

    ## 制御後- 制御前の降水強度
    plt.figure(figsize=(10, 8))
    plt.grid(True)
    plt.imshow(TEST_CTRL_prec_matrix.T, cmap='viridis', aspect='auto',origin='lower')  # カラーマップは好みに応じて変更可能
    plt.colorbar(label="precipitation (mm/h)")  # カラーバーのラベル
    plt.title(f"Time={Objective_T}h_X={Objective_X_low}-{Objective_X_high-1}_Y={Objective_Y_low}-{Objective_Y_high-1} ")
    plt.xlabel("X")
    plt.ylabel("Y")
    filename = os.path.join(base_dir, "PREC-heatmap", f"Grid_diff_C-noC.png")
    plt.savefig(filename ,dpi=300)
    return

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
    config_file_path = os.path.join(base_dir, filename)  # 修正ポイント

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


        for pe in range(nofpe):
            org = org_file.replace('######', str(pe).zfill(6))
            init = init_file.replace('######', str(pe).zfill(6))
            subprocess.run(["cp", org, init])

        TEST_prec_sum, CTRL_prec_sum, TEST_prec_matrix, CTRL_prec_matrix , TEST_CTRL_prec_matrix=sim(optimal_inputs[0],optimal_inputs[1],optimal_inputs[2],optimal_inputs[3])#  change


            # # 結果のクイック描画
            # subprocess.run(["mpirun", "-n", str(nofpe), "./sno", "sno.vgridope.d01.conf"])
            # subprocess.run(["mpirun", "-n", "1","./sno", "sno.hgridope.d01.conf"])
            # subprocess.run(["grads", "-blc", "checkfig_real.gs"])

    # シミュレーション結果
    print(f"CTRL_prec_matrix={CTRL_prec_matrix}")
    print(f"TEST_prec_matrix={TEST_prec_matrix}")

    print(f"{TEST_prec_sum=}")
    print(f"{CTRL_prec_sum=}")
    print(f"%={TEST_prec_sum/CTRL_prec_sum*100}%")

    plot_PREC(TEST_prec_matrix, TEST_CTRL_prec_matrix)
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
    main()

    webhook_url =os.getenv("SLACK_WEBHOOK_URL") # export SLACK_WEBHOOK_URL="OOOO"したらOK
    # 送信するメッセージを設定
    message = f"✅ {get_script_name()}の処理が完了しました。"
    notify_slack(webhook_url, message, channel="webhook")



