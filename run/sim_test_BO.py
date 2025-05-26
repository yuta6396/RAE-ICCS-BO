import os
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import subprocess

import cartopy.crs as ccrs
import gradsio2 as gio
import numpy.ma as ma
import matplotlib.ticker as mticker
import matplotlib.colors as mcolors
import imageio  # GIF作成用
from datetime import datetime

# 時刻を計測するライブラリ
import time
import pytz
from datetime import datetime
from zoneinfo import ZoneInfo

import warnings
# BO用
from skopt.space import Real,Integer
from skopt import Optimizer

import random
import requests
import json
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


# 関数評価回数
update_batch_times = 20
batch_size=5 #WSでは<=6
opt_num=batch_size*update_batch_times  

trial_num = 10#  試行回数
trial_base = 0
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
Control_Z_high = 2#MAX =36 36層存在 2の時429mが中心 4:1000m付近　2:248m 風車による制御効果が～200ｍ程度

# 介入領域の大きさ
Control_X_size = 5 
Control_Y_size = 5
Control_Z_size = Control_Z_high+1
#　介入幅
Control_Var = 0.3

optimal_inputs= [71,57]




#BOの獲得関数
Base_estimator="GP"
Acq_func="EI"

base_dir = f"../test_result/reduce_MOMXY/t=0-{Objective_T_period}_CVar={Control_Var}_{Control_X_size}-{Control_Y_size}-{Control_Z_size}grids_FET={opt_num}_trials={trial_base}-{trial_base+trial_num-1}_{current_time}"


def predict(inputs, f):
    """
    与えられた制御変数の値の配列(=input)を用いて、並列的に目的関数の値を評価し、配列の形で目的関数の値を返す。
    """
    #inputs= (batch_size,input_space_size)の構造
    
    global sub_history_file,sub_init_file

    for i in range(batch_size):
        print(f"predict_input{i}:{inputs[i]}\n")
        for pe in range(nofpe):
            sub_init_file = f"000{i}/init_d01_20070714-180000.000.pe######.nc"
            if i>=10:
                sub_init_file = f"00{i}/init_d01_20070714-180000.000.pe######.nc"
            sub_init = sub_init_file.replace('######', str(pe).zfill(6))
            init = init_file.replace('######', str(pe).zfill(6))
            subprocess.run(["cp", init, sub_init])

        init_val_intervation(i,inputs[i][0],inputs[i][1], f) # ここで初期値を書き換える

    # 目的関数の値の計算
    result = [0]*batch_size 
    for loop in range(0, Objective_T_period):
        # runの前に毎回histoyを消さないとERRORが起こる(history_fileの中身が消える)
        for i in range(batch_size):
            for pe in range(nofpe):
                sub_history_file = f"000{i}/history_d01.pe######.nc"
                sub_history = sub_history_file.replace('######', str(pe).zfill(6))
                history_path = file_path+'/'+sub_history
                if (os.path.isfile(history_path)):
                    subprocess.run(["rm", sub_history])
        # historyの時刻から runを回す run.d03.confで1h先まで
        result_mpi = subprocess.run(
            ["mpirun", "-n", str(nofpe*batch_size), "./scale-rm", "run.launch_prec_ave.conf"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(result_mpi.stdout.decode())
        print(result_mpi.stderr.decode())
        result_step = calc_objective_func(f, loop)
        for i in range(0, batch_size):
            result[i] += result_step[i]
    print(result)
    print(f"t=0-{Objective_T_period}:Sum of PREC in X[{Objective_X_low}, {Objective_X_high-1}], Y[{Objective_Y_low}, {Objective_Y_high-1}] = {result} [mm/h]\n")
    f.write(f"t=0-{Objective_T_period}:Sum of PREC in X[{Objective_X_low}, {Objective_X_high-1}], Y[{Objective_Y_low}, {Objective_Y_high-1}] = {result} [mm/h]\n")
    return result

def calc_objective_func(f, loop):
    result = [0]*batch_size # ステップごとの結果を初期化
    for i in range(batch_size):
        for pe in range(nofpe):
            sub_history_file = f"000{i}/history_d01.pe######.nc"
            fiy, fix = np.unravel_index(pe, (fny, fnx))
            nc = netCDF4.Dataset(sub_history_file.replace('######', str(pe).zfill(6)))
            nt = nc.dimensions['time'].size
            nx = nc.dimensions['x'].size
            ny = nc.dimensions['y'].size
            nz = nc.dimensions['z'].size
            gx1 = nx * fix
            gx2 = nx * (fix + 1)
            gy1 = ny * fiy
            gy2 = ny * (fiy + 1)
            if pe == 0:
                dat = np.zeros((nt, nz, fny*ny, fnx*nx)) # (シミュレーション時間, 36, 90, 90)
            dat[:, 0, gy1:gy2, gx1:gx2] = nc[varname][:]
        for j in range(Objective_X_low,Objective_X_high):
            for k in range(Objective_Y_low,Objective_Y_high):
                # for l in range(1, Objective_T_period+1):
                    result[i] += dat[1, 0, k, j]*TIME_INTERVAL
    print(f"t={loop}-{loop+1}:Sum of PREC in X[{Objective_X_low}, {Objective_X_high-1}], Y[{Objective_Y_low}, {Objective_Y_high-1}] = {result} [mm/h]\n")
    f.write(f"t={loop}-{loop+1}:Sum of PREC in X[{Objective_X_low}, {Objective_X_high-1}], Y[{Objective_Y_low}, {Objective_Y_high-1}] = {result} [mm/h]\n")
    return result

def init_val_intervation(num,LB_Control_X,LB_Control_Y, f):
    """
    与えられた制御変数の値（LB_Control_X~4）を用いて初期値を変更する。
    """
    
    global org_file
    f.write(f"LB_Control_X={LB_Control_X},LB_Control_Y={LB_Control_Y}\n")
    print(f"LB_Control_X={LB_Control_X},LB_Control_Y={LB_Control_Y}\n")
    for pe in range(nofpe):
        output_file = f"000{num}/outpe######.nc"
        sub_init_file = f"000{num}/init_d01_20070714-180000.000.pe######.nc"
        if num>=10:
            output_file = f"00{num}/out.pe######.nc"
            sub_init_file = f"00{num}/init_d01_20070714-180000.000.pe######.nc"
        sub_init = sub_init_file.replace('######', str(pe).zfill(6))
        output = output_file.replace('######', str(pe).zfill(6))

        with netCDF4.Dataset(sub_init) as src, netCDF4.Dataset(output, "w") as dst:
            dst.setncatts(src.__dict__)
            for name, dimension in src.dimensions.items():
                dst.createDimension(
                    name, (len(dimension) if not dimension.isunlimited() else None))
            for name, variable in src.variables.items():
                x = dst.createVariable(
                    name, variable.datatype, variable.dimensions)
                dst[name].setncatts(src[name].__dict__)
                if name == "MOMX":
                    var = src[name][:]
                    dst[name][:] = control_var_operation(var, pe, LB_Control_X, LB_Control_Y, Control_Y_size, Control_X_size, 0)
                elif name == "MOMY":
                    var = src[name][:]
                    dst[name][:] = control_var_operation(var, pe, LB_Control_X, LB_Control_Y, Control_Y_size, Control_X_size, 0)
                else:
                    dst[name][:] = src[name][:]
        subprocess.run(["cp", output, sub_init ])
    return

def control_var_operation(var, pe, LB_Control_X, LB_Control_Y, Control_Y_size, Control_X_size, B_Control_Z):
    if 0 <= LB_Control_X <= Control_X_low - Control_X_size and 0 <= LB_Control_Y <= Control_X_low - Control_X_size:
        if pe == 0:
            var[LB_Control_Y + 2 : LB_Control_Y + Control_Y_size + 2, LB_Control_X + 2 : LB_Control_X + Control_X_size + 2, B_Control_Z: B_Control_Z + Control_Z_size]*= Control_Var # (y, x, z) 3:5なら3, 4のみ
    elif  45 <= LB_Control_X <= Control_X_high- Control_X_size and 0 <= LB_Control_Y <= Control_X_low - Control_X_size: # 右下
        if pe == 1:
            var[LB_Control_Y + 2 : LB_Control_Y + Control_Y_size + 2, LB_Control_X - 45 : LB_Control_X + Control_X_size - 45, B_Control_Z: B_Control_Z + Control_Z_size] *= Control_Var
    elif   0 <= LB_Control_X <= Control_X_low - Control_X_size and 45 <= LB_Control_Y <= Control_X_high- Control_X_size: # 左上
        if pe == 2:
            var[LB_Control_Y - 45 : LB_Control_Y + Control_Y_size - 45, LB_Control_X + 2 : LB_Control_X + Control_X_size + 2, B_Control_Z: B_Control_Z + Control_Z_size] *= Control_Var
    elif   45 <= LB_Control_X <= Control_X_high- Control_X_size and 45 <= LB_Control_Y <= Control_X_high- Control_X_size: # 右上
        if pe == 3:
            var[LB_Control_Y - 45 : LB_Control_Y + Control_Y_size - 45, LB_Control_X -45 : LB_Control_X + Control_X_size -45, B_Control_Z: B_Control_Z + Control_Z_size] *= Control_Var

    elif 0 <= LB_Control_X <= Control_X_low - Control_X_size: # Y=41~44
        for Y_i in range(LB_Control_Y, 45): # pe == 0
            if pe == 0:
                var[Y_i + 2, LB_Control_X + 2 : LB_Control_X + Control_X_size + 2, B_Control_Z: B_Control_Z + Control_Z_size]*= Control_Var
        for Y_i in range(45, LB_Control_Y + Control_Y_size): # pe == 2
            if pe == 2:
                var[Y_i - 45, LB_Control_X + 2 : LB_Control_X + Control_X_size + 2, B_Control_Z: B_Control_Z + Control_Z_size]*= Control_Var

    elif 0 <= LB_Control_Y <= Control_X_low - Control_X_size: # X=41~44
        for X_i in range(LB_Control_X, 45): # pe == 0
            if pe == 0:
                var[LB_Control_Y + 2 : LB_Control_Y + Control_Y_size + 2, X_i + 2, B_Control_Z: B_Control_Z + Control_Z_size]*= Control_Var
        for X_i in range(45, LB_Control_X + Control_X_size): # pe == 1
            if pe == 1:
                var[LB_Control_Y + 2 : LB_Control_Y + Control_Y_size + 2, X_i - 45, B_Control_Z: B_Control_Z + Control_Z_size]*= Control_Var

    elif 45 <= LB_Control_X <= Control_X_high- Control_X_size: # Y=41~44
        for Y_i in range(LB_Control_Y, 45): # pe == 1
            if pe == 1:
                var[Y_i + 2, LB_Control_X - 45 : LB_Control_X + Control_X_size - 45, B_Control_Z: B_Control_Z + Control_Z_size]*= Control_Var
        for Y_i in range(45, LB_Control_Y + Control_Y_size): # pe == 3
            if pe == 3:
                var[Y_i - 45, LB_Control_X - 45 : LB_Control_X + Control_X_size - 45, B_Control_Z: B_Control_Z + Control_Z_size]*= Control_Var

    elif 45 <= LB_Control_Y <= Control_X_high- Control_X_size: # X=41~44
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

def sim(trial_i, LB_Control_X,LB_Control_Y, f):
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
                    dst[name][:] = control_var_operation(var, pe, LB_Control_X, LB_Control_Y, Control_Y_size, Control_X_size, 0)
                elif name == "MOMY":
                    var = src[name][:]
                    dst[name][:] = control_var_operation(var, pe, LB_Control_X, LB_Control_Y, Control_Y_size, Control_X_size, 0)
                else:
                    dst[name][:] = src[name][:]
        subprocess.run(["cp", output, init])

    for loop in range(0, Objective_T_period):
        TEST_prec_sum, CTRL_prec_sum, TEST_prec_matrix, CTRL_prec_matrix, TEST_CTRL_prec_matrix = state_update(trial_i, loop, TEST_prec_matrix, CTRL_prec_matrix, TEST_CTRL_prec_matrix, TEST_prec_sum, CTRL_prec_sum)

    return TEST_prec_sum, CTRL_prec_sum, TEST_prec_matrix, CTRL_prec_matrix ,TEST_CTRL_prec_matrix
   
def state_update(trial_i, loop, TEST_prec_matrix, CTRL_prec_matrix, TEST_CTRL_prec_matrix, TEST_prec_sum, CTRL_prec_sum):
    orgfile = f"no-control_t={loop+11}.pe######.nc"
    gpyoptfile=f"RS_remake_Init_MOMXY-pe=1,3_FET={opt_num}_t={loop}_seed={trial_i}.pe######.nc" #history
    gpyoptfile=f"{base_dir}/t={loop}_seed={trial_i}.pe######.nc" 

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
            TEST_prec_matrix[loop, j-Objective_X_low, k-Objective_Y_low] += dat[1,0,k,j]*TIME_INTERVAL #(y,x,z)=(time,z,y,x)　j-Objective_X_low,k-Objective_Y_lowは[0,0]->[5,5]とか
            CTRL_prec_matrix[loop,j-Objective_X_low, k-Objective_Y_low] += odat[1,0,k,j]*TIME_INTERVAL
            TEST_CTRL_prec_matrix[loop, j-Objective_X_low, k-Objective_Y_low] = TEST_prec_matrix[loop, j-Objective_X_low, k-Objective_Y_low] - CTRL_prec_matrix[loop, j-Objective_X_low, k-Objective_Y_low]
            TEST_prec_sum+=TEST_prec_matrix[loop, j-Objective_X_low, k-Objective_Y_low]
            CTRL_prec_sum+=CTRL_prec_matrix[loop, j-Objective_X_low, k-Objective_Y_low]
    print(f"{TEST_prec_sum=}")
    print(f"{CTRL_prec_sum=}")
    return TEST_prec_sum, CTRL_prec_sum, TEST_prec_matrix, CTRL_prec_matrix ,TEST_CTRL_prec_matrix



def plot_PREC(trial_i, TEST_prec_matrix, TEST_CTRL_prec_matrix):
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
    plt_extent_sub = [139, 143, 36, 40]     # 詳細表示用の範囲
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
        history_name      = f"{base_dir}/t={loop}_seed={trial_i}"
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
    plt.savefig(os.path.join(out_dir, f"LonLat_seed={trial_i}.pdf"), dpi=1200)
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
    plt.savefig(os.path.join(out_dir, f"Sub_LonLat_seed={trial_i}.pdf"), dpi=1200, bbox_inches = 'tight')
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
    plt.savefig(os.path.join(out_dir, f"Loop_by_loop_seed={trial_i}.png"), dpi=1200)
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

def convert_np_types(obj):
    # BO_search_data 内の NumPy の数値を、通常の Python の int に変換
    if isinstance(obj, dict):
        return {k: convert_np_types(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_np_types(i) for i in obj]
    elif isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    else:
        return obj

### 実行
def main():
    make_directory(base_dir)
    filename = f"config.txt"
    config_file_path = os.path.join(base_dir, filename)  # 修正ポイント

    # ファイルに書き込む
    with open(config_file_path,  "w") as f:
        # 関数評価回数
        f.write(f"update_batch_times = {update_batch_times}\n")
        f.write(f"batch_size = {batch_size}\n")
        f.write(f"opt_num = {opt_num}  # WSでは<=24\n")
        
        # 試行回数
        f.write(f"\ntrial_num = {trial_num}  # 試行回数\n")
        f.write(f"\n{trial_base=}  # seedの最初の値\n")
        
        # 降水強度を最小化したい領域
        f.write("\n# 降水強度を最小化したい領域\n")
        f.write(f"Objective_X_low = {Objective_X_low}\n")
        f.write(f"Objective_X_high = {Objective_X_high}\n")
        f.write(f"Objective_Y_low = {Objective_Y_low}\n")
        f.write(f"Objective_Y_high = {Objective_Y_high}\n")
        f.write(f"Area_size = {Area_size}\n")
        
        # 降水強度を最小化したい時刻
        f.write("\n# 降水強度を最小化したい時間\n")
        f.write(f"{Objective_T_period=}\n")
        
        # 制御対象範囲
        f.write("\n# 制御対象範囲\n")
        f.write(f"Control_X_low = {Control_X_low}\n")
        f.write(f"Control_X_high = {Control_X_high}\n")
        f.write(f"Control_Y_low = {Control_Y_low}\n")
        f.write(f"Control_Y_high = {Control_Y_high}  # MAX =90\n")
        f.write(f"Control_Z_high = {Control_Z_high}  # MAX =35?\n")
        
        # 介入領域の大きさ
        f.write("\n# 介入領域の大きさ\n")
        f.write(f"Control_X_size = {Control_X_size}\n")
        f.write(f"Control_Y_size = {Control_Y_size}\n")
        f.write(f"Control_Z_size = {Control_Z_size}\n")
        
        # BOの獲得関数
        f.write("\n# BOの獲得関数\n")
        f.write(f"Base_estimator = \"{Base_estimator}\"\n")
        f.write(f"Acq_func = \"{Acq_func}\"\n")

    log_file = os.path.join(base_dir, "summary", f"BO_log.txt")
    summary_file = os.path.join(base_dir, "summary", f"BO_summary.txt")
    analysis_file = os.path.join(base_dir, "summary", f"BO_analysis.txt")

    # 入力次元と最小値・最大値の定義
    bounds = [ Integer(Control_X_low, Control_X_high-Control_X_size),Integer(Control_Y_low, Control_Y_high-Control_Y_size)] # Int(2, 4) なら2, 3, 4からランダム

    with open(log_file, 'w') as f, open(summary_file, 'w') as f_s, open(analysis_file, 'w') as f_a:
        BO_data = []
        BO_search_data = []
        for trial_i in range(trial_base, trial_base+trial_num):
            f.write(f"{trial_i=}\n")
            f_a.write(f"\n{trial_i=}\n")
            trial_best_values = []
            trial_best_points = []
            random.seed(trial_i)

            start = time.time()  # 現在時刻（処理開始前）を取得
            opt = Optimizer(bounds, base_estimator=Base_estimator, acq_func=Acq_func, random_state=trial_i)

            for pe in range(nofpe):
                org = org_file.replace('######', str(pe).zfill(6))
                init = init_file.replace('######', str(pe).zfill(6))
                subprocess.run(["cp", org, init])

            best_values = []
            current_best = float('inf')

            for i in range(2):  # 最初の10点はRS
                # 各範囲でランダム値を生成
                random_samples1 = [[random.randint(Control_X_low, Control_X_high-Control_X_size),] for _ in range(batch_size)]
                random_samples2 = [[random.randint(Control_Y_low, Control_Y_high-Control_Y_size) ] for _ in range(batch_size)]

                # リストを結合してサンプルを作成
                combined_samples = [sample1 + sample2  for sample1, sample2 in zip(random_samples1, random_samples2)]

                Y = predict(combined_samples, f)
                opt.tell(combined_samples,Y)
                f_a.write(f"Batch {i}: Best value so far: {min(opt.yi)}\n")
                trial_best_values.append(min(opt.yi))
                best_index = np.argmin(opt.yi)
                trial_best_points.append(opt.Xi[best_index])

            for j in range(update_batch_times-2):
                # アクイジション関数を用いて次の探索点を決定
                next_points = opt.ask(n_points=batch_size)

                # 並列で評価関数を計算
                values = predict(next_points, f)
                print(f"values{values}")

                # 評価結果をモデルに反映
                opt.tell(next_points, values)
                best_index = np.argmin(opt.yi)
                print(f"Batch {j+2}: Best value so far: {min(opt.yi)}:best_X={opt.Xi[best_index]}\n")
                f.write(f"Batch {j+2}: Best value so far: {min(opt.yi)}\n\n")
                f_a.write(f"Batch {j+2}: Best value so far: {min(opt.yi)}:best_X={opt.Xi[best_index]}\n")
                trial_best_values.append(min(opt.yi))
                trial_best_points.append(opt.Xi[best_index])

            # 結果の取得
            BO_data.append(trial_best_values)
            BO_search_data.append(trial_best_points)
            best_value = min(opt.yi)
            # 最小値のインデックスを取得
            min_index = opt.yi.index(min(opt.yi))

            # 対応するベストポイントを取得
            best_point = opt.Xi[min_index]
            # best_point = opt.Xi[np.argmin(opt.yi)]
            print(f"Best value: {best_value} at point {best_point}")

            optimal_inputs = best_point
            # ここまでコメントアウトして下のoptimalで探索地点手入力
            # optimal_inputs=[70, 50, 20, 0]

            # 得られた最適解でもう一度シミュレーションする
            for pe in range(nofpe):
                org = org_file.replace('######', str(pe).zfill(6))
                init = init_file.replace('######', str(pe).zfill(6))
                subprocess.run(["cp", org, init])

            TEST_prec_sum, CTRL_prec_sum, TEST_prec_matrix, CTRL_prec_matrix , TEST_CTRL_prec_matrix=sim(trial_i, optimal_inputs[0],optimal_inputs[1], f)#  change

            end=time.time()
            time_diff = end - start
            print(f'\n\n実行時間:{time_diff}\n')

            # シミュレーション結果
            print(f"{TEST_prec_sum=}")
            print(f"{CTRL_prec_sum=}")
            print(f"%={TEST_prec_sum/CTRL_prec_sum*100}%")
            f_s.write(f'{trial_i=}')

            f_s.write(f'実行時間:{time_diff}\n')
            f_s.write(f'{opt_num}回の評価で得られた最適解：{optimal_inputs}\n')
            f_s.write(f"CTRL_prec_matrix={CTRL_prec_matrix}\n")
            f_s.write(f"TEST_prec_matrix={TEST_prec_matrix}\n")

            f_s.write(f"{TEST_prec_sum=}")
            f_s.write(f"{CTRL_prec_sum=}")
            f_s.write(f"%={TEST_prec_sum/CTRL_prec_sum*100}%\n\n")
            plot_PREC(trial_i, TEST_prec_matrix, TEST_CTRL_prec_matrix)

    print(base_dir)
    with open(f'{base_dir}/summary/BO_data.json', 'w') as f, open(f'{base_dir}/summary/BO_search_data.json', 'w') as f_s:
        BO_search_data_converted = convert_np_types(BO_search_data)
        json.dump(BO_search_data_converted, f_s, indent=4)
        json.dump(BO_data, f, indent=4)

def notify_slack(webhook_url, message, channel=None, username=None, icon_emoji=None):
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