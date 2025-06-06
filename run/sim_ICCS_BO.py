import os
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import subprocess

import cartopy.crs as ccrs
import gradsio2 as gio
import matplotlib.ticker as mticker
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
update_batch_times = 60
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
Control_Z_high = 4#MAX =36 36層存在 2の時429mが中心

# 介入領域の大きさ
Control_X_size = 5 #あんま変えられない　いくつか同時に変更する地点ありrandom_samples1 とか
Control_Y_size = 5
Control_Z_size = Control_Z_high+1
#　介入幅
Control_Var_low = -20   #30~30でエラー確認済み(MOMY)
Control_Var_high = 20




#BOの獲得関数
"""
ガウス過程：GaussianProcessRegressor – カーネル例：Matern(nu=1.5) / RBF(length_scale=1.0) / RationalQuadratic	・連続変数の探索に最適・不確かさ（予測分散）を自然に扱える
ランダムフォレスト	RandomForestRegressor(n_estimators=100)	・離散変数やカテゴリ変数を含む場合に強い・ノイズ耐性が高い
Extra-Trees	ExtraTreesRegressor(n_estimators=100)	・ランダムフォレストよりもランダム性を強めたい場合
勾配ブースティング	HistGradientBoostingRegressor()	・非常に大規模なデータや高次元の探索にも対応しやすい
カスタム回帰器	scikit-learn 準拠の任意の Regressor	・ニューラルネットワークや他ライブラリのモデルも可（要 predict_std 実装）

"""
Base_estimator="RF" #RF', 'ET', 'GP', 'GBRT' or 'DUMMY'
"""
"EI"	Expected Improvement
（期待改善幅）	xi（改善探索のバランス調整）
"PI"	Probability of Improvement
（改善確率）	xi（確率探索の閾値）
"LCB"	(Lower) Confidence Bound
（下側信頼境界）	kappa（探索度合いの重み）
"gp_hedge"	複数の取得関数を動的に組み合わせる手法	内部で EI/PI/LCB を自動選択
"""
Acq_func="gp_hedge"
num_rand =2 # *batch_sizeだけRS

base_dir = f"../test_result/Acq_MOMXY/{Base_estimator}_{Acq_func}_RS={num_rand*batch_size}_t=0-{Objective_T_period}_{Control_X_size}-{Control_Y_size}-{Control_Z_size}grids_FET={opt_num}_trials={trial_base}-{trial_base+trial_num-1}_{current_time}"


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

        init_val_intervation(i,inputs[i][0],inputs[i][1],inputs[i][2], inputs[i][3], f) # ここで初期値を書き換える

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

def init_val_intervation(num,LB_Control_X,LB_Control_Y, Control_MOMX, Control_MOMY, f):
    """
    与えられた制御変数の値（LB_Control_X~4）を用いて初期値を変更する。
    """
    
    global org_file
    f.write(f"LB_Control_X={LB_Control_X},LB_Control_Y={LB_Control_Y},{ Control_MOMX=},{Control_MOMY=}\n")
    print(f"LB_Control_X={LB_Control_X},LB_Control_Y={LB_Control_Y},{ Control_MOMX=},{Control_MOMY=}\n")
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
                    dst[name][:] = control_var_operation(var, pe, LB_Control_X, LB_Control_Y, Control_Y_size, Control_X_size, 0, Control_MOMX)
                elif name == "MOMY":
                    var = src[name][:]
                    dst[name][:] = control_var_operation(var, pe, LB_Control_X, LB_Control_Y, Control_Y_size, Control_X_size, 0, Control_MOMY)
                else:
                    dst[name][:] = src[name][:]
        subprocess.run(["cp", output, sub_init ])
    return

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

def sim(trial_i, LB_Control_X,LB_Control_Y,Control_MOMX, Control_MOMY):
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
                    dst[name][:] = control_var_operation(var, pe, LB_Control_X, LB_Control_Y, Control_Y_size, Control_X_size, 0, Control_MOMX)
                elif name == "MOMY":
                    var = src[name][:]
                    dst[name][:] = control_var_operation(var, pe, LB_Control_X, LB_Control_Y, Control_Y_size, Control_X_size, 0, Control_MOMY)
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
    得られた最適解（制御入力）でシミュレーションをした結果を、"可視化する"関数。
    目的時刻 までのの降水強度の総和を可視化する
    """
    outtype = 'fcst'                    # outtype in convert_letkfout.py
    ncsuf = 'pe______.nc'

    # params for plot
    cmap = plt.cm.jet
    cmap.set_under('lightgray')
    m_levels = [0,10,20,30,40]
    s_levels = [25, 50, 75, 100, 125, 150, 175] 
    plt_extent = [126,144,26,44]

    for loop in range(0, Objective_T_period) :
        history_name = f"{base_dir}/t={loop}_seed={trial_i}" 
        ncfilebasetmp = '{:s}.{:}'.format(history_name, ncsuf)
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
        var *=  TIME_INTERVAL
        if loop == 0:
                # get lon lat
            slon = gio.readvar2d('lon',ncfilebasetmp,outtype,fnx,fny)
            slat = gio.readvar2d('lat',ncfilebasetmp,outtype,fnx,fny)  
            prec_sum_matrix = var[1]
        else:
            prec_sum_matrix += var[1]

    # plot
    fig = plt.figure()

    # scale
    ax = fig.add_subplot(1,1,1,projection=ccrs.PlateCarree())
    ax.set_extent(plt_extent)
    ax.coastlines(resolution='50m',linewidth=0.5)

    # grid line
    gl = ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=True,linewidth=0.5,linestyle='--',color='gray')
    gl.xlocator = mticker.FixedLocator(range(126,146)) #146
    gl.ylocator = mticker.FixedLocator(range(26,44))
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size':6}
    gl.ylabel_style = {'size':6}

    im = ax.contourf(slon,slat,prec_sum_matrix,cmap=cmap,levels=s_levels, extend='both',transform=ccrs.PlateCarree())

    pos = ax.get_position()
    fig.subplots_adjust(right=0.85)
    cax = fig.add_axes([pos.x1,pos.y0,0.02,pos.y1-pos.y0])
    clb = fig.colorbar(im, orientation='vertical', cax=cax)

    filename = os.path.join(base_dir, "PREC-heatmap", f"LonLat_t=0-{Objective_T_period}_seed={trial_i}.png")
    plt.savefig(filename ,dpi=300)
    plt.close()

    # ヒートマップの作成 緯度経度でなくグリッド
    ## 制御後の降水強度
    vmin = 0  # 最小値
    vmax = 160  # 最大値
    TEST_prec_sum_matrix = TEST_prec_matrix[0]
    for t in range(1, Objective_T_period):
        TEST_prec_sum_matrix += TEST_prec_matrix[t]

    plt.figure(figsize=(10, 8))
    plt.grid(True)
    plt.imshow(TEST_prec_sum_matrix.T, cmap='viridis', aspect='auto',origin='lower',vmin = vmin, vmax=vmax) #       カラーマップは好みに応じて変更可能
    plt.colorbar(label="Accumulated PREC (mm)")  # カラーバーのラベル
    plt.title(f"Time=0-{Objective_T_period}h_X={Objective_X_low}-{Objective_X_high-1}_Y={Objective_Y_low}-{Objective_Y_high-1} ")
    plt.xlabel("X")
    plt.ylabel("Y")
    filename = os.path.join(base_dir, "PREC-heatmap", f"Grid_seed={trial_i}.png")
    plt.savefig(filename ,dpi=300)

    ## 制御後- 制御前の降水強度
    TEST_CTRL_prec_sum_matrix = TEST_CTRL_prec_matrix[0]
    for t in range(1, Objective_T_period):
        TEST_CTRL_prec_sum_matrix += TEST_CTRL_prec_matrix[t]

    plt.figure(figsize=(10, 8))
    plt.grid(True)
    plt.imshow(TEST_CTRL_prec_sum_matrix.T, cmap='viridis', aspect='auto',origin='lower') #       カラーマップは好みに応じて変更可能
    plt.colorbar(label="Accumulated PREC (mm)")  # カラーバーのラベル
    plt.title(f"Time=0-{Objective_T_period}h_X={Objective_X_low}-{Objective_X_high-1}_Y={Objective_Y_low}-{Objective_Y_high-1} ")
    plt.xlabel("X")
    plt.ylabel("Y")
    filename = os.path.join(base_dir, "PREC-heatmap", f"Grid_diff_C-noC_seed={trial_i}.png")
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
        # 介入幅
        f.write("\n# 介入幅\n")
        f.write(f"Control_Var_low = {Control_Var_low}  # -30~30でエラー確認済み\n")
        f.write(f"Control_Var_high = {Control_Var_high}\n")
        
        # BOの獲得関数
        f.write("\n# BOの獲得関数\n")
        f.write(f"Base_estimator = \"{Base_estimator}\"\n")
        f.write(f"Acq_func = \"{Acq_func}\"\n")

    log_file = os.path.join(base_dir, "summary", f"BO_log.txt")
    summary_file = os.path.join(base_dir, "summary", f"BO_summary.txt")
    analysis_file = os.path.join(base_dir, "summary", f"BO_analysis.txt")

    # 入力次元と最小値・最大値の定義
    bounds = [ Integer(Control_X_low, Control_X_high-Control_X_size),Integer(Control_Y_low, Control_Y_high-Control_Y_size) ,Real(Control_Var_low,Control_Var_high), Real(Control_Var_low,Control_Var_high)] # Int(2, 4) なら2, 3, 4からランダム

    with open(log_file, 'w') as f, open(summary_file, 'w') as f_s, open(analysis_file, 'w') as f_a:
        BO_data = []
        for trial_i in range(trial_base, trial_base+trial_num):
            f.write(f"{trial_i=}\n")
            f_a.write(f"\n{trial_i=}\n")
            trial_best_values = []
            random.seed(trial_i)

            start = time.time()  # 現在時刻（処理開始前）を取得
            opt = Optimizer(bounds, base_estimator=Base_estimator, acq_func=Acq_func, random_state=trial_i)

            for pe in range(nofpe):
                org = org_file.replace('######', str(pe).zfill(6))
                init = init_file.replace('######', str(pe).zfill(6))
                subprocess.run(["cp", org, init])

            best_values = []
            current_best = float('inf')


            for i in range(num_rand):  # 最初の10点はRS
                # 各範囲でランダム値を生成
                random_samples1 = [[random.randint(Control_X_low, Control_X_high-Control_X_size),] for _ in range(batch_size)]
                random_samples2 = [[random.randint(Control_Y_low, Control_Y_high-Control_Y_size) ] for _ in range(batch_size)]
                random_samples3 = [[random.uniform(Control_Var_low,Control_Var_high)] for _ in range(batch_size)]
                random_samples4 = [[random.uniform(Control_Var_low,Control_Var_high)] for _ in range(batch_size)]

                # リストを結合してサンプルを作成
                combined_samples = [sample1 + sample2 + sample3 + sample4 for sample1, sample2, sample3, sample4 in zip(random_samples1, random_samples2, random_samples3, random_samples4)]

                Y = predict(combined_samples, f)
                opt.tell(combined_samples,Y)
                f_a.write(f"Batch {i}: Best value so far: {min(opt.yi)}\n")
                trial_best_values.append(min(opt.yi))

            for j in range(update_batch_times-num_rand):
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

            # 結果の取得
            BO_data.append(trial_best_values)
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

            TEST_prec_sum, CTRL_prec_sum, TEST_prec_matrix, CTRL_prec_matrix , TEST_CTRL_prec_matrix=sim(trial_i, optimal_inputs[0],optimal_inputs[1],optimal_inputs[2], optimal_inputs[3])#  change

            end=time.time()
            time_diff = end - start
            print(f'\n\n実行時間:{time_diff}\n')


            # 結果のクイック描画
            subprocess.run(["mpirun", "-n", str(nofpe), "./sno", "sno.vgridope.d01.conf"])
            subprocess.run(["mpirun", "-n", "1","./sno", "sno.hgridope.d01.conf"])
            subprocess.run(["grads", "-blc", "checkfig_real.gs"])

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
    with open(f'{base_dir}/summary/BO_data.json', 'w') as f:
        json.dump(BO_data, f, indent=4)

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