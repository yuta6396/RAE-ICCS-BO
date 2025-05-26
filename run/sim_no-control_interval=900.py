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
from skopt import gp_minimize
from skopt.space import Real,Integer
import warnings
matplotlib.use('Agg')

# 警告を抑制
warnings.filterwarnings('ignore', category=DeprecationWarning)


nofpe = 4
fny = 2
fnx = 2

interval=900
num_t = 24*int(3600/interval)


varname = 'PREC'

init_file = "../init/init_d01_20070714-180000.000.pe######.nc"
org_file = "../init/init_d01_20070714-180000.000.pe######.org.nc"
history_file = "history_interval=900.pe######.nc"

no_control = 'no-control_t=0.pe######.nc'
file_path = '/home/yuta/scale-5.5.3/scale-rm/test/tutorial/real/experiment_init_BO'
gpyoptfile=f"no-control_t=0.pe######.nc"

sum_gpy=np.zeros((90,90)) #ベイズ最適化の累積降水量
sum_no=np.zeros((num_t,90,90)) #制御前の累積降水量
# sum_prec=np.zeros((6,90,90))

cnt=0
# opt_num=100

for pe in range(nofpe):
    init = init_file.replace('######', str(pe).zfill(6))
    org = org_file.replace('######', str(pe).zfill(6))
    subprocess.run(["cp", org, init]) #初期化

for i in range(num_t):
    
    subprocess.run(["mpirun", "-n", "4", "./scale-rm", "run_interval=900.conf"])
    for pe in range(nofpe):
        no_control=f"no-control_interval={interval}_t={i}.pe######.nc"
        no_control = no_control.replace('######', str(pe).zfill(6))
        history = history_file.replace('######', str(pe).zfill(6))
        subprocess.run(["cp", history,no_control])

        restart=f"restart_interval={interval}_t={i}.pe######.nc"
        restart= restart.replace('######', str(pe).zfill(6))
        init = init_file.replace('######', str(pe).zfill(6))
        subprocess.run(["cp", init, restart])
    print(f"loop={i}")

gpy=np.zeros(num_t)
no=np.zeros(num_t)

outtype = 'fcst'                    # outtype in convert_letkfout.py
member = '0046'                     # member in convert_letkfout.py
ncsuf = 'pe______.nc'
# params for plot
cmap = plt.cm.jet
cmap.set_under('lightgray')
m_levels = [0,10,20,30,40]
s_levels = [50,100,150,200] 
plt_extent = [126,144,26,44]

for i in range(num_t):
    for pe in range(nofpe):  # history処理
            fiy, fix = np.unravel_index(pe, (fny, fnx))
            
            no_control=f"no-control_interval={interval}_t={i}.pe######.nc"
            no_control= no_control.replace('######', str(pe).zfill(6))

            nc = netCDF4.Dataset(no_control.replace('######', str(pe).zfill(6)))
            
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
                # odat = np.zeros((nt, nz, fny*ny, fnx*nx))
            dat[:, 0, gy1:gy2, gx1:gx2] = nc[varname][:]
            for j in range(0,90):
                for k in range(0,90):
                    # sum_gpy[j,k]+=dat[5,0,k,j]*3600 #(y,x,z)=(time,z,y,x)
                    
                    sum_no[i,j,k]+=dat[1,0,k,j]*3600
                    # sum+=sum_gpy[j,k]
                    # no+=sum_no[j,k]
    
        

for i in range(num_t):
    history_name = f'no-control_interval={interval}_t={i}'
    ncfilebasetmp = '{:s}.{:}'.format(history_name, ncsuf)

    # get lon lat
    slon = gio.readvar2d('lon',ncfilebasetmp,outtype,fnx,fny)
    slat = gio.readvar2d('lat',ncfilebasetmp,outtype,fnx,fny)
    var = gio.readvar2d('PREC',ncfilebasetmp,outtype,fnx,fny)
    # unit conversion (kg/m/s --> mm/d)
    var *=  60 * 60


    # plot
    # base
    fig = plt.figure()

    # scale
    ax = fig.add_subplot(1,1,1,projection=ccrs.PlateCarree())
    ax.set_extent(plt_extent)
    ax.coastlines(resolution='50m',linewidth=0.5)

    # grid line
    gl = ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=True,linewidth=0.5,linestyle='--',color='gray')
    gl.xlocator = mticker.FixedLocator(range(126,146))
    gl.ylocator = mticker.FixedLocator(range(26,44))
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size':6}
    gl.ylabel_style = {'size':6}

    im = ax.contourf(slon,slat,var[1],cmap=cmap,levels=m_levels, extend='both',transform=ccrs.PlateCarree())

    pos = ax.get_position()
    fig.subplots_adjust(right=0.85)
    cax = fig.add_axes([pos.x1,pos.y0,0.02,pos.y1-pos.y0])
    clb = fig.colorbar(im, orientation='vertical', cax=cax)

    plt.savefig(f'result_noControl/no-control_interval={interval}_t={i}.png',dpi=300)
    plt.close()




    # ヒートマップの作成
    plt.figure(figsize=(10, 8))
    plt.grid(True)
    plt.imshow(sum_no[i,:,:].T, cmap='viridis', aspect='auto',origin='lower')  # カラーマップは好みに応じて変更可能
    plt.colorbar(label="precipitation (mm/h)")  # カラーバーのラベル
    plt.title("precipitation")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.savefig(f"result_noControl/PREC_visualize_no-control_interval={interval}_t={i}.png", dpi= 300)
    plt.close()



