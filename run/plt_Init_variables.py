import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import gradsio2 as gio
import matplotlib.ticker as mticker
import matplotlib.colors as mcolors
import cartopy.feature as cfeature
"""
.nc file list
"../init/init_d01_20070714-180000.000.pe00000{pe}.nc" : init_file 
"history_d01"
"""
Z_low = 0#どの高さにするか
Z_high = 36
name = "init_file"
var = "QC"
dpi = 300
for Z in range(Z_low, Z_high):
    for pe in range(4):
        # NetCDFファイルのパスを指定
        nc_file = f"../init/init_d01_20070714-180000.000.pe00000{pe}.nc"  # ここを実際のファイル名に変更してください
        # NetCDFファイルを読み込み
        if pe == 0:
            dataset0 = netCDF4.Dataset(nc_file, "r")
        elif pe == 1:
            dataset1 = netCDF4.Dataset(nc_file, "r")
        elif pe == 2:
            dataset2 = netCDF4.Dataset(nc_file, "r")
        else:
            dataset3 = netCDF4.Dataset(nc_file, "r")

    # これらを辞書にまとめると扱いやすいです。
    datasets = {
        0: dataset0,
        1: dataset1,
        2: dataset2,
        3: dataset3
    }
    # それぞれの変数について、pe ごとに45×45の部分を切り出す
    # ※ 元のデータサイズが (47,47) と仮定し、インデックス [2:, 2:] 等で切り出します

    # --- var ---
    Y_data0 = datasets[0].variables[var][2:, 2:,Z]
    Y_data1 = datasets[1].variables[var][2:, 0:-2,Z]
    Y_data2 = datasets[2].variables[var][0:-2, 2:,Z]
    Y_data3 = datasets[3].variables[var][0:-2, 0:-2,Z]

    # --- lon ---
    lon_data0 = datasets[0].variables["lon"][2:, 2:]
    lon_data1 = datasets[1].variables["lon"][2:, 0:-2]
    lon_data2 = datasets[2].variables["lon"][0:-2, 2:]
    lon_data3 = datasets[3].variables["lon"][0:-2, 0:-2]

    # --- lat ---
    lat_data0 = datasets[0].variables["lat"][2:, 2:]
    lat_data1 = datasets[1].variables["lat"][2:, 0:-2]
    lat_data2 = datasets[2].variables["lat"][0:-2, 2:]
    lat_data3 = datasets[3].variables["lat"][0:-2, 0:-2]

    # ここで、各切り出し済み配列のサイズは (45,45) となっているはずです

    # 配置順については、例として以下のように配置します：
    # 上段：pe=0 (右下部分) を左、pe=1 (左下部分) を右
    # 下段：pe=2 (右上部分) を左、pe=3 (左上部分) を右
    #
    # ※ 注意：配置順は用途に合わせて変更してください

    # var の連結
    upper_Y = np.hstack((Y_data0, Y_data1))
    lower_Y = np.hstack((Y_data2, Y_data3))
    combined_Var = np.vstack((upper_Y, lower_Y))

    # lon の連結
    upper_lon = np.hstack((lon_data0, lon_data1))
    lower_lon = np.hstack((lon_data2, lon_data3))
    combined_lon = np.vstack((upper_lon, lower_lon))

    # lat の連結
    upper_lat = np.hstack((lat_data0, lat_data1))
    lower_lat = np.hstack((lat_data2, lat_data3))
    combined_lat = np.vstack((upper_lat, lower_lat))

    # 各連結後のサイズを確認
    print(f"combined_{var}.shape:", combined_Var.shape)  # (90, 90)
    print("combined_lon.shape: ", combined_lon.shape)   # (90, 90)
    print("combined_lat.shape: ", combined_lat.shape)   # (90, 90)



    # Var_sub の最大値を取得
    max_Var = np.max(combined_Var)
    min_Var = np.min(combined_Var)
    print(f"{var} の最大値:", max_Var)
    print(f"{var} の最小値:", min_Var)

    # プロット作成
    fig = plt.figure(figsize=(12, 8))
    plt_extent = [126,144,26,43]
    cmap = plt.cm.jet
    cmap.set_under('lightgray')
    norm = mcolors.Normalize(vmin=min_Var, vmax=max_Var)
    # カラープロット
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_extent(plt_extent)
    ax.coastlines(resolution='50m', linewidth=0.5)
    #ax.tick_params(axis='both', which='major', labelsize=9)
    ax.set_title(f"{var} Z={Z}", fontsize=20, fontweight='bold', loc='center')
    # グリッド線
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, linestyle='--', color='gray')
    gl.xlocator = mticker.FixedLocator(range(125, 146))
    gl.ylocator = mticker.FixedLocator(range(26, 43)) # 44
    #gl.ylocator = mticker.FixedLocator(np.arange(26, 43 + 1, 3))    # 26°から46°まで3度ごと
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 15}
    gl.ylabel_style = {'size': 15}

    # `pcolormesh`で連続的な色分布を描画
    im = ax.pcolormesh(combined_lon, combined_lat, combined_Var, cmap=cmap, norm=norm,shading='auto', transform=ccrs.PlateCarree())

    # カラーバーを追加
    cbar = fig.colorbar(im, ax=ax, orientation='vertical', pad=0.05, aspect=30)
    cbar.set_label('kg/kg', fontsize=18)
    cbar.ax.tick_params(labelsize=15)
    plt.savefig(f"/home/yuta/scale-5.5.3/scale-rm/test/tutorial/real/data_values/{var}/{name}_Z={Z}.png", dpi=dpi)
    plt.close()
    # ファイルクローズ
    dataset0.close()
    dataset1.close()
    dataset2.close()
    dataset3.close()