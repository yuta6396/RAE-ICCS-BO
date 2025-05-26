import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
"""
.nc file list
"../init/init_d01_20070714-180000.000.pe00000{pe}.nc" : init_file 
"history_d01"
"""
Z_low= 0 #どの高さにするか
Z_high = 0
name = "init_file"

# 71〜75, 57〜61 の部分（5×5 = 25個）の始点座標を取得
L_X = 71
B_Y = 57
d_X = 20
d_Y = -20
Input_X_size = 5

# グリッドが密すぎる場合は、サブサンプリング（例：5点ごとに表示）
step = 1
for Z in range(Z_low, Z_high+1):
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

    # --- MOMX ---
    X_data0 = datasets[0].variables["MOMX"][2:, 2:,Z]    # pe=0: 右下部分
    X_data1 = datasets[1].variables["MOMX"][2:, 0:-2,Z]   # pe=1: 左下部分
    X_data2 = datasets[2].variables["MOMX"][0:-2, 2:,Z]   # pe=2: 右上部分
    X_data3 = datasets[3].variables["MOMX"][0:-2, 0:-2, Z]  # pe=3: 左上部分

    # --- MOMY ---
    Y_data0 = datasets[0].variables["MOMY"][2:, 2:,Z]
    Y_data1 = datasets[1].variables["MOMY"][2:, 0:-2,Z]
    Y_data2 = datasets[2].variables["MOMY"][0:-2, 2:,Z]
    Y_data3 = datasets[3].variables["MOMY"][0:-2, 0:-2,Z]

    # # MOMX・MOMYの値を介入後の値に変更する
    # X_data3[B_Y-45:B_Y-45+Input_X_size,L_X-45:L_X-45+Input_X_size] += d_X
    # Y_data3[B_Y-45:B_Y-45+Input_X_size,L_X-45:L_X-45+Input_X_size] += d_Y
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

    # MOMX の連結
    upper_X = np.hstack((X_data0, X_data1))
    lower_X = np.hstack((X_data2, X_data3))
    combined_MOMX = np.vstack((upper_X, lower_X))

    # MOMY の連結
    upper_Y = np.hstack((Y_data0, Y_data1))
    lower_Y = np.hstack((Y_data2, Y_data3))
    combined_MOMY = np.vstack((upper_Y, lower_Y))

    # lon の連結
    upper_lon = np.hstack((lon_data0, lon_data1))
    lower_lon = np.hstack((lon_data2, lon_data3))
    combined_lon = np.vstack((upper_lon, lower_lon))

    # lat の連結
    upper_lat = np.hstack((lat_data0, lat_data1))
    lower_lat = np.hstack((lat_data2, lat_data3))
    combined_lat = np.vstack((upper_lat, lower_lat))

    # 各連結後のサイズを確認
    print("combined_MOMX.shape:", combined_MOMX.shape)  # (90, 90)
    print("combined_MOMY.shape:", combined_MOMY.shape)  # (90, 90)
    print("combined_lon.shape: ", combined_lon.shape)   # (90, 90)
    print("combined_lat.shape: ", combined_lat.shape)   # (90, 90)



    # サブサンプリング：元の lon, lat, MOMX, MOMY の shape が一致するようにする
    Lon_sub = combined_lon[::step, ::step]   # 例: (47,47) -> (ceil(47/5), ceil(47/5)) ≒ (10,10)
    Lat_sub = combined_lat[::step, ::step]
    MOMX_sub = combined_MOMX[::step, ::step]
    MOMY_sub = combined_MOMY[::step, ::step]

    print("lon shape:", Lon_sub.shape)
    print("lat shape:", Lat_sub.shape)
    print("MOMX_sub shape:", MOMX_sub.shape)
    print("Lon_sub shape:", Lon_sub.shape)
    # MOMX_sub の最大値を取得
    max_MOMX = np.max(MOMX_sub)
    min_MOMX = np.min(MOMX_sub)
    print("MOMX_sub の最大値:", max_MOMX)
    print("MOMX_sub の最小値:", min_MOMX)
    # MOMY_sub の最大値を取得
    max_MOMY = np.max(MOMY_sub)
    min_MOMY = np.min(MOMY_sub)
    print("MOMY_sub の最大値:", max_MOMY)
    print("MOMY_sub の最小値:", min_MOMY)
    # ベクトルの大きさ（ノルム）を計算（必要に応じて、scale オプションに利用できます）
    speed = np.sqrt(MOMX_sub**2 + MOMY_sub**2)

    # プロットの作成
    # 投影法の設定（例: PlateCarree は経緯度座標系）
    projection = ccrs.PlateCarree()

    # 図と軸を作成
    fig = plt.figure(figsize=(12, 8))
    ax = plt.axes(projection=projection)
    plt_extent = [125,146,26,43]
    # 日本全域の範囲を設定 (経度: 122〜154, 緯度: 24〜46)
    ax.set_extent(plt_extent, crs=projection)

    # 地図の背景情報を追加
    ax.add_feature(cfeature.LAND, facecolor='lightgray')
    ax.add_feature(cfeature.OCEAN, facecolor='lightblue')
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    # グリッド線を追加（オプション）
    gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
    gl.top_labels = gl.right_labels = False

    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    # quiver プロット
    # # scale: 矢印の長さのスケーリング（データに合わせて調整してください）
    q = plt.quiver(Lon_sub, Lat_sub, MOMX_sub, MOMY_sub, speed, 
                cmap="jet", scale=500, headwidth=3, headlength=4, width=0.005)

    # # 介入を矢印で描画
    # arrow_lon_arr = combined_lon[B_Y:B_Y+Input_X_size, L_X:L_X+Input_X_size]  # 71〜75行, 57〜61列
    # arrow_lat_arr = combined_lat[B_Y:B_Y+Input_X_size, L_X:L_X+Input_X_size]
    # u_arr = np.full(arrow_lon_arr.shape, d_X)
    # v_arr = np.full(arrow_lat_arr.shape, d_Y)
    # ax.quiver(arrow_lon_arr, arrow_lat_arr, u_arr, v_arr,
    #       color=(0, 0, 0, 0.5), scale=500, headwidth=3, headlength=4,width=0.005,
    #       transform=ccrs.PlateCarree())


    cbar = fig.colorbar(q, pad=0.05, aspect=30)
    cbar.set_label(r"$(\text{MOMX}^2+\text{MOMY}^2)^{1/2}$", fontsize=17) # こんな色だけどうまくいっている
    # plt.savefig(f"/home/yuta/scale-5.5.3/scale-rm/test/tutorial/real/data_MOMXY_values/test1_{name}_step={step}_input=({L_X},{B_Y},{d_X},{d_Y}_MOMXY_Z={Z}.png",  dpi=300)
    plt.savefig(f"/home/yuta/scale-5.5.3/scale-rm/test/tutorial/real/data_MOMXY_values/test1_{name}_step={step}_input=({L_X},{B_Y},{d_X},{d_Y}_MOMXY_Z={Z}.pdf" ,dpi=1200, bbox_inches = 'tight')
    plt.close()
    # ファイルクローズ
    dataset0.close()
    dataset1.close()
    dataset2.close()
    dataset3.close()