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
Z= 0 #どの高さにするか
name = "init_file"
for Z in range(5):
    print(f"{Z=}")
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

    # --- DENS ---
    DENS_data0 = datasets[0].variables["DENS"][2:, 2:,Z]
    DENS_data1 = datasets[1].variables["DENS"][2:, 0:-2,Z]
    DENS_data2 = datasets[2].variables["DENS"][0:-2, 2:,Z]
    DENS_data3 = datasets[3].variables["DENS"][0:-2, 0:-2,Z]

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


    # DENS の連結
    upper_DENS = np.hstack((DENS_data0, DENS_data1))
    lower_DENS = np.hstack((DENS_data2, DENS_data3))
    combined_DENS = np.vstack((upper_DENS, lower_DENS))
    # lon の連結
    upper_lon = np.hstack((lon_data0, lon_data1))
    lower_lon = np.hstack((lon_data2, lon_data3))
    combined_lon = np.vstack((upper_lon, lower_lon))

    Wind_X = combined_MOMX/combined_DENS
    Wind_Y = combined_MOMY/combined_DENS
    # lat の連結
    upper_lat = np.hstack((lat_data0, lat_data1))
    lower_lat = np.hstack((lat_data2, lat_data3))
    combined_lat = np.vstack((upper_lat, lower_lat))

    # 各連結後のサイズを確認
    # print("Wind_X.shape:", Wind_X.shape)  # (90, 90)
    # print("Wind_Y.shape:", Wind_Y.shape)  # (90, 90)
    # print("combined_lon.shape: ", combined_lon.shape)   # (90, 90)
    # print("combined_lat.shape: ", combined_lat.shape)   # (90, 90)


    # グリッドが密すぎる場合は、サブサンプリング（例：5点ごとに表示）
    step = 1
    # サブサンプリング：元の lon, lat, MOMX, MOMY の shape が一致するようにする
    Lon_sub = combined_lon[::step, ::step]   # 例: (47,47) -> (ceil(47/5), ceil(47/5)) ≒ (10,10)
    Lat_sub = combined_lat[::step, ::step]
    Wind_X_sub = Wind_X[::step, ::step]
    Wind_Y_sub = Wind_Y[::step, ::step]

    # print("lon shape:", Lon_sub.shape)
    # print("lat shape:", Lat_sub.shape)
    # print("Wind_X_sub shape:", Wind_X_sub.shape)
    # print("Lon_sub shape:", Lon_sub.shape)
    # Wind_X_sub の最大値を取得
    max_Wind_X = np.max(Wind_X_sub)
    min_Wind_X = np.min(Wind_X_sub)
    print("Wind_X_sub の最大値:", max_Wind_X)
    print("Wind_X_sub の最小値:", min_Wind_X)
    # Wind_Y_sub の最大値を取得
    max_Wind_Y = np.max(Wind_Y_sub)
    min_Wind_Y = np.min(Wind_Y_sub)
    print("Wind_Y_sub の最大値:", max_Wind_Y)
    print("Wind_Y_sub の最小値:", min_Wind_Y)
    # ベクトルの大きさ（ノルム）を計算（必要に応じて、scale オプションに利用できます）
    speed = np.sqrt(Wind_X_sub**2 + Wind_Y_sub**2)

    # プロットの作成
    # 投影法の設定（例: PlateCarree は経緯度座標系）
    projection = ccrs.PlateCarree()

    # 図と軸を作成
    fig = plt.figure(figsize=(8, 8))
    ax = plt.axes(projection=projection)

    # 日本全域の範囲を設定 (経度: 122〜154, 緯度: 24〜46)
    ax.set_extent([124, 146, 26, 44], crs=projection)

    # 地図の背景情報を追加
    ax.add_feature(cfeature.LAND, facecolor='lightgray')
    ax.add_feature(cfeature.OCEAN, facecolor='lightblue')
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linestyle=':')

    # グリッド線を追加（オプション）
    gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
    gl.top_labels = gl.right_labels = False
    # plt.figure(figsize=(10, 8))
    plt.title("Initial Time Vector Field (Wind_X, Wind_Y)")
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    # quiver プロット
    # scale: 矢印の長さのスケーリング（データに合わせて調整してください）
    q = plt.quiver(Lon_sub, Lat_sub, Wind_X_sub, Wind_Y_sub, speed, 
                cmap="jet", scale=500, headwidth=3, headlength=4)
    # plt.scatter(Lon_sub, Lat_sub, color='k', s=10, zorder=3)
    plt.colorbar(q, label="Vector Magnitude")
    plt.savefig(f"/home/yuta/scale-5.5.3/scale-rm/test/tutorial/real/data_Wind_XY_values/{name}_step={step}_Wind_XY_Z={Z}.png", dpi=300)
    # ファイルクローズ
    dataset0.close()
    dataset1.close()
    dataset2.close()
    dataset3.close()