import re

# ログファイルのパス
log_file_path = '../test_result/accumlated_4var/MOMY_t=0-6_5*5*1grids_FET=200_trials=8-9_12-21-09-52/summary/BO_analysis.txt'
# BO_data を格納するリスト
BO_data = []

# 現在の trial_i のデータを格納するリスト
current_trial = []

# trial_i を追跡する変数
current_trial_i = None

# ログファイルを読み込む
with open(log_file_path, 'r') as file:
    for line in file:
        # trial_i の行を検出
        trial_match = re.match(r'trial_i=(\d+)', line)
        if trial_match:
            # 既存の trial のデータを BO_data に追加
            if current_trial:
                BO_data.append(current_trial)
                current_trial = []
            current_trial_i = int(trial_match.group(1))
            continue  # 次の行へ

        # Batch の行を検出
        batch_match = re.match(r'Batch\s+\d+:\s+Best value so far:\s+([\d\.]+)', line)
        if batch_match and current_trial_i is not None:
            best_value = float(batch_match.group(1))
            current_trial.append(best_value)

# 最後の trial のデータを追加
if current_trial:
    BO_data.append(current_trial)

# 結果を表示
for idx, trial in enumerate(BO_data):
    print(f"trial_i={idx}: {trial}")

# 必要に応じてファイルに保存
import json

with open('../test_result/accumlated_4var/MOMY_t=0-6_5*5*1grids_FET=200_trials=8-9_12-21-09-52/summary/BO_data_parsed.json', 'w') as f:
    json.dump(BO_data, f, indent=4)
