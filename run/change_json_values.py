import json

def subtract_value_from_element(element, value):
    """
    要素がリストの場合は再帰的に処理し、数値の場合は値を引く。
    """
    if isinstance(element, list):
        return [subtract_value_from_element(item, value) for item in element]
    elif isinstance(element, dict):
        return {key: subtract_value_from_element(val, value) for key, val in element.items()}
    elif isinstance(element, (int, float)):
        return value - element
    else:
        # 数値以外のデータはそのまま返す
        return element

def process_json(input_file, output_file, subtract_value):
    # JSONファイルを読み込む
    with open(input_file, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    # 処理を実行
    processed_data = subtract_value_from_element(data, subtract_value)
    
    # 結果を新しいJSONファイルに保存
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(processed_data, f, ensure_ascii=False, indent=4)

if __name__ == "__main__":
    input_file = '../test_result/accumlated_xy_MOMXY/t=0-6_5*5*5grids_FET=200_trials=0-9_01-07-18-13/summary/BO_data.json'       # 入力JSONファイルのパス
    output_file = '../test_result/accumlated_xy_MOMXY/t=0-6_5*5*5grids_FET=200_trials=0-9_01-07-18-13/summary/BO_data_diff.json'     # 出力JSONファイルのパス
    subtract_value = 7162.777587601443              # 引きたい値

    process_json(input_file, output_file, subtract_value)
    #print(f"{input_file} の各数値から {subtract_value} を引き、{output_file} に保存しました。")
