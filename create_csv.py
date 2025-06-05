from obspy import read
import pandas as pd
import os
import csv  # ← これを忘れずに追加

# 出力フォルダ（必要に応じて変更）
output_dir = "converted_csv"
os.makedirs(output_dir, exist_ok=True)

# 0.mseed から 9.mseed まで処理
for i in range(10):
    filename = f"{i}.mseed"
    st = read(filename)
    
    for trace in st:
        # トレース情報からファイル名を生成
        net = trace.stats.network
        sta = trace.stats.station
        loc = trace.stats.location or "00"  # locationが空の場合は"00"にする
        cha = trace.stats.channel
        start_time = trace.stats.starttime
        end_time = trace.stats.endtime
        sampling_rate = trace.stats.sampling_rate
        delta = trace.stats.delta
        npts = trace.stats.npts
        calib = trace.stats.calib
        _format = trace.stats._format
        mseed = trace.stats.mseed  # mseed情報

        # CSVファイル名
        csv_filename = f"{net}_{sta}_{loc}_{cha}_{i}.csv"
        csv_path = os.path.join(output_dir, csv_filename)

        # 時系列データをDataFrameに変換
        times = trace.times("utcdatetime")  # UTCDateTime型
        values = trace.data  # NumPy array
        
        # DataFrameにすべての情報を追加
        df = pd.DataFrame({
            "time": times,
            "amplitude": values,
            "network": net,
            "station": sta,
            "location": loc,
            "channel": cha,
            "starttime": start_time,
            "endtime": end_time,
            "sampling_rate": sampling_rate,
            "delta": delta,
            "npts": npts,
            "calib": calib,
            "_format": _format,
            "mseed_dataquality": mseed["dataquality"],
            "mseed_number_of_records": mseed["number_of_records"],
            "mseed_encoding": mseed["encoding"],
            "mseed_byteorder": mseed["byteorder"],
            "mseed_record_length": mseed["record_length"],
            "mseed_filesize": mseed["filesize"]
        })

        # 書き出し
        df.to_csv(csv_path, index=False)

print("全ファイルの変換が完了しました。")
