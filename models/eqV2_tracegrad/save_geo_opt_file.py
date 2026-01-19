import gzip
import json

import pandas as pd


def detect_file_format(file_path):
    """自动检测JSON文件的格式类型"""
    with gzip.open(file_path, "rt") as f:
        first_char = f.read(1)
        # JSON数组格式以 [ 开头
        if first_char == "[":
            return "array"
        # 行分隔JSON通常以 { 开头
        try:
            f.seek(0)
            line = f.readline().strip()
            json.loads(line)  # 验证是否为完整JSON对象
            return "lines"
        except json.JSONDecodeError:
            return "unknown"


# 原始文件路径
ml_relaxed_structs_path = "/data/fywang/code/matbench-discovery/models/eqV2_tracegrad/eqV2-small-dens-tracegrad.json.gz"
geo_opt_path = "/data/fywang/code/matbench-discovery/models/eqV2_tracegrad/eqV2-small-dens-tracegrad-geo-opt.json.gz"

# 自动检测文件格式
# file_format = detect_file_format(ml_relaxed_structs_path)

# # 根据检测结果设置保存参数
# save_params = {
#     'compression': 'gzip',
#     'index': False,
#     'orient': 'records',
#     'lines': False
# }

# if file_format == 'array':
#     save_params.update(lines=False)
# elif file_format == 'lines':
#     save_params.update(lines=True)
# else:
#     raise ValueError(f"无法识别的文件格式: {file_format}")

df_ml_structs = pd.read_json(ml_relaxed_structs_path)

# del df_ml_structs["eqV2-31M-dens-MP-p5_energy"]
df_ml_structs["index"] = range(len(df_ml_structs))
# 保存清洗后的数据
df_ml_structs.to_json(geo_opt_path)
