#!/usr/bin/env python3
"""
脚本 1 (多目录版): 原始数据特征提取器 (Raw Data to Classification Dataset)
功能: 读取多个目录中的 DEM 模拟原始 CSV，逐帧计算速度 PDF，生成用于机器学习的 dataset.csv

改动: 支持多个 BASE_PATH，适用于不同材料类型放在不同目录的情况
"""

import os
import glob
import re
import pandas as pd
import numpy as np
from pathlib import Path
from tqdm import tqdm

# ================= 配置区域 =================
# 只读取这些 SE 类型（文件夹名中的 se 字段，如 se000 → '000'）
# 设为 None 则读取全部
VALID_SE = {'000', '010', '020', '040'}

GLOBAL_BINS = 50
V_RANGE = (0, 2.0)
BIN_EDGES = np.linspace(V_RANGE[0], V_RANGE[1], GLOBAL_BINS + 1)

def extract_features_from_frame(file_path):
    """读取单个CSV文件，计算速度PDF特征"""
    try:
        df = pd.read_csv(file_path, usecols=['v_x', 'v_y', 'v_z'])
        v_total = np.sqrt(df['v_x']**2 + df['v_y']**2 + df['v_z']**2)
        hist, _ = np.histogram(v_total, bins=BIN_EDGES, density=True)
        return hist
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return None

def process_all_data(base_path):
    dataset = []
    base_path = Path(base_path)
    
    target_dirs = [
        ('Standard', r'SUPMixerOutput_f(\d+)se(\d+)$', 0),
        ('SizeDiff', r'SUPMixerOutput_SizeDiff_f(\d+)se(\d+)$', 1),
        ('DualDensity', r'SUPMixerOutput_DualDensity_f(\d+)se(\d+)$', 2)
    ]

    print(f"正在扫描目录: {base_path} ...")
    
    for folder_path in sorted(base_path.iterdir()):
        if not folder_path.is_dir():
            continue
            
        folder_name = folder_path.name
        matched_type = None
        matched_info = None

        for type_name, pattern, label_code in target_dirs:
            match = re.match(pattern, folder_name)
            if match:
                matched_type = type_name
                matched_info = match
                current_label = label_code
                break
        
        if not matched_type:
            continue

        scale_factor = int(matched_info.group(1))
        surface_energy = matched_info.group(2)

        # 跳过不在白名单中的 SE 类型
        if VALID_SE is not None and surface_energy not in VALID_SE:
            continue

        csv_pattern = folder_path / 'mixer_output_*.csv'
        files = sorted(glob.glob(str(csv_pattern)))
        
        if not files:
            print(f"  [跳过] 空文件夹: {folder_name}")
            continue

        frames = []
        for f in files:
            m = re.search(r'_(\d+)\.csv', f)
            if m:
                frames.append((f, int(m.group(1))))
        
        frames.sort(key=lambda x: x[1])
        
        cutoff_index = int(len(frames) * 0.5)
        steady_frames = frames[cutoff_index:]
        
        step = 1 
        selected_frames = steady_frames[::step]

        print(f"  处理: {folder_name} | 类型: {matched_type} | 样本数: {len(selected_frames)}")

        for file_path, frame_num in tqdm(selected_frames, desc="    提取特征", leave=False):
            features = extract_features_from_frame(file_path)
            
            if features is not None:
                sample = {
                    'material_type': matched_type,
                    'label': current_label,
                    'surface_energy': surface_energy,
                    'scale_factor': scale_factor,
                    'frame': frame_num,
                    'source_folder': folder_name
                }
                
                for i, val in enumerate(features):
                    sample[f'bin_{i}'] = val
                
                dataset.append(sample)

    return dataset


# ================= 主程序 =================
if __name__ == "__main__":
    # ============================================================
    #  在这里填写你的目录路径，每种材料一个或多个目录都行
    #  脚本会依次扫描所有目录，自动识别 Standard/SizeDiff/DualDensity
    # ============================================================
    BASE_PATHS = [
        "/media/huyuze/Fanxiang1/DEME_Data/backup251115",
        "/media/huyuze/Fanxiang/backup251025",
        "/media/huyuze/Fanxiang/backup251109"
        # 可以继续添加更多目录...
    ]

    OUTPUT_FILENAME = "Midterm_Raw_PDF_Dataset.csv"

    # 逐目录扫描，合并结果
    all_samples = []
    for path in BASE_PATHS:
        if not os.path.exists(path):
            print(f"[警告] 路径不存在，跳过: {path}")
            continue
        samples = process_all_data(path)
        print(f"  -> 从 {path} 提取了 {len(samples)} 个样本")
        all_samples.extend(samples)

    df_final = pd.DataFrame(all_samples)

    if not df_final.empty:
        df_final.to_csv(OUTPUT_FILENAME, index=False)
        print("\n" + "=" * 60)
        print(f" 成功! 数据集已保存: {OUTPUT_FILENAME}")
        print("=" * 60)
        print(f"总样本数: {len(df_final)}")
        print(f"特征维度: {GLOBAL_BINS}")
        print(f"\n材料类型分布:")
        print(df_final['material_type'].value_counts().to_string())
        print(f"\nSE分布:")
        print(df_final['surface_energy'].value_counts().to_string())

        # 检查是否集齐 3 种
        found = set(df_final['material_type'].unique())
        required = {'Standard', 'SizeDiff', 'DualDensity'}
        missing = required - found
        if missing:
            print(f"\n[注意] 缺少材料类型: {missing}")
            print(f"请在 BASE_PATHS 中添加包含这些材料的目录路径")
        else:
            print(f"\n 3 种材料类型齐全，可以进行下一步:")
            print(f"  python generate_midterm_dataset.py")
    else:
        print("未提取到任何数据，请检查路径或文件名格式。")
