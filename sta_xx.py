#!/usr/bin/env python3
"""
脚本 1: 原始数据特征提取器 (Raw Data to Classification Dataset)
功能: 读取 DEM 模拟的原始 CSV，逐帧计算速度 PDF，生成用于机器学习的 dataset.csv
"""

import os
import glob
import re
import pandas as pd
import numpy as np
from pathlib import Path
from tqdm import tqdm  # 进度条库，如果没有请运行 pip install tqdm

# ================= 配置区域 =================
# 全局直方图设置 (必须保证所有样本使用相同的 Bin 定义)
# 根据你的物理经验，最大速度大约是 1.5m/s 左右，我们设为 2.0 以防万一
GLOBAL_BINS = 50
V_RANGE = (0, 2.0)
BIN_EDGES = np.linspace(V_RANGE[0], V_RANGE[1], GLOBAL_BINS + 1)

def extract_features_from_frame(file_path):
    """读取单个CSV文件，计算速度PDF特征"""
    try:
        # 只读取我们需要的列，加快速度
        df = pd.read_csv(file_path, usecols=['v_x', 'v_y', 'v_z'])
        
        # 计算合速度 v_total
        v_total = np.sqrt(df['v_x']**2 + df['v_y']**2 + df['v_z']**2)
        
        # 计算直方图 (PDF)
        # density=True 保证面积为1，消除了粒子数量不同带来的影响
        hist, _ = np.histogram(v_total, bins=BIN_EDGES, density=True)
        
        return hist
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return None

def process_all_data(base_path):
    dataset = []
    base_path = Path(base_path)
    
    # 定义我们要寻找的三种材料类型的文件夹特征
    # 格式: (类型名称, 文件夹匹配正则, 对应的分类标签Label)
    # Label: 0=Standard, 1=SizeDiff, 2=DualDensity
    target_dirs = [
        ('Standard', r'SUPMixerOutput_f(\d+)se(\d+)$', 0),
        ('SizeDiff', r'SUPMixerOutput_SizeDiff_f(\d+)se(\d+)$', 1),
        ('DualDensity', r'SUPMixerOutput_DualDensity_f(\d+)se(\d+)$', 2)
    ]

    print(f"正在扫描目录: {base_path} ...")
    
    # 遍历主目录下的所有文件夹
    for folder_path in sorted(base_path.iterdir()):
        if not folder_path.is_dir():
            continue
            
        folder_name = folder_path.name
        matched_type = None
        matched_info = None

        # 检查文件夹是否匹配我们的目标类型
        for type_name, pattern, label_code in target_dirs:
            match = re.match(pattern, folder_name)
            if match:
                matched_type = type_name
                matched_info = match
                current_label = label_code
                break
        
        if not matched_type:
            continue

        # 解析文件夹参数
        scale_factor = int(matched_info.group(1))
        surface_energy = matched_info.group(2) # 保持字符串格式 "000", "010" 等

        # 查找该文件夹下的所有 output csv
        csv_pattern = folder_path / 'mixer_output_*.csv'
        files = sorted(glob.glob(str(csv_pattern)))
        
        if not files:
            print(f"  [跳过] 空文件夹: {folder_name}")
            continue

        # --- 稳态筛选逻辑 ---
        # 我们通过文件名提取帧号，只取后 50% 的数据作为稳态
        # 例如: mixer_output_1000.csv -> frame 1000
        frames = []
        for f in files:
            m = re.search(r'_(\d+)\.csv', f)
            if m:
                frames.append((f, int(m.group(1))))
        
        frames.sort(key=lambda x: x[1]) # 按帧号排序
        
        # 截取后 50% (Steady State)
        cutoff_index = int(len(frames) * 0.5)
        steady_frames = frames[cutoff_index:]
        
        # 为了防止数据量过大且相邻帧太相似，我们可以每隔几帧取一个
        # 如果你数据不够，把 step 改成 1
        step = 1 
        selected_frames = steady_frames[::step]

        print(f"  处理: {folder_name} | 类型: {matched_type} | 样本数: {len(selected_frames)}")

        # --- 逐帧处理 ---
        for file_path, frame_num in tqdm(selected_frames, desc="    提取特征", leave=False):
            features = extract_features_from_frame(file_path)
            
            if features is not None:
                # 构建单行样本数据
                sample = {
                    'material_type': matched_type,      # 原始字符串标签 (用于看)
                    'label': current_label,             # 数字标签 (0,1,2 - 用于训练)
                    'surface_energy': surface_energy,   # SE 值 (也可以作为分类目标)
                    'scale_factor': scale_factor,
                    'frame': frame_num,
                    'source_folder': folder_name
                }
                
                # 将 50 个 bin 的值放入字典
                for i, val in enumerate(features):
                    sample[f'bin_{i}'] = val
                
                dataset.append(sample)

    return pd.DataFrame(dataset)

# ================= 主程序 =================
if __name__ == "__main__":
    # --- 请修改这里的路径 ---
    BASE_PATH = "/media/huyuze/Fanxiang1/backup251025"
    OUTPUT_FILENAME = "Assignment2_PDF_Dataset.csv"

    # 1. 处理数据
    if not os.path.exists(BASE_PATH):
        print(f"错误: 路径不存在 {BASE_PATH}")
    else:
        df_final = process_all_data(BASE_PATH)

        # 2. 保存结果
        if not df_final.empty:
            df_final.to_csv(OUTPUT_FILENAME, index=False)
            print("\n" + "="*50)
            print(f"成功! 数据集已保存为: {OUTPUT_FILENAME}")
            print("="*50)
            print(f"总样本数: {len(df_final)}")
            print(f"特征维度: {GLOBAL_BINS}")
            print("\n类别分布:")
            print(df_final['material_type'].value_counts())
            print("\nSE分布:")
            print(df_final['surface_energy'].value_counts())
        else:
            print("未提取到任何数据，请检查路径或文件名格式。")