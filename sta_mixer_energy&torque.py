#!/usr/bin/env python3
"""
SUP混合器模拟数据处理脚本
用于提取和比较不同scale_factor下的动能和扭矩数据
"""

import os
import re
import pandas as pd
import numpy as np
from pathlib import Path

def find_sup_directories(base_path="build"):
    """
    查找所有符合SUPMixerOutput_f{scale_factor}格式的目录
    
    参数:
        base_path: 基础搜索路径，默认为"build"
    
    返回:
        包含(目录路径, scale_factor)元组的列表
    """
    pattern = re.compile(r'SUPMixerOutput_f(\d+(?:\.\d+)?)')
    directories = []
    
    # 确保基础路径存在
    if not os.path.exists(base_path):
        print(f"警告：路径 {base_path} 不存在")
        return directories
    
    # 遍历目录查找匹配的文件夹
    for item in os.listdir(base_path):
        item_path = os.path.join(base_path, item)
        if os.path.isdir(item_path):
            match = pattern.match(item)
            if match:
                scale_factor = float(match.group(1))
                directories.append((item_path, scale_factor))
    
    # 按scale_factor排序
    directories.sort(key=lambda x: x[1])
    
    print(f"找到 {len(directories)} 个SUP输出目录:")
    for dir_path, sf in directories:
        print(f"  - {dir_path} (scale_factor={sf})")
    
    return directories

def read_kinetic_energy_file(directory_path):
    """
    读取指定目录下的kinetic_energy_history.csv文件
    
    参数:
        directory_path: 包含CSV文件的目录路径
    
    返回:
        DataFrame或None（如果文件不存在）
    """
    csv_path = os.path.join(directory_path, "kinetic_energy_history.csv")
    
    if not os.path.exists(csv_path):
        print(f"警告：文件 {csv_path} 不存在")
        return None
    
    try:
        df = pd.read_csv(csv_path)
        # 检查必需的列是否存在
        required_columns = ['Time(s)', 'KineticEnergy(J)', 'MixerTorque(Nm)']
        if not all(col in df.columns for col in required_columns):
            print(f"警告：文件 {csv_path} 缺少必需的列")
            return None
        return df
    except Exception as e:
        print(f"错误：读取文件 {csv_path} 时出错: {e}")
        return None

def create_comparison_dataframes(directories):
    """
    创建动能和扭矩的比较数据框
    
    参数:
        directories: 包含(目录路径, scale_factor)元组的列表
    
    返回:
        (动能DataFrame, 扭矩DataFrame)元组
    """
    kinetic_energy_data = {}
    mixer_torque_data = {}
    
    for dir_path, scale_factor in directories:
        df = read_kinetic_energy_file(dir_path)
        if df is not None:
            # 使用scale_factor作为列名
            col_name = f"f{scale_factor}"
            
            # 提取时间、动能和扭矩数据
            time_col = df['Time(s)'].values
            ke_col = df['KineticEnergy(J)'].values
            torque_col = df['MixerTorque(Nm)'].values
            
            # 存储数据
            if 'Time(s)' not in kinetic_energy_data:
                kinetic_energy_data['Time(s)'] = time_col
                mixer_torque_data['Time(s)'] = time_col
            
            kinetic_energy_data[col_name] = ke_col
            mixer_torque_data[col_name] = torque_col
    
    # 创建DataFrame
    ke_df = pd.DataFrame(kinetic_energy_data) if kinetic_energy_data else pd.DataFrame()
    torque_df = pd.DataFrame(mixer_torque_data) if mixer_torque_data else pd.DataFrame()
    
    return ke_df, torque_df

def interpolate_data(dataframes, directories):
    """
    插值数据以处理不同时间步长的情况
    
    参数:
        dataframes: (动能DataFrame, 扭矩DataFrame)元组
        directories: 包含(目录路径, scale_factor)元组的列表
    
    返回:
        插值后的(动能DataFrame, 扭矩DataFrame)元组
    """
    ke_df, torque_df = dataframes
    
    if ke_df.empty or torque_df.empty:
        return ke_df, torque_df
    
    # 找到公共时间范围
    time_min = ke_df['Time(s)'].min()
    time_max = ke_df['Time(s)'].max()
    
    # 创建统一的时间网格（使用最细的时间步长）
    n_points = 1000  # 可以根据需要调整
    common_time = np.linspace(time_min, time_max, n_points)
    
    # 插值动能数据
    ke_interpolated = {'Time(s)': common_time}
    torque_interpolated = {'Time(s)': common_time}
    
    for col in ke_df.columns[1:]:  # 跳过Time列
        ke_interpolated[col] = np.interp(common_time, ke_df['Time(s)'], ke_df[col])
        torque_interpolated[col] = np.interp(common_time, torque_df['Time(s)'], torque_df[col])
    
    return pd.DataFrame(ke_interpolated), pd.DataFrame(torque_interpolated)

def save_comparison_files(ke_df, torque_df):
    """
    保存比较数据到CSV文件（直接保存在当前目录）
    
    参数:
        ke_df: 动能比较DataFrame
        torque_df: 扭矩比较DataFrame
    """
    # 保存文件到当前目录
    ke_output_path = "kinetic_energy_comparison.csv"
    torque_output_path = "mixer_torque_comparison.csv"
    
    ke_df.to_csv(ke_output_path, index=False, float_format='%.6e')
    torque_df.to_csv(torque_output_path, index=False, float_format='%.6e')
    
    print(f"\n文件已保存:")
    print(f"  - 动能比较: {ke_output_path}")
    print(f"  - 扭矩比较: {torque_output_path}")

def main():
    """
    主函数：执行完整的数据处理流程
    """
    print("=" * 60)
    print("SUP混合器数据处理脚本")
    print("=" * 60)
    
    # 1. 查找所有SUP输出目录
    directories = find_sup_directories()
    
    if not directories:
        print("\n错误：未找到任何SUP输出目录")
        return
    
    # 2. 创建比较数据框
    print("\n正在读取和处理数据...")
    ke_df, torque_df = create_comparison_dataframes(directories)
    
    if ke_df.empty or torque_df.empty:
        print("\n错误：无法创建比较数据")
        return
    
    # 3. 可选：插值数据（如果不同模拟有不同的时间步长）
    # 如果需要插值，取消下面这行的注释
    # ke_df, torque_df = interpolate_data((ke_df, torque_df), directories)
    
    # 4. 保存比较文件
    save_comparison_files(ke_df, torque_df)
    
    # 5. 输出统计信息
    print("\n" + "=" * 60)
    print("数据处理完成！")
    print(f"  - 处理了 {len(directories)} 个模拟结果")
    print(f"  - 时间范围: {ke_df['Time(s)'].min():.3f} - {ke_df['Time(s)'].max():.3f} 秒")
    print(f"  - 数据点数: {len(ke_df)}")
    
    # 输出各scale_factor的最大值统计
    print("\n最大值统计:")
    for col in ke_df.columns[1:]:
        max_ke = ke_df[col].max()
        max_torque = torque_df[col].max()
        print(f"  {col}: 最大动能={max_ke:.3e} J, 最大扭矩={max_torque:.3e} Nm")

if __name__ == "__main__":
    main()