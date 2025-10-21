#!/usr/bin/env python3
"""
SUP混合器模拟数据处理脚本
用于提取和比较不同scale_factor下的动能和扭矩数据
适配目录格式：SUPMixerOutput_f{scale_factor}se{cohesion_code}
"""

import os
import re
import pandas as pd
import numpy as np
from pathlib import Path

def find_sup_directories(base_path=".", pattern_str=None):
    """
    查找所有符合SUPMixerOutput_f{scale_factor}se{cohesion}格式的目录
    
    参数:
        base_path: 基础搜索路径，默认为当前目录"."
        pattern_str: 自定义正则表达式模式（可选）
    
    返回:
        包含(目录路径, scale_factor, cohesion_code)元组的列表
    """
    # 默认模式：匹配 SUPMixerOutput_f{数字}se{三位数字}
    if pattern_str is None:
        pattern = re.compile(r'SUPMixerOutput_f(\d+(?:\.\d+)?)se(\d{3})')
    else:
        pattern = re.compile(pattern_str)
    
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
                # 如果有cohesion代码，也提取出来
                cohesion_code = match.group(2) if len(match.groups()) > 1 else "000"
                directories.append((item_path, scale_factor, cohesion_code))
    
    # 按scale_factor排序
    directories.sort(key=lambda x: x[1])
    
    print(f"找到 {len(directories)} 个SUP输出目录:")
    for dir_path, sf, cc in directories:
        print(f"  - {dir_path} (scale_factor={sf}, cohesion_code={cc})")
    
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

def create_comparison_dataframes(directories, filter_cohesion="000"):
    """
    创建动能和扭矩的比较数据框
    
    参数:
        directories: 包含(目录路径, scale_factor, cohesion_code)元组的列表
        filter_cohesion: 只处理特定cohesion代码的数据，默认为"000"
    
    返回:
        (动能DataFrame, 扭矩DataFrame)元组
    """
    kinetic_energy_data = {}
    mixer_torque_data = {}
    
    for dir_path, scale_factor, cohesion_code in directories:
        # 如果指定了过滤条件，只处理符合条件的数据
        if filter_cohesion and cohesion_code != filter_cohesion:
            continue
            
        df = read_kinetic_energy_file(dir_path)
        if df is not None:
            # 使用scale_factor作为列名
            col_name = f"f{scale_factor}"
            
            # 如果有多个cohesion值，可以在列名中包含
            if filter_cohesion is None:
                col_name = f"f{scale_factor}_se{cohesion_code}"
            
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

def interpolate_data(dataframes, n_points=1000):
    """
    插值数据以处理不同时间步长的情况
    
    参数:
        dataframes: (动能DataFrame, 扭矩DataFrame)元组
        n_points: 插值点数，默认1000
    
    返回:
        插值后的(动能DataFrame, 扭矩DataFrame)元组
    """
    ke_df, torque_df = dataframes
    
    if ke_df.empty or torque_df.empty:
        return ke_df, torque_df
    
    # 找到公共时间范围
    time_min = ke_df['Time(s)'].min()
    time_max = ke_df['Time(s)'].max()
    
    # 创建统一的时间网格
    common_time = np.linspace(time_min, time_max, n_points)
    
    # 插值数据
    ke_interpolated = {'Time(s)': common_time}
    torque_interpolated = {'Time(s)': common_time}
    
    for col in ke_df.columns[1:]:  # 跳过Time列
        ke_interpolated[col] = np.interp(common_time, ke_df['Time(s)'], ke_df[col])
        torque_interpolated[col] = np.interp(common_time, torque_df['Time(s)'], torque_df[col])
    
    return pd.DataFrame(ke_interpolated), pd.DataFrame(torque_interpolated)

def save_comparison_files(ke_df, torque_df, output_dir="sta_results/torque&energy"):
    """
    保存比较数据到CSV文件
    
    参数:
        ke_df: 动能比较DataFrame
        torque_df: 扭矩比较DataFrame
        output_dir: 输出目录
    """
    # 创建输出目录（如果不存在）
    os.makedirs(output_dir, exist_ok=True)
    print(f"\n输出目录: {output_dir}")
    
    # 构建输出文件路径
    ke_output_path = os.path.join(output_dir, "kinetic_energy_comparison.csv")
    torque_output_path = os.path.join(output_dir, "mixer_torque_comparison.csv")
    
    # 保存文件
    ke_df.to_csv(ke_output_path, index=False, float_format='%.6e')
    torque_df.to_csv(torque_output_path, index=False, float_format='%.6e')
    
    print(f"文件已保存:")
    print(f"  - 动能比较: {ke_output_path}")
    print(f"  - 扭矩比较: {torque_output_path}")

def analyze_scale_factor_groups(base_path=".", cohesion_codes=["000"]):
    """
    分析不同cohesion代码下的scale factor组
    
    参数:
        base_path: 基础搜索路径
        cohesion_codes: 要分析的cohesion代码列表
    """
    all_dirs = find_sup_directories(base_path)
    
    # 按cohesion代码分组
    groups = {}
    for dir_path, sf, cc in all_dirs:
        if cc not in groups:
            groups[cc] = []
        groups[cc].append((dir_path, sf, cc))
    
    print(f"\n发现的cohesion组:")
    for cc, dirs in groups.items():
        print(f"  Cohesion {cc}: {len(dirs)} 个scale factors")
        for _, sf, _ in sorted(dirs, key=lambda x: x[1]):
            print(f"    - f{sf}")
    
    return groups

def main():
    """
    主函数：执行完整的数据处理流程
    """
    print("=" * 60)
    print("SUP混合器数据处理脚本")
    print("=" * 60)
    
    # 可配置参数
    base_path = "build"  # 搜索目录，可以改为"build"或其他路径
    filter_cohesion = "010"  # 只处理cohesion=0的数据，设为None处理所有
    use_interpolation = False  # 是否使用插值
    
    # 1. 查找所有SUP输出目录
    directories = find_sup_directories(base_path)
    
    if not directories:
        print("\n错误：未找到任何SUP输出目录")
        print("请确保目录格式为: SUPMixerOutput_f{数字}se{三位数字}")
        return
    
    # 2. 分析发现的目录
    print(f"\n分析配置:")
    print(f"  - 基础路径: {base_path}")
    print(f"  - Cohesion过滤: {filter_cohesion if filter_cohesion else '无(处理所有)'}")
    print(f"  - 使用插值: {'是' if use_interpolation else '否'}")
    
    # 3. 创建比较数据框
    print("\n正在读取和处理数据...")
    ke_df, torque_df = create_comparison_dataframes(directories, filter_cohesion)
    
    if ke_df.empty or torque_df.empty:
        print("\n错误：无法创建比较数据")
        return
    
    # 4. 可选：插值数据
    if use_interpolation:
        print("正在插值数据...")
        ke_df, torque_df = interpolate_data((ke_df, torque_df))
    
    # 5. 保存比较文件
    save_comparison_files(ke_df, torque_df)
    
    # 6. 输出统计信息
    print("\n" + "=" * 60)
    print("数据处理完成！")
    print(f"  - 处理了 {len(ke_df.columns) - 1} 个模拟结果")
    print(f"  - 时间范围: {ke_df['Time(s)'].min():.3f} - {ke_df['Time(s)'].max():.3f} 秒")
    print(f"  - 数据点数: {len(ke_df)}")
    
    # 输出各scale_factor的统计
    print("\n最大值统计:")
    for col in ke_df.columns[1:]:
        max_ke = ke_df[col].max()
        max_torque = torque_df[col].abs().max()  # 使用绝对值
        avg_ke = ke_df[col].mean()
        avg_torque = torque_df[col].abs().mean()
        print(f"  {col}:")
        print(f"    - 动能: 最大={max_ke:.3e} J, 平均={avg_ke:.3e} J")
        print(f"    - 扭矩: 最大={max_torque:.3e} Nm, 平均={avg_torque:.3e} Nm")

if __name__ == "__main__":
    main()