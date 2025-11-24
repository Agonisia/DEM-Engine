#!/usr/bin/env python3
"""
SUP混合器模拟数据处理脚本
用于提取和比较不同scale_factor下的动能和扭矩数据
适配三种目录格式：
  - SUPMixerOutput_f{scale_factor}se{cohesion_code}
  - SUPMixerOutput_SizeDiff_f{scale_factor}se{cohesion_code}
  - SUPMixerOutput_DualDensity_f{scale_factor}se{cohesion_code}
"""

import os
import re
import pandas as pd
import numpy as np
from pathlib import Path

def find_sup_directories(base_path=".", pattern_str=None):
    """
    查找所有符合SUPMixerOutput格式的目录（支持三种格式）
    
    参数:
        base_path: 基础搜索路径，默认为当前目录"."
        pattern_str: 自定义正则表达式模式（可选）
    
    返回:
        包含(目录路径, scale_factor, cohesion_code, experiment_type)元组的列表
    """
    directories = []
    
    # 确保基础路径存在
    if not os.path.exists(base_path):
        print(f"警告：路径 {base_path} 不存在")
        return directories
    
    # 定义三种模式
    patterns = [
        ('Standard', re.compile(r'SUPMixerOutput_f(\d+(?:\.\d+)?)se(\d{3})$')),
        ('SizeDiff', re.compile(r'SUPMixerOutput_SizeDiff_f(\d+(?:\.\d+)?)se(\d{3})$')),
        ('DualDensity', re.compile(r'SUPMixerOutput_DualDensity_f(\d+(?:\.\d+)?)se(\d{3})$'))
    ]
    
    # 遍历目录查找匹配的文件夹
    for item in os.listdir(base_path):
        item_path = os.path.join(base_path, item)
        if os.path.isdir(item_path):
            for exp_type, pattern in patterns:
                match = pattern.match(item)
                if match:
                    scale_factor = float(match.group(1))
                    cohesion_code = match.group(2)
                    directories.append((item_path, scale_factor, cohesion_code, exp_type))
                    break  # 找到匹配就跳出内层循环
    
    # 按experiment_type，然后scale_factor排序
    directories.sort(key=lambda x: (x[3], x[1]))
    
    # 按实验类型分组显示
    print(f"找到 {len(directories)} 个SUP输出目录:")
    exp_groups = {}
    for dir_path, sf, cc, exp_type in directories:
        if exp_type not in exp_groups:
            exp_groups[exp_type] = []
        exp_groups[exp_type].append((dir_path, sf, cc))
    
    for exp_type in ['Standard', 'SizeDiff', 'DualDensity']:
        if exp_type in exp_groups:
            print(f"\n  {exp_type}:")
            for dir_path, sf, cc in exp_groups[exp_type]:
                dir_name = os.path.basename(dir_path)
                print(f"    - {dir_name} (f={sf}, se={cc})")
    
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
        
        # 检查必需的列是否存在（支持新旧两种格式）
        # 新格式使用 TotalKE(J)，旧格式使用 KineticEnergy(J)
        has_new_format = all(col in df.columns for col in ['Time(s)', 'TotalKE(J)', 'MixerTorque(Nm)'])
        has_old_format = all(col in df.columns for col in ['Time(s)', 'KineticEnergy(J)', 'MixerTorque(Nm)'])
        
        if has_new_format:
            # 将新格式列名重命名为统一格式
            df = df.rename(columns={'TotalKE(J)': 'KineticEnergy(J)'})
        elif not has_old_format:
            print(f"警告：文件 {csv_path} 缺少必需的列")
            print(f"  需要: Time(s), TotalKE(J)/KineticEnergy(J), MixerTorque(Nm)")
            print(f"  实际列: {list(df.columns)}")
            return None
            
        return df
    except Exception as e:
        print(f"错误：读取文件 {csv_path} 时出错: {e}")
        return None

def create_comparison_dataframes(directories, filter_cohesion="000", filter_exp_type=None):
    """
    创建动能和扭矩的比较数据框
    
    参数:
        directories: 包含(目录路径, scale_factor, cohesion_code, exp_type)元组的列表
        filter_cohesion: 只处理特定cohesion代码的数据，默认为"000"
        filter_exp_type: 只处理特定实验类型的数据，None表示处理所有
    
    返回:
        (动能DataFrame, 扭矩DataFrame)元组
    """
    kinetic_energy_data = {}
    mixer_torque_data = {}
    
    for dir_path, scale_factor, cohesion_code, exp_type in directories:
        # 如果指定了过滤条件，只处理符合条件的数据
        if filter_cohesion and cohesion_code != filter_cohesion:
            continue
        if filter_exp_type and exp_type != filter_exp_type:
            continue
            
        df = read_kinetic_energy_file(dir_path)
        if df is not None:
            # 使用scale_factor和实验类型作为列名
            if exp_type == 'Standard':
                col_name = f"f{scale_factor}"
            else:
                col_name = f"{exp_type}_f{scale_factor}"
            
            # 如果有多个cohesion值，可以在列名中包含
            if filter_cohesion is None:
                col_name = f"{col_name}_se{cohesion_code}"
            
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

def save_comparison_files(ke_df, torque_df, output_dir="sta_results/torque&energy", suffix=""):
    """
    保存比较数据到CSV文件
    
    参数:
        ke_df: 动能比较DataFrame
        torque_df: 扭矩比较DataFrame
        output_dir: 输出目录
        suffix: 文件名后缀
    """
    # 创建输出目录（如果不存在）
    os.makedirs(output_dir, exist_ok=True)
    print(f"\n输出目录: {output_dir}")
    
    # 构建输出文件路径
    ke_filename = f"kinetic_energy_comparison{suffix}.csv"
    torque_filename = f"mixer_torque_comparison{suffix}.csv"
    ke_output_path = os.path.join(output_dir, ke_filename)
    torque_output_path = os.path.join(output_dir, torque_filename)
    
    # 保存文件
    ke_df.to_csv(ke_output_path, index=False, float_format='%.6e')
    torque_df.to_csv(torque_output_path, index=False, float_format='%.6e')
    
    print(f"文件已保存:")
    print(f"  - 动能比较: {ke_output_path}")
    print(f"  - 扭矩比较: {torque_output_path}")

def analyze_scale_factor_groups(base_path="."):
    """
    分析不同实验类型和cohesion代码下的scale factor组
    
    参数:
        base_path: 基础搜索路径
    """
    all_dirs = find_sup_directories(base_path)
    
    # 按实验类型和cohesion代码分组
    groups = {}
    for dir_path, sf, cc, exp_type in all_dirs:
        key = (exp_type, cc)
        if key not in groups:
            groups[key] = []
        groups[key].append((dir_path, sf, cc))
    
    print(f"\n发现的实验组:")
    for (exp_type, cc), dirs in sorted(groups.items()):
        print(f"  {exp_type} - Cohesion {cc}: {len(dirs)} 个scale factors")
        for _, sf, _ in sorted(dirs, key=lambda x: x[1]):
            print(f"    - f{sf}")
    
    return groups

def main():
    """
    主函数：执行完整的数据处理流程
    """
    print("=" * 60)
    print("SUP混合器数据处理脚本 - 支持多种目录格式")
    print("=" * 60)
    
    # 可配置参数
    base_path = "/media/huyuze/Fanxiang/DEME_Data/backup251115"  # 搜索目录
    filter_cohesion = "020"  # 只处理特定cohesion的数据，设为None处理所有
    filter_exp_type = None  # 只处理特定实验类型，None处理所有，可选值: 'Standard', 'SizeDiff', 'DualDensity'
    use_interpolation = False  # 是否使用插值
    separate_by_type = False  # 是否按实验类型分别保存文件
    
    # 1. 查找所有SUP输出目录
    directories = find_sup_directories(base_path)
    
    if not directories:
        print("\n错误：未找到任何SUP输出目录")
        print("请确保目录格式为以下之一:")
        print("  - SUPMixerOutput_f{数字}se{三位数字}")
        print("  - SUPMixerOutput_SizeDiff_f{数字}se{三位数字}")
        print("  - SUPMixerOutput_DualDensity_f{数字}se{三位数字}")
        return
    
    # 2. 分析发现的目录
    print(f"\n分析配置:")
    print(f"  - 基础路径: {base_path}")
    print(f"  - Cohesion过滤: {filter_cohesion if filter_cohesion else '无(处理所有)'}")
    print(f"  - 实验类型过滤: {filter_exp_type if filter_exp_type else '无(处理所有)'}")
    print(f"  - 使用插值: {'是' if use_interpolation else '否'}")
    print(f"  - 按类型分别保存: {'是' if separate_by_type else '否'}")
    
    if separate_by_type:
        # 按实验类型分别处理和保存
        exp_types = set(d[3] for d in directories)
        for exp_type in sorted(exp_types):
            print(f"\n{'='*40}")
            print(f"处理 {exp_type} 类型数据...")
            print('='*40)
            
            # 创建比较数据框
            ke_df, torque_df = create_comparison_dataframes(
                directories, filter_cohesion, filter_exp_type=exp_type
            )
            
            if ke_df.empty or torque_df.empty:
                print(f"  警告：{exp_type} 类型无有效数据")
                continue
            
            # 可选：插值数据
            if use_interpolation:
                print(f"  正在插值 {exp_type} 数据...")
                ke_df, torque_df = interpolate_data((ke_df, torque_df))
            
            # 保存比较文件
            output_dir = f"sta_results/torque&energy/{exp_type}"
            save_comparison_files(ke_df, torque_df, output_dir)
            
            # 输出统计信息
            print(f"\n  {exp_type} 统计:")
            for col in ke_df.columns[1:]:
                max_ke = ke_df[col].max()
                max_torque = torque_df[col].abs().max()
                print(f"    {col}: KE_max={max_ke:.3e}J, T_max={max_torque:.3e}Nm")
    else:
        # 处理所有数据并合并保存
        print("\n正在读取和处理数据...")
        ke_df, torque_df = create_comparison_dataframes(
            directories, filter_cohesion, filter_exp_type
        )
        
        if ke_df.empty or torque_df.empty:
            print("\n错误：无法创建比较数据")
            return
        
        # 可选：插值数据
        if use_interpolation:
            print("正在插值数据...")
            ke_df, torque_df = interpolate_data((ke_df, torque_df))
        
        # 保存比较文件
        save_comparison_files(ke_df, torque_df)
        
        # 输出统计信息
        print("\n" + "=" * 60)
        print("数据处理完成！")
        print(f"  - 处理了 {len(ke_df.columns) - 1} 个模拟结果")
        print(f"  - 时间范围: {ke_df['Time(s)'].min():.3f} - {ke_df['Time(s)'].max():.3f} 秒")
        print(f"  - 数据点数: {len(ke_df)}")
        
        # 输出各配置的统计（按实验类型分组显示）
        print("\n最大值统计:")
        
        # 分组显示
        standard_cols = [col for col in ke_df.columns[1:] if not col.startswith(('SizeDiff', 'DualDensity'))]
        sizediff_cols = [col for col in ke_df.columns[1:] if col.startswith('SizeDiff')]
        dualdensity_cols = [col for col in ke_df.columns[1:] if col.startswith('DualDensity')]
        
        for cols, label in [(standard_cols, "Standard"), 
                           (sizediff_cols, "SizeDiff"), 
                           (dualdensity_cols, "DualDensity")]:
            if cols:
                print(f"\n  {label}:")
                for col in sorted(cols):
                    max_ke = ke_df[col].max()
                    max_torque = torque_df[col].abs().max()
                    avg_ke = ke_df[col].mean()
                    avg_torque = torque_df[col].abs().mean()
                    print(f"    {col}:")
                    print(f"      - 动能: 最大={max_ke:.3e} J, 平均={avg_ke:.3e} J")
                    print(f"      - 扭矩: 最大={max_torque:.3e} Nm, 平均={avg_torque:.3e} Nm")

if __name__ == "__main__":
    main()
