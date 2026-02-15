#!/usr/bin/env python3
"""
SUP混合器模拟数据处理脚本
用于提取和比较不同scale_factor下的动能和扭矩数据
适配三种目录格式：
  - SUPMixerOutput_f{scale_factor}se{cohesion_code}
  - SUPMixerOutput_SizeDiff_f{scale_factor}se{cohesion_code}
  - SUPMixerOutput_DualDensity_f{scale_factor}se{cohesion_code}
Fixed: 添加对f1帧时间间隔的特殊处理
Updated: 输出结构与PDF分析脚本保持一致
"""

import os
import re
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Optional

def get_frame_dt(scale_factor):
    """根据scale_factor获取帧时间间隔"""
    if scale_factor == 1:
        return 5e-3  # f1 using 5ms interval
    else:
        return 1e-3  # others using 1ms interval

def find_sup_directories(base_path="."):
    """
    查找所有符合SUPMixerOutput格式的目录（支持三种格式）
    
    返回:
        包含目录信息字典的列表
    """
    directories = []
    
    if not os.path.exists(base_path):
        print(f"警告：路径 {base_path} 不存在")
        return directories
    
    # 定义三种模式
    patterns = [
        ('standard', re.compile(r'SUPMixerOutput_f(\d+(?:\.\d+)?)se(\d{3})$')),
        ('SizeDiff', re.compile(r'SUPMixerOutput_SizeDiff_f(\d+(?:\.\d+)?)se(\d{3})$')),
        ('DualDensity', re.compile(r'SUPMixerOutput_DualDensity_f(\d+(?:\.\d+)?)se(\d{3})$'))
    ]
    
    for item in os.listdir(base_path):
        item_path = os.path.join(base_path, item)
        if os.path.isdir(item_path):
            for exp_type, pattern in patterns:
                match = pattern.match(item)
                if match:
                    scale_factor = float(match.group(1))
                    cohesion_code = match.group(2)
                    frame_dt = get_frame_dt(scale_factor)
                    directories.append({
                        'path': item_path,
                        'name': item,
                        'scale_factor': scale_factor,
                        'surface_energy': cohesion_code,
                        'dir_type': exp_type,
                        'frame_dt': frame_dt
                    })
                    break
    
    # 按dir_type, surface_energy, scale_factor排序
    type_order = {'standard': 0, 'SizeDiff': 1, 'DualDensity': 2}
    directories.sort(key=lambda x: (type_order.get(x['dir_type'], 99), 
                                    x['surface_energy'], 
                                    x['scale_factor']))
    
    return directories

def group_directories(directories: List[Dict]) -> Dict[Tuple[str, str], List[Dict]]:
    """
    按 (dir_type, surface_energy) 分组目录
    
    返回:
        {(dir_type, surface_energy): [dir_info, ...], ...}
    """
    groups = {}
    for d in directories:
        key = (d['dir_type'], d['surface_energy'])
        if key not in groups:
            groups[key] = []
        groups[key].append(d)
    return groups

def read_kinetic_energy_file(directory_path):
    """读取指定目录下的kinetic_energy_history.csv文件"""
    csv_path = os.path.join(directory_path, "kinetic_energy_history.csv")
    
    if not os.path.exists(csv_path):
        print(f"  警告：文件 {csv_path} 不存在")
        return None
    
    try:
        df = pd.read_csv(csv_path)
        
        # 检查必需的列是否存在（支持新旧两种格式）
        has_new_format = all(col in df.columns for col in ['Time(s)', 'TotalKE(J)', 'MixerTorque(Nm)'])
        has_old_format = all(col in df.columns for col in ['Time(s)', 'KineticEnergy(J)', 'MixerTorque(Nm)'])
        
        if has_new_format:
            df = df.rename(columns={'TotalKE(J)': 'KineticEnergy(J)'})
        elif not has_old_format:
            print(f"  警告：文件 {csv_path} 缺少必需的列")
            return None
            
        return df
    except Exception as e:
        print(f"  错误：读取文件 {csv_path} 时出错: {e}")
        return None

def process_group(group_dirs: List[Dict], auto_interpolate: bool = True, n_points: int = 1000) -> Tuple[pd.DataFrame, pd.DataFrame, Dict]:
    """
    处理一个分组内的所有目录，生成比较数据
    
    参数:
        group_dirs: 目录信息列表
        auto_interpolate: 当数据长度不一致时自动插值
        n_points: 插值点数
    
    返回:
        (动能DataFrame, 扭矩DataFrame, 元数据字典)
    """
    # 先读取所有数据
    raw_data = {}
    metadata = {}
    
    for dir_info in group_dirs:
        df = read_kinetic_energy_file(dir_info['path'])
        if df is None:
            continue
        
        # 列名只用 f{scale_factor}
        sf = dir_info['scale_factor']
        col_name = f"f{int(sf) if sf == int(sf) else sf}"
        
        raw_data[col_name] = {
            'time': df['Time(s)'].values,
            'ke': df['KineticEnergy(J)'].values,
            'torque': df['MixerTorque(Nm)'].values
        }
        
        metadata[col_name] = {
            'scale_factor': dir_info['scale_factor'],
            'surface_energy': dir_info['surface_energy'],
            'dir_type': dir_info['dir_type'],
            'frame_dt': dir_info['frame_dt'],
            'n_points': len(df)
        }
    
    if not raw_data:
        return pd.DataFrame(), pd.DataFrame(), metadata
    
    # 检查数据长度是否一致
    lengths = [len(d['time']) for d in raw_data.values()]
    lengths_consistent = len(set(lengths)) == 1
    
    if lengths_consistent:
        # 长度一致，直接构建 DataFrame
        print(f"  ✓ 数据长度一致 ({lengths[0]} 点)，直接合并")
        first_key = list(raw_data.keys())[0]
        kinetic_energy_data = {'Time(s)': raw_data[first_key]['time']}
        mixer_torque_data = {'Time(s)': raw_data[first_key]['time']}
        
        for col_name, data in raw_data.items():
            kinetic_energy_data[col_name] = data['ke']
            mixer_torque_data[col_name] = data['torque']
    else:
        # 长度不一致
        print(f"  ⚠ 数据长度不一致: {dict(zip(raw_data.keys(), lengths))}")
        
        if auto_interpolate:
            print(f"  ✓ 自动插值到统一时间网格 ({n_points} 点)")
            
            # 找到公共时间范围
            time_min = max(d['time'].min() for d in raw_data.values())
            time_max = min(d['time'].max() for d in raw_data.values())
            common_time = np.linspace(time_min, time_max, n_points)
            
            kinetic_energy_data = {'Time(s)': common_time}
            mixer_torque_data = {'Time(s)': common_time}
            
            for col_name, data in raw_data.items():
                kinetic_energy_data[col_name] = np.interp(common_time, data['time'], data['ke'])
                mixer_torque_data[col_name] = np.interp(common_time, data['time'], data['torque'])
        else:
            # 不插值，用 NaN 填充
            print(f"  ✓ 使用最长时间序列，短数据用 NaN 填充")
            max_len = max(lengths)
            max_key = [k for k, d in raw_data.items() if len(d['time']) == max_len][0]
            
            kinetic_energy_data = {'Time(s)': raw_data[max_key]['time']}
            mixer_torque_data = {'Time(s)': raw_data[max_key]['time']}
            
            for col_name, data in raw_data.items():
                if len(data['ke']) == max_len:
                    kinetic_energy_data[col_name] = data['ke']
                    mixer_torque_data[col_name] = data['torque']
                else:
                    # 用 NaN 填充
                    ke_padded = np.full(max_len, np.nan)
                    torque_padded = np.full(max_len, np.nan)
                    ke_padded[:len(data['ke'])] = data['ke']
                    torque_padded[:len(data['torque'])] = data['torque']
                    kinetic_energy_data[col_name] = ke_padded
                    mixer_torque_data[col_name] = torque_padded
    
    ke_df = pd.DataFrame(kinetic_energy_data)
    torque_df = pd.DataFrame(mixer_torque_data)
    
    return ke_df, torque_df, metadata

def interpolate_data(ke_df, torque_df, metadata, n_points=1000):
    """插值数据以处理不同时间步长的情况"""
    if ke_df.empty or torque_df.empty:
        return ke_df, torque_df
    
    time_min = ke_df['Time(s)'].min()
    time_max = ke_df['Time(s)'].max()
    common_time = np.linspace(time_min, time_max, n_points)
    
    ke_interpolated = {'Time(s)': common_time}
    torque_interpolated = {'Time(s)': common_time}
    
    for col in ke_df.columns[1:]:
        ke_interpolated[col] = np.interp(common_time, ke_df['Time(s)'], ke_df[col])
        torque_interpolated[col] = np.interp(common_time, torque_df['Time(s)'], torque_df[col])
    
    return pd.DataFrame(ke_interpolated), pd.DataFrame(torque_interpolated)

def save_group_results(ke_df, torque_df, metadata, dir_type, surface_energy, base_output="sta_results"):
    """
    保存分组结果到对应目录
    
    输出结构: sta_results/{dir_type}_se{surface_energy}/
    """
    # 构建输出路径 - 与PDF脚本保持一致
    output_dir = Path(base_output) / f"{dir_type}_se{surface_energy}"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 保存数据文件
    ke_df.to_csv(output_dir / "kinetic_energy_comparison.csv", index=False, float_format='%.6e')
    torque_df.to_csv(output_dir / "mixer_torque_comparison.csv", index=False, float_format='%.6e')
    
    # 保存元数据
    with open(output_dir / "torque_energy_metadata.txt", 'w') as f:
        f.write(f"SUP Mixer Torque & Energy Comparison\n")
        f.write(f"Directory Type: {dir_type}\n")
        f.write(f"Surface Energy: se{surface_energy}\n")
        f.write("=" * 50 + "\n\n")
        
        for col_name, meta in sorted(metadata.items(), key=lambda x: x[1]['scale_factor']):
            f.write(f"{col_name}:\n")
            f.write(f"  Scale Factor: f{meta['scale_factor']}\n")
            f.write(f"  Frame Interval: {meta['frame_dt']*1000:.1f} ms\n")
            f.write(f"  Data Points: {meta['n_points']}\n")
            
            # 添加统计信息
            if col_name in ke_df.columns:
                f.write(f"  KE Max: {ke_df[col_name].max():.4e} J\n")
                f.write(f"  KE Mean: {ke_df[col_name].mean():.4e} J\n")
                f.write(f"  Torque Max: {torque_df[col_name].abs().max():.4e} Nm\n")
                f.write(f"  Torque Mean: {torque_df[col_name].abs().mean():.4e} Nm\n")
            f.write("\n")
    
    return output_dir

def print_directories_info(directories):
    """打印目录信息"""
    print(f"\n找到 {len(directories)} 个SUP输出目录:")
    
    # 按类型分组显示
    groups = group_directories(directories)
    for (dir_type, se), dirs in sorted(groups.items()):
        print(f"\n  {dir_type}_se{se}:")
        for d in sorted(dirs, key=lambda x: x['scale_factor']):
            print(f"    - f{d['scale_factor']} (dt={d['frame_dt']*1000:.1f}ms)")

def main():
    """主函数"""
    print("=" * 60)
    print("SUP混合器数据处理脚本 - 统一输出结构版本")
    print("输出结构: sta_results/{dir_type}_se{surface_energy}/")
    print("=" * 60)
    
    # ==================== 可配置参数 ====================
    base_path = "/media/huyuze/Fanxiang/DEME_Data/backup251115"
    filter_cohesion = "000"       # 只处理特定cohesion，None表示处理所有
    filter_type = None            # 只处理特定类型，None表示处理所有
    auto_interpolate = True       # 数据长度不一致时自动插值
    interpolate_points = 1000     # 插值点数
    base_output = "sta_results"   # 基础输出目录
    # ====================================================
    
    # 1. 查找所有SUP输出目录
    directories = find_sup_directories(base_path)
    
    if not directories:
        print("\n错误：未找到任何SUP输出目录")
        return
    
    print_directories_info(directories)
    
    # 2. 显示配置
    print(f"\n分析配置:")
    print(f"  - 基础路径: {base_path}")
    print(f"  - Cohesion过滤: {filter_cohesion if filter_cohesion else '无(处理所有)'}")
    print(f"  - 类型过滤: {filter_type if filter_type else '无(处理所有)'}")
    print(f"  - 自动插值: {'是' if auto_interpolate else '否'}")
    print(f"  - 插值点数: {interpolate_points}")
    print(f"  - 输出目录: {base_output}/")
    
    # 3. 按 (dir_type, surface_energy) 分组处理
    groups = group_directories(directories)
    processed_count = 0
    
    for (dir_type, surface_energy), group_dirs in sorted(groups.items()):
        # 应用过滤条件
        if filter_cohesion and surface_energy != filter_cohesion:
            continue
        if filter_type and dir_type != filter_type:
            continue
        
        print(f"\n{'='*60}")
        print(f"处理: {dir_type}_se{surface_energy}")
        print(f"  包含 {len(group_dirs)} 个scale factors")
        print('='*60)
        
        # 处理该分组
        ke_df, torque_df, metadata = process_group(
            group_dirs, 
            auto_interpolate=auto_interpolate,
            n_points=interpolate_points
        )
        
        if ke_df.empty or torque_df.empty:
            print(f"  ✗ 无有效数据，跳过")
            continue
        
        # 保存结果
        output_dir = save_group_results(
            ke_df, torque_df, metadata, 
            dir_type, surface_energy, 
            base_output
        )
        
        print(f"  ✓ 结果已保存到: {output_dir}")
        
        # 显示统计
        print(f"\n  统计信息:")
        print(f"    时间范围: {ke_df['Time(s)'].min():.3f} - {ke_df['Time(s)'].max():.3f} s")
        print(f"    数据点数: {len(ke_df)}")
        for col in sorted(ke_df.columns[1:], key=lambda x: metadata.get(x, {}).get('scale_factor', 0)):
            if col in metadata:
                meta = metadata[col]
                print(f"    {col} (dt={meta['frame_dt']*1000:.1f}ms):")
                print(f"      KE: max={ke_df[col].max():.3e}J, mean={ke_df[col].mean():.3e}J")
                print(f"      Torque: max={torque_df[col].abs().max():.3e}Nm")
        
        processed_count += 1
    
    # 4. 总结
    print("\n" + "=" * 60)
    print("处理完成!")
    print(f"  - 处理了 {processed_count} 个分组")
    print(f"  - 输出目录: {base_output}/")
    print("=" * 60)

if __name__ == "__main__":
    main()