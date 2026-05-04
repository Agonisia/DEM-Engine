#!/usr/bin/env python3
"""
SUP混合器模拟数据处理脚本 v2
==============================
提取和比较不同scale_factor下的动能、扭矩和配位数数据

更新内容 (v2):
  - 新增: 配位数 (AvgContacts) 提取和比较
  - 新增: 支持多个基础目录路径
  - 保持: 原有动能和扭矩提取功能

适配三种目录格式：
  - SUPMixerOutput_f{scale_factor}se{cohesion_code}
  - SUPMixerOutput_SizeDiff_f{scale_factor}se{cohesion_code}
  - SUPMixerOutput_DualDensity_f{scale_factor}se{cohesion_code}

输出结构: sta_results/{dir_type}_se{surface_energy}/
  - kinetic_energy_comparison.csv
  - mixer_torque_comparison.csv
  - avg_contacts_comparison.csv     [NEW]
  - torque_energy_metadata.txt
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
        return 5e-3
    else:
        return 1e-3


def find_sup_directories(base_paths):
    """
    从多个基础路径中查找所有SUP输出目录
    
    参数:
        base_paths: str 或 List[str]，一个或多个基础路径
    """
    if isinstance(base_paths, str):
        base_paths = [base_paths]
    
    directories = []
    
    patterns = [
        ('standard', re.compile(r'SUPMixerOutput_f(\d+(?:\.\d+)?)se(\d{3})$')),
        ('SizeDiff', re.compile(r'SUPMixerOutput_SizeDiff_f(\d+(?:\.\d+)?)se(\d{3})$')),
        ('DualDensity', re.compile(r'SUPMixerOutput_DualDensity_f(\d+(?:\.\d+)?)se(\d{3})$'))
    ]
    
    for base_path in base_paths:
        if not os.path.exists(base_path):
            print(f"警告：路径 {base_path} 不存在，跳过")
            continue
        
        for item in os.listdir(base_path):
            item_path = os.path.join(base_path, item)
            if os.path.isdir(item_path):
                for exp_type, pattern in patterns:
                    match = pattern.match(item)
                    if match:
                        scale_factor = float(match.group(1))
                        cohesion_code = match.group(2)
                        directories.append({
                            'path': item_path,
                            'name': item,
                            'scale_factor': scale_factor,
                            'surface_energy': cohesion_code,
                            'dir_type': exp_type,
                            'frame_dt': get_frame_dt(scale_factor),
                            'base_path': base_path
                        })
                        break
    
    # 去重：同一 (dir_type, surface_energy, scale_factor) 只保留第一个
    seen = set()
    unique_dirs = []
    for d in directories:
        key = (d['dir_type'], d['surface_energy'], d['scale_factor'])
        if key not in seen:
            seen.add(key)
            unique_dirs.append(d)
    
    # 排序
    type_order = {'standard': 0, 'SizeDiff': 1, 'DualDensity': 2}
    unique_dirs.sort(key=lambda x: (type_order.get(x['dir_type'], 99),
                                     x['surface_energy'],
                                     x['scale_factor']))
    
    return unique_dirs


def group_directories(directories: List[Dict]) -> Dict[Tuple[str, str], List[Dict]]:
    """按 (dir_type, surface_energy) 分组"""
    groups = {}
    for d in directories:
        key = (d['dir_type'], d['surface_energy'])
        if key not in groups:
            groups[key] = []
        groups[key].append(d)
    return groups


def read_history_file(directory_path):
    """
    读取 kinetic_energy_history.csv，提取动能、扭矩和配位数
    
    返回: DataFrame with columns [Time(s), KineticEnergy(J), MixerTorque(Nm), AvgContacts]
    """
    csv_path = os.path.join(directory_path, "kinetic_energy_history.csv")
    
    if not os.path.exists(csv_path):
        print(f"  警告：文件 {csv_path} 不存在")
        return None
    
    try:
        df = pd.read_csv(csv_path)
        
        # 统一列名 (支持新旧格式)
        if 'TotalKE(J)' in df.columns:
            df = df.rename(columns={'TotalKE(J)': 'KineticEnergy(J)'})
        
        # 检查必需列
        required = ['Time(s)', 'KineticEnergy(J)', 'MixerTorque(Nm)']
        if not all(col in df.columns for col in required):
            print(f"  警告：{csv_path} 缺少必需列: {required}")
            return None
        
        # AvgContacts 是可选的
        if 'AvgContacts' not in df.columns:
            print(f"  注意：{csv_path} 无 AvgContacts 列")
        
        return df
    except Exception as e:
        print(f"  错误：读取 {csv_path} 失败: {e}")
        return None


def process_group(group_dirs: List[Dict], auto_interpolate: bool = True,
                  n_points: int = 1000) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, Dict]:
    """
    处理一个分组，生成动能、扭矩、配位数的比较数据
    
    返回: (ke_df, torque_df, contacts_df, metadata)
    """
    raw_data = {}
    metadata = {}
    
    for dir_info in group_dirs:
        df = read_history_file(dir_info['path'])
        if df is None:
            continue
        
        sf = dir_info['scale_factor']
        col_name = f"f{int(sf) if sf == int(sf) else sf}"
        
        data_entry = {
            'time': df['Time(s)'].values,
            'ke': df['KineticEnergy(J)'].values,
            'torque': df['MixerTorque(Nm)'].values,
        }
        # AvgContacts 可选
        if 'AvgContacts' in df.columns:
            data_entry['contacts'] = df['AvgContacts'].values
        
        raw_data[col_name] = data_entry
        metadata[col_name] = {
            'scale_factor': dir_info['scale_factor'],
            'surface_energy': dir_info['surface_energy'],
            'dir_type': dir_info['dir_type'],
            'frame_dt': dir_info['frame_dt'],
            'n_points': len(df),
            'has_contacts': 'AvgContacts' in df.columns
        }
    
    if not raw_data:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), metadata
    
    # 构建统一时间网格
    lengths = [len(d['time']) for d in raw_data.values()]
    needs_interp = len(set(lengths)) > 1
    
    if not needs_interp:
        print(f"  ✓ 数据长度一致 ({lengths[0]} 点)")
        first_key = list(raw_data.keys())[0]
        common_time = raw_data[first_key]['time']
    elif auto_interpolate:
        print(f"  ⚠ 数据长度不一致: {dict(zip(raw_data.keys(), lengths))}")
        print(f"  ✓ 自动插值到 {n_points} 点")
        time_min = max(d['time'].min() for d in raw_data.values())
        time_max = min(d['time'].max() for d in raw_data.values())
        common_time = np.linspace(time_min, time_max, n_points)
    else:
        # 用最长序列
        max_len = max(lengths)
        max_key = [k for k, d in raw_data.items() if len(d['time']) == max_len][0]
        common_time = raw_data[max_key]['time']
    
    # 构建各 DataFrame
    ke_data = {'Time(s)': common_time}
    torque_data = {'Time(s)': common_time}
    contacts_data = {'Time(s)': common_time}
    has_any_contacts = False
    
    for col_name, data in sorted(raw_data.items()):
        if needs_interp and auto_interpolate:
            ke_data[col_name] = np.interp(common_time, data['time'], data['ke'])
            torque_data[col_name] = np.interp(common_time, data['time'], data['torque'])
            if 'contacts' in data:
                contacts_data[col_name] = np.interp(common_time, data['time'], data['contacts'])
                has_any_contacts = True
        elif not needs_interp:
            ke_data[col_name] = data['ke']
            torque_data[col_name] = data['torque']
            if 'contacts' in data:
                contacts_data[col_name] = data['contacts']
                has_any_contacts = True
        else:
            n = len(common_time)
            dn = len(data['ke'])
            if dn == n:
                ke_data[col_name] = data['ke']
                torque_data[col_name] = data['torque']
                if 'contacts' in data:
                    contacts_data[col_name] = data['contacts']
                    has_any_contacts = True
            else:
                pad_ke = np.full(n, np.nan)
                pad_torque = np.full(n, np.nan)
                pad_ke[:dn] = data['ke']
                pad_torque[:dn] = data['torque']
                ke_data[col_name] = pad_ke
                torque_data[col_name] = pad_torque
                if 'contacts' in data:
                    pad_contacts = np.full(n, np.nan)
                    pad_contacts[:dn] = data['contacts']
                    contacts_data[col_name] = pad_contacts
                    has_any_contacts = True
    
    ke_df = pd.DataFrame(ke_data)
    torque_df = pd.DataFrame(torque_data)
    contacts_df = pd.DataFrame(contacts_data) if has_any_contacts else pd.DataFrame()
    
    return ke_df, torque_df, contacts_df, metadata


def save_group_results(ke_df, torque_df, contacts_df, metadata,
                       dir_type, surface_energy, base_output="sta_results"):
    """保存分组结果"""
    output_dir = Path(base_output) / f"{dir_type}_se{surface_energy}"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    ke_df.to_csv(output_dir / "kinetic_energy_comparison.csv", index=False, float_format='%.6e')
    torque_df.to_csv(output_dir / "mixer_torque_comparison.csv", index=False, float_format='%.6e')
    
    if not contacts_df.empty:
        contacts_df.to_csv(output_dir / "avg_contacts_comparison.csv", index=False, float_format='%.6e')
    
    # 元数据
    with open(output_dir / "torque_energy_metadata.txt", 'w') as f:
        f.write(f"SUP Mixer Data Comparison\n")
        f.write(f"Directory Type: {dir_type}\n")
        f.write(f"Surface Energy: se{surface_energy}\n")
        f.write("=" * 50 + "\n\n")
        
        for col_name, meta in sorted(metadata.items(), key=lambda x: x[1]['scale_factor']):
            f.write(f"{col_name}:\n")
            f.write(f"  Scale Factor: f{meta['scale_factor']}\n")
            f.write(f"  Frame Interval: {meta['frame_dt']*1000:.1f} ms\n")
            f.write(f"  Data Points: {meta['n_points']}\n")
            f.write(f"  Has AvgContacts: {meta.get('has_contacts', False)}\n")
            
            if col_name in ke_df.columns:
                f.write(f"  KE Max: {ke_df[col_name].max():.4e} J\n")
                f.write(f"  KE Mean: {ke_df[col_name].mean():.4e} J\n")
                f.write(f"  Torque Max: {torque_df[col_name].abs().max():.4e} Nm\n")
                f.write(f"  Torque Mean: {torque_df[col_name].abs().mean():.4e} Nm\n")
            if not contacts_df.empty and col_name in contacts_df.columns:
                f.write(f"  AvgContacts Mean: {contacts_df[col_name].mean():.4f}\n")
            f.write("\n")
    
    return output_dir


def print_directories_info(directories):
    """打印目录信息"""
    print(f"\n找到 {len(directories)} 个SUP输出目录:")
    groups = group_directories(directories)
    for (dir_type, se), dirs in sorted(groups.items()):
        print(f"\n  {dir_type}_se{se}:")
        for d in sorted(dirs, key=lambda x: x['scale_factor']):
            print(f"    - f{d['scale_factor']} (dt={d['frame_dt']*1000:.1f}ms) [{d['base_path']}]")


def main():
    """主函数"""
    print("=" * 60)
    print("SUP混合器数据处理脚本 v2")
    print("提取: 动能 + 扭矩 + 配位数")
    print("输出: sta_results/{dir_type}_se{surface_energy}/")
    print("=" * 60)
    
    # ==================== 可配置参数 ====================
    # 支持多个目录路径，会自动合并去重
    base_paths = [
        "/media/huyuze/Fanxiang1/DEME_Data/backup251115",
        # "/media/huyuze/Fanxiang/backup251025",
        # "/media/huyuze/Fanxiang/backup251109"
        # "/path/to/another/data/directory",    # 取消注释添加更多路径
        # "/path/to/yet/another/directory",
    ]
    
    filter_cohesion = "015"        # 只处理特定cohesion (如 "000")，None=处理所有
    filter_type = None            # 只处理特定类型 (如 "SizeDiff")，None=处理所有
    auto_interpolate = True       # 数据长度不一致时自动插值
    interpolate_points = 1000     # 插值点数
    base_output = "sta_results"   # 基础输出目录
    # ====================================================
    
    # 1. 查找目录
    directories = find_sup_directories(base_paths)
    
    if not directories:
        print("\n错误：未找到任何SUP输出目录")
        return
    
    print_directories_info(directories)
    
    # 2. 配置
    print(f"\n分析配置:")
    print(f"  - 基础路径: {base_paths}")
    print(f"  - Cohesion过滤: {filter_cohesion or '无(处理所有)'}")
    print(f"  - 类型过滤: {filter_type or '无(处理所有)'}")
    print(f"  - 自动插值: {'是' if auto_interpolate else '否'}")
    print(f"  - 输出目录: {base_output}/")
    
    # 3. 分组处理
    groups = group_directories(directories)
    processed_count = 0
    
    for (dir_type, surface_energy), group_dirs in sorted(groups.items()):
        if filter_cohesion and surface_energy != filter_cohesion:
            continue
        if filter_type and dir_type != filter_type:
            continue
        
        print(f"\n{'='*60}")
        print(f"处理: {dir_type}_se{surface_energy} ({len(group_dirs)} scale factors)")
        print('='*60)
        
        ke_df, torque_df, contacts_df, metadata = process_group(
            group_dirs, auto_interpolate=auto_interpolate, n_points=interpolate_points
        )
        
        if ke_df.empty:
            print(f"  ✗ 无有效数据，跳过")
            continue
        
        output_dir = save_group_results(
            ke_df, torque_df, contacts_df, metadata,
            dir_type, surface_energy, base_output
        )
        
        print(f"  ✓ 已保存到: {output_dir}")
        print(f"    时间范围: {ke_df['Time(s)'].min():.3f} - {ke_df['Time(s)'].max():.3f} s")
        print(f"    数据点数: {len(ke_df)}")
        
        for col in sorted(ke_df.columns[1:], key=lambda x: metadata.get(x, {}).get('scale_factor', 0)):
            if col in metadata:
                meta = metadata[col]
                contacts_info = ""
                if not contacts_df.empty and col in contacts_df.columns:
                    contacts_info = f", Contacts={contacts_df[col].mean():.2f}"
                print(f"    {col}: KE_mean={ke_df[col].mean():.3e}J{contacts_info}")
        
        processed_count += 1
    
    print("\n" + "=" * 60)
    print(f"完成! 处理了 {processed_count} 个分组 → {base_output}/")
    print("=" * 60)


if __name__ == "__main__":
    main()
    