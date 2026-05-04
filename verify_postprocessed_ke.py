#!/usr/bin/env python3
"""
SUP 后处理动能一致性核算脚本
用于验证: sta_results 中的总动能 == extended_features_results 中的 (平动动能 + 转动动能)
"""

import os
import pandas as pd
import numpy as np
from pathlib import Path

def verify_processed_data(sta_base="sta_results", ext_base="extended_features_results"):
    sta_path = Path(sta_base)
    ext_path = Path(ext_base)

    if not sta_path.exists() or not ext_path.exists():
        print(f"❌ 找不到数据目录。请确保 {sta_base}/ 和 {ext_base}/ 存在于当前目录。")
        return

    # 遍历所有的实验组 (例如 SizeDiff_se000, SizeDiff_se015 等)
    for group_dir in sta_path.iterdir():
        if not group_dir.is_dir():
            continue
            
        group_name = group_dir.name
        print(f"\n{'='*70}")
        print(f"正在核算分组: {group_name}")
        print(f"{'='*70}")

        # 对应的宏观和微观文件路径
        sta_ke_file = group_dir / "kinetic_energy_comparison.csv"
        ext_trans_file = ext_path / group_name / "comparison_KE_trans.csv"
        ext_rot_file = ext_path / group_name / "comparison_KE_rot.csv"

        if not (sta_ke_file.exists() and ext_trans_file.exists() and ext_rot_file.exists()):
            print(f"⚠️ 分组 {group_name} 缺少必要的比对文件，跳过。")
            continue

        # 读取数据
        df_sta = pd.read_csv(sta_ke_file)
        df_trans = pd.read_csv(ext_trans_file)
        df_rot = pd.read_csv(ext_rot_file)

        # 找出共同存在的 scale factor 列 (例如 'f1.0', 'f2.0', 'f4.0')
        cols_sta = set(df_sta.columns) - {'Time(s)'}
        cols_ext = set(df_trans.columns) - {'Time(s)'}
        common_cols = sorted(list(cols_sta & cols_ext))

        if not common_cols:
            print("⚠️ 没有找到共同的 scale factor 列，跳过。")
            continue

        for col in common_cols:
            print(f"\n>>> 校验 Scale Factor: {col}")
            
            # 以 sta_results 的时间轴为基准
            time_ref = df_sta['Time(s)'].values
            ke_total_sta = df_sta[col].values

            # 微观特征的数据时间轴可能与 sta 不一致，需要插值对齐
            ke_trans_ext = np.interp(time_ref, df_trans['Time(s)'].values, df_trans[col].values)
            ke_rot_ext = np.interp(time_ref, df_rot['Time(s)'].values, df_rot[col].values)

            # 计算总和
            ke_total_ext_sum = ke_trans_ext + ke_rot_ext

            # 计算相对误差 (排除分母为0的情况)
            valid_mask = ke_total_sta > 1e-10
            if not np.any(valid_mask):
                print("  ⚠️ 该列数据全为0或极小值，无法计算相对误差。")
                continue

            abs_error = np.abs(ke_total_ext_sum[valid_mask] - ke_total_sta[valid_mask])
            rel_error = abs_error / ke_total_sta[valid_mask]

            max_rel_err = np.max(rel_error)
            mean_rel_err = np.mean(rel_error)

            # 打印几个采样点展示
            print(f"  {'Time(s)':<10} | {'sta_KE_total':<15} | {'ext_KE_trans':<15} | {'ext_KE_rot':<15} | {'ext_Sum':<15} | {'Rel_Error':<10}")
            print("  " + "-"*90)
            
            # 均匀抽取5个点展示
            sample_indices = np.linspace(0, len(time_ref)-1, 5, dtype=int)
            for idx in sample_indices:
                t = time_ref[idx]
                sta_val = ke_total_sta[idx]
                trans_val = ke_trans_ext[idx]
                rot_val = ke_rot_ext[idx]
                sum_val = ke_total_ext_sum[idx]
                err_val = abs(sum_val - sta_val)/sta_val if sta_val > 1e-10 else 0.0
                print(f"  {t:<10.3f} | {sta_val:<15.6e} | {trans_val:<15.6e} | {rot_val:<15.6e} | {sum_val:<15.6e} | {err_val:<10.2%}")

            print("  " + "-"*90)
            if max_rel_err < 0.01: # 允许 1% 的插值误差
                print(f"  ✅ 完美对账! 最大相对误差: {max_rel_err:.4%} (平均: {mean_rel_err:.4%})")
            else:
                print(f"  ❌ 存在较大偏差! 最大相对误差: {max_rel_err:.4%} (平均: {mean_rel_err:.4%})")

if __name__ == "__main__":
    verify_processed_data()