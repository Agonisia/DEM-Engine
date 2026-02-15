#!/usr/bin/env python3
"""
SUP Mixer Extended Features Analysis
=====================================
扩展特征分析脚本，计算完整的动力学特征集用于SUP模型训练

特征列表:
---------
A. 不需要柱坐标分解（直接计算）:
   - KE_trans_sum: 平动动能代理 Σ(vx²+vy²+vz²)
   - KE_rot_sum: 转动动能代理 Σ(ωx²+ωy²+ωz²)
   - vx_mean, vy_mean, vz_mean: 各方向平均速度
   - speed_mean, speed_std: 平均速率和标准差
   - PE_sum: 重力势能代理 g·Σz
   - com_z: 质心高度 mean(Z)

B. 需要柱坐标分解:
   - v_r_mean, v_theta_mean, v_z_mean: 柱坐标平均速度
   - v_r_std, v_theta_std, v_z_std: 柱坐标速度标准差
   - omega_mean: 平均角速度 median(v_theta/r)
   - T_r, T_theta, T_z, T_total: 粒子温度(granular temperature)

输出结构: sta_results/{dir_type}_se{surface_energy}/extended_features.csv
"""

import os
import re
import glob
import pandas as pd
import numpy as np
from numba import jit, prange
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
import warnings
warnings.filterwarnings('ignore')

# ======================== Constants ========================
G = 9.81  # m/s², 重力加速度

# ======================== Configuration ========================
@dataclass
class AnalysisConfig:
    """分析配置参数"""
    base_path: str = "."
    filter_cohesion: Optional[str] = None  # None = 处理所有
    filter_type: Optional[str] = None      # None = 处理所有
    steady_state_start_time: float = 0.5   # 从0.5秒开始作为稳态
    frame_skip: int = 100                  # 每100帧取1帧
    output_base: str = "sta_results"
    
    def get_frame_dt(self, scale_factor: int) -> float:
        """根据scale_factor获取帧时间间隔"""
        return 5e-3 if scale_factor == 1 else 1e-3


# ======================== JIT Accelerated Functions ========================
@jit(nopython=True, parallel=True)
def compute_cylindrical_velocities(x, y, vx, vy):
    """
    计算柱坐标速度分量
    v_r = (x·vx + y·vy) / r
    v_theta = (-y·vx + x·vy) / r
    """
    n = len(x)
    v_r = np.zeros(n)
    v_theta = np.zeros(n)
    r_dist = np.zeros(n)
    
    for i in prange(n):
        r = np.sqrt(x[i]**2 + y[i]**2)
        r_dist[i] = r
        if r > 1e-10:  # 避免除零
            v_r[i] = (x[i] * vx[i] + y[i] * vy[i]) / r
            v_theta[i] = (-y[i] * vx[i] + x[i] * vy[i]) / r
        else:
            v_r[i] = 0.0
            v_theta[i] = 0.0
    
    return v_r, v_theta, r_dist


@jit(nopython=True)
def compute_granular_temperature(v_r, v_theta, v_z, r_dist):
    """
    计算粒子温度 (granular temperature)
    T_r = mean((v_r - v_r_mean)²)
    T_theta = mean((v_theta - ω·r)²)  # 扣除整体旋转
    T_z = mean((v_z - v_z_mean)²)
    """
    n = len(v_r)
    
    # 计算平均值
    v_r_mean = np.mean(v_r)
    v_z_mean = np.mean(v_z)
    
    # 估计整体角速度 (中位数更稳健)
    omega_vals = np.zeros(n)
    valid_count = 0
    for i in range(n):
        if r_dist[i] > 1e-10:
            omega_vals[valid_count] = v_theta[i] / r_dist[i]
            valid_count += 1
    
    omega_mean = np.median(omega_vals[:valid_count]) if valid_count > 0 else 0.0
    
    # 计算各方向温度
    T_r = 0.0
    T_theta = 0.0
    T_z = 0.0
    
    for i in range(n):
        T_r += (v_r[i] - v_r_mean) ** 2
        T_theta += (v_theta[i] - omega_mean * r_dist[i]) ** 2
        T_z += (v_z[i] - v_z_mean) ** 2
    
    T_r /= n
    T_theta /= n
    T_z /= n
    T_total = (T_r + T_theta + T_z) / 3.0
    
    return T_r, T_theta, T_z, T_total, omega_mean


@jit(nopython=True)
def compute_basic_stats(vx, vy, vz, wx, wy, wz, z):
    """计算基本统计量（数量加权）"""
    n = len(vx)
    
    # 速度平方和 (动能代理)
    KE_trans_sum = 0.0
    KE_rot_sum = 0.0
    PE_sum = 0.0
    
    speed_arr = np.zeros(n)
    
    for i in range(n):
        KE_trans_sum += vx[i]**2 + vy[i]**2 + vz[i]**2
        KE_rot_sum += wx[i]**2 + wy[i]**2 + wz[i]**2
        PE_sum += z[i]
        speed_arr[i] = np.sqrt(vx[i]**2 + vy[i]**2 + vz[i]**2)
    
    PE_sum *= G  # 乘以重力加速度
    
    # 速度统计
    vx_mean = np.mean(vx)
    vy_mean = np.mean(vy)
    vz_mean = np.mean(vz)
    speed_mean = np.mean(speed_arr)
    speed_std = np.std(speed_arr)
    com_z = np.mean(z)
    
    return (KE_trans_sum, KE_rot_sum, PE_sum, 
            vx_mean, vy_mean, vz_mean, 
            speed_mean, speed_std, com_z)


# ======================== Directory Discovery ========================
def find_sup_directories(base_path: str) -> List[Dict]:
    """查找所有SUP输出目录"""
    directories = []
    base = Path(base_path)
    
    if not base.exists():
        print(f"警告：路径 {base_path} 不存在")
        return directories
    
    patterns = [
        ('standard', re.compile(r'SUPMixerOutput_f(\d+)se(\d{3})$')),
        ('SizeDiff', re.compile(r'SUPMixerOutput_SizeDiff_f(\d+)se(\d{3})$')),
        ('DualDensity', re.compile(r'SUPMixerOutput_DualDensity_f(\d+)se(\d{3})$'))
    ]
    
    for item in base.iterdir():
        if not item.is_dir():
            continue
        for dir_type, pattern in patterns:
            match = pattern.match(item.name)
            if match:
                directories.append({
                    'path': item,
                    'name': item.name,
                    'scale_factor': int(match.group(1)),
                    'surface_energy': match.group(2),
                    'dir_type': dir_type
                })
                break
    
    # 排序：类型 -> 表面能 -> 缩放因子
    type_order = {'standard': 0, 'SizeDiff': 1, 'DualDensity': 2}
    directories.sort(key=lambda x: (
        type_order.get(x['dir_type'], 99),
        x['surface_energy'],
        x['scale_factor']
    ))
    
    return directories


def group_directories(directories: List[Dict]) -> Dict[Tuple[str, str], List[Dict]]:
    """按 (dir_type, surface_energy) 分组"""
    groups = {}
    for d in directories:
        key = (d['dir_type'], d['surface_energy'])
        if key not in groups:
            groups[key] = []
        groups[key].append(d)
    return groups


# ======================== Data Processing ========================
class FrameFeatureExtractor:
    """单帧特征提取器"""
    
    @staticmethod
    def extract_features(df: pd.DataFrame) -> Dict[str, float]:
        """从单帧数据提取所有特征"""
        
        # 提取数组
        x = df['X'].values
        y = df['Y'].values
        z = df['Z'].values
        vx = df['v_x'].values
        vy = df['v_y'].values
        vz = df['v_z'].values
        wx = df['w_x'].values
        wy = df['w_y'].values
        wz = df['w_z'].values
        
        # A. 基本统计量
        (KE_trans_sum, KE_rot_sum, PE_sum,
         vx_mean, vy_mean, vz_mean,
         speed_mean, speed_std, com_z) = compute_basic_stats(
            vx, vy, vz, wx, wy, wz, z
        )
        
        # B. 柱坐标速度
        v_r, v_theta, r_dist = compute_cylindrical_velocities(x, y, vx, vy)
        
        v_r_mean = np.mean(v_r)
        v_theta_mean = np.mean(v_theta)
        v_z_mean = np.mean(vz)
        v_r_std = np.std(v_r)
        v_theta_std = np.std(v_theta)
        v_z_std = np.std(vz)
        
        # C. 粒子温度
        T_r, T_theta, T_z, T_total, omega_mean = compute_granular_temperature(
            v_r, v_theta, vz, r_dist
        )
        
        # 粒子数
        n_particles = len(df)
        
        return {
            # 基本统计
            'n_particles': n_particles,
            'KE_trans_sum': KE_trans_sum,
            'KE_rot_sum': KE_rot_sum,
            'PE_sum': PE_sum,
            'vx_mean': vx_mean,
            'vy_mean': vy_mean,
            'vz_mean': vz_mean,
            'speed_mean': speed_mean,
            'speed_std': speed_std,
            'com_z': com_z,
            # 柱坐标速度
            'v_r_mean': v_r_mean,
            'v_theta_mean': v_theta_mean,
            'v_z_cyl_mean': v_z_mean,  # 与vz_mean相同，但明确是柱坐标
            'v_r_std': v_r_std,
            'v_theta_std': v_theta_std,
            'v_z_cyl_std': v_z_std,
            'omega_mean': omega_mean,
            # 粒子温度
            'T_r': T_r,
            'T_theta': T_theta,
            'T_z': T_z,
            'T_total': T_total,
        }


class DirectoryProcessor:
    """单目录处理器"""
    
    def __init__(self, dir_info: Dict, config: AnalysisConfig):
        self.dir_info = dir_info
        self.config = config
        self.frame_dt = config.get_frame_dt(dir_info['scale_factor'])
        
    def load_steady_state_frames(self) -> List[Tuple[str, int]]:
        """加载稳态帧文件列表"""
        pattern = self.dir_info['path'] / 'mixer_output_*.csv'
        files = sorted(glob.glob(str(pattern)))
        
        if not files:
            raise FileNotFoundError(f"No CSV files in {self.dir_info['path']}")
        
        # 提取帧号并排序
        all_frames = []
        for f in files:
            match = re.search(r'_(\d+)\.csv', f)
            if match:
                all_frames.append((f, int(match.group(1))))
        all_frames.sort(key=lambda x: x[1])
        
        # 计算起始帧号（从 steady_state_start_time 开始）
        start_frame = int(self.config.steady_state_start_time / self.frame_dt)
        
        steady_frames = [
            (f, n) for f, n in all_frames
            if n >= start_frame and n % self.config.frame_skip == 0
        ]
        
        return steady_frames
    
    def process(self) -> pd.DataFrame:
        """处理目录，返回时间序列特征DataFrame"""
        frames = self.load_steady_state_frames()
        
        if not frames:
            raise ValueError(f"No valid frames found in {self.dir_info['path']}")
        
        print(f"  处理 {len(frames)} 帧 (dt={self.frame_dt*1000:.1f}ms)")
        
        results = []
        for file_path, frame_num in frames:
            df = pd.read_csv(file_path)
            features = FrameFeatureExtractor.extract_features(df)
            features['frame'] = frame_num
            features['time'] = frame_num * self.frame_dt
            results.append(features)
        
        return pd.DataFrame(results)


# ======================== Multi-Directory Comparison ========================
class ExtendedFeatureAnalyzer:
    """扩展特征分析器 - 处理多目录比较"""
    
    def __init__(self, config: AnalysisConfig):
        self.config = config
        self.directories = find_sup_directories(config.base_path)
        
    def process_all(self):
        """处理所有匹配的目录"""
        if not self.directories:
            print("未找到任何SUP输出目录")
            return
        
        # 显示找到的目录
        print(f"\n找到 {len(self.directories)} 个目录:")
        groups = group_directories(self.directories)
        for (dir_type, se), dirs in sorted(groups.items()):
            sfs = [d['scale_factor'] for d in dirs]
            print(f"  {dir_type}_se{se}: f{sfs}")
        
        # 按组处理
        for (dir_type, surface_energy), group_dirs in sorted(groups.items()):
            # 应用过滤
            if self.config.filter_cohesion and surface_energy != self.config.filter_cohesion:
                continue
            if self.config.filter_type and dir_type != self.config.filter_type:
                continue
            
            print(f"\n{'='*60}")
            print(f"处理组: {dir_type}_se{surface_energy}")
            print(f"{'='*60}")
            
            self._process_group(dir_type, surface_energy, group_dirs)
    
    def _process_group(self, dir_type: str, surface_energy: str, group_dirs: List[Dict]):
        """处理一个分组"""
        
        # 创建输出目录
        output_dir = Path(self.config.output_base) / f"{dir_type}_se{surface_energy}"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        all_time_series = {}
        all_steady_stats = []
        
        for dir_info in sorted(group_dirs, key=lambda x: x['scale_factor']):
            sf = dir_info['scale_factor']
            print(f"\n>>> f{sf}: {dir_info['name']}")
            
            try:
                processor = DirectoryProcessor(dir_info, self.config)
                df_features = processor.process()
                
                # 保存时间序列
                col_prefix = f"f{sf}"
                all_time_series[col_prefix] = df_features
                
                # 计算稳态统计量
                steady_stats = self._compute_steady_stats(df_features, sf, dir_type)
                all_steady_stats.append(steady_stats)
                
                # 保存单独的时间序列文件
                ts_file = output_dir / f"timeseries_f{sf}.csv"
                df_features.to_csv(ts_file, index=False, float_format='%.6e')
                print(f"  ✓ 时间序列保存: {ts_file.name}")
                
            except Exception as e:
                print(f"  ✗ 错误: {e}")
                continue
        
        if not all_steady_stats:
            print("  无有效数据")
            return
        
        # 合并稳态统计
        df_steady = pd.DataFrame(all_steady_stats)
        df_steady = df_steady.sort_values('scale_factor').reset_index(drop=True)
        
        # 保存稳态统计
        steady_file = output_dir / "extended_features_steady.csv"
        df_steady.to_csv(steady_file, index=False, float_format='%.6e')
        print(f"\n✓ 稳态特征保存: {steady_file}")
        
        # 保存合并的时间序列比较文件
        self._save_comparison_files(all_time_series, output_dir)
        
        # 打印摘要
        self._print_summary(df_steady)
    
    def _compute_steady_stats(self, df: pd.DataFrame, scale_factor: int, dir_type: str) -> Dict:
        """计算稳态统计量（时间平均）"""
        
        # 需要求平均的特征
        mean_features = [
            'n_particles', 'KE_trans_sum', 'KE_rot_sum', 'PE_sum',
            'vx_mean', 'vy_mean', 'vz_mean', 'speed_mean', 'speed_std', 'com_z',
            'v_r_mean', 'v_theta_mean', 'v_z_cyl_mean',
            'v_r_std', 'v_theta_std', 'v_z_cyl_std', 'omega_mean',
            'T_r', 'T_theta', 'T_z', 'T_total'
        ]
        
        stats = {
            'scale_factor': scale_factor,
            'dir_type': dir_type,
            'n_frames': len(df),
            'time_start': df['time'].min(),
            'time_end': df['time'].max(),
        }
        
        for feat in mean_features:
            if feat in df.columns:
                stats[f'{feat}_mean'] = df[feat].mean()
                stats[f'{feat}_std'] = df[feat].std()
        
        return stats
    
    def _save_comparison_files(self, all_time_series: Dict[str, pd.DataFrame], output_dir: Path):
        """保存跨scale_factor的比较文件"""
        
        if len(all_time_series) < 2:
            return
        
        # 选择关键特征做比较
        key_features = [
            'KE_trans_sum', 'KE_rot_sum', 'speed_mean', 
            'T_total', 'omega_mean', 'com_z'
        ]
        
        for feat in key_features:
            comparison_data = {'Time(s)': None}
            
            for prefix, df in sorted(all_time_series.items()):
                if comparison_data['Time(s)'] is None:
                    comparison_data['Time(s)'] = df['time'].values
                
                if feat in df.columns:
                    comparison_data[prefix] = df[feat].values
            
            if len(comparison_data) > 1:
                # 处理不同长度的时间序列 - 使用插值
                time_ref = comparison_data['Time(s)']
                df_comp = pd.DataFrame({'Time(s)': time_ref})
                
                for prefix, df in sorted(all_time_series.items()):
                    if feat in df.columns:
                        # 插值到参考时间网格
                        interp_vals = np.interp(
                            time_ref, 
                            df['time'].values, 
                            df[feat].values,
                            left=np.nan, right=np.nan
                        )
                        df_comp[prefix] = interp_vals
                
                comp_file = output_dir / f"comparison_{feat}.csv"
                df_comp.to_csv(comp_file, index=False, float_format='%.6e')
        
        print(f"✓ 比较文件保存: comparison_*.csv")
    
    def _print_summary(self, df_steady: pd.DataFrame):
        """打印摘要统计"""
        print(f"\n{'─'*40}")
        print("稳态特征摘要 (时间平均):")
        print(f"{'─'*40}")
        
        key_cols = ['scale_factor', 'n_frames', 
                    'speed_mean_mean', 'T_total_mean', 'omega_mean_mean']
        
        available_cols = [c for c in key_cols if c in df_steady.columns]
        if available_cols:
            print(df_steady[available_cols].to_string(index=False))


# ======================== Main ========================
def main():
    """主函数"""
    print("=" * 60)
    print("SUP Mixer Extended Features Analysis")
    print("数量加权版本 - 适用于所有材料配置")
    print("=" * 60)
    
    # ==================== 配置区域 ====================
    config = AnalysisConfig(
        base_path="/media/huyuze/Fanxiang/DEME_Data/backup251115",
        filter_cohesion=None,    # 只处理特定cohesion，None=处理所有
        filter_type=None,         # 只处理特定类型，None=处理所有
        steady_state_start_time=0.5,  # 从0.5秒开始作为稳态
        frame_skip=10,
        output_base="extended_features_results"
    )
    # =================================================
    
    print(f"\n配置:")
    print(f"  基础路径: {config.base_path}")
    print(f"  Cohesion过滤: {config.filter_cohesion or '全部'}")
    print(f"  类型过滤: {config.filter_type or '全部'}")
    print(f"  稳态起始时间: {config.steady_state_start_time} s")
    print(f"  帧采样间隔: {config.frame_skip}")
    
    # 运行分析
    analyzer = ExtendedFeatureAnalyzer(config)
    analyzer.process_all()
    
    print("\n" + "=" * 60)
    print("分析完成!")
    print("=" * 60)


if __name__ == "__main__":
    main()