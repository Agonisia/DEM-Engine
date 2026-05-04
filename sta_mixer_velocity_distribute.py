#!/usr/bin/env python3
"""
SUP Mixer Velocity Distribution Analysis v2 - Multi-Directory Comparison Version
=================================================================================
Analyzes particle velocity distributions in DEM mixer simulations with SUP models

更新内容 (v2):
  - 支持多个基础目录路径 (base_paths: List[str])
  - 自动跨路径搜索、去重 (同 sta_other_info_v2 模式)
  - 保持: 原有全部分析功能不变

Fixed: Flattened MultiIndex headers in CSV outputs for better compatibility.
"""

import os
import re
import datetime
import pandas as pd
import numpy as np
from numba import jit, prange
import glob
from pathlib import Path
from typing import Tuple, Dict, Optional, List
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import warnings
warnings.filterwarnings('ignore')

# ======================== Configuration Parameters ========================
class Config:
    """Configuration parameter management for single-directory analysis"""
    def __init__(self, scale_factor: int = 4, surface_energy: str = "000",
                 input_path: Optional[Path] = None, base_path: str = ".",
                 enable_plotting: bool = False,
                 steady_state_fraction: float = 0.5, dir_type: str = "standard"):
        self.scale_factor = scale_factor
        self.surface_energy = surface_energy
        self.dir_type = dir_type

        # 支持直接传入已发现的路径 (多目录模式)
        if input_path is not None:
            self.input_path = Path(input_path)
        else:
            self.base_path = Path(base_path)
            if dir_type == "standard":
                self.input_path = self.base_path / f'SUPMixerOutput_f{scale_factor}se{surface_energy}'
            elif dir_type == "SizeDiff":
                self.input_path = self.base_path / f'SUPMixerOutput_SizeDiff_f{scale_factor}se{surface_energy}'
            elif dir_type == "DualDensity":
                self.input_path = self.base_path / f'SUPMixerOutput_DualDensity_f{scale_factor}se{surface_energy}'
            else:
                raise ValueError(f"Unknown directory type: {dir_type}")

        self.output_base = Path(f'sta_results/{dir_type}_se{surface_energy}/')
        self.output_path = self.output_base / f'f{scale_factor}/'

        if scale_factor == 1:
            self.frame_dt = 5e-3
        else:
            self.frame_dt = 1e-3

        self.RPM = 300
        self.frame_skip = 100
        self.steady_state_fraction = steady_state_fraction
        self.enable_plotting = enable_plotting
        self.plot_dpi = 300

        if enable_plotting:
            plt.rcParams['font.family'] = 'serif'
            plt.rcParams['axes.unicode_minus'] = False
        else:
            matplotlib.use('Agg')

        self.output_path.mkdir(parents=True, exist_ok=True)


class MultiDirectoryConfig:
    """Configuration for multi-directory comparison analysis — supports multiple base paths"""
    def __init__(self, base_paths=None, filter_cohesion: str = "000",
                 use_interpolation: bool = False, enable_plotting: bool = True,
                 steady_state_fraction: float = 0.5, filter_type: Optional[str] = None):
        # 兼容旧接口：接受 str 或 List[str]
        if base_paths is None:
            base_paths = ["."]
        elif isinstance(base_paths, str):
            base_paths = [base_paths]
        self.base_paths = base_paths

        self.filter_cohesion = filter_cohesion
        self.filter_type = filter_type
        self.enable_plotting = enable_plotting
        self.steady_state_fraction = steady_state_fraction
        self.directories = self._find_matching_directories()

    def _find_matching_directories(self) -> List[Dict]:
        """从多个基础路径查找所有匹配目录，并去重"""
        dirs = []
        patterns = [
            ('standard', r'SUPMixerOutput_f(\d+)se(\d+)$'),
            ('SizeDiff', r'SUPMixerOutput_SizeDiff_f(\d+)se(\d+)$'),
            ('DualDensity', r'SUPMixerOutput_DualDensity_f(\d+)se(\d+)$')
        ]

        for base_path in self.base_paths:
            bp = Path(base_path)
            if not bp.exists():
                print(f"警告：路径 {base_path} 不存在，跳过")
                continue

            for path in bp.iterdir():
                if not path.is_dir():
                    continue
                for dir_type, pattern in patterns:
                    match = re.match(pattern, path.name)
                    if match:
                        scale_factor = int(match.group(1))
                        surface_energy = match.group(2)
                        if self.filter_cohesion and surface_energy != self.filter_cohesion:
                            continue
                        if self.filter_type and dir_type != self.filter_type:
                            continue
                        dirs.append({
                            'path': path,
                            'scale_factor': scale_factor,
                            'surface_energy': surface_energy,
                            'dir_type': dir_type,
                            'name': path.name,
                            'base_path': str(base_path)
                        })
                        break

        # 去重：(dir_type, surface_energy, scale_factor) 只保留第一个
        seen = set()
        unique = []
        for d in dirs:
            key = (d['dir_type'], d['surface_energy'], d['scale_factor'])
            if key not in seen:
                seen.add(key)
                unique.append(d)

        type_order = {'standard': 0, 'SizeDiff': 1, 'DualDensity': 2}
        return sorted(unique, key=lambda x: (
            type_order.get(x['dir_type'], 99),
            x['surface_energy'],
            x['scale_factor']
        ))


# ======================== JIT Accelerated Functions ========================
@jit(nopython=True, parallel=True)
def compute_angles_vectorized(x, y, frame_numbers, frame_dt, RPM):
    n = len(x)
    polar_angles = np.zeros(n); rotor_angles = np.zeros(n); relative_angles = np.zeros(n)
    RPS = RPM / 60.0
    for i in prange(n):
        polar_angles[i] = np.arctan2(y[i], x[i])
        time = frame_numbers[i] * frame_dt
        rotor_angles[i] = (RPS * time * 2 * np.pi) % (2 * np.pi)
        rel_angle = polar_angles[i] - rotor_angles[i]
        relative_angles[i] = np.mod(rel_angle + np.pi, 2 * np.pi) - np.pi
    return polar_angles, rotor_angles, relative_angles

@jit(nopython=True, parallel=True)
def classify_angles_vectorized(angles):
    n = len(angles)
    intervals = np.zeros(n, dtype=np.int32)
    for i in prange(n):
        angle_deg = np.mod(np.rad2deg(angles[i]) + 360, 360)
        intervals[i] = 1 + int(angle_deg // 10)
    return intervals

@jit(nopython=True, parallel=True)
def calculate_all_velocities(v_x, v_y, v_z, w_x, w_y, w_z, x, y):
    n = len(v_x)
    v_total = np.sqrt(v_x**2 + v_y**2 + v_z**2)
    v_xy = np.sqrt(v_x**2 + v_y**2)
    v_z_abs = np.abs(v_z)
    w_total = np.sqrt(w_x**2 + w_y**2 + w_z**2)
    r_dist = np.sqrt(x**2 + y**2)
    return v_total, v_xy, v_z_abs, w_total, r_dist


# ======================== Data Processing & Analysis ========================
class MixerDataProcessor:
    def __init__(self, config: Config):
        self.config = config
        self.df = None

    def load_data(self) -> pd.DataFrame:
        pattern = self.config.input_path / 'mixer_output_*.csv'
        files = sorted(glob.glob(str(pattern)))
        if not files:
            raise FileNotFoundError(f"No CSV files in {self.config.input_path}")

        all_frames = sorted(
            [(f, int(re.search(r'_(\d+)\.csv', f).group(1))) for f in files],
            key=lambda x: x[1]
        )
        start_idx = int(len(all_frames) * (1 - self.config.steady_state_fraction))
        steady_frames = [f for f in all_frames[start_idx:] if f[1] % self.config.frame_skip == 0]

        all_data = []
        for file_path, frame_num in steady_frames:
            df = pd.read_csv(file_path)
            df['frame'] = frame_num
            df['time'] = frame_num * self.config.frame_dt
            all_data.append(df)
        self.df = pd.concat(all_data, ignore_index=True)
        return self.df

    def process_velocities(self):
        v_t, v_xy, v_za, w_t, r_d = calculate_all_velocities(
            self.df['v_x'].values, self.df['v_y'].values, self.df['v_z'].values,
            self.df['w_x'].values, self.df['w_y'].values, self.df['w_z'].values,
            self.df['X'].values, self.df['Y'].values
        )
        self.df['v_total'], self.df['v_xy'], self.df['v_z_abs'] = v_t, v_xy, v_za
        self.df['w_total'], self.df['r_dist'] = w_t, r_d

    def process_angles(self):
        pa, ra, rela = compute_angles_vectorized(
            self.df['X'].values, self.df['Y'].values,
            self.df['frame'].values, self.config.frame_dt, self.config.RPM
        )
        self.df['relative_angle'] = rela
        self.df['angle_interval'] = classify_angles_vectorized(rela)


class StatisticalAnalyzer:
    """核心修复：处理多级表头和索引"""
    def __init__(self, df: pd.DataFrame):
        self.df = df

    def _flatten_cols(self, df_stats):
        """内部工具：将 (变量, 统计量) 转为 '变量_统计量' 格式"""
        df_stats.columns = [f"{c[0]}_{c[1]}" if c[1] else c[0] for c in df_stats.columns]
        return df_stats.reset_index()

    def compute_angle_statistics(self) -> pd.DataFrame:
        stats = self.df.groupby('angle_interval').agg({
            'v_z': ['mean', 'std'],
            'v_z_abs': ['mean'],
            'v_total': ['mean', 'std'],
            'Z': ['mean'],
            'r': ['mean']
        }).round(4)
        return self._flatten_cols(stats)

    def compute_radial_statistics(self, bin_size_mm: float = 1.0) -> pd.DataFrame:
        self.df['r_group'] = (self.df['r_dist'] * 1000 / bin_size_mm).astype(int) * bin_size_mm / 1000
        stats = self.df.groupby('r_group').agg({
            'v_total': ['mean', 'std'],
            'v_z': ['mean', 'std'],
            'Z': ['mean']
        }).round(4)
        return self._flatten_cols(stats)

    def compute_particle_statistics(self) -> pd.DataFrame:
        stats = self.df.groupby('r').agg({
            'v_total': ['mean', 'std'],
            'v_z': ['mean', 'std'],
            'w_total': ['mean', 'std']
        }).round(4)
        return self._flatten_cols(stats)

    def compute_velocity_pdf(self, bins: int = 76, v_range: Tuple[float, float] = (0, 1.52)) -> pd.DataFrame:
        hist, edges = np.histogram(
            self.df['v_total'].values,
            bins=np.linspace(v_range[0], v_range[1], bins + 1), density=True
        )
        return pd.DataFrame({
            'velocity_center': (edges[:-1] + edges[1:]) / 2,
            'probability_density': hist
        })

    def get_summary_statistics(self) -> dict:
        return {
            'total_particles': len(self.df),
            'num_frames': self.df['frame'].nunique(),
            'velocity_stats': self.df['v_total'].describe().to_dict(),
            'time_range': (self.df['time'].min(), self.df['time'].max())
        }


# ======================== Multi-Directory Comparison Class ========================
class MultiDirectoryComparison:
    def __init__(self, multi_config: MultiDirectoryConfig):
        self.multi_config = multi_config
        self.results = {}

    def process_all_directories(self):
        for dir_info in self.multi_config.directories:
            print(f"\n>>> Processing: {dir_info['name']}  [{dir_info['base_path']}]")
            try:
                config = Config(
                    scale_factor=dir_info['scale_factor'],
                    surface_energy=dir_info['surface_energy'],
                    input_path=dir_info['path'],  # 直接传入已发现的路径
                    enable_plotting=self.multi_config.enable_plotting,
                    steady_state_fraction=self.multi_config.steady_state_fraction,
                    dir_type=dir_info['dir_type']
                )
                processor = MixerDataProcessor(config)
                df = processor.load_data()
                processor.process_velocities()
                processor.process_angles()

                analyzer = StatisticalAnalyzer(df)
                self.results[dir_info['name']] = {
                    'config': config,
                    'dir_info': dir_info,
                    'analyzer': analyzer,
                    'angle_stats': analyzer.compute_angle_statistics(),
                    'radial_stats': analyzer.compute_radial_statistics(),
                    'particle_stats': analyzer.compute_particle_statistics(),
                    'pdf_data': analyzer.compute_velocity_pdf(),
                    'summary': analyzer.get_summary_statistics()
                }
                self.save_individual_results(dir_info['name'])
            except Exception as e:
                print(f"  Error: {e}")
                import traceback
                traceback.print_exc()
                continue

    def save_individual_results(self, dir_name: str):
        res = self.results[dir_name]
        path = res['config'].output_path
        sf = res['config'].scale_factor
        dt = res['dir_info']['dir_type']

        res['angle_stats'].to_csv(path / f'AngleStats_{dt}_f{sf}.csv', index=False)
        res['radial_stats'].to_csv(path / f'RadialStats_{dt}_f{sf}.csv', index=False)
        res['particle_stats'].to_csv(path / f'ParticleStats_{dt}_f{sf}.csv', index=False)
        res['pdf_data'].to_csv(path / f'PDF_{dt}_f{sf}.csv', index=False)

        with open(path / f'Summary_{dt}_f{sf}.txt', 'w') as f:
            f.write(f"Analysis Summary for {dir_name}\n" + "=" * 30 + "\n")
            for k, v in res['summary'].items():
                f.write(f"{k}: {v}\n")

        print(f"  ✓ 结果已保存到: {path}")


# ======================== Main Program ========================
def main():
    # ==================== 用户配置区 ====================
    base_paths = [
        "/media/huyuze/Fanxiang1/DEME_Data/backup251115",
        "/media/huyuze/Fanxiang/backup251025",
        "/media/huyuze/Fanxiang/backup251109",
    ]
    filter_cohesion = "015"       # None = 处理所有
    filter_type = None            # None = 处理所有
    # ===================================================

    print("=" * 60)
    print("SUP Mixer Velocity Distribution Analysis v2")
    print("支持多基础路径搜索")
    print("=" * 60)
    print(f"\n配置:")
    print(f"  基础路径: {base_paths}")
    print(f"  Cohesion过滤: {filter_cohesion or '全部'}")
    print(f"  类型过滤: {filter_type or '全部'}")

    multi_config = MultiDirectoryConfig(
        base_paths=base_paths,
        filter_cohesion=filter_cohesion,
        filter_type=filter_type,
    )

    if not multi_config.directories:
        print("No matching directories found.")
        return

    print(f"\n找到 {len(multi_config.directories)} 个目录:")
    for d in multi_config.directories:
        print(f"  - {d['name']}  [{d['base_path']}]")

    comparison = MultiDirectoryComparison(multi_config)
    comparison.process_all_directories()
    print("\n[DONE] All analysis completed. Check 'sta_results/' folder.")


if __name__ == "__main__":
    main()