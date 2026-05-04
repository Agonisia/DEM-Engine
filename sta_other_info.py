#!/usr/bin/env python3
"""
SUP Mixer Extended Features Analysis v3.5 (Aligned & Dense Version)
================================================================
扩展特征分析脚本，计算完整的动力学特征集用于SUP模型训练

更新内容 (v3.5 — 强制高密度对齐版):
  - 核心改动: 修改了数据采样间隔和插值逻辑，强制输出与 energy&torque 脚本同密度的结果。
  - frame_skip 默认改为 1，提取所有可用帧文件。
  - steady_state_start_time 默认改为 0.0，对齐能量数据的时间轴起点。
  - 插值时间轴直接从 sta_results 的 kinetic_energy_comparison.csv 读取，
    确保 comparison CSV 与 KE/Torque 的 Time(s) 列完全一致。
"""

import os
import re
import glob
import pandas as pd
import numpy as np
from numba import jit, prange
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass, field
import warnings
warnings.filterwarnings('ignore')

# ======================== Constants ========================
G = 9.81  # m/s²
PI = np.pi


# ======================== Configuration ========================
@dataclass
class AnalysisConfig:
    """分析配置参数"""
    base_paths: List[str] = field(default_factory=lambda: ["."])

    filter_cohesion: Optional[str] = None
    filter_type: Optional[str] = None

    # [v3.5] 确保提取全量数据且起点一致
    steady_state_start_time: float = 0.0
    frame_skip: int = 1

    # [v3.5] 从 sta_results 读取 KE 时间轴用于对齐
    auto_interpolate: bool = True
    sta_base: str = "sta_results"

    output_base: str = "extended_features_results"

    # 密度配置
    default_density: float = 1000.0  # kg/m³
    dual_density_map: Dict[int, float] = field(default_factory=lambda: {
        2: 1000.0,
        3: 2000.0,
    })

    # Cylindrical Binning 参数
    nr_bins: int = 8
    nz_bins: int = 6
    min_particles_per_bin: int = 5

    def get_frame_dt(self, scale_factor: int) -> float:
        return 5e-3 if scale_factor == 1 else 1e-3


# ======================== JIT Accelerated Functions ========================
@jit(nopython=True, parallel=True)
def compute_cylindrical_velocities(x, y, vx, vy):
    n = len(x)
    v_r = np.zeros(n)
    v_theta = np.zeros(n)
    r_dist = np.zeros(n)

    for i in prange(n):
        r = np.sqrt(x[i]**2 + y[i]**2)
        r_dist[i] = r
        if r > 1e-10:
            v_r[i] = (x[i] * vx[i] + y[i] * vy[i]) / r
            v_theta[i] = (-y[i] * vx[i] + x[i] * vy[i]) / r

    return v_r, v_theta, r_dist


@jit(nopython=True)
def compute_global_omega_mean(v_theta, r_dist, mass):
    num = 0.0
    den = 0.0
    for i in range(len(v_theta)):
        m = mass[i]
        num += m * v_theta[i] * r_dist[i]
        den += m * r_dist[i] * r_dist[i]
    return num / den if den > 1e-10 else 0.0


@jit(nopython=True)
def compute_granular_temperature_binned(v_r, v_theta, v_z, r_dist, z, mass,
                                         nr_bins, nz_bins, r_max, z_min, z_max,
                                         min_count):
    n = len(v_r)
    n_bins_total = nr_bins * nz_bins

    dr = r_max / nr_bins
    dz = (z_max - z_min) / nz_bins

    bin_mass = np.zeros(n_bins_total)
    bin_mvr = np.zeros(n_bins_total)
    bin_mvt = np.zeros(n_bins_total)
    bin_mvz = np.zeros(n_bins_total)
    bin_count = np.zeros(n_bins_total, dtype=np.int64)

    particle_bin = np.full(n, -1, dtype=np.int64)

    for i in range(n):
        ir = int(r_dist[i] / dr)
        if ir >= nr_bins:
            ir = nr_bins - 1
        if ir < 0:
            ir = 0

        iz = int((z[i] - z_min) / dz)
        if iz >= nz_bins:
            iz = nz_bins - 1
        if iz < 0:
            iz = 0

        b = ir * nz_bins + iz
        particle_bin[i] = b

        m = mass[i]
        bin_mass[b] += m
        bin_mvr[b] += m * v_r[i]
        bin_mvt[b] += m * v_theta[i]
        bin_mvz[b] += m * v_z[i]
        bin_count[b] += 1

    bin_vr_mean = np.zeros(n_bins_total)
    bin_vt_mean = np.zeros(n_bins_total)
    bin_vz_mean = np.zeros(n_bins_total)

    for b in range(n_bins_total):
        if bin_mass[b] > 1e-30 and bin_count[b] >= min_count:
            bin_vr_mean[b] = bin_mvr[b] / bin_mass[b]
            bin_vt_mean[b] = bin_mvt[b] / bin_mass[b]
            bin_vz_mean[b] = bin_mvz[b] / bin_mass[b]

    m_included = 0.0
    T_r = 0.0
    T_theta = 0.0
    T_z = 0.0
    n_included = 0

    for i in range(n):
        b = particle_bin[i]
        if b < 0 or bin_count[b] < min_count:
            continue

        m = mass[i]
        dvr = v_r[i] - bin_vr_mean[b]
        dvt = v_theta[i] - bin_vt_mean[b]
        dvz = v_z[i] - bin_vz_mean[b]

        T_r += m * dvr * dvr
        T_theta += m * dvt * dvt
        T_z += m * dvz * dvz
        m_included += m
        n_included += 1

    if m_included > 1e-30:
        T_r /= m_included
        T_theta /= m_included
        T_z /= m_included

    T_total = (T_r + T_theta + T_z) / 3.0

    return T_r, T_theta, T_z, T_total, n_included


@jit(nopython=True)
def compute_all_physics(x, y, z, vx, vy, vz, wx, wy, wz, mass, inertia):
    n = len(x)
    result = np.zeros(28)

    for i in prange(n):
        m = mass[i]
        I = inertia[i]

        v_sq = vx[i]**2 + vy[i]**2 + vz[i]**2
        w_sq = wx[i]**2 + wy[i]**2 + wz[i]**2
        speed = np.sqrt(v_sq)
        omega_mag = np.sqrt(w_sq)

        result[0] += 0.5 * m * v_sq
        result[1] += 0.5 * I * w_sq
        result[2] += m * G * z[i]
        result[3] += v_sq
        result[4] += w_sq

        result[5] += m * (y[i] * vz[i] - z[i] * vy[i])
        result[6] += m * (z[i] * vx[i] - x[i] * vz[i])
        result[7] += m * (x[i] * vy[i] - y[i] * vx[i])

        result[8] += I * wx[i]
        result[9] += I * wy[i]
        result[10] += I * wz[i]

        result[11] += vx[i]
        result[12] += vy[i]
        result[13] += vz[i]
        result[14] += speed
        result[15] += speed * speed
        result[16] += z[i]
        result[18] += m * z[i]
        result[19] += m

        result[20] += m * vx[i]
        result[21] += m * vy[i]
        result[22] += m * vz[i]

        result[23] += wz[i]
        result[24] += wz[i] * wz[i]
        result[25] += omega_mag
        result[26] += omega_mag * omega_mag

    result[17] = n
    return result


# ======================== Directory Discovery ========================
def find_sup_directories(base_paths) -> List[Dict]:
    if isinstance(base_paths, str):
        base_paths = [base_paths]

    directories = []
    patterns = [
        ('standard', re.compile(r'SUPMixerOutput_f(\d+)se(\d{3})$')),
        ('SizeDiff', re.compile(r'SUPMixerOutput_SizeDiff_f(\d+)se(\d{3})$')),
        ('DualDensity', re.compile(r'SUPMixerOutput_DualDensity_f(\d+)se(\d{3})$'))
    ]

    for base_path in base_paths:
        base = Path(base_path)
        if not base.exists():
            print(f"警告：路径 {base_path} 不存在，跳过")
            continue

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
                        'dir_type': dir_type,
                        'base_path': str(base_path)
                    })
                    break

    seen = set()
    unique = []
    for d in directories:
        key = (d['dir_type'], d['surface_energy'], d['scale_factor'])
        if key not in seen:
            seen.add(key)
            unique.append(d)

    type_order = {'standard': 0, 'SizeDiff': 1, 'DualDensity': 2}
    unique.sort(key=lambda x: (
        type_order.get(x['dir_type'], 99),
        x['surface_energy'],
        x['scale_factor']
    ))
    return unique

def group_directories(directories: List[Dict]) -> Dict[Tuple[str, str], List[Dict]]:
    groups = {}
    for d in directories:
        key = (d['dir_type'], d['surface_energy'])
        if key not in groups:
            groups[key] = []
        groups[key].append(d)
    return groups


# ======================== Data Processing ========================
class FrameFeatureExtractor:
    def __init__(self, config: AnalysisConfig):
        self.config = config

    def _compute_mass_inertia(self, r: np.ndarray, family: np.ndarray,
                               dir_type: str) -> Tuple[np.ndarray, np.ndarray]:
        n = len(r)
        mass = np.zeros(n)
        inertia = np.zeros(n)

        if dir_type == 'DualDensity':
            for i in range(n):
                fam = int(family[i])
                rho = self.config.dual_density_map.get(fam, self.config.default_density)
                m = (4.0 / 3.0) * PI * r[i]**3 * rho
                mass[i] = m
                inertia[i] = 0.4 * m * r[i]**2
        else:
            rho = self.config.default_density
            for i in range(n):
                m = (4.0 / 3.0) * PI * r[i]**3 * rho
                mass[i] = m
                inertia[i] = 0.4 * m * r[i]**2
        return mass, inertia

    def extract_features(self, df: pd.DataFrame, dir_type: str) -> Dict[str, float]:
        x = df['X'].values.astype(np.float64)
        y = df['Y'].values.astype(np.float64)
        z = df['Z'].values.astype(np.float64)
        vx = df['v_x'].values.astype(np.float64)
        vy = df['v_y'].values.astype(np.float64)
        vz = df['v_z'].values.astype(np.float64)
        wx = df['w_x'].values.astype(np.float64)
        wy = df['w_y'].values.astype(np.float64)
        wz = df['w_z'].values.astype(np.float64)
        r = df['r'].values.astype(np.float64)
        family = df['family'].values.astype(np.float64) if 'family' in df.columns else np.zeros(len(df))

        n = len(df)
        mass, inertia = self._compute_mass_inertia(r, family, dir_type)
        phys = compute_all_physics(x, y, z, vx, vy, vz, wx, wy, wz, mass, inertia)

        KE_trans = phys[0]
        KE_rot = phys[1]
        PE = phys[2]
        KE_trans_proxy = phys[3]
        KE_rot_proxy = phys[4]
        L_orb_x, L_orb_y, L_orb_z = phys[5], phys[6], phys[7]
        L_spin_x, L_spin_y, L_spin_z = phys[8], phys[9], phys[10]

        vx_mean = phys[11] / n
        vy_mean = phys[12] / n
        vz_mean = phys[13] / n
        com_z = phys[16] / n

        L_orbital_mag = np.sqrt(L_orb_x**2 + L_orb_y**2 + L_orb_z**2)
        L_spin_mag = np.sqrt(L_spin_x**2 + L_spin_y**2 + L_spin_z**2)
        M_total = phys[19]

        wz_mean = phys[23] / n
        wz_sq_mean = phys[24] / n
        wz_std = np.sqrt(max(wz_sq_mean - wz_mean**2, 0.0))

        omega_mag_mean = phys[25] / n
        omega_mag_sq_mean = phys[26] / n
        omega_mag_std = np.sqrt(max(omega_mag_sq_mean - omega_mag_mean**2, 0.0))

        v_r, v_theta, r_dist = compute_cylindrical_velocities(x, y, vx, vy)
        omega_mean = compute_global_omega_mean(v_theta, r_dist, mass)

        r_max = np.max(r_dist) * 1.001
        z_min = np.min(z) - 1e-10
        z_max = np.max(z) + 1e-10

        T_r, T_theta, T_z, T_total, n_included = compute_granular_temperature_binned(
            v_r, v_theta, vz, r_dist, z, mass,
            self.config.nr_bins, self.config.nz_bins,
            r_max, z_min, z_max, self.config.min_particles_per_bin
        )

        KE_trans_fluct = 1.5 * M_total * T_total
        T_components = [T_r, T_theta, T_z]
        T_min = min(T_components)
        T_max_val = max(T_components)
        GT_aniso = T_max_val / T_min if T_min > 1e-20 else 0.0

        return {
            'n_particles': n,
            'n_GT_included': n_included,
            'KE_trans': KE_trans, 'KE_rot': KE_rot, 'KE_total': KE_trans + KE_rot,
            'PE': PE, 'KE_trans_fluct': KE_trans_fluct,
            'L_orbital_z': L_orb_z, 'L_spin_z': L_spin_z, 'L_total_z': L_orb_z + L_spin_z,
            'L_orbital_mag': L_orbital_mag, 'L_spin_mag': L_spin_mag,
            'KE_trans_proxy': KE_trans_proxy, 'KE_rot_proxy': KE_rot_proxy,
            'KE_trans_sum': KE_trans_proxy, 'KE_rot_sum': KE_rot_proxy,
            'vx_mean': vx_mean, 'vy_mean': vy_mean, 'vz_mean': vz_mean, 'com_z': com_z,
            'v_r_mean': np.mean(v_r), 'v_theta_mean': np.mean(v_theta), 'v_z_cyl_mean': np.mean(vz),
            'omega_mean': omega_mean, 'omega_z_mean': wz_mean, 'omega_z_std': wz_std,
            'omega_mag_mean': omega_mag_mean, 'omega_mag_std': omega_mag_std,
            'T_r': T_r, 'T_theta': T_theta, 'T_z': T_z, 'T_total': T_total, 'GT_aniso': GT_aniso,
            'v_r_std': np.sqrt(T_r), 'v_theta_std': np.sqrt(T_theta), 'v_z_cyl_std': np.sqrt(T_z),
            'speed_mean': phys[14] / n, 'speed_std': np.sqrt(T_total * 3.0),
        }


class DirectoryProcessor:
    def __init__(self, dir_info: Dict, config: AnalysisConfig):
        self.dir_info = dir_info
        self.config = config
        self.frame_dt = config.get_frame_dt(dir_info['scale_factor'])
        self.extractor = FrameFeatureExtractor(config)

    def load_steady_state_frames(self) -> List[Tuple[str, int]]:
        pattern = self.dir_info['path'] / 'mixer_output_*.csv'
        files = sorted(glob.glob(str(pattern)))

        if not files:
            raise FileNotFoundError(f"No CSV files in {self.dir_info['path']}")

        all_frames = []
        for f in files:
            match = re.search(r'_(\d+)\.csv', f)
            if match:
                all_frames.append((f, int(match.group(1))))
        all_frames.sort(key=lambda x: x[1])

        start_frame = int(self.config.steady_state_start_time / self.frame_dt)
        return [
            (f, n) for f, n in all_frames
            if n >= start_frame and n % self.config.frame_skip == 0
        ]

    def process(self) -> pd.DataFrame:
        frames = self.load_steady_state_frames()
        if not frames:
            raise ValueError(f"No valid frames in {self.dir_info['path']}")

        print(f"  处理 {len(frames)} 帧 (dt={self.frame_dt*1000:.1f}ms, skip={self.config.frame_skip})")

        results = []
        for file_path, frame_num in frames:
            df = pd.read_csv(file_path)
            features = self.extractor.extract_features(df, self.dir_info['dir_type'])
            features['frame'] = frame_num
            features['time'] = frame_num * self.frame_dt
            results.append(features)
        return pd.DataFrame(results)


# ======================== Multi-Directory Comparison ========================
class ExtendedFeatureAnalyzer:
    def __init__(self, config: AnalysisConfig):
        self.config = config
        self.directories = find_sup_directories(config.base_paths)

    def process_all(self):
        if not self.directories:
            print("未找到任何SUP输出目录")
            return

        print(f"\n找到 {len(self.directories)} 个目录:")
        groups = group_directories(self.directories)
        for (dir_type, se), dirs in sorted(groups.items()):
            sfs = [d['scale_factor'] for d in dirs]
            print(f"  {dir_type}_se{se}: f{sfs}")

        for (dir_type, surface_energy), group_dirs in sorted(groups.items()):
            if self.config.filter_cohesion and surface_energy != self.config.filter_cohesion:
                continue
            if self.config.filter_type and dir_type != self.config.filter_type:
                continue

            print(f"\n{'='*60}")
            print(f"处理组: {dir_type}_se{surface_energy}")
            print(f"{'='*60}")

            self._process_group(dir_type, surface_energy, group_dirs)

    def _get_common_time(self, dir_type: str, surface_energy: str,
                          raw_time_series: Dict[str, pd.DataFrame]) -> np.ndarray:
        """
        [v3.5 核心] 获取公共时间轴。
        优先从 sta_results 的 kinetic_energy_comparison.csv 读取，
        确保与 KE/Torque 输出完全一致。
        """
        ke_file = (Path(self.config.sta_base)
                   / f"{dir_type}_se{surface_energy}"
                   / "kinetic_energy_comparison.csv")

        if ke_file.exists():
            ke_df = pd.read_csv(ke_file)
            common_time = ke_df['Time(s)'].values
            print(f"  ✓ 时间轴来源: {ke_file.name} ({len(common_time)} 点)")
            return common_time

        # Fallback: 用时间点最多的 scale factor
        densest_prefix = max(raw_time_series.keys(),
                             key=lambda k: len(raw_time_series[k]))
        common_time = raw_time_series[densest_prefix]['time'].values
        print(f"  ⚠ KE文件不存在，使用 {densest_prefix} 的时间轴 ({len(common_time)} 点)")
        return common_time

    def _process_group(self, dir_type: str, surface_energy: str, group_dirs: List[Dict]):
        output_dir = Path(self.config.output_base) / f"{dir_type}_se{surface_energy}"
        output_dir.mkdir(parents=True, exist_ok=True)

        raw_time_series = {}
        all_steady_stats = []

        # 1. 提取原始数据
        for dir_info in sorted(group_dirs, key=lambda x: x['scale_factor']):
            sf = dir_info['scale_factor']
            print(f"\n>>> f{sf}: {dir_info['name']}")
            try:
                processor = DirectoryProcessor(dir_info, self.config)
                df_features = processor.process()
                col_prefix = f"f{sf}"
                raw_time_series[col_prefix] = df_features
            except Exception as e:
                print(f"  ✗ 错误: {e}")
                import traceback
                traceback.print_exc()
                continue

        if not raw_time_series:
            print("  无有效数据")
            return

        # 2. [v3.5] 对齐到 KE 时间轴
        aligned_time_series = {}
        if self.config.auto_interpolate:
            common_time = self._get_common_time(dir_type, surface_energy, raw_time_series)

            for prefix, df in raw_time_series.items():
                interp_data = {'time': common_time}
                for col in df.columns:
                    if col not in ['time', 'frame']:
                        interp_data[col] = np.interp(
                            common_time, df['time'].values, df[col].values,
                            left=np.nan, right=np.nan
                        )
                aligned_time_series[prefix] = pd.DataFrame(interp_data)
        else:
            aligned_time_series = raw_time_series

        # 3. 保存时序和稳态统计
        for prefix, df_features in aligned_time_series.items():
            sf = int(prefix[1:])

            steady_stats = self._compute_steady_stats(df_features, sf, dir_type)
            all_steady_stats.append(steady_stats)

            ts_file = output_dir / f"timeseries_{prefix}.csv"
            df_features.to_csv(ts_file, index=False, float_format='%.6e')
            print(f"  ✓ 时间序列: {ts_file.name} ({len(df_features)} 点)")

        if not all_steady_stats:
            print("  无有效数据")
            return

        df_steady = pd.DataFrame(all_steady_stats)
        df_steady = df_steady.sort_values('scale_factor').reset_index(drop=True)
        steady_file = output_dir / "extended_features_steady.csv"
        df_steady.to_csv(steady_file, index=False, float_format='%.6e')
        print(f"\n✓ 稳态特征: {steady_file}")

        self._save_comparison_files(aligned_time_series, output_dir)
        self._print_summary(df_steady)

    def _compute_steady_stats(self, df: pd.DataFrame, scale_factor: int, dir_type: str) -> Dict:
        feature_cols = [c for c in df.columns if c not in ('frame', 'time')]
        stats = {
            'scale_factor': scale_factor,
            'dir_type': dir_type,
            'n_frames': len(df),
            'time_start': df['time'].min(),
            'time_end': df['time'].max(),
        }
        for feat in feature_cols:
            stats[f'{feat}_mean'] = df[feat].mean()
            stats[f'{feat}_std'] = df[feat].std()
        return stats

    def _save_comparison_files(self, aligned_time_series: Dict[str, pd.DataFrame], output_dir: Path):
        if len(aligned_time_series) < 2:
            return

        key_features = [
            'KE_trans', 'KE_rot', 'KE_total', 'PE',
            'KE_trans_fluct', 'L_orbital_z', 'L_spin_z', 'L_total_z',
            'L_orbital_mag', 'L_spin_mag', 'KE_trans_proxy', 'KE_rot_proxy',
            'KE_trans_sum', 'KE_rot_sum', 'speed_mean', 'speed_std', 'com_z',
            'vx_mean', 'vy_mean', 'vz_mean', 'v_r_mean', 'v_r_std',
            'v_theta_mean', 'v_theta_std', 'v_z_cyl_mean', 'v_z_cyl_std',
            'omega_mean', 'omega_z_mean', 'omega_z_std', 'omega_mag_mean', 'omega_mag_std',
            'T_total', 'T_r', 'T_theta', 'T_z', 'GT_aniso',
        ]

        # 所有 scale factor 已对齐到同一时间轴，取任意一个的 time 列
        first_prefix = sorted(aligned_time_series.keys())[0]
        time_ref = aligned_time_series[first_prefix]['time'].values

        n_saved = 0
        for feat in key_features:
            df_comp = pd.DataFrame({'Time(s)': time_ref})
            has_feat = False
            for prefix, df in sorted(aligned_time_series.items()):
                if feat in df.columns:
                    df_comp[prefix] = df[feat].values
                    has_feat = True

            if has_feat:
                comp_file = output_dir / f"comparison_{feat}.csv"
                df_comp.to_csv(comp_file, index=False, float_format='%.6e')
                n_saved += 1

        print(f"✓ 比较文件: {n_saved} comparison_*.csv ({len(time_ref)} 点, 与 KE 时间轴一致)")

    def _print_summary(self, df_steady: pd.DataFrame):
        print(f"\n{'─'*50}")
        print("稳态特征摘要 (时间平均):")
        print(f"{'─'*50}")

        key_cols = [
            'scale_factor', 'n_frames',
            'KE_trans_mean', 'KE_rot_mean', 'KE_trans_fluct_mean',
            'L_orbital_z_mean', 'L_spin_z_mean', 'T_total_mean',
            'omega_mag_std_mean', 'GT_aniso_mean',
        ]
        available = [c for c in key_cols if c in df_steady.columns]
        if available:
            display = df_steady[available].copy()
            for col in available:
                if col not in ('scale_factor', 'n_frames'):
                    display[col] = display[col].apply(lambda x: f"{x:.4e}")
            print(display.to_string(index=False))


# ======================== Main ========================
def main():
    print("=" * 60)
    print("SUP Mixer Extended Features Analysis v3.5")
    print("Cylindrical Binning 颗粒温度 | 强制与 KE 时间轴对齐")
    print("=" * 60)

    # ==================== 配置区域 ====================
    config = AnalysisConfig(
        base_paths=[
            "/media/huyuze/Fanxiang1/DEME_Data/backup251115",
            # "/media/huyuze/Fanxiang/backup251025",
            # "/media/huyuze/Fanxiang/backup251109"
        ],
        filter_cohesion=None,
        filter_type=None,

        # [v3.5 核心配置]
        steady_state_start_time=0.0,   # 从 t=0 开始，对齐 KE
        frame_skip=1,                   # 每帧都处理
        auto_interpolate=True,          # 开启对齐插值
        sta_base="sta_results",         # KE 时间轴来源

        output_base="extended_features_results",
        default_density=1000.0,
        dual_density_map={2: 1000.0, 3: 2000.0},
        nr_bins=8,
        nz_bins=6,
        min_particles_per_bin=5,
    )
    # =================================================

    print(f"\n配置:")
    print(f"  基础路径: {config.base_paths}")
    print(f"  采样: frame_skip={config.frame_skip}, start_time={config.steady_state_start_time}s")
    print(f"  时间轴对齐: sta_results/{'{type}_se{coh}'}/kinetic_energy_comparison.csv")
    print(f"  Cylindrical Binning: {config.nr_bins}(r) × {config.nz_bins}(z)")

    analyzer = ExtendedFeatureAnalyzer(config)
    analyzer.process_all()

    print("\n" + "=" * 60)
    print("分析完成! comparison CSV 时间轴已与 KE 完全一致")
    print("=" * 60)


if __name__ == "__main__":
    main()
