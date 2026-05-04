#!/usr/bin/env python3
"""
SUP Mixer Lacey Index Analysis v2 - Multi-Directory Comparison Version
=======================================================================
Analyzes mixing quality in DEM mixer simulations using Lacey mixing index
Supports multiple scale factors and surface energy values comparison

更新内容 (v2):
  - 支持多个基础目录路径 (base_paths: List[str])
  - 自动跨路径搜索、去重 (同 sta_other_info_v2 模式)
  - 保持: 原有全部分析功能不变

Supports two directory formats:
  - SUPMixerOutput_SizeDiff_f{scale_factor}se{surface_energy}
  - SUPMixerOutput_DualDensity_f{scale_factor}se{surface_energy}
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
from matplotlib.animation import FuncAnimation
import warnings
warnings.filterwarnings('ignore')

# ======================== Configuration Parameters ========================
class Config:
    """Configuration parameter management for Lacey index analysis"""
    def __init__(self, scale_factor: int = 4, surface_energy: str = "000",
                 input_path: Optional[Path] = None, base_path: str = ".",
                 enable_plotting: bool = False,
                 steady_state_fraction: float = 0.5, experiment_type: str = "SizeDiff"):
        self.scale_factor = scale_factor
        self.surface_energy = surface_energy
        self.experiment_type = experiment_type

        # 支持直接传入已发现的路径 (多目录模式)
        if input_path is not None:
            self.input_path = Path(input_path)
        else:
            self.base_path = Path(base_path)
            self.input_path = self.base_path / f'SUPMixerOutput_{experiment_type}_f{scale_factor}se{surface_energy}'

        # Output path includes experiment type
        self.output_base = Path(f'sta_results/lacey_results/{experiment_type}/se{surface_energy}/')
        self.output_path = self.output_base / f'f{scale_factor}/'

        if scale_factor == 1:
            self.frame_dt = 5e-3
        else:
            self.frame_dt = 1e-3

        self.RPM = 300
        self.frame_skip = 10
        self.max_sample_size = 100000
        self.chunk_size = 10000
        self.steady_state_fraction = steady_state_fraction

        # Particle densities (kg/m³)
        self.density_light = 1000.0  # Family 2
        self.density_heavy = 2000.0  # Family 3

        # Grid parameters
        self.grid_divisions = 5
        self.min_particles_per_cell = 10

        # Plotting
        self.enable_plotting = enable_plotting
        self.plot_dpi = 300
        self.plot_show = False

        if enable_plotting:
            plt.rcParams['font.family'] = 'serif'
            plt.rcParams['font.serif'] = ['DejaVu Serif', 'Times New Roman']
            plt.rcParams['font.size'] = 12
            plt.rcParams['axes.unicode_minus'] = False
            plt.rcParams['figure.dpi'] = 100
            plt.rcParams['savefig.dpi'] = self.plot_dpi
            plt.rcParams['axes.linewidth'] = 1.8
            plt.rcParams['xtick.major.width'] = 1.2
            plt.rcParams['ytick.major.width'] = 1.2
        else:
            matplotlib.use('Agg')

        self.output_path.mkdir(parents=True, exist_ok=True)

    def __str__(self):
        return f"Config({self.experiment_type}_f{self.scale_factor}se{self.surface_energy}, frame_dt={self.frame_dt*1000:.1f}ms)"


class MultiDirectoryConfig:
    """Configuration for multi-directory comparison — supports multiple base paths"""
    def __init__(self, base_paths=None, filter_cohesion: str = "000",
                 enable_plotting: bool = True, steady_state_fraction: float = 0.5):
        # 兼容旧接口
        if base_paths is None:
            base_paths = ["."]
        elif isinstance(base_paths, str):
            base_paths = [base_paths]
        self.base_paths = base_paths

        self.filter_cohesion = filter_cohesion
        self.enable_plotting = enable_plotting
        self.steady_state_fraction = steady_state_fraction
        self.directories = self._find_matching_directories()

    def _find_matching_directories(self) -> List[Dict]:
        """从多个基础路径查找所有匹配目录 (SizeDiff & DualDensity)，并去重"""
        dirs = []
        patterns = [
            ('SizeDiff', r'SUPMixerOutput_SizeDiff_f(\d+)se(\d+)$'),
            ('DualDensity', r'SUPMixerOutput_DualDensity_f(\d+)se(\d+)$'),
        ]

        for base_path in self.base_paths:
            bp = Path(base_path)
            if not bp.exists():
                print(f"警告：路径 {base_path} 不存在，跳过")
                continue

            for path in bp.iterdir():
                if not path.is_dir():
                    continue
                for exp_type, pattern in patterns:
                    match = re.match(pattern, path.name)
                    if match:
                        scale_factor = int(match.group(1))
                        surface_energy = match.group(2)
                        if self.filter_cohesion and surface_energy != self.filter_cohesion:
                            continue
                        dirs.append({
                            'path': path,
                            'scale_factor': scale_factor,
                            'surface_energy': surface_energy,
                            'experiment_type': exp_type,
                            'name': path.name,
                            'base_path': str(base_path)
                        })
                        break

        # 去重
        seen = set()
        unique = []
        for d in dirs:
            key = (d['experiment_type'], d['surface_energy'], d['scale_factor'])
            if key not in seen:
                seen.add(key)
                unique.append(d)

        return sorted(unique, key=lambda x: (x['experiment_type'], x['scale_factor'], x['surface_energy']))


# ======================== JIT Accelerated Functions ========================
@jit(nopython=True, parallel=True)
def calculate_particle_masses(radii: np.ndarray, families: np.ndarray,
                             density_light: float, density_heavy: float) -> np.ndarray:
    n = len(radii)
    masses = np.zeros(n)
    for i in prange(n):
        volume = (4.0/3.0) * np.pi * radii[i]**3
        if families[i] == 2:
            masses[i] = volume * density_light
        else:
            masses[i] = volume * density_heavy
    return masses

@jit(nopython=True)
def assign_particles_to_grid(x, y, z, x_min, x_max, y_min, y_max,
                            z_min, z_max, n_divisions):
    n = len(x)
    cell_indices = np.zeros(n, dtype=np.int32)
    dx = (x_max - x_min) / n_divisions
    dy = (y_max - y_min) / n_divisions
    dz = (z_max - z_min) / n_divisions
    for i in range(n):
        ix = min(int((x[i] - x_min) / dx), n_divisions - 1)
        iy = min(int((y[i] - y_min) / dy), n_divisions - 1)
        iz = min(int((z[i] - z_min) / dz), n_divisions - 1)
        ix = max(0, min(ix, n_divisions - 1))
        iy = max(0, min(iy, n_divisions - 1))
        iz = max(0, min(iz, n_divisions - 1))
        cell_indices[i] = ix * n_divisions * n_divisions + iy * n_divisions + iz
    return cell_indices

@jit(nopython=True)
def calculate_lacey_index(mass_fractions, particles_per_cell, min_particles):
    valid_mask = particles_per_cell >= min_particles
    valid_fractions = mass_fractions[valid_mask]
    if len(valid_fractions) < 2:
        return np.nan, np.nan, np.nan, np.nan, 0
    x_bar = np.mean(valid_fractions)
    n = len(valid_fractions)
    sum_sq_diff = 0.0
    for i in range(n):
        sum_sq_diff += (valid_fractions[i] - x_bar) ** 2
    S2 = sum_sq_diff / (n - 1)
    S02 = x_bar * (1.0 - x_bar)
    avg_particles = np.mean(particles_per_cell[valid_mask])
    SR2 = x_bar * (1.0 - x_bar) / avg_particles
    if S02 - SR2 > 1e-10:
        M = (S02 - S2) / (S02 - SR2)
    else:
        M = np.nan
    if not np.isnan(M):
        M = max(0.0, min(1.0, M))
    return M, S2, S02, SR2, len(valid_fractions)


# ======================== Data Processing Class ========================
class LaceyDataProcessor:
    def __init__(self, config: Config):
        self.config = config
        self.df = None
        self.lacey_history = []

    def extract_frame_number(self, filename: str) -> int:
        match = re.search(r'mixer_output_(\d+)\.csv', filename)
        return int(match.group(1)) if match else 0

    def get_files_to_read(self, for_time_series: bool = False) -> list:
        pattern = self.config.input_path / 'mixer_output_*.csv'
        files = sorted(glob.glob(str(pattern)))
        if not files:
            raise FileNotFoundError(f"No CSV files found in {self.config.input_path}")

        all_frames = []
        for file in files:
            frame_num = self.extract_frame_number(os.path.basename(file))
            all_frames.append((file, frame_num))
        all_frames.sort(key=lambda x: x[1])

        if for_time_series:
            files_to_read = [(f, n) for f, n in all_frames if n % self.config.frame_skip == 0]
        else:
            total_frames = len(all_frames)
            start_idx = int(total_frames * (1 - self.config.steady_state_fraction))
            steady_state_frames = all_frames[start_idx:]
            print(f"  ✔ Total frames found: {total_frames}")
            print(f"  ✔ Using last {self.config.steady_state_fraction:.0%} ({len(steady_state_frames)} frames) for steady state")
            files_to_read = [(f, n) for f, n in steady_state_frames if n % self.config.frame_skip == 0]

        return files_to_read

    def calculate_lacey_for_frame(self, df: pd.DataFrame) -> Dict:
        masses = calculate_particle_masses(
            df['r'].values, df['family'].values,
            self.config.density_light, self.config.density_heavy
        )
        df['mass'] = masses

        margin = 0.001
        x_min, x_max = df['X'].min() - margin, df['X'].max() + margin
        y_min, y_max = df['Y'].min() - margin, df['Y'].max() + margin
        z_min, z_max = df['Z'].min() - margin, df['Z'].max() + margin

        cell_indices = assign_particles_to_grid(
            df['X'].values, df['Y'].values, df['Z'].values,
            x_min, x_max, y_min, y_max, z_min, z_max,
            self.config.grid_divisions
        )
        df['cell'] = cell_indices

        n_cells = self.config.grid_divisions ** 3
        mass_fractions = np.zeros(n_cells)
        particles_per_cell = np.zeros(n_cells)

        for cell_id in range(n_cells):
            cell_particles = df[df['cell'] == cell_id]
            if len(cell_particles) > 0:
                particles_per_cell[cell_id] = len(cell_particles)
                heavy_mass = cell_particles[cell_particles['family'] == 3]['mass'].sum()
                total_mass = cell_particles['mass'].sum()
                if total_mass > 0:
                    mass_fractions[cell_id] = heavy_mass / total_mass

        M, S2, S02, SR2, valid_cells = calculate_lacey_index(
            mass_fractions, particles_per_cell, self.config.min_particles_per_cell
        )

        non_empty_cells = particles_per_cell[particles_per_cell > 0]

        return {
            'lacey_index': M,
            'actual_variance': S2,
            'segregated_variance': S02,
            'random_variance': SR2,
            'valid_cells': valid_cells,
            'total_cells': n_cells,
            'mean_mass_fraction': np.mean(mass_fractions[particles_per_cell > 0]) if len(non_empty_cells) > 0 else np.nan,
            'mean_particles_per_cell': np.mean(non_empty_cells) if len(non_empty_cells) > 0 else 0,
            'std_particles_per_cell': np.std(non_empty_cells) if len(non_empty_cells) > 0 else 0,
            'total_particles': int(np.sum(particles_per_cell))
        }

    def process_time_series(self) -> pd.DataFrame:
        files_to_read = self.get_files_to_read(for_time_series=True)
        print(f"  ✔ Processing {len(files_to_read)} frames for time series")

        lacey_data = []
        for i, (file_path, frame_num) in enumerate(files_to_read, 1):
            try:
                df = pd.read_csv(file_path)
                result = self.calculate_lacey_for_frame(df)
                result['frame'] = frame_num
                result['time'] = frame_num * self.config.frame_dt
                result['revolutions'] = result['time'] * self.config.RPM / 60.0
                lacey_data.append(result)
                if i % 10 == 0 or i == len(files_to_read):
                    print(f"    Processing: {(i / len(files_to_read)) * 100:5.1f}% ({i}/{len(files_to_read)})")
            except Exception as e:
                print(f"  ✗ Error processing frame {frame_num}: {e}")
                continue

        if not lacey_data:
            raise ValueError("No data successfully processed")

        self.lacey_history = pd.DataFrame(lacey_data)
        return self.lacey_history

    def process_steady_state(self) -> Dict:
        files_to_read = self.get_files_to_read(for_time_series=False)
        print(f"  ✔ Processing {len(files_to_read)} frames for steady state")

        lacey_values = []
        for i, (file_path, frame_num) in enumerate(files_to_read, 1):
            try:
                df = pd.read_csv(file_path)
                result = self.calculate_lacey_for_frame(df)
                lacey_values.append(result['lacey_index'])
                if i % 5 == 0 or i == len(files_to_read):
                    print(f"    Processing: {(i / len(files_to_read)) * 100:5.1f}% ({i}/{len(files_to_read)})")
            except Exception as e:
                print(f"  ✗ Error processing frame {frame_num}: {e}")
                continue

        lacey_values = np.array(lacey_values)
        lacey_values = lacey_values[~np.isnan(lacey_values)]

        if len(lacey_values) == 0:
            raise ValueError("No valid Lacey index values calculated")

        steady_state_result = {
            'mean_lacey': np.mean(lacey_values),
            'std_lacey': np.std(lacey_values),
            'min_lacey': np.min(lacey_values),
            'max_lacey': np.max(lacey_values),
            'median_lacey': np.median(lacey_values),
            'n_samples': len(lacey_values),
            'all_values': lacey_values
        }
        print(f"  ✔ Steady state Lacey index: {steady_state_result['mean_lacey']:.4f} ± {steady_state_result['std_lacey']:.4f}")
        return steady_state_result


# ======================== Visualization Class ========================
class LaceyVisualizer:
    def __init__(self, config: Config):
        self.config = config

    def plot_time_series(self, lacey_history: pd.DataFrame, output_path: Path):
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

        ax1.plot(lacey_history['time'], lacey_history['lacey_index'],
                'b-', linewidth=2, label='Lacey Index')
        ax1.axhline(y=0.95, color='g', linestyle='--', linewidth=1.5, label='Good Mixing (M=0.95)')
        ax1.axhline(y=0.75, color='orange', linestyle='--', linewidth=1.5, label='Moderate Mixing (M=0.75)')
        ax1.set_xlabel('Time [s]', fontsize=12, fontweight='bold')
        ax1.set_ylabel('Lacey Index', fontsize=12, fontweight='bold')
        ax1.set_title(f'Mixing Quality Evolution - {self.config.experiment_type} f{self.config.scale_factor} se{self.config.surface_energy}',
                     fontsize=14, fontweight='bold')
        ax1.grid(True, alpha=0.3)
        ax1.legend(loc='lower right')
        ax1.set_ylim([0, 1.05])

        ax2.plot(lacey_history['revolutions'], lacey_history['lacey_index'], 'r-', linewidth=2)
        ax2.axhline(y=0.95, color='g', linestyle='--', linewidth=1.5)
        ax2.axhline(y=0.75, color='orange', linestyle='--', linewidth=1.5)
        ax2.set_xlabel('Revolutions', fontsize=12, fontweight='bold')
        ax2.set_ylabel('Lacey Index', fontsize=12, fontweight='bold')
        ax2.set_title('Mixing Quality vs Rotor Revolutions', fontsize=14, fontweight='bold')
        ax2.grid(True, alpha=0.3)
        ax2.set_ylim([0, 1.05])

        steady_data = lacey_history[lacey_history['time'] > lacey_history['time'].max() * 0.5]
        if len(steady_data) > 0:
            mean_l = steady_data['lacey_index'].mean()
            std_l = steady_data['lacey_index'].std()
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            ax1.text(0.02, 0.98, f'Steady State:\nM = {mean_l:.3f} ± {std_l:.3f}',
                    transform=ax1.transAxes, fontsize=10, verticalalignment='top', bbox=props)

        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()

    def plot_comparison_bar_chart(self, results: Dict, output_path: Path):
        fig, ax = plt.subplots(figsize=(12, 8))
        configs, mean_values, std_values = [], [], []

        for name, result in results.items():
            if 'steady_state' in result:
                exp_type = result['dir_info']['experiment_type']
                scale_f = result['dir_info']['scale_factor']
                configs.append(f"{exp_type}_f{scale_f}")
                mean_values.append(result['steady_state']['mean_lacey'])
                std_values.append(result['steady_state']['std_lacey'])

        x_pos = np.arange(len(configs))
        colors = plt.cm.viridis(np.linspace(0, 0.9, len(configs)))
        bars = ax.bar(x_pos, mean_values, yerr=std_values, capsize=5,
                      color=colors, edgecolor='black', linewidth=1.5, alpha=0.8)

        ax.set_xlabel('Configuration', fontsize=14, fontweight='bold')
        ax.set_ylabel('Lacey Index', fontsize=14, fontweight='bold')
        ax.set_title('Steady State Lacey Index Comparison', fontsize=16, fontweight='bold')
        ax.set_xticks(x_pos)
        ax.set_xticklabels(configs, rotation=45, ha='right')
        ax.set_ylim([0, 1.1])
        ax.axhline(y=0.95, color='g', linestyle='--', linewidth=1.5, alpha=0.5, label='Excellent (M=0.95)')
        ax.axhline(y=0.85, color='orange', linestyle='--', linewidth=1.5, alpha=0.5, label='Good (M=0.85)')
        ax.axhline(y=0.75, color='red', linestyle='--', linewidth=1.5, alpha=0.5, label='Moderate (M=0.75)')
        ax.grid(True, axis='y', alpha=0.3)
        ax.legend(loc='upper right', frameon=True, fancybox=False, edgecolor='black', fontsize=10)

        for bar, mean, std in zip(bars, mean_values, std_values):
            ax.text(bar.get_x() + bar.get_width()/2., bar.get_height() + std + 0.02,
                    f'{mean:.3f}', ha='center', va='bottom', fontsize=10, fontweight='bold')

        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()

    def plot_comparison_time_series(self, results: Dict, output_path: Path):
        fig, ax = plt.subplots(figsize=(14, 8))
        if any('lacey_history' in r for r in results.values()):
            colors = plt.cm.viridis(np.linspace(0, 0.9, len(results)))
            for (name, result), color in zip(results.items(), colors):
                if 'lacey_history' in result:
                    history = result['lacey_history']
                    label = f"{result['dir_info']['experiment_type']}_f{result['dir_info']['scale_factor']}"
                    ax.plot(history['time'], history['lacey_index'],
                           '-', linewidth=2.4, color=color, label=label, alpha=0.8)
            ax.set_xlabel('Time [s]', fontsize=14, fontweight='bold')
            ax.set_ylabel('Lacey Index', fontsize=14, fontweight='bold')
            ax.set_title('Mixing Evolution Comparison', fontsize=16, fontweight='bold')
            ax.axhline(y=0.95, color='g', linestyle='--', linewidth=1.5, alpha=0.5)
            ax.axhline(y=0.85, color='orange', linestyle='--', linewidth=1.5, alpha=0.5)
            ax.axhline(y=0.75, color='red', linestyle='--', linewidth=1.5, alpha=0.5)
            ax.legend(loc='lower right', frameon=True, fancybox=False, edgecolor='black', fontsize=10)
            ax.grid(True, alpha=0.3)
            ax.set_ylim([0, 1.05])
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()


# ======================== Multi-Directory Comparison Class ========================
class MultiDirectoryLaceyAnalysis:
    def __init__(self, multi_config: MultiDirectoryConfig):
        self.multi_config = multi_config
        self.results = {}

    def process_all_directories(self):
        print("\n" + "=" * 60)
        print("MULTI-DIRECTORY LACEY INDEX ANALYSIS")
        print("(DualDensity & SizeDiff experiments only)")
        print("=" * 60)
        print(f"Found {len(self.multi_config.directories)} directories to process")

        for dir_info in self.multi_config.directories:
            print(f"\n{'='*60}")
            print(f"Processing: {dir_info['name']}  [{dir_info['base_path']}]")
            print(f"  Experiment Type: {dir_info['experiment_type']}")
            print(f"  Scale Factor: f{dir_info['scale_factor']}")
            print(f"  Surface Energy: se{dir_info['surface_energy']}")
            print('='*60)

            try:
                config = Config(
                    scale_factor=dir_info['scale_factor'],
                    surface_energy=dir_info['surface_energy'],
                    input_path=dir_info['path'],  # 直接传入已发现的路径
                    enable_plotting=self.multi_config.enable_plotting,
                    steady_state_fraction=self.multi_config.steady_state_fraction,
                    experiment_type=dir_info['experiment_type']
                )

                processor = LaceyDataProcessor(config)

                print("\n  Processing time series...")
                lacey_history = processor.process_time_series()

                print("\n  Processing steady state...")
                steady_state = processor.process_steady_state()

                self.results[dir_info['name']] = {
                    'config': config,
                    'dir_info': dir_info,
                    'lacey_history': lacey_history,
                    'steady_state': steady_state
                }

                self.save_individual_results(dir_info['name'])

                if self.multi_config.enable_plotting:
                    visualizer = LaceyVisualizer(config)
                    plot_path = config.output_path / f'lacey_time_series_{config.experiment_type}_f{config.scale_factor}c{config.surface_energy}.png'
                    visualizer.plot_time_series(lacey_history, plot_path)
                    print(f"  ✔ Plot saved to: {plot_path}")

            except Exception as e:
                print(f"  ✗ Error processing {dir_info['name']}: {e}")
                import traceback
                traceback.print_exc()
                continue

        if len(self.results) > 1:
            self.generate_comparison()

    def save_individual_results(self, dir_name: str):
        result = self.results[dir_name]
        config = result['config']
        exp_type = config.experiment_type

        result['lacey_history'].to_csv(
            config.output_path / f'lacey_time_series_{exp_type}_f{config.scale_factor}c{config.surface_energy}.csv',
            index=False
        )

        steady_state = result['steady_state']
        with open(config.output_path / f'lacey_summary_{exp_type}_f{config.scale_factor}c{config.surface_energy}.txt', 'w') as f:
            f.write(f"Lacey Index Analysis Summary\n")
            f.write(f"Experiment Type: {exp_type}\n")
            f.write(f"Scale Factor: f{config.scale_factor}\n")
            f.write(f"Surface Energy: se{config.surface_energy}\n")
            f.write(f"Frame time interval: {config.frame_dt*1000:.1f} ms\n")
            f.write("=" * 60 + "\n\n")
            f.write(f"Steady State Analysis (last {config.steady_state_fraction:.0%} of data):\n")
            f.write(f"  Mean Lacey Index:   {steady_state['mean_lacey']:.4f}\n")
            f.write(f"  Std Deviation:      {steady_state['std_lacey']:.4f}\n")
            f.write(f"  Min Lacey Index:    {steady_state['min_lacey']:.4f}\n")
            f.write(f"  Max Lacey Index:    {steady_state['max_lacey']:.4f}\n")
            f.write(f"  Median Lacey Index: {steady_state['median_lacey']:.4f}\n")
            f.write(f"  Number of samples:  {steady_state['n_samples']}\n")
            f.write("\nGrid Configuration:\n")
            f.write(f"  Grid divisions: {config.grid_divisions}x{config.grid_divisions}x{config.grid_divisions}\n")
            f.write(f"  Min particles per cell: {config.min_particles_per_cell}\n")
            f.write("\nParticle Properties:\n")
            f.write(f"  Light particle density: {config.density_light} kg/m³\n")
            f.write(f"  Heavy particle density: {config.density_heavy} kg/m³\n")
            f.write("\nMixing Quality Assessment:\n")
            ml = steady_state['mean_lacey']
            if ml >= 0.95:
                f.write("  ★★★★★ Excellent mixing (M >= 0.95)\n")
            elif ml >= 0.85:
                f.write("  ★★★★☆ Good mixing (0.85 <= M < 0.95)\n")
            elif ml >= 0.75:
                f.write("  ★★★☆☆ Moderate mixing (0.75 <= M < 0.85)\n")
            elif ml >= 0.65:
                f.write("  ★★☆☆☆ Poor mixing (0.65 <= M < 0.75)\n")
            else:
                f.write("  ★☆☆☆☆ Very poor mixing (M < 0.65)\n")

        print(f"  ✔ Results saved to: {config.output_path}")

    def generate_comparison(self):
        print("\n" + "=" * 60)
        print("GENERATING COMPARISON ANALYSIS")
        print("=" * 60)

        if self.multi_config.filter_cohesion:
            comparison_dir = Path(f'sta_results/lacey_results/comparison_se{self.multi_config.filter_cohesion}/')
        else:
            comparison_dir = Path('sta_results/lacey_results/comparison/')
        comparison_dir.mkdir(parents=True, exist_ok=True)

        comparison_data = []
        for name, result in self.results.items():
            steady_state = result['steady_state']
            comparison_data.append({
                'Directory': name,
                'Experiment_Type': result['dir_info']['experiment_type'],
                'Scale_Factor': result['dir_info']['scale_factor'],
                'Surface_Energy': result['dir_info']['surface_energy'],
                'Frame_Interval_ms': result['config'].frame_dt * 1000,
                'Mean_Lacey': steady_state['mean_lacey'],
                'Std_Lacey': steady_state['std_lacey'],
                'Min_Lacey': steady_state['min_lacey'],
                'Max_Lacey': steady_state['max_lacey'],
                'Median_Lacey': steady_state['median_lacey'],
                'N_Samples': steady_state['n_samples']
            })

        comparison_df = pd.DataFrame(comparison_data)
        comparison_df = comparison_df.sort_values(['Experiment_Type', 'Scale_Factor', 'Surface_Energy'])
        comparison_df.to_csv(comparison_dir / 'lacey_comparison_DualDensity_SizeDiff.csv', index=False)
        print(f"  ✔ Comparison data saved")

        if self.multi_config.enable_plotting:
            visualizer = LaceyVisualizer(Config(experiment_type='SizeDiff'))
            visualizer.plot_comparison_bar_chart(
                self.results, comparison_dir / 'lacey_comparison_bar.png'
            )
            visualizer.plot_comparison_time_series(
                self.results, comparison_dir / 'lacey_comparison_time_series.png'
            )
            print(f"  ✔ Comparison plots saved to: {comparison_dir}")


# ======================== Main Program ========================
def main():
    # ==================== 用户配置区 ====================
    base_paths = [
        "/media/huyuze/Fanxiang1/DEME_Data/backup251115",
        "/media/huyuze/Fanxiang/backup251025",
        "/media/huyuze/Fanxiang/backup251109",
    ]
    filter_cohesion = "015"       # None = 处理所有
    enable_plotting = True
    steady_state_fraction = 1
    # ===================================================

    print("=" * 60)
    print("SUP Mixer Lacey Index Analysis v2")
    print("支持多基础路径搜索")
    print("=" * 60)
    print(f"\n配置:")
    print(f"  基础路径: {base_paths}")
    print(f"  Cohesion过滤: {filter_cohesion or '全部'}")
    print(f"  稳态比例: {steady_state_fraction}")

    start_time = datetime.datetime.now()

    try:
        multi_config = MultiDirectoryConfig(
            base_paths=base_paths,
            filter_cohesion=filter_cohesion,
            enable_plotting=enable_plotting,
            steady_state_fraction=steady_state_fraction
        )

        if not multi_config.directories:
            print("✗ No matching directories found!")
            return

        print(f"\n找到 {len(multi_config.directories)} 个目录:")
        for d in multi_config.directories:
            frame_dt_ms = 5.0 if d['scale_factor'] == 1 else 1.0
            print(f"  - {d['name']} (dt={frame_dt_ms}ms)  [{d['base_path']}]")

        analysis = MultiDirectoryLaceyAnalysis(multi_config)
        analysis.process_all_directories()

        elapsed = (datetime.datetime.now() - start_time).total_seconds()
        print(f"\n{'='*60}")
        print(f"完成! 处理 {len(analysis.results)} 个目录, 耗时 {elapsed:.1f}s")

        print("\n  Mixing Quality Summary:")
        for name, result in analysis.results.items():
            ss = result['steady_state']
            quality = "★" * min(5, max(1, int(ss['mean_lacey'] * 5 + 0.5)))
            print(f"    {name}: M = {ss['mean_lacey']:.3f} {quality}")

        print("=" * 60)

    except Exception as e:
        print(f"\n✗ ERROR: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
    