#!/usr/bin/env python3
"""
SUP Mixer Lacey Index Analysis - Multi-Directory Comparison Version
Analyzes mixing quality in DEM mixer simulations using Lacey mixing index
Supports multiple scale factors and surface energy values comparison
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
                 base_path: str = ".", enable_plotting: bool = False,
                 steady_state_fraction: float = 0.5):
        self.scale_factor = scale_factor
        self.surface_energy = surface_energy
        self.base_path = Path(base_path)
        
        # Construct input and output paths
        self.input_path = self.base_path / f'SUPMixerOutput_f{scale_factor}se{surface_energy}'
        self.output_base = Path(f'lacey_results/se{surface_energy}/')
        self.output_path = self.output_base / f'f{scale_factor}/'
        
        # Frame time interval based on scale factor
        if scale_factor == 1:
            self.frame_dt = 5e-3  # f1 using 5ms interval
        else:
            self.frame_dt = 1e-3  # others using 1ms interval

        # Simulation parameters
        self.RPM = 300  # Mixer rotation speed
        self.frame_skip = 10  # Read every nth frame for Lacey calculation
        self.max_sample_size = 100000  # Maximum sampling rows
        self.chunk_size = 10000  # Chunk reading size
        self.steady_state_fraction = steady_state_fraction
        
        # Particle densities (kg/m³)
        self.density_light = 1000.0  # Family 2
        self.density_heavy = 2000.0  # Family 3
        
        # Grid parameters for Lacey index calculation
        self.grid_divisions = 5  # 5x5x5 grid by default
        self.min_particles_per_cell = 10  # Minimum particles for valid cell
        
        # Plotting control parameters
        self.enable_plotting = enable_plotting
        self.plot_dpi = 300
        self.plot_show = False
        
        # Set matplotlib parameters for scientific style
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
        
        # Create output directories
        self.output_path.mkdir(parents=True, exist_ok=True)
        
    def __str__(self):
        return f"Config(f{self.scale_factor}se{self.surface_energy}, frame_dt={self.frame_dt*1000:.1f}ms)"

class MultiDirectoryConfig:
    """Configuration for multi-directory comparison analysis"""
    def __init__(self, base_path: str = ".", filter_cohesion: str = "000",
                 enable_plotting: bool = True, steady_state_fraction: float = 0.5):
        self.base_path = Path(base_path)
        self.filter_cohesion = filter_cohesion
        self.enable_plotting = enable_plotting
        self.steady_state_fraction = steady_state_fraction
        
        # Find all matching directories
        self.directories = self.find_matching_directories()
        
    def find_matching_directories(self) -> List[Dict]:
        """Find all directories matching the pattern"""
        pattern = f'SUPMixerOutput_f*se*'
        if self.filter_cohesion:
            pattern = f'SUPMixerOutput_f*se{self.filter_cohesion}'
        
        dirs = []
        for path in self.base_path.glob(pattern):
            if path.is_dir():
                match = re.match(r'SUPMixerOutput_f(\d+)se(\d+)', path.name)
                if match:
                    dirs.append({
                        'path': path,
                        'scale_factor': int(match.group(1)),
                        'surface_energy': match.group(2),
                        'name': path.name
                    })
        
        return sorted(dirs, key=lambda x: (x['scale_factor'], x['surface_energy']))

# ======================== JIT Accelerated Functions ========================
@jit(nopython=True, parallel=True)
def calculate_particle_masses(radii: np.ndarray, families: np.ndarray,
                             density_light: float, density_heavy: float) -> np.ndarray:
    """Calculate particle masses based on radius and family"""
    n = len(radii)
    masses = np.zeros(n)
    
    for i in prange(n):
        volume = (4.0/3.0) * np.pi * radii[i]**3
        if families[i] == 2:  # Light particle
            masses[i] = volume * density_light
        else:  # Heavy particle (family == 3)
            masses[i] = volume * density_heavy
    
    return masses

@jit(nopython=True)
def assign_particles_to_grid(x: np.ndarray, y: np.ndarray, z: np.ndarray,
                            x_min: float, x_max: float,
                            y_min: float, y_max: float,
                            z_min: float, z_max: float,
                            n_divisions: int) -> np.ndarray:
    """Assign particles to grid cells"""
    n = len(x)
    cell_indices = np.zeros(n, dtype=np.int32)
    
    dx = (x_max - x_min) / n_divisions
    dy = (y_max - y_min) / n_divisions
    dz = (z_max - z_min) / n_divisions
    
    for i in range(n):
        ix = min(int((x[i] - x_min) / dx), n_divisions - 1)
        iy = min(int((y[i] - y_min) / dy), n_divisions - 1)
        iz = min(int((z[i] - z_min) / dz), n_divisions - 1)
        
        # Ensure indices are within bounds
        ix = max(0, min(ix, n_divisions - 1))
        iy = max(0, min(iy, n_divisions - 1))
        iz = max(0, min(iz, n_divisions - 1))
        
        cell_indices[i] = ix * n_divisions * n_divisions + iy * n_divisions + iz
    
    return cell_indices

@jit(nopython=True)
def calculate_lacey_index(mass_fractions: np.ndarray, 
                         particles_per_cell: np.ndarray,
                         min_particles: int) -> Tuple[float, float, float, float, int]:
    """Calculate Lacey mixing index from mass fractions"""
    # Filter out cells with too few particles
    valid_mask = particles_per_cell >= min_particles
    valid_fractions = mass_fractions[valid_mask]
    
    if len(valid_fractions) < 2:
        return np.nan, np.nan, np.nan, np.nan, 0
    
    # Calculate mean mass fraction
    x_bar = np.mean(valid_fractions)
    
    # Calculate actual variance S² with ddof=1 (sample variance)
    # Manual calculation since numba doesn't support ddof parameter
    n = len(valid_fractions)
    sum_sq_diff = 0.0
    for i in range(n):
        sum_sq_diff += (valid_fractions[i] - x_bar) ** 2
    S2 = sum_sq_diff / (n - 1)  # Using n-1 for sample variance (ddof=1)
    
    # Calculate variance for complete segregation S₀²
    S02 = x_bar * (1.0 - x_bar)
    
    # Calculate variance for random mixing S_R²
    # Assuming average particles per cell for simplification
    avg_particles = np.mean(particles_per_cell[valid_mask])
    SR2 = x_bar * (1.0 - x_bar) / avg_particles
    
    # Calculate Lacey index
    if S02 - SR2 > 1e-10:  # Avoid division by zero
        M = (S02 - S2) / (S02 - SR2)
    else:
        M = np.nan
    
    # Ensure M is in [0, 1]
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
        """Extract frame number from filename"""
        match = re.search(r'mixer_output_(\d+)\.csv', filename)
        return int(match.group(1)) if match else 0
    
    def get_files_to_read(self, for_time_series: bool = False) -> list:
        """Get list of files to read"""
        pattern = self.config.input_path / 'mixer_output_*.csv'
        files = sorted(glob.glob(str(pattern)))
        
        if not files:
            raise FileNotFoundError(f"No CSV files found in {self.config.input_path}")
        
        # Extract all frame numbers
        all_frames = []
        for file in files:
            frame_num = self.extract_frame_number(os.path.basename(file))
            all_frames.append((file, frame_num))
        
        # Sort by frame number
        all_frames.sort(key=lambda x: x[1])
        
        if for_time_series:
            # For time series, use all frames with frame_skip
            files_to_read = [(f, n) for f, n in all_frames if n % self.config.frame_skip == 0]
        else:
            # For steady state analysis
            total_frames = len(all_frames)
            start_idx = int(total_frames * (1 - self.config.steady_state_fraction))
            steady_state_frames = all_frames[start_idx:]
            
            print(f"  ✔ Total frames found: {total_frames}")
            print(f"  ✔ Using last {self.config.steady_state_fraction:.0%} ({len(steady_state_frames)} frames) for steady state")
            
            files_to_read = [(f, n) for f, n in steady_state_frames if n % self.config.frame_skip == 0]
        
        return files_to_read
    
    def calculate_lacey_for_frame(self, df: pd.DataFrame) -> Dict:
        """Calculate Lacey index for a single frame"""
        # Calculate particle masses
        masses = calculate_particle_masses(
            df['r'].values,
            df['family'].values,
            self.config.density_light,
            self.config.density_heavy
        )
        df['mass'] = masses
        
        # Determine grid boundaries
        x_min, x_max = df['X'].min(), df['X'].max()
        y_min, y_max = df['Y'].min(), df['Y'].max()
        z_min, z_max = df['Z'].min(), df['Z'].max()
        
        # Add small margin to avoid edge issues
        margin = 0.001
        x_min -= margin
        x_max += margin
        y_min -= margin
        y_max += margin
        z_min -= margin
        z_max += margin
        
        # Assign particles to grid cells
        cell_indices = assign_particles_to_grid(
            df['X'].values, df['Y'].values, df['Z'].values,
            x_min, x_max, y_min, y_max, z_min, z_max,
            self.config.grid_divisions
        )
        df['cell'] = cell_indices
        
        # Calculate mass fractions for each cell
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
        
        # Calculate Lacey index
        M, S2, S02, SR2, valid_cells = calculate_lacey_index(
            mass_fractions, particles_per_cell, self.config.min_particles_per_cell
        )
        
        return {
            'lacey_index': M,
            'actual_variance': S2,
            'segregated_variance': S02,
            'random_variance': SR2,
            'valid_cells': valid_cells,
            'total_cells': n_cells,
            'mean_mass_fraction': np.mean(mass_fractions[particles_per_cell > 0]),
            'particles_per_cell': particles_per_cell
        }
    
    def process_time_series(self) -> pd.DataFrame:
        """Process Lacey index over time"""
        files_to_read = self.get_files_to_read(for_time_series=True)
        
        print(f"  ✔ Processing {len(files_to_read)} frames for time series")
        print(f"  ✔ Reading every {self.config.frame_skip} frames")
        
        lacey_data = []
        
        for i, (file_path, frame_num) in enumerate(files_to_read, 1):
            try:
                df = pd.read_csv(file_path)
                
                # Calculate Lacey index for this frame
                result = self.calculate_lacey_for_frame(df)
                
                # Add time information
                result['frame'] = frame_num
                result['time'] = frame_num * self.config.frame_dt
                result['revolutions'] = result['time'] * self.config.RPM / 60.0
                
                lacey_data.append(result)
                
                # Progress indicator
                if i % 10 == 0 or i == len(files_to_read):
                    progress = (i / len(files_to_read)) * 100
                    print(f"    Processing: {progress:5.1f}% complete ({i}/{len(files_to_read)})")
                    
            except Exception as e:
                print(f"  ✗ Error processing frame {frame_num}: {e}")
                continue
        
        if not lacey_data:
            raise ValueError("No data successfully processed")
        
        # Convert to DataFrame
        self.lacey_history = pd.DataFrame(lacey_data)
        print(f"  ✔ Time series analysis complete")
        
        return self.lacey_history
    
    def process_steady_state(self) -> Dict:
        """Process steady state Lacey index"""
        files_to_read = self.get_files_to_read(for_time_series=False)
        
        print(f"  ✔ Processing {len(files_to_read)} frames for steady state")
        
        lacey_values = []
        
        for i, (file_path, frame_num) in enumerate(files_to_read, 1):
            try:
                df = pd.read_csv(file_path)
                result = self.calculate_lacey_for_frame(df)
                lacey_values.append(result['lacey_index'])
                
                # Progress indicator
                if i % 5 == 0 or i == len(files_to_read):
                    progress = (i / len(files_to_read)) * 100
                    print(f"    Processing: {progress:5.1f}% complete ({i}/{len(files_to_read)})")
                    
            except Exception as e:
                print(f"  ✗ Error processing frame {frame_num}: {e}")
                continue
        
        # Calculate steady state statistics
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
        """Plot Lacey index time series"""
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
        
        # Plot 1: Lacey Index vs Time
        ax1.plot(lacey_history['time'], lacey_history['lacey_index'], 
                'b-', linewidth=2, label='Lacey Index')
        ax1.axhline(y=0.95, color='g', linestyle='--', linewidth=1.5, 
                   label='Good Mixing (M=0.95)')
        ax1.axhline(y=0.75, color='orange', linestyle='--', linewidth=1.5,
                   label='Moderate Mixing (M=0.75)')
        ax1.set_xlabel('Time [s]', fontsize=12, fontweight='bold')
        ax1.set_ylabel('Lacey Index', fontsize=12, fontweight='bold')
        ax1.set_title(f'Mixing Quality Evolution - f{self.config.scale_factor} se{self.config.surface_energy}',
                     fontsize=14, fontweight='bold')
        ax1.grid(True, alpha=0.3)
        ax1.legend(loc='lower right')
        ax1.set_ylim([0, 1.05])
        
        # Plot 2: Lacey Index vs Revolutions
        ax2.plot(lacey_history['revolutions'], lacey_history['lacey_index'],
                'r-', linewidth=2)
        ax2.axhline(y=0.95, color='g', linestyle='--', linewidth=1.5)
        ax2.axhline(y=0.75, color='orange', linestyle='--', linewidth=1.5)
        ax2.set_xlabel('Revolutions', fontsize=12, fontweight='bold')
        ax2.set_ylabel('Lacey Index', fontsize=12, fontweight='bold')
        ax2.set_title('Mixing Quality vs Rotor Revolutions', fontsize=14, fontweight='bold')
        ax2.grid(True, alpha=0.3)
        ax2.set_ylim([0, 1.05])
        
        # Add text box with statistics
        steady_state_data = lacey_history[lacey_history['time'] > lacey_history['time'].max() * 0.5]
        if len(steady_state_data) > 0:
            mean_lacey = steady_state_data['lacey_index'].mean()
            std_lacey = steady_state_data['lacey_index'].std()
            textstr = f'Steady State:\nM = {mean_lacey:.3f} ± {std_lacey:.3f}'
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            ax1.text(0.02, 0.98, textstr, transform=ax1.transAxes, fontsize=10,
                    verticalalignment='top', bbox=props)
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
    def plot_comparison_bar_chart(self, results: Dict, output_path: Path):
        """Plot bar chart comparison of steady state Lacey indices"""
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Prepare data
        configs = []
        mean_values = []
        std_values = []
        
        for name, result in results.items():
            if 'steady_state' in result:
                configs.append(name.replace('SUPMixerOutput_', ''))
                mean_values.append(result['steady_state']['mean_lacey'])
                std_values.append(result['steady_state']['std_lacey'])
        
        x_pos = np.arange(len(configs))
        colors = plt.cm.viridis(np.linspace(0, 0.9, len(configs)))
        
        # Bar chart with error bars
        bars = ax.bar(x_pos, mean_values, yerr=std_values, capsize=5,
                      color=colors, edgecolor='black', linewidth=1.5, alpha=0.8)
        
        ax.set_xlabel('Configuration', fontsize=14, fontweight='bold')
        ax.set_ylabel('Lacey Index', fontsize=14, fontweight='bold')
        ax.set_title('Steady State Lacey Index Comparison', fontsize=16, fontweight='bold')
        ax.set_xticks(x_pos)
        ax.set_xticklabels(configs, rotation=45, ha='right')
        ax.set_ylim([0, 1.1])
        
        # Add reference lines
        ax.axhline(y=0.95, color='g', linestyle='--', linewidth=1.5, alpha=0.5, label='Excellent (M=0.95)')
        ax.axhline(y=0.85, color='orange', linestyle='--', linewidth=1.5, alpha=0.5, label='Good (M=0.85)')
        ax.axhline(y=0.75, color='red', linestyle='--', linewidth=1.5, alpha=0.5, label='Moderate (M=0.75)')
        
        ax.grid(True, axis='y', alpha=0.3, linewidth=0.8)
        ax.legend(loc='upper right', frameon=True, fancybox=False, edgecolor='black', fontsize=10)
        
        # Add value labels on bars
        for bar, mean, std in zip(bars, mean_values, std_values):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + std + 0.02,
                    f'{mean:.3f}', ha='center', va='bottom', fontsize=10, fontweight='bold')
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
    def plot_comparison_time_series(self, results: Dict, output_path: Path):
        """Plot time series comparison of Lacey indices"""
        fig, ax = plt.subplots(figsize=(14, 8))
        
        if any('lacey_history' in result for result in results.values()):
            colors = plt.cm.viridis(np.linspace(0, 0.9, len(results)))
            
            for (name, result), color in zip(results.items(), colors):
                if 'lacey_history' in result:
                    history = result['lacey_history']
                    label = name.replace('SUPMixerOutput_', '')
                    ax.plot(history['time'], history['lacey_index'],
                           '-', linewidth=2.4, color=color, label=label, alpha=0.8)
            
            ax.set_xlabel('Time [s]', fontsize=14, fontweight='bold')
            ax.set_ylabel('Lacey Index', fontsize=14, fontweight='bold')
            ax.set_title('Mixing Evolution Comparison - Multiple Configurations', 
                        fontsize=16, fontweight='bold')
            
            # Add reference lines
            ax.axhline(y=0.95, color='g', linestyle='--', linewidth=1.5, alpha=0.5)
            ax.axhline(y=0.85, color='orange', linestyle='--', linewidth=1.5, alpha=0.5)
            ax.axhline(y=0.75, color='red', linestyle='--', linewidth=1.5, alpha=0.5)
            
            ax.legend(loc='lower right', frameon=True, fancybox=False,
                     edgecolor='black', fontsize=10)
            ax.grid(True, alpha=0.3, linewidth=0.8)
            ax.set_ylim([0, 1.05])
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()

# ======================== Multi-Directory Comparison Class ========================
class MultiDirectoryLaceyAnalysis:
    """Handle Lacey index analysis across multiple directories"""
    def __init__(self, multi_config: MultiDirectoryConfig):
        self.multi_config = multi_config
        self.results = {}
        
    def process_all_directories(self):
        """Process all matching directories"""
        print("\n" + "=" * 60)
        print("MULTI-DIRECTORY LACEY INDEX ANALYSIS")
        print("=" * 60)
        print(f"Found {len(self.multi_config.directories)} directories to process")
        
        for dir_info in self.multi_config.directories:
            print(f"\n{'='*60}")
            print(f"Processing: {dir_info['name']}")
            print(f"  Scale Factor: f{dir_info['scale_factor']}")
            print(f"  Surface Energy: se{dir_info['surface_energy']}")
            print('='*60)
            
            try:
                # Create config for this directory
                config = Config(
                    scale_factor=dir_info['scale_factor'],
                    surface_energy=dir_info['surface_energy'],
                    base_path=self.multi_config.base_path,
                    enable_plotting=self.multi_config.enable_plotting,
                    steady_state_fraction=self.multi_config.steady_state_fraction
                )
                
                # Process directory
                processor = LaceyDataProcessor(config)
                
                # Calculate time series
                print("\n  Processing time series...")
                lacey_history = processor.process_time_series()
                
                # Calculate steady state
                print("\n  Processing steady state...")
                steady_state = processor.process_steady_state()
                
                # Store results
                self.results[dir_info['name']] = {
                    'config': config,
                    'dir_info': dir_info,
                    'lacey_history': lacey_history,
                    'steady_state': steady_state
                }
                
                # Save individual results
                self.save_individual_results(dir_info['name'])
                
                # Generate plots if enabled
                if self.multi_config.enable_plotting:
                    visualizer = LaceyVisualizer(config)
                    plot_path = config.output_path / f'lacey_time_series_f{config.scale_factor}.png'
                    visualizer.plot_time_series(lacey_history, plot_path)
                    print(f"  ✔ Plot saved to: {plot_path}")
                
            except Exception as e:
                print(f"  ✗ Error processing {dir_info['name']}: {e}")
                import traceback
                traceback.print_exc()
                continue
        
        # Generate comparison results
        if len(self.results) > 1:
            self.generate_comparison()
    
    def save_individual_results(self, dir_name: str):
        """Save results for individual directory"""
        result = self.results[dir_name]
        config = result['config']
        
        # Save time series data
        result['lacey_history'].to_csv(
            config.output_path / f'lacey_time_series_f{config.scale_factor}.csv',
            index=False
        )
        
        # Save steady state summary
        steady_state = result['steady_state']
        with open(config.output_path / f'lacey_summary_f{config.scale_factor}.txt', 'w') as f:
            f.write(f"Lacey Index Analysis Summary - Scale Factor f{config.scale_factor}\n")
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
            f.write(f"  Grid divisions: {config.grid_divisions}×{config.grid_divisions}×{config.grid_divisions}\n")
            f.write(f"  Min particles per cell: {config.min_particles_per_cell}\n")
            f.write("\nParticle Properties:\n")
            f.write(f"  Light particle density: {config.density_light} kg/m³\n")
            f.write(f"  Heavy particle density: {config.density_heavy} kg/m³\n")
            
            # Add mixing quality assessment
            f.write("\nMixing Quality Assessment:\n")
            if steady_state['mean_lacey'] >= 0.95:
                f.write("  ★★★★★ Excellent mixing (M ≥ 0.95)\n")
            elif steady_state['mean_lacey'] >= 0.85:
                f.write("  ★★★★☆ Good mixing (0.85 ≤ M < 0.95)\n")
            elif steady_state['mean_lacey'] >= 0.75:
                f.write("  ★★★☆☆ Moderate mixing (0.75 ≤ M < 0.85)\n")
            elif steady_state['mean_lacey'] >= 0.65:
                f.write("  ★★☆☆☆ Poor mixing (0.65 ≤ M < 0.75)\n")
            else:
                f.write("  ★☆☆☆☆ Very poor mixing (M < 0.65)\n")
        
        print(f"  ✔ Results saved to: {config.output_path}")
    
    def generate_comparison(self):
        """Generate comparison analysis"""
        print("\n" + "=" * 60)
        print("GENERATING COMPARISON ANALYSIS")
        print("=" * 60)
        
        # Create comparison directory
        if self.multi_config.filter_cohesion:
            comparison_dir = Path(f'lacey_results/se{self.multi_config.filter_cohesion}/comparison/')
        else:
            comparison_dir = Path('lacey_results/comparison/')
        comparison_dir.mkdir(parents=True, exist_ok=True)
        
        # Collect comparison data
        comparison_data = []
        for name, result in self.results.items():
            steady_state = result['steady_state']
            comparison_data.append({
                'Directory': name,
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
        
        # Save comparison CSV
        comparison_df = pd.DataFrame(comparison_data)
        comparison_df = comparison_df.sort_values(['Scale_Factor', 'Surface_Energy'])
        comparison_df.to_csv(comparison_dir / 'lacey_comparison.csv', index=False)
        print(f"  ✔ Comparison data saved to: {comparison_dir / 'lacey_comparison.csv'}")
        
        # Generate comparison plots
        if self.multi_config.enable_plotting:
            visualizer = LaceyVisualizer(Config())  # Use default config for plotting
            
            # Generate separate plots
            visualizer.plot_comparison_bar_chart(
                self.results, 
                comparison_dir / 'lacey_comparison_bar.png'
            )
            print(f"  ✔ Bar chart saved to: {comparison_dir / 'lacey_comparison_bar.png'}")
            
            visualizer.plot_comparison_time_series(
                self.results,
                comparison_dir / 'lacey_comparison_time_series.png'
            )
            print(f"  ✔ Time series comparison saved to: {comparison_dir / 'lacey_comparison_time_series.png'}")

# ======================== Print Functions ========================
def print_header(title: str, width: int = 60):
    """Print formatted header"""
    print("\n" + "=" * width)
    print(title)
    print("=" * width)

def print_section(title: str, width: int = 40):
    """Print formatted section"""
    print("\n" + title)
    print("-" * width)

# ======================== Main Program ========================
def main():
    """
    Main program for Lacey index analysis
    """
    # ==================== Configurable Parameters ====================
    base_path = "build"  # Search directory
    filter_cohesion = "020"  # Only process cohesion=0 data, set to None for all
    enable_plotting = True  # Generate plots
    steady_state_fraction = 1  # Use last 50% of data for steady state
    
    # ==================== Run Analysis ====================
    print_header("SUP MIXER LACEY INDEX ANALYSIS")
    print("Version: Multi-directory comparison with time series")
    start_time = datetime.datetime.now()
    
    try:
        # Initialize multi-directory configuration
        multi_config = MultiDirectoryConfig(
            base_path=base_path,
            filter_cohesion=filter_cohesion,
            enable_plotting=enable_plotting,
            steady_state_fraction=steady_state_fraction
        )
        
        if not multi_config.directories:
            print("✗ No matching directories found!")
            print(f"  Searched for: SUPMixerOutput_f*se{filter_cohesion if filter_cohesion else '*'}")
            print(f"  In path: {base_path}")
            return
        
        print(f"Configuration:")
        print(f"  Base path: {base_path}")
        print(f"  Filter cohesion: {filter_cohesion if filter_cohesion else 'None (process all)'}")
        print(f"  Steady state fraction: {steady_state_fraction:.0%}")
        print(f"  Enable plotting: {enable_plotting}")
        print(f"\nFound {len(multi_config.directories)} directories to process:")
        for d in multi_config.directories:
            scale_f = d['scale_factor']
            frame_dt_ms = 5.0 if scale_f == 1 else 1.0
            print(f"  - {d['name']} (frame interval: {frame_dt_ms} ms)")
        
        # Process all directories
        analysis = MultiDirectoryLaceyAnalysis(multi_config)
        analysis.process_all_directories()
        
        # Performance summary
        print_header("ANALYSIS COMPLETE")
        elapsed = (datetime.datetime.now() - start_time).total_seconds()
        print(f"  Total Processing Time: {elapsed:.2f} seconds")
        print(f"  Directories Processed: {len(analysis.results)}")
        
        if filter_cohesion:
            print(f"  Results saved to: lacey_results/se{filter_cohesion}/")
        else:
            print(f"  Results saved to: lacey_results/")
        
        # Print summary of mixing quality
        print("\n  Mixing Quality Summary:")
        for name, result in analysis.results.items():
            steady_state = result['steady_state']
            quality = "★" * min(5, max(1, int(steady_state['mean_lacey'] * 5 + 0.5)))
            print(f"    {name.replace('SUPMixerOutput_', '')}: M = {steady_state['mean_lacey']:.3f} {quality}")
        
        print("\n✓ All analyses complete!")
        print("=" * 60)
        
    except Exception as e:
        print(f"\n✗ ERROR: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
