#!/usr/bin/env python3
"""
SUP Mixer Velocity Distribution Analysis - Multi-Directory Comparison Version
Analyzes particle velocity distributions in DEM mixer simulations with SUP models
Supports multiple surface energy values comparison
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
    """Configuration parameter management for multi-directory analysis"""
    def __init__(self, scale_factor: int = 4, surface_energy: str = "000", 
                 base_path: str = ".", enable_plotting: bool = False,
                 steady_state_fraction: float = 0.5):
        self.scale_factor = scale_factor
        self.surface_energy = surface_energy
        self.base_path = Path(base_path)
        
        # Construct input and output paths
        self.input_path = self.base_path / f'SUPMixerOutput_f{scale_factor}se{surface_energy}'
        self.output_base = Path(f'sta_results/se{surface_energy}/')
        self.output_path = self.output_base / f'f{scale_factor}/'
        
        # Simulation parameters
        self.RPM = 300  # Mixer rotation speed
        self.dt = 1e-6  # Time step
        self.frame_skip = 100  # Read every nth frame
        self.max_sample_size = 100000  # Maximum sampling rows
        self.chunk_size = 10000  # Chunk reading size
        self.steady_state_fraction = steady_state_fraction  # Use last fraction of data (0.5 = last half)
        
        # Plotting control parameters
        self.enable_plotting = enable_plotting
        self.plot_dpi = 300  # Publication quality DPI
        self.plot_show = False  # Show plots or just save
        
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
            matplotlib.use('Agg')  # Non-interactive backend
        
        # Create output directories
        self.output_path.mkdir(parents=True, exist_ok=True)
        
    def __str__(self):
        return f"Config(f{self.scale_factor}se{self.surface_energy})"

class MultiDirectoryConfig:
    """Configuration for multi-directory comparison analysis"""
    def __init__(self, base_path: str = ".", filter_cohesion: str = "000",
                 use_interpolation: bool = False, enable_plotting: bool = True,
                 steady_state_fraction: float = 0.5):
        self.base_path = Path(base_path)
        self.filter_cohesion = filter_cohesion
        self.use_interpolation = use_interpolation
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
                # Extract scale_factor and surface_energy from directory name
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
def compute_angles_vectorized(x: np.ndarray, y: np.ndarray, 
                             frame_numbers: np.ndarray, 
                             dt: float, RPM: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Vectorized calculation of all angle-related values"""
    n = len(x)
    polar_angles = np.zeros(n)
    rotor_angles = np.zeros(n)
    relative_angles = np.zeros(n)
    
    RPS = RPM / 60.0
    
    for i in prange(n):
        # Calculate polar angle
        polar_angles[i] = np.arctan2(y[i], x[i])
        
        # Calculate rotor angle
        time = frame_numbers[i] * dt
        rotor_angles[i] = (RPS * time * 2 * np.pi) % (2 * np.pi)
        
        # Calculate relative angle
        rel_angle = polar_angles[i] - rotor_angles[i]
        relative_angles[i] = np.mod(rel_angle + np.pi, 2 * np.pi) - np.pi
    
    return polar_angles, rotor_angles, relative_angles

@jit(nopython=True, parallel=True)
def classify_angles_vectorized(angles: np.ndarray) -> np.ndarray:
    """Vectorized angle classification"""
    n = len(angles)
    intervals = np.zeros(n, dtype=np.int32)
    
    for i in prange(n):
        angle_deg = np.rad2deg(angles[i])
        angle_deg = np.mod(angle_deg + 360, 360)
        intervals[i] = 1 + int(angle_deg // 10)
    
    return intervals

@jit(nopython=True, parallel=True)
def calculate_all_velocities(v_x: np.ndarray, v_y: np.ndarray, v_z: np.ndarray,
                            w_x: np.ndarray, w_y: np.ndarray, w_z: np.ndarray,
                            x: np.ndarray, y: np.ndarray) -> Dict:
    """Calculate all velocity and distance related values at once"""
    n = len(v_x)
    v_total = np.zeros(n)
    v_xy = np.zeros(n)
    v_z_abs = np.zeros(n)
    w_total = np.zeros(n)
    r_dist = np.zeros(n)
    
    for i in prange(n):
        v_total[i] = np.sqrt(v_x[i]**2 + v_y[i]**2 + v_z[i]**2)
        v_xy[i] = np.sqrt(v_x[i]**2 + v_y[i]**2)
        v_z_abs[i] = np.abs(v_z[i])
        w_total[i] = np.sqrt(w_x[i]**2 + w_y[i]**2 + w_z[i]**2)
        r_dist[i] = np.sqrt(x[i]**2 + y[i]**2)
    
    return {
        'v_total': v_total,
        'v_xy': v_xy,
        'v_z_abs': v_z_abs,
        'w_total': w_total,
        'r_dist': r_dist
    }

# ======================== Data Processing Class ========================
class MixerDataProcessor:
    def __init__(self, config: Config):
        self.config = config
        self.df = None
        
    def extract_frame_number(self, filename: str) -> int:
        """Extract frame number from filename"""
        match = re.search(r'mixer_output_(\d+)\.csv', filename)
        return int(match.group(1)) if match else 0
    
    def get_files_to_read(self, max_files: Optional[int] = None) -> list:
        """Get list of files to read (only steady state portion)"""
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
        
        # Calculate steady state portion
        total_frames = len(all_frames)
        start_idx = int(total_frames * (1 - self.config.steady_state_fraction))
        steady_state_frames = all_frames[start_idx:]
        
        print(f"  ✓ Total frames found: {total_frames}")
        print(f"  ✓ Using last {self.config.steady_state_fraction:.0%} ({len(steady_state_frames)} frames) for steady state")
        
        # Filter by frame_skip
        files_to_read = []
        for file, frame_num in steady_state_frames:
            if frame_num % self.config.frame_skip == 0:
                files_to_read.append((file, frame_num))
        
        if max_files:
            files_to_read = files_to_read[:max_files]
            
        return files_to_read
    
    def read_data_chunked(self, file_path: str) -> pd.DataFrame:
        """Read large files in chunks"""
        chunks = []
        for chunk in pd.read_csv(file_path, chunksize=self.config.chunk_size):
            chunks.append(chunk)
        return pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame()
    
    def load_data(self, max_files: Optional[int] = None) -> pd.DataFrame:
        """Load mixer data (steady state only)"""
        files_to_read = self.get_files_to_read(max_files)
        
        print(f"  ✓ Will process {len(files_to_read)} files")
        print(f"  ✓ Reading every {self.config.frame_skip} frames")
        print("-" * 40)
        
        # Check first file structure
        first_file_path, _ = files_to_read[0]
        first_df = pd.read_csv(first_file_path, nrows=1)
        has_velocity = 'velocity' in first_df.columns
        
        if has_velocity:
            print("  ✓ Detected pre-processed files with velocity column")
        else:
            print("  ⚠ No velocity column found - will calculate from v_x, v_y, v_z")
        
        # Read data files
        all_data = []
        for i, (file_path, frame_num) in enumerate(files_to_read, 1):
            try:
                # Use chunked reading for large files
                file_size = os.path.getsize(file_path)
                if file_size > 50 * 1024 * 1024:  # > 50MB
                    df = self.read_data_chunked(file_path)
                else:
                    df = pd.read_csv(file_path)
                
                df['frame'] = frame_num
                df['time'] = frame_num * self.config.dt
                all_data.append(df)
                
                # Progress indicator
                if i % 10 == 0 or i == len(files_to_read):
                    progress = (i / len(files_to_read)) * 100
                    print(f"    Reading files: {progress:5.1f}% complete ({i}/{len(files_to_read)})")
                    
            except Exception as e:
                print(f"  ✗ Error reading {file_path}: {e}")
                continue
        
        if not all_data:
            raise ValueError("No data successfully loaded")
        
        # Merge data
        print("-" * 40)
        print("  ✓ Merging data frames...")
        self.df = pd.concat(all_data, ignore_index=True)
        print(f"  ✓ Total data points: {len(self.df):,}")
        
        return self.df
    
    def process_velocities(self) -> None:
        """Process velocity data"""
        if self.df is None:
            raise ValueError("No data loaded. Call load_data() first.")
        
        # Check for pre-calculated velocity
        if 'velocity' in self.df.columns:
            print("  ✓ Using pre-calculated velocity column")
            results = {
                'v_total': self.df['velocity'].values,
                'v_xy': np.sqrt(self.df['v_x'].values**2 + self.df['v_y'].values**2),
                'v_z_abs': np.abs(self.df['v_z'].values),
                'w_total': np.sqrt(self.df['w_x'].values**2 + 
                                 self.df['w_y'].values**2 + 
                                 self.df['w_z'].values**2),
                'r_dist': np.sqrt(self.df['X'].values**2 + self.df['Y'].values**2)
            }
        else:
            print("  ✓ Calculating velocities using JIT acceleration...")
            results = calculate_all_velocities(
                self.df['v_x'].values, self.df['v_y'].values, self.df['v_z'].values,
                self.df['w_x'].values, self.df['w_y'].values, self.df['w_z'].values,
                self.df['X'].values, self.df['Y'].values
            )
        
        # Batch assignment
        for key, value in results.items():
            self.df[key] = value
        
        print(f"  ✓ Velocity processing complete")
    
    def process_angles(self) -> None:
        """Process angle data"""
        if self.df is None:
            raise ValueError("No data loaded. Call load_data() first.")
        
        print("  ✓ Calculating angles using JIT acceleration...")
        
        # Use vectorized JIT functions
        polar_angles, rotor_angles, relative_angles = compute_angles_vectorized(
            self.df['X'].values,
            self.df['Y'].values,
            self.df['frame'].values,
            self.config.dt,
            self.config.RPM
        )
        
        # Classify angles
        intervals = classify_angles_vectorized(relative_angles)
        
        # Batch assignment
        self.df['polar_angle'] = polar_angles
        self.df['rotor_angle'] = rotor_angles
        self.df['relative_angle'] = relative_angles
        self.df['angle_interval'] = intervals
        
        print(f"  ✓ Angle processing complete")

# ======================== Statistical Analysis Class ========================
class StatisticalAnalyzer:
    def __init__(self, df: pd.DataFrame):
        self.df = df
        
    def compute_angle_statistics(self) -> pd.DataFrame:
        """Compute statistics by angle interval"""
        return self.df.groupby('angle_interval').agg({
            'v_z': ['mean', 'std'],
            'v_z_abs': 'mean',
            'v_total': ['mean', 'std'],
            'Z': 'mean',
            'r': 'mean'
        }).round(4)
    
    def compute_radial_statistics(self, bin_size_mm: float = 1.0) -> pd.DataFrame:
        """Compute statistics by radial distance"""
        self.df['r_group'] = (self.df['r_dist'] * 1000 / bin_size_mm).astype(int) * bin_size_mm / 1000
        return self.df.groupby('r_group').agg({
            'v_total': ['mean', 'std'],
            'v_z': ['mean', 'std'],
            'Z': 'mean'
        }).round(4)
    
    def compute_particle_statistics(self) -> pd.DataFrame:
        """Compute statistics by particle radius"""
        return self.df.groupby('r').agg({
            'v_total': ['mean', 'std'],
            'v_z': ['mean', 'std'],
            'w_total': ['mean', 'std']
        }).round(4)
    
    def compute_velocity_pdf(self, bins: int = 76, v_range: Tuple[float, float] = (0, 1.52)) -> pd.DataFrame:
        """Calculate velocity probability density distribution"""
        bin_edges = np.linspace(v_range[0], v_range[1], bins + 1)
        hist, edges = np.histogram(self.df['v_total'].values, bins=bin_edges, density=True)
        bin_centers = (edges[:-1] + edges[1:]) / 2
        
        return pd.DataFrame({
            'velocity_center': bin_centers,
            'probability_density': hist
        })
    
    def get_summary_statistics(self) -> dict:
        """Get summary statistics"""
        return {
            'total_particles': len(self.df),
            'num_frames': self.df['frame'].nunique(),
            'time_range': (self.df['time'].min(), self.df['time'].max()),
            'velocity_stats': {
                'mean': self.df['v_total'].mean(),
                'std': self.df['v_total'].std(),
                'max': self.df['v_total'].max(),
                'min': self.df['v_total'].min(),
                'median': self.df['v_total'].median()
            },
            'particle_sizes': sorted(self.df['r'].unique()) if 'r' in self.df.columns else [],
            'spatial_extent': {
                'x_range': (self.df['X'].min(), self.df['X'].max()),
                'y_range': (self.df['Y'].min(), self.df['Y'].max()),
                'z_range': (self.df['Z'].min(), self.df['Z'].max())
            }
        }

# ======================== Multi-Directory Comparison Class ========================
class MultiDirectoryComparison:
    """Handle comparison across multiple directories"""
    def __init__(self, multi_config: MultiDirectoryConfig):
        self.multi_config = multi_config
        self.results = {}
        
    def process_all_directories(self):
        """Process all matching directories"""
        print("\n" + "=" * 60)
        print("MULTI-DIRECTORY ANALYSIS")
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
                processor = MixerDataProcessor(config)
                df = processor.load_data()
                processor.process_velocities()
                processor.process_angles()
                
                # Statistical analysis
                analyzer = StatisticalAnalyzer(processor.df)
                
                # Store results
                self.results[dir_info['name']] = {
                    'config': config,
                    'dir_info': dir_info,
                    'df': processor.df,
                    'analyzer': analyzer,
                    'angle_stats': analyzer.compute_angle_statistics(),
                    'radial_stats': analyzer.compute_radial_statistics(),
                    'particle_stats': analyzer.compute_particle_statistics(),
                    'pdf_data': analyzer.compute_velocity_pdf(),
                    'summary': analyzer.get_summary_statistics()
                }
                
                # Save individual results
                self.save_individual_results(dir_info['name'])
                
            except Exception as e:
                print(f"  ✗ Error processing {dir_info['name']}: {e}")
                continue
        
        # Generate comparison results
        if len(self.results) > 1:
            self.generate_comparison()
    
    def save_individual_results(self, dir_name: str):
        """Save results for individual directory"""
        result = self.results[dir_name]
        config = result['config']
        
        # Save statistical results
        result['angle_stats'].to_csv(config.output_path / f'AngleStats_f{config.scale_factor}.csv')
        result['radial_stats'].to_csv(config.output_path / f'RadialStats_f{config.scale_factor}.csv')
        result['particle_stats'].to_csv(config.output_path / f'ParticleStats_f{config.scale_factor}.csv')
        result['pdf_data'].to_csv(config.output_path / f'PDF_f{config.scale_factor}.csv', index=False)
        
        # Save summary
        with open(config.output_path / f'Summary_f{config.scale_factor}.txt', 'w') as f:
            summary = result['summary']
            f.write(f"SUP Mixer Analysis Summary - Scale Factor f{config.scale_factor}\n")
            f.write(f"Surface Energy: se{config.surface_energy}\n")
            f.write("=" * 60 + "\n\n")
            f.write(f"Total particles analyzed: {summary['total_particles']}\n")
            f.write(f"Number of frames: {summary['num_frames']}\n")
            f.write(f"Time range: {summary['time_range'][0]:.6f} - {summary['time_range'][1]:.6f} seconds\n")
            f.write(f"\nVelocity Statistics:\n")
            f.write(f"  Mean:   {summary['velocity_stats']['mean']:.4f} m/s\n")
            f.write(f"  Std:    {summary['velocity_stats']['std']:.4f} m/s\n")
            f.write(f"  Max:    {summary['velocity_stats']['max']:.4f} m/s\n")
            f.write(f"  Min:    {summary['velocity_stats']['min']:.4f} m/s\n")
            f.write(f"  Median: {summary['velocity_stats']['median']:.4f} m/s\n")
        
        print(f"  ✓ Results saved to: {config.output_path}")
    
    def generate_comparison(self):
        """Generate comparison plots and summary"""
        print("\n" + "=" * 60)
        print("GENERATING COMPARISON ANALYSIS")
        print("=" * 60)
        
        # Create comparison directory
        if self.multi_config.filter_cohesion:
            comparison_dir = Path(f'sta_results/se{self.multi_config.filter_cohesion}/comparison/')
        else:
            comparison_dir = Path('sta_results/comparison/')
        comparison_dir.mkdir(parents=True, exist_ok=True)
        
        # Collect comparison data
        comparison_data = []
        for name, result in self.results.items():
            summary = result['summary']
            comparison_data.append({
                'Directory': name,
                'Scale_Factor': result['dir_info']['scale_factor'],
                'Surface_Energy': result['dir_info']['surface_energy'],
                'Mean_Velocity': summary['velocity_stats']['mean'],
                'Std_Velocity': summary['velocity_stats']['std'],
                'Max_Velocity': summary['velocity_stats']['max'],
                'Min_Velocity': summary['velocity_stats']['min'],
                'Median_Velocity': summary['velocity_stats']['median'],
                'Total_Particles': summary['total_particles'],
                'Num_Frames': summary['num_frames']
            })
        
        # Save comparison CSV
        comparison_df = pd.DataFrame(comparison_data)
        comparison_df = comparison_df.sort_values(['Scale_Factor', 'Surface_Energy'])
        comparison_df.to_csv(comparison_dir / 'velocity_comparison.csv', index=False)
        print(f"  ✓ Comparison data saved to: {comparison_dir / 'velocity_comparison.csv'}")
        
        # Generate comparison plots if enabled
        if self.multi_config.enable_plotting:
            self.plot_comparison(comparison_dir)
    
    def plot_comparison(self, output_dir: Path):
        """Generate comparison plots"""
        # Set up scientific style
        plt.rcParams.update({
            'font.family': 'serif',
            'font.size': 12,
            'axes.linewidth': 1.8,
            'xtick.major.width': 1.2,
            'ytick.major.width': 1.2
        })
        
        # Plot 1: PDF Comparison
        fig, ax = plt.subplots(figsize=(14, 8))
        
        colors = plt.cm.viridis(np.linspace(0, 0.9, len(self.results)))
        
        for (name, result), color in zip(self.results.items(), colors):
            pdf_data = result['pdf_data']
            dir_info = result['dir_info']
            label = f"f{dir_info['scale_factor']} (se{dir_info['surface_energy']})"
            
            ax.plot(pdf_data['velocity_center'], pdf_data['probability_density'],
                   '-', linewidth=2.4, color=color, label=label, alpha=0.8)
            
            # Fill area with transparency
            ax.fill_between(pdf_data['velocity_center'], pdf_data['probability_density'],
                          alpha=0.1, color=color)
        
        ax.set_xlabel('Velocity [m/s]', fontsize=14, fontweight='bold')
        ax.set_ylabel('Probability Density', fontsize=14, fontweight='bold')
        ax.set_title('Velocity PDF Comparison - Multiple Configurations', 
                    fontsize=16, fontweight='bold')
        ax.legend(loc='upper right', frameon=True, fancybox=False,
                 edgecolor='black', fontsize=10)
        ax.grid(True, alpha=0.3, linewidth=0.8)
        
        plt.tight_layout()
        plt.savefig(output_dir / 'PDF_comparison.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # Plot 2: Mean Velocity Bar Chart
        fig, ax = plt.subplots(figsize=(12, 8))
        
        comparison_data = []
        for name, result in self.results.items():
            summary = result['summary']
            dir_info = result['dir_info']
            comparison_data.append({
                'label': f"f{dir_info['scale_factor']}\nse{dir_info['surface_energy']}",
                'mean': summary['velocity_stats']['mean'],
                'std': summary['velocity_stats']['std']
            })
        
        x_pos = np.arange(len(comparison_data))
        means = [d['mean'] for d in comparison_data]
        stds = [d['std'] for d in comparison_data]
        labels = [d['label'] for d in comparison_data]
        
        bars = ax.bar(x_pos, means, yerr=stds, capsize=5, 
                     color=colors, edgecolor='black', linewidth=1.5, alpha=0.8)
        
        ax.set_xlabel('Configuration', fontsize=14, fontweight='bold')
        ax.set_ylabel('Mean Velocity [m/s]', fontsize=14, fontweight='bold')
        ax.set_title('Mean Velocity Comparison', fontsize=16, fontweight='bold')
        ax.set_xticks(x_pos)
        ax.set_xticklabels(labels)
        ax.grid(True, axis='y', alpha=0.3, linewidth=0.8)
        
        # Add value labels on bars
        for bar, mean, std in zip(bars, means, stds):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + std,
                   f'{mean:.3f}', ha='center', va='bottom', fontsize=10)
        
        plt.tight_layout()
        plt.savefig(output_dir / 'mean_velocity_comparison.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"  ✓ Comparison plots saved to: {output_dir}")

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
    Main program for multi-directory analysis
    """
    # ==================== Configurable Parameters ====================
    base_path = "build"  # Search directory, can be "build" or other path
    filter_cohesion = "001"  # Only process cohesion=0 data, set to None for all
    use_interpolation = False  # Whether to use interpolation
    enable_plotting = True  # Generate plots
    steady_state_fraction = 0.5  # Use last 50% of data (can be 0.33 for last 1/3)
    
    # ==================== Run Analysis ====================
    print_header("SUP MIXER MULTI-DIRECTORY VELOCITY ANALYSIS")
    start_time = datetime.datetime.now()
    
    try:
        # Initialize multi-directory configuration
        multi_config = MultiDirectoryConfig(
            base_path=base_path,
            filter_cohesion=filter_cohesion,
            use_interpolation=use_interpolation,
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
        print(f"  Steady state fraction: {steady_state_fraction:.0%} (using last {steady_state_fraction:.0%} of data)")
        print(f"  Enable plotting: {enable_plotting}")
        print(f"  Use interpolation: {use_interpolation}")
        print(f"\nFound {len(multi_config.directories)} directories to process:")
        for d in multi_config.directories:
            print(f"  - {d['name']}")
        
        # Process all directories
        comparison = MultiDirectoryComparison(multi_config)
        comparison.process_all_directories()
        
        # Performance summary
        print_header("ANALYSIS COMPLETE")
        elapsed = (datetime.datetime.now() - start_time).total_seconds()
        print(f"  Total Processing Time: {elapsed:.2f} seconds")
        print(f"  Directories Processed: {len(comparison.results)}")
        
        if filter_cohesion:
            print(f"  Results saved to: sta_results/se{filter_cohesion}/")
        else:
            print(f"  Results saved to: sta_results/")
        
        print("\n✓ All analyses complete!")
        print("=" * 60)
        
    except Exception as e:
        print(f"\n✗ ERROR: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
