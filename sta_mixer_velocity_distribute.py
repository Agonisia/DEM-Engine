#!/usr/bin/env python3
"""
SUP Mixer Velocity Distribution Analysis - Scientific Publication Style
Analyzes particle velocity distributions in DEM mixer simulations with SUP models
"""

import os
import re
import datetime
import pandas as pd
import numpy as np
from numba import jit, prange
import glob
from pathlib import Path
from typing import Tuple, Dict, Optional
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ======================== Configuration Parameters ========================
class Config:
    """Configuration parameter management"""
    def __init__(self, scale_factor: int = 4, enable_plotting: bool = False):
        self.scale_factor = scale_factor
        self.input_path = Path(f'build/SUPMixerOutput_f{scale_factor}/')
        self.output_path = Path(f'analysis_results/f{scale_factor}/')
        self.RPM = 300  # Mixer rotation speed
        self.dt = 1e-6  # Time step
        self.frame_skip = 100  # Read every nth frame
        self.max_sample_size = 100000  # Maximum sampling rows
        self.chunk_size = 10000  # Chunk reading size
        
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
        
        # Create output directory
        self.output_path.mkdir(parents=True, exist_ok=True)

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
        """Get list of files to read"""
        pattern = self.config.input_path / 'mixer_output_*.csv'
        files = sorted(glob.glob(str(pattern)))
        
        if not files:
            raise FileNotFoundError(f"No CSV files found in {self.config.input_path}")
        
        # Filter files
        files_to_read = []
        for file in files:
            frame_num = self.extract_frame_number(os.path.basename(file))
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
        """Load mixer data"""
        files_to_read = self.get_files_to_read(max_files)
        
        print(f"  ✓ Found {len(files_to_read)} files to process")
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
            'particle_sizes': sorted(self.df['r'].unique()),
            'spatial_extent': {
                'x_range': (self.df['X'].min(), self.df['X'].max()),
                'y_range': (self.df['Y'].min(), self.df['Y'].max()),
                'z_range': (self.df['Z'].min(), self.df['Z'].max())
            }
        }

# ======================== Scientific Style Plotting Class ========================
class PlotGenerator:
    """Generate scientific publication style analysis plots"""
    def __init__(self, config: Config):
        self.config = config
    
    def apply_scientific_style(self, ax):
        """Apply scientific paper style to axes"""
        # Axis frame settings
        for spine in ax.spines.values():
            spine.set_linewidth(1.8)
            spine.set_color('black')
        
        # Major tick settings
        ax.tick_params(axis='both', which='major', 
                      labelsize=14, width=1.2, length=5,
                      direction='in', top=True, right=True)
        
        # Minor tick settings
        ax.tick_params(axis='both', which='minor', 
                      width=0.6, length=2.4,
                      direction='in', top=True, right=True)
        
        # Grid settings
        ax.grid(True, which='major', alpha=0.15, linewidth=0.8, 
               linestyle='-', color='gray', zorder=0)
        ax.grid(True, which='minor', alpha=0.08, linewidth=0.5, 
               linestyle='--', color='gray', zorder=0)
        ax.minorticks_on()
        
        # Ensure data is above grid
        ax.set_axisbelow(True)
    
    def plot_velocity_pdf(self, pdf_data: pd.DataFrame) -> None:
        """Plot velocity probability density distribution - Scientific style"""
        fig, ax = plt.subplots(figsize=(14, 8))
        
        # Plot PDF curve
        line = ax.plot(pdf_data['velocity_center'], 
                      pdf_data['probability_density'], 
                      '-', linewidth=2.4, color='#0066CC',
                      label=f'SUP Scale Factor f{self.config.scale_factor}',
                      zorder=5)[0]
        
        # Fill area under curve
        ax.fill_between(pdf_data['velocity_center'], 
                       pdf_data['probability_density'],
                       alpha=0.15, color='#0066CC', zorder=3)
        
        # Mark peak
        max_idx = pdf_data['probability_density'].idxmax()
        max_v = pdf_data.loc[max_idx, 'velocity_center']
        max_p = pdf_data.loc[max_idx, 'probability_density']
        
        ax.scatter([max_v], [max_p], s=120, color='red', 
                  marker='o', zorder=6, edgecolors='darkred', linewidth=1.5)
        
        # Peak annotation
        ax.annotate(f'Peak: v={max_v:.3f} m/s\nP={max_p:.3f}',
                   xy=(max_v, max_p), xytext=(max_v + 0.15, max_p),
                   fontsize=11, fontweight='bold',
                   bbox=dict(boxstyle='round,pad=0.3', 
                            facecolor='white', edgecolor='black', linewidth=1),
                   arrowprops=dict(arrowstyle='->', color='black', lw=1.2),
                   zorder=7)
        
        # Statistics box
        mean_v = np.sum(pdf_data['velocity_center'] * pdf_data['probability_density'] * 
                       np.diff(np.append(pdf_data['velocity_center'], pdf_data['velocity_center'].iloc[-1])))
        
        stats_text = (f'Velocity PDF Statistics:\n'
                     f'{"─" * 30}\n'
                     f'Peak Velocity:  {max_v:.3f} m/s\n'
                     f'Peak Density:   {max_p:.3f}\n'
                     f'Mean Velocity:  {mean_v:.3f} m/s\n'
                     f'Scale Factor:   f{self.config.scale_factor}')
        
        ax.text(0.98, 0.02, stats_text, transform=ax.transAxes,
               verticalalignment='bottom', horizontalalignment='right',
               fontsize=10, fontweight='normal', fontfamily='monospace',
               bbox=dict(boxstyle='square,pad=0.5', facecolor='white',
                        edgecolor='black', linewidth=1.2),
               zorder=10)
        
        # Axis labels
        ax.set_xlabel('Velocity [m/s]', fontsize=17, fontweight='bold')
        ax.set_ylabel('Probability Density', fontsize=17, fontweight='bold')
        ax.set_title(f'(a) Velocity Probability Density Function - SUP f{self.config.scale_factor}', 
                    fontsize=19, fontweight='bold', pad=15)
        
        # Set axis ranges
        ax.set_xlim(0, pdf_data['velocity_center'].max() * 1.05)
        ax.set_ylim(0, pdf_data['probability_density'].max() * 1.1)
        
        # Apply scientific style
        self.apply_scientific_style(ax)
        
        # Legend
        legend = ax.legend(loc='upper right', frameon=True, fancybox=False,
                          shadow=False, framealpha=1.0, edgecolor='black',
                          facecolor='white', fontsize=12)
        legend.get_frame().set_linewidth(1.2)
        legend.set_zorder(10)
        
        plt.tight_layout()
        plt.savefig(self.config.output_path / f'PDF_f{self.config.scale_factor}.png', 
                   dpi=self.config.plot_dpi, bbox_inches='tight')
        if self.config.plot_show:
            plt.show()
        plt.close()
    
    def plot_radial_velocity(self, radial_stats: pd.DataFrame) -> None:
        """Plot radial velocity distribution - Scientific style"""
        fig, ax = plt.subplots(figsize=(14, 8))
        
        radial_means = radial_stats[('v_total', 'mean')]
        radial_stds = radial_stats[('v_total', 'std')]
        radial_distances = radial_means.index * 1000  # Convert to mm
        
        # Main curve with markers
        line = ax.plot(radial_distances, radial_means.values, 
                      '-', linewidth=2.4, color='#CC0000',
                      label='Mean Velocity', marker='o', markersize=6,
                      markerfacecolor='white', markeredgecolor='#CC0000',
                      markeredgewidth=1.5, zorder=5)[0]
        
        # Error band
        ax.fill_between(radial_distances, 
                       radial_means.values - radial_stds.values,
                       radial_means.values + radial_stds.values,
                       alpha=0.2, color='#CC0000', 
                       label='±1σ', zorder=3)
        
        # Mark key points
        max_idx = radial_means.idxmax()
        min_idx = radial_means.idxmin()
        
        # Annotate maximum
        ax.annotate(f'Max: {radial_means[max_idx]:.3f} m/s\n@ {max_idx*1000:.1f} mm',
                   xy=(max_idx*1000, radial_means[max_idx]),
                   xytext=(max_idx*1000, radial_means[max_idx] + 0.05),
                   fontsize=10, fontweight='bold',
                   bbox=dict(boxstyle='round,pad=0.3', 
                            facecolor='yellow', alpha=0.8, edgecolor='black'),
                   ha='center', va='bottom', zorder=6)
        
        # Reference line (mixer radius)
        mixer_radius = 25  # mm
        ax.axvline(x=mixer_radius, color='green', linestyle='--', 
                  linewidth=1.5, alpha=0.7, label=f'Mixer Radius ({mixer_radius} mm)')
        
        # Statistics box
        stats_text = (f'Radial Velocity Statistics:\n'
                     f'{"─" * 30}\n'
                     f'Maximum:  {radial_means.max():.3f} m/s @ {radial_means.idxmax()*1000:.1f} mm\n'
                     f'Minimum:  {radial_means.min():.3f} m/s @ {radial_means.idxmin()*1000:.1f} mm\n'
                     f'Range:    {radial_means.max() - radial_means.min():.3f} m/s\n'
                     f'Samples:  {len(radial_means)} bins')
        
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
               verticalalignment='top', horizontalalignment='left',
               fontsize=10, fontweight='normal', fontfamily='monospace',
               bbox=dict(boxstyle='square,pad=0.5', facecolor='white',
                        edgecolor='black', linewidth=1.2),
               zorder=10)
        
        # Axis labels
        ax.set_xlabel('Radial Distance [mm]', fontsize=17, fontweight='bold')
        ax.set_ylabel('Mean Velocity [m/s]', fontsize=17, fontweight='bold')
        ax.set_title(f'(b) Radial Velocity Distribution - SUP f{self.config.scale_factor}', 
                    fontsize=19, fontweight='bold', pad=15)
        
        # Apply scientific style
        self.apply_scientific_style(ax)
        
        # Set axis ranges
        ax.set_xlim(radial_distances.min() * 0.95, radial_distances.max() * 1.05)
        ax.set_ylim(0, radial_means.max() * 1.15)
        
        # Legend
        legend = ax.legend(loc='upper right', frameon=True, fancybox=False,
                          shadow=False, framealpha=1.0, edgecolor='black',
                          facecolor='white', fontsize=12)
        legend.get_frame().set_linewidth(1.2)
        legend.set_zorder(10)
        
        plt.tight_layout()
        plt.savefig(self.config.output_path / f'RadialVelocity_f{self.config.scale_factor}.png',
                   dpi=self.config.plot_dpi, bbox_inches='tight')
        if self.config.plot_show:
            plt.show()
        plt.close()
    
    def plot_combined_analysis(self, pdf_data, radial_stats, angle_stats):
        """Create combined analysis plot - 2x2 scientific layout"""
        fig = plt.figure(figsize=(18, 14))
        
        # Create 2x2 subplot layout
        gs = fig.add_gridspec(2, 2, hspace=0.25, wspace=0.25)
        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[0, 1])
        ax3 = fig.add_subplot(gs[1, 0])
        ax4 = fig.add_subplot(gs[1, 1])
        
        # ====== Panel 1: PDF ======
        ax1.plot(pdf_data['velocity_center'], pdf_data['probability_density'],
                '-', linewidth=2.4, color='#0066CC')
        ax1.fill_between(pdf_data['velocity_center'], pdf_data['probability_density'],
                        alpha=0.15, color='#0066CC')
        ax1.set_xlabel('Velocity [m/s]', fontsize=14, fontweight='bold')
        ax1.set_ylabel('Probability Density', fontsize=14, fontweight='bold')
        ax1.set_title('(1) Velocity PDF', fontsize=15, fontweight='bold')
        self.apply_scientific_style(ax1)
        
        # ====== Panel 2: Radial ======
        radial_means = radial_stats[('v_total', 'mean')]
        ax2.plot(radial_means.index * 1000, radial_means.values,
                '-o', linewidth=2.4, color='#CC0000', markersize=4)
        ax2.set_xlabel('Radial Distance [mm]', fontsize=14, fontweight='bold')
        ax2.set_ylabel('Mean Velocity [m/s]', fontsize=14, fontweight='bold')
        ax2.set_title('(2) Radial Distribution', fontsize=15, fontweight='bold')
        self.apply_scientific_style(ax2)
        
        # ====== Panel 3: Angular ======
        angles = (angle_stats.index - 1) * 10 + 5
        ax3.plot(angles, angle_stats[('v_total', 'mean')],
                '-s', linewidth=2.4, color='#0066CC', markersize=4)
        ax3.set_xlabel('Angle [degrees]', fontsize=14, fontweight='bold')
        ax3.set_ylabel('Mean Velocity [m/s]', fontsize=14, fontweight='bold')
        ax3.set_title('(3) Angular Distribution', fontsize=15, fontweight='bold')
        ax3.set_xlim(0, 360)
        self.apply_scientific_style(ax3)
        
        # ====== Panel 4: Vertical Velocity ======
        ax4.plot(angles, angle_stats[('v_z', 'mean')],
                '-o', linewidth=2.4, color='#00AA44', markersize=4)
        ax4.axhline(y=0, color='black', linestyle='--', linewidth=1.5, alpha=0.7)
        ax4.set_xlabel('Angle [degrees]', fontsize=14, fontweight='bold')
        ax4.set_ylabel('Vertical Velocity [m/s]', fontsize=14, fontweight='bold')
        ax4.set_title('(4) Vertical Flow Pattern', fontsize=15, fontweight='bold')
        ax4.set_xlim(0, 360)
        self.apply_scientific_style(ax4)
        
        # Overall title
        fig.suptitle(f'DEM Mixer Comprehensive Analysis - SUP Scale Factor f{self.config.scale_factor}',
                    fontsize=20, fontweight='bold')
        
        plt.tight_layout(rect=[0, 0.02, 1, 0.96])
        plt.savefig(self.config.output_path / f'Combined_Analysis_f{self.config.scale_factor}.png',
                   dpi=300, bbox_inches='tight')
        if self.config.plot_show:
            plt.show()
        plt.close()

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

def print_statistics_summary(summary: dict, config: Config):
    """Print detailed statistics summary"""
    print_header("DATA STATISTICS SUMMARY")
    
    # Basic information
    print_section("Basic Information")
    print(f"  Scale Factor:        f{config.scale_factor}")
    print(f"  Total Particles:     {summary['total_particles']:,}")
    print(f"  Number of Frames:    {summary['num_frames']}")
    print(f"  Time Range:          {summary['time_range'][0]:.6f} - {summary['time_range'][1]:.6f} seconds")
    print(f"  Particle Sizes:      {summary['particle_sizes']}")
    
    # Velocity statistics
    print_section("Velocity Statistics")
    v_stats = summary['velocity_stats']
    print(f"  Minimum:             {v_stats['min']:.6f} m/s")
    print(f"  Maximum:             {v_stats['max']:.6f} m/s")
    print(f"  Mean:                {v_stats['mean']:.6f} m/s")
    print(f"  Standard Deviation:  {v_stats['std']:.6f} m/s")
    print(f"  Median:              {v_stats['median']:.6f} m/s")
    
    # Spatial extent
    print_section("Spatial Extent")
    s_ext = summary['spatial_extent']
    print(f"  X Range: [{s_ext['x_range'][0]:.4f}, {s_ext['x_range'][1]:.4f}] m")
    print(f"  Y Range: [{s_ext['y_range'][0]:.4f}, {s_ext['y_range'][1]:.4f}] m")
    print(f"  Z Range: [{s_ext['z_range'][0]:.4f}, {s_ext['z_range'][1]:.4f}] m")

# ======================== Main Program ========================
def main(enable_plotting: bool = True, plot_show: bool = False, 
         scale_factor: int = 4, interactive: bool = True):
    """
    Main program
    
    Parameters:
        enable_plotting: Whether to generate plots
        plot_show: Whether to show plot windows
        scale_factor: SUP scaling factor
        interactive: Whether to run in interactive mode
    """
    # Initialize configuration
    config = Config(scale_factor=scale_factor, enable_plotting=enable_plotting)
    config.plot_show = plot_show
    
    print_header(f"SUP MIXER VELOCITY DISTRIBUTION ANALYSIS - SCALE FACTOR f{scale_factor}")
    
    start_time = datetime.datetime.now()
    
    try:
        # 1. Data loading and processing
        print_section("1. Loading and Processing Data")
        processor = MixerDataProcessor(config)
        df = processor.load_data()
        
        # 2. Calculate velocities
        print_section("2. Processing Velocities")
        velocity_start = datetime.datetime.now()
        processor.process_velocities()
        velocity_time = (datetime.datetime.now() - velocity_start).total_seconds()
        print(f"  ⚡ Velocity processing time: {velocity_time:.3f} seconds")
        
        # 3. Calculate angles
        print_section("3. Processing Angles")
        processor.process_angles()
        
        # 4. Statistical analysis
        print_section("4. Computing Statistics")
        analyzer = StatisticalAnalyzer(processor.df)
        
        print("  ✓ Computing angle statistics...")
        angle_stats = analyzer.compute_angle_statistics()
        
        print("  ✓ Computing radial statistics...")
        radial_stats = analyzer.compute_radial_statistics()
        
        print("  ✓ Computing particle statistics...")
        particle_stats = analyzer.compute_particle_statistics()
        
        print("  ✓ Computing velocity PDF...")
        pdf_data = analyzer.compute_velocity_pdf()
        
        print("  ✓ Generating summary statistics...")
        summary = analyzer.get_summary_statistics()
        
        # 5. Save results
        print_section("5. Saving Results")
        
        # Save processed data (sampled)
        sample_size = min(config.max_sample_size, len(processor.df))
        df_sample = processor.df.sample(n=sample_size) if len(processor.df) > sample_size else processor.df
        
        output_file = config.output_path / f'RAW_f{config.scale_factor}_processed.csv.gz'
        df_sample.to_csv(output_file, index=False, compression='gzip')
        print(f"  ✓ Processed data saved: {output_file}")
        
        # Save statistical results
        angle_stats.to_csv(config.output_path / f'AngleStats_f{config.scale_factor}.csv')
        print(f"  ✓ Angle statistics saved")
        
        radial_stats.to_csv(config.output_path / f'RadialStats_f{config.scale_factor}.csv')
        print(f"  ✓ Radial statistics saved")
        
        particle_stats.to_csv(config.output_path / f'ParticleStats_f{config.scale_factor}.csv')
        print(f"  ✓ Particle statistics saved")
        
        pdf_data.to_csv(config.output_path / f'PDF_f{config.scale_factor}.csv', index=False)
        print(f"  ✓ PDF data saved")
        
        # Save summary
        with open(config.output_path / f'Summary_f{config.scale_factor}.txt', 'w') as f:
            f.write(f"SUP Mixer Analysis Summary - Scale Factor f{config.scale_factor}\n")
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
            f.write(f"\nParticle sizes: {summary['particle_sizes']}\n")
        print(f"  ✓ Summary report saved")
        
        # 6. Generate plots (optional)
        if config.enable_plotting:
            print_section("6. Generating Plots")
            
            if interactive:
                print("\nPlot Options:")
                print("1. Individual plots (PDF, Radial, Angular)")
                print("2. Combined analysis plot (2x2 layout)")
                print("3. All plots")
                print("=" * 40)
                choice = input("Select option (1/2/3) [default: 3]: ").strip() or "3"
            else:
                choice = "3"  # Generate all plots in non-interactive mode
            
            plotter = PlotGenerator(config)
            
            if choice in ["1", "3"]:
                print("  ✓ Generating velocity PDF plot...")
                plotter.plot_velocity_pdf(pdf_data)
                
                print("  ✓ Generating radial velocity plot...")
                plotter.plot_radial_velocity(radial_stats)
            
            if choice in ["2", "3"]:
                print("  ✓ Generating combined analysis plot...")
                plotter.plot_combined_analysis(pdf_data, radial_stats, angle_stats)
            
            print(f"  ✓ All plots saved to: {config.output_path}")
        else:
            print_section("6. Plotting")
            print("  ⚠ Plotting disabled (set ENABLE_PLOTTING=True to generate plots)")
        
        # 7. Print analysis summary
        print_statistics_summary(summary, config)
        
        # Performance summary
        print_header("PERFORMANCE SUMMARY")
        elapsed = (datetime.datetime.now() - start_time).total_seconds()
        print(f"  Total Processing Time: {elapsed:.2f} seconds")
        
        if 'velocity' in processor.df.columns:
            print(f"  ✓ Pre-calculated velocity used (saved time)")
        
        memory_usage_mb = processor.df.memory_usage(deep=True).sum() / 1024 / 1024
        print(f"  Peak Memory Usage:     {memory_usage_mb:.2f} MB")
        print(f"  Output Directory:      {config.output_path}")
        
        print("\n" + "=" * 60)
        print("✓ ANALYSIS COMPLETE!")
        print("=" * 60)
        
        # Ask if user wants to save plots
        if interactive and config.enable_plotting and config.plot_show:
            save_option = input("\nSave additional high-resolution copies? (y/n) [default: n]: ").strip().lower()
            if save_option == 'y':
                print("Saving high-resolution plots (600 DPI)...")
                config.plot_dpi = 600
                plotter = PlotGenerator(config)
                plotter.plot_velocity_pdf(pdf_data)
                plotter.plot_radial_velocity(radial_stats)
                plotter.plot_combined_analysis(pdf_data, radial_stats, angle_stats)
                print("✓ High-resolution plots saved")
        
    except FileNotFoundError as e:
        print(f"\n✗ ERROR: {e}")
        print("Please check that the input directory exists and contains CSV files.")
    except Exception as e:
        print(f"\n✗ ERROR: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    # ==================== Configuration Options ====================
    # Set running parameters here
    ENABLE_PLOTTING = True   # Generate plots (True/False)
    SHOW_PLOTS = True      # Show plot windows (True/False), False only saves files
    SCALE_FACTOR = 4        # SUP scaling factor
    INTERACTIVE = True      # Run in interactive mode (True/False)
    
    # Run main program
    main(enable_plotting=ENABLE_PLOTTING, 
         plot_show=SHOW_PLOTS,
         scale_factor=SCALE_FACTOR,
         interactive=INTERACTIVE)