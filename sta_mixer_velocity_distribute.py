import os
import re
import datetime
import pandas as pd
import numpy as np
from numba import jit, prange
import glob
from pathlib import Path
from typing import Tuple, Dict, Optional

# ======================== 配置参数 ========================
class Config:
    """配置参数集中管理"""
    def __init__(self, scale_factor: int = 4):
        self.scale_factor = scale_factor
        self.input_path = Path(f'build/SUPMixerOutput_f{scale_factor}/')
        self.output_path = Path(f'sta_results/f{scale_factor}/')
        self.RPM = 300  # 搅拌器转速
        self.dt = 1e-6  # 时间步长
        self.frame_skip = 100  # 每隔多少帧读取一次数据
        self.max_sample_size = 100000  # 最大采样行数
        self.chunk_size = 10000  # 分块读取大小
        
        # 创建输出目录
        self.output_path.mkdir(parents=True, exist_ok=True)

# ======================== JIT加速函数 ========================
@jit(nopython=True, parallel=True)
def compute_angles_vectorized(x: np.ndarray, y: np.ndarray, 
                             frame_numbers: np.ndarray, 
                             dt: float, RPM: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """向量化计算所有角度相关值"""
    n = len(x)
    polar_angles = np.zeros(n)
    rotor_angles = np.zeros(n)
    relative_angles = np.zeros(n)
    
    RPS = RPM / 60.0
    
    for i in prange(n):
        # 计算极角
        polar_angles[i] = np.arctan2(y[i], x[i])
        
        # 计算搅拌器角度
        time = frame_numbers[i] * dt
        rotor_angles[i] = (RPS * time * 2 * np.pi) % (2 * np.pi)
        
        # 计算相对角度
        rel_angle = polar_angles[i] - rotor_angles[i]
        relative_angles[i] = np.mod(rel_angle + np.pi, 2 * np.pi) - np.pi
    
    return polar_angles, rotor_angles, relative_angles

@jit(nopython=True, parallel=True)
def classify_angles_vectorized(angles: np.ndarray) -> np.ndarray:
    """向量化角度分类"""
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
    """一次性计算所有速度和距离相关值"""
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

# ======================== 数据处理类 ========================
class MixerDataProcessor:
    def __init__(self, config: Config):
        self.config = config
        self.df = None
        
    def extract_frame_number(self, filename: str) -> int:
        """从文件名中提取帧数"""
        match = re.search(r'mixer_output_(\d+)\.csv', filename)
        return int(match.group(1)) if match else 0
    
    def get_files_to_read(self, max_files: Optional[int] = None) -> list:
        """获取要读取的文件列表"""
        pattern = self.config.input_path / 'mixer_output_*.csv'
        files = sorted(glob.glob(str(pattern)))
        
        if not files:
            raise FileNotFoundError(f"No CSV files found in {self.config.input_path}")
        
        # 筛选文件
        files_to_read = []
        for file in files:
            frame_num = self.extract_frame_number(os.path.basename(file))
            if frame_num % self.config.frame_skip == 0:
                files_to_read.append((file, frame_num))
        
        if max_files:
            files_to_read = files_to_read[:max_files]
            
        return files_to_read
    
    def read_data_chunked(self, file_path: str) -> pd.DataFrame:
        """分块读取大文件"""
        chunks = []
        for chunk in pd.read_csv(file_path, chunksize=self.config.chunk_size):
            chunks.append(chunk)
        return pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame()
    
    def load_data(self, max_files: Optional[int] = None) -> pd.DataFrame:
        """加载混合器数据"""
        files_to_read = self.get_files_to_read(max_files)
        
        print(f"Reading {len(files_to_read)} files (every {self.config.frame_skip} frames)")
        
        # 检查第一个文件的列结构
        first_file_path, _ = files_to_read[0]
        first_df = pd.read_csv(first_file_path, nrows=1)
        has_velocity = 'velocity' in first_df.columns
        
        if has_velocity:
            print("✅ Detected pre-processed files with velocity column")
        else:
            print("⚠️  No velocity column found - will calculate from v_x, v_y, v_z")
        
        # 并行读取数据（使用列表推导式更高效）
        all_data = []
        for file_path, frame_num in files_to_read:
            try:
                # 对大文件使用分块读取
                file_size = os.path.getsize(file_path)
                if file_size > 50 * 1024 * 1024:  # 大于50MB
                    df = self.read_data_chunked(file_path)
                else:
                    df = pd.read_csv(file_path)
                
                df['frame'] = frame_num
                df['time'] = frame_num * self.config.dt
                all_data.append(df)
                print(f"  Read frame {frame_num}: {len(df)} particles")
            except Exception as e:
                print(f"  Error reading {file_path}: {e}")
                continue
        
        if not all_data:
            raise ValueError("No data successfully loaded")
        
        # 合并数据
        self.df = pd.concat(all_data, ignore_index=True)
        print(f"\nTotal data points: {len(self.df)}")
        
        return self.df
    
    def process_velocities(self) -> None:
        """处理速度数据"""
        if self.df is None:
            raise ValueError("No data loaded. Call load_data() first.")
        
        # 如果已有velocity列，只需计算其他值
        if 'velocity' in self.df.columns:
            print("   Using pre-calculated velocity column")
            # 仍需计算其他速度分量
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
            print("   Calculating velocities...")
            # 使用JIT加速计算
            results = calculate_all_velocities(
                self.df['v_x'].values, self.df['v_y'].values, self.df['v_z'].values,
                self.df['w_x'].values, self.df['w_y'].values, self.df['w_z'].values,
                self.df['X'].values, self.df['Y'].values
            )
        
        # 批量赋值（更高效）
        for key, value in results.items():
            self.df[key] = value
    
    def process_angles(self) -> None:
        """处理角度数据"""
        if self.df is None:
            raise ValueError("No data loaded. Call load_data() first.")
        
        print("   Calculating angles...")
        
        # 使用向量化JIT函数
        polar_angles, rotor_angles, relative_angles = compute_angles_vectorized(
            self.df['X'].values,
            self.df['Y'].values,
            self.df['frame'].values,
            self.config.dt,
            self.config.RPM
        )
        
        # 分类角度
        intervals = classify_angles_vectorized(relative_angles)
        
        # 批量赋值
        self.df['polar_angle'] = polar_angles
        self.df['rotor_angle'] = rotor_angles
        self.df['relative_angle'] = relative_angles
        self.df['angle_interval'] = intervals

# ======================== 统计分析类 ========================
class StatisticalAnalyzer:
    def __init__(self, df: pd.DataFrame):
        self.df = df
        
    def compute_angle_statistics(self) -> pd.DataFrame:
        """按角度区间统计"""
        return self.df.groupby('angle_interval').agg({
            'v_z': ['mean', 'std'],
            'v_z_abs': 'mean',
            'v_total': ['mean', 'std'],
            'Z': 'mean',
            'r': 'mean'
        }).round(4)
    
    def compute_radial_statistics(self, bin_size_mm: float = 1.0) -> pd.DataFrame:
        """按径向距离统计"""
        self.df['r_group'] = (self.df['r_dist'] * 1000 / bin_size_mm).astype(int) * bin_size_mm / 1000
        return self.df.groupby('r_group').agg({
            'v_total': ['mean', 'std'],
            'v_z': ['mean', 'std'],
            'Z': 'mean'
        }).round(4)
    
    def compute_particle_statistics(self) -> pd.DataFrame:
        """按颗粒半径统计"""
        return self.df.groupby('r').agg({
            'v_total': ['mean', 'std'],
            'v_z': ['mean', 'std'],
            'w_total': ['mean', 'std']
        }).round(4)
    
    def compute_velocity_pdf(self, bins: int = 75, v_range: Tuple[float, float] = (0, 1.5)) -> pd.DataFrame:
        """计算速度概率密度分布"""
        bin_edges = np.linspace(v_range[0], v_range[1], bins + 1)
        hist, edges = np.histogram(self.df['v_total'].values, bins=bin_edges, density=True)
        bin_centers = (edges[:-1] + edges[1:]) / 2
        
        return pd.DataFrame({
            'velocity_center': bin_centers,
            'probability_density': hist
        })
    
    def get_summary_statistics(self) -> dict:
        """获取汇总统计信息"""
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

# ======================== 主程序 ========================
def main():
    # 初始化配置
    config = Config(scale_factor=4)  # 可以通过参数调整缩放因子
    
    print(f"{'='*60}")
    print(f"DEME Mixer Analysis - Scale Factor: f{config.scale_factor}")
    print(f"{'='*60}")
    
    start_time = datetime.datetime.now()
    
    try:
        # 1. 数据加载和处理
        print("\n1. Loading and processing data...")
        processor = MixerDataProcessor(config)
        df = processor.load_data()
        
        # 2. 计算速度
        print("\n2. Processing velocities...")
        velocity_start = datetime.datetime.now()
        processor.process_velocities()
        velocity_time = (datetime.datetime.now() - velocity_start).total_seconds()
        
        # 3. 计算角度
        print("\n3. Processing angles...")
        processor.process_angles()
        
        # 4. 统计分析
        print("\n4. Computing statistics...")
        analyzer = StatisticalAnalyzer(processor.df)
        
        angle_stats = analyzer.compute_angle_statistics()
        radial_stats = analyzer.compute_radial_statistics()
        particle_stats = analyzer.compute_particle_statistics()
        pdf_data = analyzer.compute_velocity_pdf()
        summary = analyzer.get_summary_statistics()
        
        # 5. 保存结果
        print("\n5. Saving results...")
        
        # 保存处理数据（采样）
        sample_size = min(config.max_sample_size, len(processor.df))
        df_sample = processor.df.sample(n=sample_size) if len(processor.df) > sample_size else processor.df
        
        # 使用压缩保存以节省空间
        df_sample.to_csv(
            config.output_path / f'RAW_f{config.scale_factor}_processed.csv.gz',
            index=False,
            compression='gzip'
        )
        
        # 保存统计结果
        angle_stats.to_csv(config.output_path / f'AngleStats_f{config.scale_factor}.csv')
        radial_stats.to_csv(config.output_path / f'RadialStats_f{config.scale_factor}.csv')
        particle_stats.to_csv(config.output_path / f'ParticleStats_f{config.scale_factor}.csv')
        pdf_data.to_csv(config.output_path / f'PDF_f{config.scale_factor}.csv', index=False)
        
        # 保存汇总统计为文本文件
        with open(config.output_path / f'Summary_f{config.scale_factor}.txt', 'w') as f:
            f.write(f"DEME Mixer Analysis Summary - Scale Factor f{config.scale_factor}\n")
            f.write("="*60 + "\n\n")
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
            f.write(f"\nSpatial Extent:\n")
            f.write(f"  X: [{summary['spatial_extent']['x_range'][0]:.4f}, {summary['spatial_extent']['x_range'][1]:.4f}]\n")
            f.write(f"  Y: [{summary['spatial_extent']['y_range'][0]:.4f}, {summary['spatial_extent']['y_range'][1]:.4f}]\n")
            f.write(f"  Z: [{summary['spatial_extent']['z_range'][0]:.4f}, {summary['spatial_extent']['z_range'][1]:.4f}]\n")
        
        # 6. 打印分析摘要
        print("\n" + "="*60)
        print("ANALYSIS SUMMARY")
        print("="*60)
        print(f"Total particles analyzed: {summary['total_particles']:,}")
        print(f"Number of frames: {summary['num_frames']}")
        print(f"Time range: {summary['time_range'][0]:.6f} - {summary['time_range'][1]:.6f} seconds")
        print(f"\nVelocity statistics:")
        print(f"  Mean velocity: {summary['velocity_stats']['mean']:.4f} m/s")
        print(f"  Std velocity:  {summary['velocity_stats']['std']:.4f} m/s")
        print(f"  Max velocity:  {summary['velocity_stats']['max']:.4f} m/s")
        print(f"\nParticle sizes found: {summary['particle_sizes']}")
        
        # 性能统计
        if 'velocity' in processor.df.columns:
            print(f"\n⚡ Performance: Pre-calculated velocity used")
        print(f"   Velocity processing time: {velocity_time:.3f} seconds")
        
        # 总运行时间
        elapsed = (datetime.datetime.now() - start_time).total_seconds()
        print(f"\nTotal processing time: {elapsed:.2f} seconds")
        print(f"Results saved to: {config.output_path}")
        
        # 内存使用估算
        memory_usage_mb = processor.df.memory_usage(deep=True).sum() / 1024 / 1024
        print(f"Peak memory usage (DataFrame): {memory_usage_mb:.2f} MB")
        
    except Exception as e:
        print(f"\n❌ Error: {e}")
        raise

if __name__ == "__main__":
    main()