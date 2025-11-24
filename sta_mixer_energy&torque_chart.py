#!/usr/bin/env python3
"""
SUP混合器数据可视化脚本 - 科学论文风格（改进版）
用于绘制不同scale_factor下的动能和扭矩对比图
包含数据平滑和改进的扭矩显示
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import os
from scipy import signal
from scipy.ndimage import uniform_filter1d

# 设置matplotlib参数以获得更好的显示效果
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['DejaVu Serif', 'Times New Roman']
plt.rcParams['font.size'] = 12
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['figure.dpi'] = 100
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['axes.linewidth'] = 1.8
plt.rcParams['xtick.major.width'] = 1.2
plt.rcParams['ytick.major.width'] = 1.2
plt.rcParams['xtick.minor.width'] = 0.6
plt.rcParams['ytick.minor.width'] = 0.6

class SUPMixerVisualizer:
    """SUP混合器数据可视化类（改进版）"""
    
    def __init__(self, data_dir="sta_results/torque&energy"):
        """
        初始化可视化器
        
        参数:
            data_dir: 包含CSV文件的目录路径
        """
        self.data_dir = data_dir
        self.ke_file = os.path.join(data_dir, "kinetic_energy_comparison.csv")
        self.torque_file = os.path.join(data_dir, "mixer_torque_comparison.csv")
        self.ke_df = None
        self.torque_df = None
        
    def load_data(self):
        """加载CSV数据文件"""
        print("=" * 60)
        print("加载SUP混合器数据")
        print("=" * 60)
        
        # 检查文件是否存在
        if not os.path.exists(self.ke_file):
            raise FileNotFoundError(f"找不到文件: {self.ke_file}")
        if not os.path.exists(self.torque_file):
            raise FileNotFoundError(f"找不到文件: {self.torque_file}")
        
        # 读取数据
        self.ke_df = pd.read_csv(self.ke_file)
        self.torque_df = pd.read_csv(self.torque_file)
        
        print(f"✓ 成功加载动能数据: {self.ke_file}")
        print(f"  - 数据形状: {self.ke_df.shape}")
        print(f"  - 列名: {self.ke_df.columns.tolist()}")
        
        print(f"\n✓ 成功加载扭矩数据: {self.torque_file}")
        print(f"  - 数据形状: {self.torque_df.shape}")
        print(f"  - 列名: {self.torque_df.columns.tolist()}")
        
        return True
    
    def smooth_data(self, data, window_size=50, method='moving_average'):
        """
        平滑数据以减少噪声
        
        参数:
            data: 原始数据数组
            window_size: 窗口大小
            method: 平滑方法 ('moving_average', 'savgol', 'gaussian')
        
        返回:
            平滑后的数据
        """
        if method == 'moving_average':
            # 移动平均
            return uniform_filter1d(data, size=window_size, mode='nearest')
        elif method == 'savgol':
            # Savitzky-Golay滤波器
            if len(data) > window_size:
                return signal.savgol_filter(data, window_size, 3)
            else:
                return data
        elif method == 'gaussian':
            # 高斯滤波
            return signal.gaussian_filter1d(data, sigma=window_size/5)
        else:
            return data
    
    def calculate_envelope(self, data, window_size=100):
        """
        计算数据的包络线（最大值和最小值）
        
        参数:
            data: 原始数据
            window_size: 窗口大小
        
        返回:
            (上包络, 下包络)
        """
        # 使用滚动窗口计算最大值和最小值
        upper_envelope = uniform_filter1d(np.maximum.reduceat(data, np.arange(0, len(data), window_size//10)), 
                                         size=10, mode='nearest')
        lower_envelope = uniform_filter1d(np.minimum.reduceat(data, np.arange(0, len(data), window_size//10)), 
                                         size=10, mode='nearest')
        
        # 插值回原始长度
        x_orig = np.arange(len(data))
        x_envelope = np.linspace(0, len(data)-1, len(upper_envelope))
        upper_envelope = np.interp(x_orig, x_envelope, upper_envelope)
        lower_envelope = np.interp(x_orig, x_envelope, lower_envelope)
        
        return upper_envelope, lower_envelope
    
    def calculate_statistics(self, df, data_cols):
        """
        计算数据统计信息
        
        参数:
            df: DataFrame
            data_cols: 数据列名列表
        
        返回:
            包含统计信息的字典
        """
        stats = {}
        for col in data_cols:
            stats[col] = {
                'min': df[col].min(),
                'max': df[col].max(),
                'mean': df[col].mean(),
                'std': df[col].std(),
                'median': df[col].median(),
                'q25': df[col].quantile(0.25),
                'q75': df[col].quantile(0.75)
            }
        return stats
    
    def plot_kinetic_energy(self, ax):
        """
        绘制动能对比图
        
        参数:
            ax: matplotlib轴对象
        """
        time_col = 'Time(s)'
        data_cols = [col for col in self.ke_df.columns if col != time_col]
        
        # 颜色映射 - 使用科学配色
        colors = plt.cm.Set1(np.linspace(0, 0.8, len(data_cols)))
        
        # 绘制每条曲线
        for i, col in enumerate(data_cols):
            time = self.ke_df[time_col]
            data = self.ke_df[col]
            
            # 绘制原始数据线条
            ax.plot(time, data, '-', color=colors[i], linewidth=2.4,
                   label=f'{col}', alpha=0.9, zorder=5)
            
            # 添加散点标记（每隔一定数量的点）
            marker_interval = max(1, len(time) // 20)
            ax.scatter(time[::marker_interval], data[::marker_interval],
                      s=60, color=colors[i], alpha=0.8, zorder=6, edgecolors='white',
                      linewidth=0.5)
        
        # 计算统计信息
        stats = self.calculate_statistics(self.ke_df, data_cols)
        
        # 创建统计信息文本框
        stats_text = "Statistics:\n" + "─" * 20 + "\n"
        for col in data_cols:
            s = stats[col]
            stats_text += (f"{col}:\n"
                          f"  Max: {s['max']:.2e} J\n"
                          f"  Mean: {s['mean']:.2e} J\n")
        
        # 添加统计信息框
        ax.text(0.98, 0.02, stats_text.strip(), transform=ax.transAxes,
               verticalalignment='bottom', horizontalalignment='right',
               fontsize=10, fontfamily='monospace',
               bbox=dict(boxstyle='round,pad=0.5', facecolor='white',
                        edgecolor='black', linewidth=1.2, alpha=0.95),
               zorder=10)
        
        # 设置标题和标签
        ax.set_title('(a) Kinetic Energy Evolution', 
                    fontsize=16, fontweight='bold', pad=15)
        ax.set_xlabel('Time [s]', fontsize=14, fontweight='bold')
        ax.set_ylabel('Kinetic Energy [J]', fontsize=14, fontweight='bold')
        
        # 科学风格坐标轴设置
        for spine in ax.spines.values():
            spine.set_linewidth(1.8)
            spine.set_color('black')
        
        # 设置刻度
        ax.tick_params(axis='both', which='major', labelsize=12, width=1.2, length=5,
                      direction='in', top=True, right=True)
        ax.tick_params(axis='both', which='minor', width=0.6, length=2.4,
                      direction='in', top=True, right=True)
        
        # 网格设置
        ax.grid(True, which='major', alpha=0.15, linewidth=0.8, linestyle='-', color='gray')
        ax.grid(True, which='minor', alpha=0.08, linewidth=0.5, linestyle='--', color='gray')
        ax.minorticks_on()
        
        # 使用科学计数法
        ax.ticklabel_format(style='scientific', axis='y', scilimits=(0,0))
        
        # 图例设置
        legend = ax.legend(loc='upper left', frameon=True, fancybox=False,
                          shadow=False, framealpha=1.0, edgecolor='black',
                          facecolor='white', fontsize=11, ncol=1,
                          bbox_to_anchor=(0.02, 0.98))
        legend.get_frame().set_linewidth(1.2)
        legend.set_zorder(12)
    
    def plot_mixer_torque_improved(self, ax, smoothing_method='moving_average', 
                                   window_size=100, show_envelope=True):
        """
        改进的扭矩绘图函数，包含平滑和包络线显示
        
        参数:
            ax: matplotlib轴对象
            smoothing_method: 平滑方法
            window_size: 窗口大小
            show_envelope: 是否显示包络线
        """
        time_col = 'Time(s)'
        data_cols = [col for col in self.torque_df.columns if col != time_col]
        
        # 颜色映射
        colors = plt.cm.Set2(np.linspace(0, 0.8, len(data_cols)))
        
        # 绘制每条曲线
        for i, col in enumerate(data_cols):
            time = self.torque_df[time_col]
            data = self.torque_df[col].values
            
            if show_envelope and len(data) > window_size:
                # 显示包络线和平滑曲线
                # 1. 绘制原始数据的阴影区域（使用包络线）
                upper, lower = self.calculate_envelope(data, window_size)
                ax.fill_between(time, lower, upper, alpha=0.2, color=colors[i], 
                               label=f'{col} (range)')
                
                # 2. 绘制平滑后的中心线
                smoothed = self.smooth_data(data, window_size, smoothing_method)
                ax.plot(time, smoothed, '-', color=colors[i], linewidth=2.4,
                       label=f'{col} (smoothed)', alpha=0.9, zorder=5)
                
            else:
                # 只绘制平滑曲线
                smoothed = self.smooth_data(data, window_size, smoothing_method)
                ax.plot(time, smoothed, '-', color=colors[i], linewidth=2.4,
                       label=f'{col}', alpha=0.9, zorder=5)
                
                # 可选：显示原始数据（半透明）
                ax.plot(time, data, '-', color=colors[i], linewidth=0.5,
                       alpha=0.2, zorder=3)
        
        # 计算统计信息（使用原始数据）
        stats = self.calculate_statistics(self.torque_df, data_cols)
        
        # 创建统计信息文本框
        stats_text = "Statistics:\n" + "─" * 20 + "\n"
        for col in data_cols:
            s = stats[col]
            stats_text += (f"{col}:\n"
                          f"  Max: {abs(s['max']):.2e} Nm\n"
                          f"  Mean: {abs(s['mean']):.2e} Nm\n")
        
        # 添加统计信息框
        ax.text(0.98, 0.02, stats_text.strip(), transform=ax.transAxes,
               verticalalignment='bottom', horizontalalignment='right',
               fontsize=10, fontfamily='monospace',
               bbox=dict(boxstyle='round,pad=0.5', facecolor='white',
                        edgecolor='black', linewidth=1.2, alpha=0.95),
               zorder=10)
        
        # 设置标题和标签
        ax.set_title('(b) Mixer Torque Evolution (Smoothed)', 
                    fontsize=16, fontweight='bold', pad=15)
        ax.set_xlabel('Time [s]', fontsize=14, fontweight='bold')
        ax.set_ylabel('Mixer Torque [Nm]', fontsize=14, fontweight='bold')
        
        # 科学风格坐标轴设置
        for spine in ax.spines.values():
            spine.set_linewidth(1.8)
            spine.set_color('black')
        
        # 设置刻度
        ax.tick_params(axis='both', which='major', labelsize=12, width=1.2, length=5,
                      direction='in', top=True, right=True)
        ax.tick_params(axis='both', which='minor', width=0.6, length=2.4,
                      direction='in', top=True, right=True)
        
        # 网格设置
        ax.grid(True, which='major', alpha=0.15, linewidth=0.8, linestyle='-', color='gray')
        ax.grid(True, which='minor', alpha=0.08, linewidth=0.5, linestyle='--', color='gray')
        ax.minorticks_on()
        
        # 使用科学计数法
        ax.ticklabel_format(style='scientific', axis='y', scilimits=(0,0))
        
        # 图例设置
        legend = ax.legend(loc='upper left', frameon=True, fancybox=False,
                          shadow=False, framealpha=1.0, edgecolor='black',
                          facecolor='white', fontsize=11, ncol=1,
                          bbox_to_anchor=(0.02, 0.98))
        legend.get_frame().set_linewidth(1.2)
        legend.set_zorder(12)
    
    def plot_torque_analysis(self):
        """创建详细的扭矩分析图（包含多种显示方式）"""
        fig = plt.figure(figsize=(16, 12))
        
        # 创建子图布局
        gs = fig.add_gridspec(3, 2, hspace=0.3, wspace=0.25)
        
        # 1. 原始数据
        ax1 = fig.add_subplot(gs[0, 0])
        self.plot_mixer_torque_original(ax1)
        ax1.set_title('(a) Original Torque Data', fontsize=14, fontweight='bold')
        
        # 2. 移动平均平滑
        ax2 = fig.add_subplot(gs[0, 1])
        self.plot_mixer_torque_improved(ax2, smoothing_method='moving_average', 
                                       window_size=100, show_envelope=False)
        ax2.set_title('(b) Moving Average Smoothed', fontsize=14, fontweight='bold')
        
        # 3. Savitzky-Golay平滑
        ax3 = fig.add_subplot(gs[1, 0])
        self.plot_mixer_torque_improved(ax3, smoothing_method='savgol', 
                                       window_size=51, show_envelope=False)
        ax3.set_title('(c) Savitzky-Golay Filtered', fontsize=14, fontweight='bold')
        
        # 4. 包络线显示
        ax4 = fig.add_subplot(gs[1, 1])
        self.plot_mixer_torque_improved(ax4, smoothing_method='moving_average', 
                                       window_size=100, show_envelope=True)
        ax4.set_title('(d) Envelope Display', fontsize=14, fontweight='bold')
        
        # 5. 频谱分析
        ax5 = fig.add_subplot(gs[2, :])
        self.plot_torque_spectrum(ax5)
        ax5.set_title('(e) Frequency Spectrum Analysis', fontsize=14, fontweight='bold')
        
        fig.suptitle('Torque Data Analysis - Multiple Views', 
                    fontsize=18, fontweight='bold', y=1.02)
        
        plt.tight_layout()
        return fig
    
    def plot_mixer_torque_original(self, ax):
        """绘制原始扭矩数据（未平滑）"""
        time_col = 'Time(s)'
        data_cols = [col for col in self.torque_df.columns if col != time_col]
        colors = plt.cm.Set2(np.linspace(0, 0.8, len(data_cols)))
        
        for i, col in enumerate(data_cols):
            time = self.torque_df[time_col]
            data = self.torque_df[col]
            ax.plot(time, data, '-', color=colors[i], linewidth=0.8,
                   label=f'{col}', alpha=0.7)
        
        ax.set_xlabel('Time [s]', fontsize=12)
        ax.set_ylabel('Torque [Nm]', fontsize=12)
        ax.legend(loc='upper right', fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.ticklabel_format(style='scientific', axis='y', scilimits=(0,0))
    
    def plot_torque_spectrum(self, ax):
        """绘制扭矩的频谱分析"""
        time_col = 'Time(s)'
        data_cols = [col for col in self.torque_df.columns if col != time_col]
        colors = plt.cm.Set2(np.linspace(0, 0.8, len(data_cols)))
        
        # 计算采样频率
        dt = self.torque_df[time_col].iloc[1] - self.torque_df[time_col].iloc[0]
        fs = 1.0 / dt
        
        for i, col in enumerate(data_cols):
            data = self.torque_df[col].values
            
            # 计算FFT
            fft = np.fft.rfft(data)
            freqs = np.fft.rfftfreq(len(data), d=dt)
            
            # 绘制功率谱
            power = np.abs(fft) ** 2
            ax.semilogy(freqs[1:], power[1:], '-', color=colors[i], 
                       linewidth=1.5, label=f'{col}', alpha=0.8)
        
        ax.set_xlabel('Frequency [Hz]', fontsize=12)
        ax.set_ylabel('Power Spectral Density', fontsize=12)
        ax.set_xlim([0, fs/2])  # Nyquist频率
        ax.legend(loc='upper right', fontsize=10)
        ax.grid(True, alpha=0.3)
    
    def plot_combined_improved(self, smoothing_window=100):
        """创建改进的组合图（动能+平滑扭矩）"""
        fig, axes = plt.subplots(2, 1, figsize=(14, 10), sharex=True)
        fig.suptitle('SUP Mixer Simulation Results (Improved)', 
                    fontsize=18, fontweight='bold', y=1.02)
        
        # 绘制动能图
        self.plot_kinetic_energy(axes[0])
        
        # 绘制改进的扭矩图
        self.plot_mixer_torque_improved(axes[1], 
                                       smoothing_method='moving_average',
                                       window_size=smoothing_window,
                                       show_envelope=True)
        
        # 调整布局
        plt.tight_layout()
        
        return fig
    
    def print_summary_statistics(self):
        """打印详细的统计摘要"""
        print("\n" + "=" * 60)
        print("数据统计摘要")
        print("=" * 60)
        
        time_col = 'Time(s)'
        
        # 动能统计
        print("\n动能数据统计:")
        print("-" * 40)
        data_cols = [col for col in self.ke_df.columns if col != time_col]
        for col in data_cols:
            print(f"\n{col}:")
            print(f"  最小值: {self.ke_df[col].min():.6e} J")
            print(f"  最大值: {self.ke_df[col].max():.6e} J")
            print(f"  平均值: {self.ke_df[col].mean():.6e} J")
            print(f"  标准差: {self.ke_df[col].std():.6e} J")
        
        # 扭矩统计
        print("\n扭矩数据统计:")
        print("-" * 40)
        data_cols = [col for col in self.torque_df.columns if col != time_col]
        for col in data_cols:
            print(f"\n{col}:")
            print(f"  最小值: {self.torque_df[col].min():.6e} Nm")
            print(f"  最大值: {self.torque_df[col].max():.6e} Nm") 
            print(f"  平均值: {self.torque_df[col].mean():.6e} Nm")
            print(f"  标准差: {self.torque_df[col].std():.6e} Nm")
            
            # 计算信噪比
            signal_power = self.torque_df[col].mean() ** 2
            noise_power = self.torque_df[col].var()
            if noise_power > 0:
                snr_db = 10 * np.log10(abs(signal_power / noise_power))
                print(f"  信噪比: {snr_db:.2f} dB")
        
        # 时间范围
        print(f"\n时间范围: {self.ke_df[time_col].min():.3f} - {self.ke_df[time_col].max():.3f} 秒")
        print(f"数据点数: {len(self.ke_df)}")

def main():
    """主函数"""
    # 创建可视化器实例
    visualizer = SUPMixerVisualizer(data_dir="sta_results/torque&energy")
    
    try:
        # 加载数据
        visualizer.load_data()
        
        # 打印统计摘要
        visualizer.print_summary_statistics()
        
        # 询问用户绘图选项
        print("\n" + "=" * 60)
        print("绘图选项:")
        print("1. 改进的组合图（动能+平滑扭矩）")
        print("2. 详细扭矩分析（多种视图）")
        print("3. 两种都生成")
        print("=" * 60)
        
        choice = input("请选择 (1/2/3) [默认: 1]: ").strip() or "1"
        
        # 询问平滑窗口大小
        window_input = input("扭矩平滑窗口大小 [默认: 100]: ").strip()
        window_size = int(window_input) if window_input else 100
        
        if choice in ["1", "3"]:
            print("\n生成改进的组合图...")
            fig_combined = visualizer.plot_combined_improved(smoothing_window=window_size)
            
        if choice in ["2", "3"]:
            print("\n生成扭矩分析图...")
            fig_analysis = visualizer.plot_torque_analysis()
        
        # 显示图形
        plt.show()
        
        # 询问是否保存
        save_option = input("\n保存图形? (y/n) [默认: n]: ").strip().lower()
        if save_option == 'y':
            output_dir = "sta_results"
            os.makedirs(output_dir, exist_ok=True)
            
            if choice in ["1", "3"]:
                output_path = os.path.join(output_dir, "sup_mixer_improved.png")
                fig_combined.savefig(output_path, dpi=300, bbox_inches='tight')
                print(f"✓ 改进图已保存: {output_path}")
            
            if choice in ["2", "3"]:
                output_path = os.path.join(output_dir, "torque_analysis.png")
                fig_analysis.savefig(output_path, dpi=300, bbox_inches='tight')
                print(f"✓ 扭矩分析图已保存: {output_path}")
        
        print("\n程序执行完成！")
        
    except FileNotFoundError as e:
        print(f"\n错误: {e}")
        print("请确保已运行数据处理脚本生成CSV文件。")
    except Exception as e:
        print(f"\n发生错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()