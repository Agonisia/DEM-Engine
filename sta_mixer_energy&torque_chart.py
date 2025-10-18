#!/usr/bin/env python3
"""
SUP混合器数据可视化脚本 - 科学论文风格
用于绘制不同scale_factor下的动能和扭矩对比图
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import os

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
    """SUP混合器数据可视化类"""
    
    def __init__(self, data_dir="sta_result/torque&energy"):
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
            # 为了清晰显示，每隔一定点数显示标记
            time = self.ke_df[time_col]
            data = self.ke_df[col]
            
            # 绘制线条
            ax.plot(time, data, '-', color=colors[i], linewidth=2.4,
                   label=f'Scale Factor: {col}', alpha=0.9, zorder=5)
            
            # 添加散点标记（每隔一定数量的点）
            marker_interval = max(1, len(time) // 20)  # 显示约20个标记点
            ax.scatter(time[::marker_interval], data[::marker_interval],
                      s=60, color=colors[i], alpha=0.8, zorder=6, edgecolors='white',
                      linewidth=0.5)
        
        # 计算统计信息
        stats = self.calculate_statistics(self.ke_df, data_cols)
        
        # 创建统计信息文本框
        stats_text = "Kinetic Energy Statistics:\n" + "─" * 30 + "\n"
        for col in data_cols:
            s = stats[col]
            stats_text += (f"{col}:\n"
                          f"  Max: {s['max']:.3e} J\n"
                          f"  Mean: {s['mean']:.3e} J\n"
                          f"  Std: {s['std']:.3e} J\n")
        
        # 添加统计信息框
        ax.text(0.98, 0.02, stats_text.strip(), transform=ax.transAxes,
               verticalalignment='bottom', horizontalalignment='right',
               fontsize=10, fontfamily='monospace',
               bbox=dict(boxstyle='round,pad=0.5', facecolor='white',
                        edgecolor='black', linewidth=1.2, alpha=0.95),
               zorder=10)
        
        # 设置标题和标签
        ax.set_title('(a) Kinetic Energy Evolution Comparison', 
                    fontsize=19, fontweight='bold', pad=15)
        ax.set_xlabel('Time [s]', fontsize=17, fontweight='bold')
        ax.set_ylabel('Kinetic Energy [J]', fontsize=17, fontweight='bold')
        
        # 科学风格坐标轴设置
        for spine in ax.spines.values():
            spine.set_linewidth(1.8)
            spine.set_color('black')
        
        # 设置刻度
        ax.tick_params(axis='both', which='major', labelsize=14, width=1.2, length=5,
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
                          facecolor='white', fontsize=12, ncol=1,
                          bbox_to_anchor=(0.02, 0.98))
        legend.get_frame().set_linewidth(1.2)
        legend.set_zorder(12)
    
    def plot_mixer_torque(self, ax):
        """
        绘制混合器扭矩对比图
        
        参数:
            ax: matplotlib轴对象
        """
        time_col = 'Time(s)'
        data_cols = [col for col in self.torque_df.columns if col != time_col]
        
        # 颜色映射 - 使用科学配色
        colors = plt.cm.Set2(np.linspace(0, 0.8, len(data_cols)))
        
        # 绘制每条曲线
        for i, col in enumerate(data_cols):
            time = self.torque_df[time_col]
            data = self.torque_df[col]
            
            # 绘制线条
            ax.plot(time, data, '-', color=colors[i], linewidth=2.4,
                   label=f'Scale Factor: {col}', alpha=0.9, zorder=5)
            
            # 添加散点标记
            marker_interval = max(1, len(time) // 20)
            ax.scatter(time[::marker_interval], data[::marker_interval],
                      s=60, color=colors[i], alpha=0.8, zorder=6, edgecolors='white',
                      linewidth=0.5, marker='s')  # 使用方形标记以区分
        
        # 计算统计信息
        stats = self.calculate_statistics(self.torque_df, data_cols)
        
        # 创建统计信息文本框
        stats_text = "Mixer Torque Statistics:\n" + "─" * 30 + "\n"
        for col in data_cols:
            s = stats[col]
            stats_text += (f"{col}:\n"
                          f"  Max: {s['max']:.3e} Nm\n"
                          f"  Mean: {s['mean']:.3e} Nm\n"
                          f"  Std: {s['std']:.3e} Nm\n")
        
        # 添加统计信息框
        ax.text(0.98, 0.02, stats_text.strip(), transform=ax.transAxes,
               verticalalignment='bottom', horizontalalignment='right',
               fontsize=10, fontfamily='monospace',
               bbox=dict(boxstyle='round,pad=0.5', facecolor='white',
                        edgecolor='black', linewidth=1.2, alpha=0.95),
               zorder=10)
        
        # 设置标题和标签
        ax.set_title('(b) Mixer Torque Evolution Comparison', 
                    fontsize=19, fontweight='bold', pad=15)
        ax.set_xlabel('Time [s]', fontsize=17, fontweight='bold')
        ax.set_ylabel('Mixer Torque [Nm]', fontsize=17, fontweight='bold')
        
        # 科学风格坐标轴设置
        for spine in ax.spines.values():
            spine.set_linewidth(1.8)
            spine.set_color('black')
        
        # 设置刻度
        ax.tick_params(axis='both', which='major', labelsize=14, width=1.2, length=5,
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
                          facecolor='white', fontsize=12, ncol=1,
                          bbox_to_anchor=(0.02, 0.98))
        legend.get_frame().set_linewidth(1.2)
        legend.set_zorder(12)
    
    def plot_combined(self):
        """创建组合图（两个子图）"""
        # 创建图形
        fig, axes = plt.subplots(2, 1, figsize=(14, 10), sharex=True)
        fig.suptitle('SUP Mixer Simulation Results Comparison', 
                    fontsize=22, fontweight='bold', y=1.02)
        
        # 绘制动能图
        self.plot_kinetic_energy(axes[0])
        
        # 绘制扭矩图
        self.plot_mixer_torque(axes[1])
        
        # 调整布局
        plt.tight_layout()
        
        return fig
    
    def plot_separate(self):
        """创建分离的图形（每个一个独立图形）"""
        # 动能图
        fig1, ax1 = plt.subplots(1, 1, figsize=(14, 6))
        self.plot_kinetic_energy(ax1)
        plt.tight_layout()
        
        # 扭矩图
        fig2, ax2 = plt.subplots(1, 1, figsize=(14, 6))
        self.plot_mixer_torque(ax2)
        plt.tight_layout()
        
        return fig1, fig2
    
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
            print(f"  中位数: {self.ke_df[col].median():.6e} J")
        
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
            print(f"  中位数: {self.torque_df[col].median():.6e} Nm")
        
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
        print("1. 组合图（两个子图在一个图形中）")
        print("2. 分离图（两个独立的图形）")
        print("3. 两种都生成")
        print("=" * 60)
        
        choice = input("请选择 (1/2/3) [默认: 1]: ").strip() or "1"
        
        if choice in ["1", "3"]:
            print("\n生成组合图...")
            fig_combined = visualizer.plot_combined()
            
        if choice in ["2", "3"]:
            print("\n生成分离图...")
            fig_ke, fig_torque = visualizer.plot_separate()
        
        # 显示图形
        plt.show()
        
        # 询问是否保存
        save_option = input("\n保存图形? (y/n) [默认: n]: ").strip().lower()
        if save_option == 'y':
            output_dir = "sta_result"
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            
            if choice in ["1", "3"]:
                output_path = os.path.join(output_dir, "sup_mixer_combined.png")
                fig_combined.savefig(output_path, dpi=300, bbox_inches='tight')
                print(f"✓ 组合图已保存: {output_path}")
            
            if choice in ["2", "3"]:
                output_path_ke = os.path.join(output_dir, "sup_mixer_kinetic_energy.png")
                output_path_torque = os.path.join(output_dir, "sup_mixer_torque.png")
                fig_ke.savefig(output_path_ke, dpi=300, bbox_inches='tight')
                fig_torque.savefig(output_path_torque, dpi=300, bbox_inches='tight')
                print(f"✓ 动能图已保存: {output_path_ke}")
                print(f"✓ 扭矩图已保存: {output_path_torque}")
        
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