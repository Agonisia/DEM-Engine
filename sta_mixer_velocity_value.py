#!/usr/bin/env python3
"""
SUP混合器数据预处理脚本（GPU加速版）
功能：检查并添加velocity列到mixer_output CSV文件
适配目录格式：SUPMixerOutput_f{scale_factor}se{cohesion_code}
作者：DEM-Engine Analysis Tools
"""

import os
import sys
import re
import glob
import shutil
import argparse
import numpy as np
from pathlib import Path
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

# 尝试导入GPU加速库
GPU_AVAILABLE = False
try:
    import cupy as cp
    import cudf
    GPU_AVAILABLE = True
    print("🚀 GPU加速已启用 (cuDF/CuPy)")
except ImportError:
    try:
        import rapids.cudf as cudf
        import cupy as cp
        GPU_AVAILABLE = True
        print("🚀 GPU加速已启用 (RAPIDS cuDF)")
    except ImportError:
        import pandas as pd
        print("💻 使用CPU模式 (提示: 安装cuDF和CuPy可启用GPU加速)")
        print("   pip install cudf-cu11 cupy-cuda11x  # 对于CUDA 11.x")
        print("   pip install cudf-cu12 cupy-cuda12x  # 对于CUDA 12.x")

# 如果没有GPU加速，使用pandas
if not GPU_AVAILABLE:
    import pandas as pd

class DEMEPreprocessor:
    """DEME数据预处理器类（支持GPU加速）"""
    
    def __init__(self, base_path='build', backup=False, verbose=True, use_gpu=None):
        self.base_path = Path(base_path)
        self.backup = backup
        self.verbose = verbose
        self.processed_files = []
        self.skipped_files = []
        self.error_files = []
        
        # 决定是否使用GPU
        if use_gpu is None:
            self.use_gpu = GPU_AVAILABLE
        else:
            self.use_gpu = use_gpu and GPU_AVAILABLE
        
        if self.use_gpu:
            print("✅ GPU加速模式")
            # 设置GPU内存池以提高性能
            try:
                mempool = cp.get_default_memory_pool()
                mempool.set_limit(size=4 * 1024**3)  # 限制4GB
            except:
                pass
        else:
            print("ℹ️ CPU模式")
        
    def find_sup_directories(self, pattern_str=None):
        """
        查找所有符合SUPMixerOutput_f{scale_factor}se{cohesion}格式的目录
        
        参数:
            pattern_str: 自定义正则表达式模式（可选）
        
        返回:
            包含(目录路径, scale_factor, cohesion_code)元组的列表
        """
        # 默认模式：匹配 SUPMixerOutput_f{数字}se{三位数字}
        if pattern_str is None:
            pattern = re.compile(r'SUPMixerOutput_f(\d+(?:\.\d+)?)se(\d{3})')
        else:
            pattern = re.compile(pattern_str)
        
        directories = []
        
        # 确保基础路径存在
        if not os.path.exists(self.base_path):
            if self.verbose:
                print(f"⚠️ 警告：路径 {self.base_path} 不存在")
            return directories
        
        # 遍历目录查找匹配的文件夹
        for item in os.listdir(self.base_path):
            item_path = os.path.join(self.base_path, item)
            if os.path.isdir(item_path):
                match = pattern.match(item)
                if match:
                    scale_factor = float(match.group(1))
                    # 如果有cohesion代码，也提取出来
                    cohesion_code = match.group(2) if len(match.groups()) > 1 else "000"
                    directories.append((item_path, scale_factor, cohesion_code))
        
        # 按scale_factor排序
        directories.sort(key=lambda x: x[1])
        
        if self.verbose:
            print(f"\n📋 找到 {len(directories)} 个SUP输出目录:")
            for dir_path, sf, cc in directories:
                print(f"   - {os.path.basename(dir_path)} (scale_factor={sf}, cohesion_code={cc})")
        
        return directories
    
    def check_csv_columns(self, file_path):
        """检查CSV文件的列（优化版）"""
        try:
            # 只读取第一行来检查列名
            if self.use_gpu:
                # GPU模式：使用cuDF读取
                df = cudf.read_csv(file_path, nrows=1)
                columns = df.columns.tolist()
            else:
                # CPU模式：使用pandas
                df = pd.read_csv(file_path, nrows=1)
                columns = df.columns.tolist()
            
            has_velocity = 'velocity' in columns
            has_components = all(col in columns for col in ['v_x', 'v_y', 'v_z'])
            
            return {
                'columns': columns,
                'has_velocity': has_velocity,
                'has_components': has_components,
                'needs_processing': (not has_velocity) and has_components
            }
        except Exception as e:
            if self.verbose:
                print(f"    ❌ Error reading {file_path}: {e}")
            return None
    
    def add_velocity_column_gpu(self, file_path, output_path=None):
        """使用GPU添加velocity列到CSV文件"""
        try:
            # 使用cuDF读取文件
            df = cudf.read_csv(file_path)
            
            # 检查是否有必要的列
            if not all(col in df.columns for col in ['v_x', 'v_y', 'v_z']):
                raise ValueError("Missing velocity component columns")
            
            # GPU上计算velocity
            df['velocity'] = cp.sqrt(df['v_x']**2 + df['v_y']**2 + df['v_z']**2)
            
            # 重新排列列顺序（将velocity放在v_z后面）
            cols = df.columns.tolist()
            cols.remove('velocity')
            v_z_index = cols.index('v_z')
            cols.insert(v_z_index + 1, 'velocity')
            df = df[cols]
            
            # 保存文件
            if output_path is None:
                output_path = file_path
            
            # cuDF直接写入CSV
            df.to_csv(output_path, index=False, float_precision='round_trip')
            return True
            
        except Exception as e:
            if self.verbose:
                print(f"    ⚠️ GPU处理失败，尝试CPU: {e}")
            # 如果GPU失败，回退到CPU处理
            return self.add_velocity_column_cpu(file_path, output_path)
    
    def add_velocity_column_cpu(self, file_path, output_path=None):
        """使用CPU添加velocity列到CSV文件"""
        try:
            # 读取完整文件
            df = pd.read_csv(file_path)
            
            # 检查是否有必要的列
            if not all(col in df.columns for col in ['v_x', 'v_y', 'v_z']):
                raise ValueError("Missing velocity component columns")
            
            # 计算velocity（向量化操作）
            df['velocity'] = np.sqrt(df['v_x']**2 + df['v_y']**2 + df['v_z']**2)
            
            # 重新排列列顺序（将velocity放在v_z后面）
            cols = df.columns.tolist()
            cols.remove('velocity')
            v_z_index = cols.index('v_z')
            cols.insert(v_z_index + 1, 'velocity')
            df = df[cols]
            
            # 保存文件
            if output_path is None:
                output_path = file_path
            
            df.to_csv(output_path, index=False, float_format='%.6f')
            return True
            
        except Exception as e:
            if self.verbose:
                print(f"    ❌ Error processing {file_path}: {e}")
            return False
    
    def add_velocity_column(self, file_path, output_path=None):
        """添加velocity列到CSV文件（自动选择GPU或CPU）"""
        if self.use_gpu:
            return self.add_velocity_column_gpu(file_path, output_path)
        else:
            return self.add_velocity_column_cpu(file_path, output_path)
    
    def process_batch_gpu(self, file_paths):
        """GPU批处理多个文件"""
        results = []
        
        # 批量处理以提高GPU利用率
        for file_path in file_paths:
            success = self.add_velocity_column(file_path)
            results.append(success)
        
        # 清理GPU内存
        if self.use_gpu:
            try:
                mempool = cp.get_default_memory_pool()
                mempool.free_all_blocks()
            except:
                pass
        
        return results
    
    def process_directory(self, dir_path, scale_factor, cohesion_code):
        """处理单个SUPMixerOutput目录"""
        dir_name = os.path.basename(dir_path)
        print(f"\n📁 Processing {dir_name}")
        print(f"   Path: {dir_path}")
        print(f"   Scale Factor: {scale_factor}")
        print(f"   Cohesion Code: {cohesion_code}")
        print(f"   Mode: {'GPU加速' if self.use_gpu else 'CPU'}")
        
        # 查找所有mixer_output_*.csv文件
        pattern = os.path.join(dir_path, "mixer_output_*.csv")
        csv_files = sorted(glob.glob(pattern))
        
        if not csv_files:
            print(f"   ⚠️ No mixer_output_*.csv files found")
            return
        
        print(f"   Found {len(csv_files)} CSV files")
        
        # 显示备份状态
        if self.backup:
            print(f"   🔒 Backup: ENABLED")
            backup_dir = Path(dir_path) / "backup_original"
            backup_dir.mkdir(exist_ok=True)
        else:
            print(f"   ⚡ Backup: DISABLED (direct modification)")
        
        # 统计
        processed_count = 0
        skipped_count = 0
        error_count = 0
        
        # 批处理准备
        files_to_process = []
        
        # 开始计时
        import time
        start_time = time.time()
        
        # 第一阶段：检查文件
        print("   🔍 检查文件...")
        for i, file_path in enumerate(csv_files):
            file_name = os.path.basename(file_path)
            
            # 显示进度
            if i % 100 == 0 and i > 0 and self.verbose:
                print(f"      已检查 {i}/{len(csv_files)} 文件...")
            
            # 检查文件
            check_result = self.check_csv_columns(file_path)
            
            if check_result is None:
                error_count += 1
                self.error_files.append(file_path)
                continue
            
            if not check_result['needs_processing']:
                if check_result['has_velocity']:
                    skipped_count += 1
                    self.skipped_files.append(file_path)
                else:
                    if self.verbose:
                        print(f"    ⚠️ {file_name}: Missing v_x, v_y, or v_z columns")
                    error_count += 1
                    self.error_files.append(file_path)
                continue
            
            # 备份原文件（如果启用）
            if self.backup:
                backup_dir = Path(dir_path) / "backup_original"
                backup_path = backup_dir / file_name
                if not backup_path.exists():
                    shutil.copy2(file_path, backup_path)
            
            files_to_process.append(file_path)
        
        # 第二阶段：批量处理
        if files_to_process:
            print(f"   ⚙️ 处理 {len(files_to_process)} 个文件...")
            
            # GPU批处理或逐个CPU处理
            if self.use_gpu and len(files_to_process) > 10:
                # GPU批处理
                batch_size = 50
                for i in range(0, len(files_to_process), batch_size):
                    batch = files_to_process[i:i+batch_size]
                    if self.verbose:
                        print(f"      处理批次 {i//batch_size + 1}/{(len(files_to_process)-1)//batch_size + 1}...")
                    
                    results = self.process_batch_gpu(batch)
                    for file_path, success in zip(batch, results):
                        if success:
                            processed_count += 1
                            self.processed_files.append(file_path)
                        else:
                            error_count += 1
                            self.error_files.append(file_path)
            else:
                # 逐个处理
                for i, file_path in enumerate(files_to_process):
                    if i % 50 == 0 and i > 0 and self.verbose:
                        print(f"      已处理 {i}/{len(files_to_process)} 文件...")
                    
                    if self.add_velocity_column(file_path):
                        processed_count += 1
                        self.processed_files.append(file_path)
                    else:
                        error_count += 1
                        self.error_files.append(file_path)
        
        # 计算耗时
        elapsed_time = time.time() - start_time
        
        # 打印统计
        print(f"\n   ✅ Summary for f{scale_factor}se{cohesion_code}:")
        print(f"      - Processed: {processed_count} files (velocity column added)")
        print(f"      - Skipped: {skipped_count} files (already have velocity)")
        print(f"      - Errors: {error_count} files")
        print(f"      - Time: {elapsed_time:.2f} seconds")
        if processed_count > 0:
            print(f"      - Speed: {processed_count/elapsed_time:.1f} files/sec")
        
        if self.backup and processed_count > 0:
            backup_dir = Path(dir_path) / "backup_original"
            print(f"      - Backup saved to: {backup_dir}")
    
    def process_all(self, scale_factors=None, cohesion_filter=None):
        """
        处理所有或指定的scale factor目录
        
        参数:
            scale_factors: 要处理的scale factor列表（可选）
            cohesion_filter: 只处理特定cohesion代码的数据（可选）
        """
        print("="*60)
        print("SUP混合器数据预处理 - 添加velocity列")
        print("="*60)
        print(f"设置:")
        print(f"  - 基础路径: {self.base_path}")
        print(f"  - 备份: {'启用' if self.backup else '禁用'}")
        print(f"  - 详细输出: {'开' if self.verbose else '关'}")
        print(f"  - GPU加速: {'启用' if self.use_gpu else '禁用'}")
        if cohesion_filter:
            print(f"  - Cohesion过滤: {cohesion_filter}")
        
        # 显示GPU信息
        if self.use_gpu:
            try:
                import cupy
                print(f"  - GPU设备: {cupy.cuda.runtime.getDeviceProperties(0)['name'].decode()}")
            except:
                pass
        
        # 查找所有目录
        directories = self.find_sup_directories()
        
        if not directories:
            return
        
        # 按cohesion过滤
        if cohesion_filter is not None:
            directories = [(d, sf, cc) for d, sf, cc in directories if cc == cohesion_filter]
            if not directories:
                print(f"⚠️ 没有找到cohesion_code={cohesion_filter}的目录")
                return
        
        # 筛选要处理的目录（按scale factor）
        if scale_factors:
            directories = [(d, sf, cc) for d, sf, cc in directories if sf in scale_factors]
            if not directories:
                print(f"⚠️ 没有找到指定scale factors的目录: {scale_factors}")
                return
        
        print(f"\n将处理 {len(directories)} 个目录")
        
        # 总计时
        import time
        total_start = time.time()
        
        # 处理每个目录
        for dir_path, scale_factor, cohesion_code in directories:
            self.process_directory(dir_path, scale_factor, cohesion_code)
        
        # 总耗时
        total_elapsed = time.time() - total_start
        
        # 总体统计
        self.print_final_summary(total_elapsed)
    
    def print_final_summary(self, total_time=None):
        """打印最终统计摘要"""
        print("\n" + "="*60)
        print("最终统计")
        print("="*60)
        print(f"✅ 总处理文件数: {len(self.processed_files)}")
        print(f"⏭️  总跳过文件数: {len(self.skipped_files)}")
        print(f"❌ 总错误文件数: {len(self.error_files)}")
        
        if total_time:
            print(f"⏱️  总耗时: {total_time:.2f} 秒")
            total_files = len(self.processed_files) + len(self.skipped_files)
            if total_files > 0:
                print(f"📊 平均速度: {total_files/total_time:.1f} 文件/秒")
        
        if not self.backup and len(self.processed_files) > 0:
            print("\n⚠️  注意: 原始文件已被直接修改（没有备份）")
            print("   要创建备份，请使用 --backup 参数")
        
        if self.error_files and self.verbose:
            print("\n⚠️ 出错的文件:")
            for file in self.error_files[:10]:  # 只显示前10个
                print(f"   - {file}")
            if len(self.error_files) > 10:
                print(f"   ... 还有 {len(self.error_files) - 10} 个")
        
        # GPU内存清理
        if self.use_gpu:
            try:
                import cupy as cp
                mempool = cp.get_default_memory_pool()
                used_bytes = mempool.used_bytes()
                total_bytes = mempool.total_bytes()
                if used_bytes > 0:
                    print(f"\n📊 GPU内存使用: {used_bytes/1024**2:.1f} MB / {total_bytes/1024**2:.1f} MB")
                mempool.free_all_blocks()
            except:
                pass

def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description='SUP混合器数据预处理 - 为CSV文件添加velocity列（GPU加速版）',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
默认行为:
  直接运行脚本将处理所有目录:
  python sta_mixer_velocity_value.py
  
其他示例:
  # 处理特定的scale factors
  python sta_mixer_velocity_value.py --scale 1 2 4
  
  # 处理特定的cohesion代码
  python sta_mixer_velocity_value.py --cohesion 000
  
  # 处理并创建备份
  python sta_mixer_velocity_value.py --backup
  
  # 使用自定义基础路径
  python sta_mixer_velocity_value.py --base-path /path/to/data
  
  # 禁用GPU加速（强制CPU模式）
  python sta_mixer_velocity_value.py --cpu
  
  # 带备份的静默模式
  python sta_mixer_velocity_value.py --backup --quiet
  
GPU加速:
  安装cuDF和CuPy以启用GPU加速:
  pip install cudf-cu11 cupy-cuda11x  # CUDA 11.x
  pip install cudf-cu12 cupy-cuda12x  # CUDA 12.x
        """
    )
    
    parser.add_argument('--all', action='store_true', default=True,
                       help='处理所有SUPMixerOutput_f*se*目录（默认）')
    parser.add_argument('--scale', nargs='+', type=float,
                       help='处理特定的scale factors（例如: --scale 1 2 4 或 --scale 1.5 2.0）')
    parser.add_argument('--cohesion', type=str, default=None,
                       help='只处理特定cohesion代码的目录（例如: --cohesion 000）')
    parser.add_argument('--base-path', default='build',
                       help='包含SUPMixerOutput目录的基础路径（默认: build）')
    parser.add_argument('--backup', action='store_true',
                       help='处理前创建原始文件的备份')
    parser.add_argument('--quiet', action='store_true',
                       help='减少输出详细程度')
    parser.add_argument('--cpu', action='store_true',
                       help='强制使用CPU模式（禁用GPU加速）')
    
    args = parser.parse_args()
    
    # 如果指定了--scale，则不使用--all
    if args.scale:
        args.all = False
    
    # 创建预处理器
    preprocessor = DEMEPreprocessor(
        base_path=args.base_path,
        backup=args.backup,
        verbose=not args.quiet,
        use_gpu=not args.cpu
    )
    
    # 运行预处理
    try:
        if args.all:
            preprocessor.process_all(cohesion_filter=args.cohesion)
        else:
            preprocessor.process_all(scale_factors=args.scale, cohesion_filter=args.cohesion)
    except KeyboardInterrupt:
        print("\n\n⚠️ 用户中断处理")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ 意外错误: {e}")
        import traceback
        if not args.quiet:
            traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()