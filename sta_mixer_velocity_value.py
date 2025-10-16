#!/usr/bin/env python3
"""
DEME混合器数据预处理脚本
功能：检查并添加velocity列到mixer_output CSV文件
作者：DEM-Engine Analysis Tools
"""

import os
import sys
import glob
import shutil
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

class DEMEPreprocessor:
    """DEME数据预处理器类"""
    
    def __init__(self, base_path='build', backup=False, verbose=True):
        self.base_path = Path(base_path)
        self.backup = backup
        self.verbose = verbose
        self.processed_files = []
        self.skipped_files = []
        self.error_files = []
        
    def find_sup_directories(self):
        """查找所有SUPMixerOutput_f*目录"""
        pattern = self.base_path / "SUPMixerOutput_f*"
        directories = glob.glob(str(pattern))
        
        if not directories:
            print(f"⚠️ No SUPMixerOutput directories found in {self.base_path}")
            return []
        
        # 提取scale factor并排序
        dir_info = []
        for dir_path in directories:
            dir_name = os.path.basename(dir_path)
            if dir_name.startswith("SUPMixerOutput_f"):
                try:
                    scale = int(dir_name.split("_f")[1])
                    dir_info.append((scale, dir_path))
                except (IndexError, ValueError):
                    continue
        
        dir_info.sort(key=lambda x: x[0])
        return dir_info
    
    def check_csv_columns(self, file_path):
        """检查CSV文件的列"""
        try:
            # 只读取第一行来检查列名
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
    
    def add_velocity_column(self, file_path, output_path=None):
        """添加velocity列到CSV文件"""
        try:
            # 读取完整文件
            df = pd.read_csv(file_path)
            
            # 检查是否有必要的列
            if not all(col in df.columns for col in ['v_x', 'v_y', 'v_z']):
                raise ValueError("Missing velocity component columns")
            
            # 计算velocity
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
    
    def process_directory(self, scale_factor, dir_path):
        """处理单个SUPMixerOutput目录"""
        print(f"\n📁 Processing SUPMixerOutput_f{scale_factor}")
        print(f"   Path: {dir_path}")
        
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
        
        # 处理每个文件
        for i, file_path in enumerate(csv_files):
            file_name = os.path.basename(file_path)
            
            # 显示进度
            if i % 100 == 0:
                print(f"   Processing file {i+1}/{len(csv_files)}...")
            
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
            
            # 添加velocity列
            if self.add_velocity_column(file_path):
                processed_count += 1
                self.processed_files.append(file_path)
            else:
                error_count += 1
                self.error_files.append(file_path)
        
        # 打印统计
        print(f"\n   ✅ Summary for f{scale_factor}:")
        print(f"      - Processed: {processed_count} files (velocity column added)")
        print(f"      - Skipped: {skipped_count} files (already have velocity)")
        print(f"      - Errors: {error_count} files")
        
        if self.backup and processed_count > 0:
            backup_dir = Path(dir_path) / "backup_original"
            print(f"      - Backup saved to: {backup_dir}")
    
    def process_all(self, scale_factors=None):
        """处理所有或指定的scale factor目录"""
        print("="*60)
        print("DEME Data Preprocessor - Adding velocity column")
        print("="*60)
        print(f"Settings:")
        print(f"  - Base path: {self.base_path}")
        print(f"  - Backup: {'ENABLED' if self.backup else 'DISABLED'}")
        print(f"  - Verbose: {'ON' if self.verbose else 'OFF'}")
        
        # 查找所有目录
        directories = self.find_sup_directories()
        
        if not directories:
            return
        
        # 筛选要处理的目录
        if scale_factors:
            directories = [(s, d) for s, d in directories if s in scale_factors]
            if not directories:
                print(f"⚠️ No directories found for scale factors: {scale_factors}")
                return
        
        print(f"\n📋 Found {len(directories)} directories to process:")
        for scale, dir_path in directories:
            print(f"   - f{scale}: {dir_path}")
        
        # 处理每个目录
        for scale, dir_path in directories:
            self.process_directory(scale, dir_path)
        
        # 总体统计
        self.print_final_summary()
    
    def print_final_summary(self):
        """打印最终统计摘要"""
        print("\n" + "="*60)
        print("FINAL SUMMARY")
        print("="*60)
        print(f"✅ Total processed: {len(self.processed_files)} files")
        print(f"⏭️ Total skipped: {len(self.skipped_files)} files")
        print(f"❌ Total errors: {len(self.error_files)} files")
        
        if not self.backup and len(self.processed_files) > 0:
            print("\n⚠️  Note: Original files were modified directly (no backup)")
            print("   To create backups, use --backup flag")
        
        if self.error_files and self.verbose:
            print("\n⚠️ Files with errors:")
            for file in self.error_files[:10]:  # 只显示前10个
                print(f"   - {file}")
            if len(self.error_files) > 10:
                print(f"   ... and {len(self.error_files) - 10} more")

def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description='DEME Data Preprocessor - Add velocity column to CSV files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process all directories (no backup)
  python preprocess_deme_data.py --all
  
  # Process specific scale factors
  python preprocess_deme_data.py --scale 1 2 4
  
  # Process with backup
  python preprocess_deme_data.py --all --backup
  
  # Process with custom base path
  python preprocess_deme_data.py --all --base-path /path/to/data
  
  # Process with backup and quiet mode
  python preprocess_deme_data.py --all --backup --quiet
        """
    )
    
    parser.add_argument('--all', action='store_true',
                       help='Process all SUPMixerOutput_f* directories')
    parser.add_argument('--scale', nargs='+', type=int,
                       help='Process specific scale factors (e.g., --scale 1 2 4)')
    parser.add_argument('--base-path', default='build',
                       help='Base path containing SUPMixerOutput directories (default: build)')
    parser.add_argument('--backup', action='store_true',
                       help='Create backup of original files before processing')
    parser.add_argument('--quiet', action='store_true',
                       help='Reduce output verbosity')
    
    args = parser.parse_args()
    
    # 验证参数
    if not args.all and not args.scale:
        parser.print_help()
        print("\n❌ Error: Must specify --all or --scale")
        sys.exit(1)
    
    # 创建预处理器
    preprocessor = DEMEPreprocessor(
        base_path=args.base_path,
        backup=args.backup,
        verbose=not args.quiet
    )
    
    # 运行预处理
    try:
        if args.all:
            preprocessor.process_all()
        else:
            preprocessor.process_all(scale_factors=args.scale)
    except KeyboardInterrupt:
        print("\n\n⚠️ Processing interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ Unexpected error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()