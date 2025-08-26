#!/usr/bin/env python3
# run_compression_batch.py
# 批量运行压缩测试程序，自动修改参数并收集结果

import os
import re
import subprocess
import shutil
import time
from pathlib import Path
import pandas as pd

class CompressionTestRunner:
    def __init__(self):
        # 获取脚本所在目录（即DEM-Engine目录）
        self.project_root = Path(__file__).parent.resolve()
        
        # 需要根据实际情况修改源文件名
        self.source_file = self.project_root / "src/demo/SUPdemo_BulkCompress.cpp"  # 请修改为实际的源文件名
        self.build_dir = self.project_root / "build"
        self.executable = self.build_dir / "bin/SUPdemo_BulkCompress"  # 请修改为实际的可执行文件名

        # 参数组合
        self.cohesion_values = [0, 0.05, 0.1, 0.2]
        self.scale_factors = [1.0, 2.0, 4.0]  # 可以添加更多缩放因子
        
        # 创建备份
        self.backup_file = str(self.source_file) + ".backup"
        
        # 保存初始工作目录
        self.initial_cwd = Path.cwd()
        
        # 结果收集
        self.results_summary = []

    def format_cohesion_string(self, cohesion):
        """将粘附力值转换为三位数字符串格式"""
        if cohesion == 0:
            return "000"
        else:
            # 0.05 -> 005, 0.1 -> 010, 0.2 -> 020
            return f"{int(cohesion * 100):03d}"
    
    def format_scale_string(self, scale_factor):
        """将缩放因子转换为字符串"""
        return str(int(scale_factor))
    
    def get_input_csv_path(self, scale_factor, cohesion):
        """根据参数生成输入CSV文件路径"""
        factor_str = self.format_scale_string(scale_factor)
        coe_str = self.format_cohesion_string(cohesion)
        
        # 构建输入文件路径
        input_path = f"/root/DEM-Engine/results/f{factor_str}c{coe_str}/SUPdemo_f{factor_str}c{coe_str}_0150.csv"
        return input_path
        
    def modify_cpp_file(self, cohesion, scale_factor):
        """修改C++源文件中的参数"""
        
        # 读取文件内容
        with open(self.source_file, 'r', encoding='utf-8') as f:
            content = f.read()
        
        # 1. 修改scale_factor
        scale_pattern = r'(const float my_scale_factor = )([\d.]+f)'
        content = re.sub(scale_pattern, rf'\g<1>{scale_factor}f', content)
        
        # 2. 修改输入CSV路径
        input_csv_path = self.get_input_csv_path(scale_factor, cohesion)
        input_pattern = r'(const std::string input_csv = ")([^"]+)(")'
        content = re.sub(input_pattern, rf'\g<1>{input_csv_path}\g<3>', content)
        
        # 3. 修改颗粒材料的Cohesion值
        lines = content.split('\n')
        in_particle_material = False
        modified_lines = []
        
        for i, line in enumerate(lines):
            # 检测是否进入颗粒材料定义部分
            if '// 定义颗粒材料' in line:
                in_particle_material = True
                modified_lines.append(line)
                continue
            
            # 如果在颗粒材料定义内，且遇到Cohesion行
            if in_particle_material and '"Cohesion",' in line:
                # 修改这一行
                new_line = re.sub(r'({"Cohesion",\s*)([\d.]+)(\})', rf'\g<1>{cohesion}\g<3>', line)
                modified_lines.append(new_line)
                # 检查是否到达材料定义结束
                if '});' in lines[i+1] if i+1 < len(lines) else False:
                    in_particle_material = False
                continue
            
            # 检测材料定义结束
            if in_particle_material and '});' in line:
                in_particle_material = False
            
            modified_lines.append(line)
        
        content = '\n'.join(modified_lines)
        
        # 4. 修改输出目录名称
        coe_str = self.format_cohesion_string(cohesion)
        factor_str = self.format_scale_string(scale_factor)
        dir_name = f"CompressionOutput_f{factor_str}c{coe_str}"
        
        # 替换输出目录
        out_dir_pattern = r'(out_dir \+= "/)([^"]+)(")'
        content = re.sub(out_dir_pattern, rf'\g<1>{dir_name}\g<3>', content)
        
        return content, dir_name, input_csv_path
    
    def check_input_file(self, input_csv_path):
        """检查输入文件是否存在"""
        if not Path(input_csv_path).exists():
            print(f"⚠️  警告：输入文件不存在: {input_csv_path}")
            print("   请确保已经运行了相应的SUP仿真生成该文件")
            return False
        return True
    
    def extract_results(self, log_file, dir_name):
        """从日志文件中提取关键结果"""
        try:
            with open(log_file, 'r') as f:
                log_content = f.read()
            
            # 提取关键信息
            results = {
                'directory': dir_name,
                'max_stress': None,
                'max_strain': None,
                'final_strain': None,
                'contact_time': None,
                'total_particles': None,
                'filtered_particles': None
            }
            
            # 使用正则表达式提取数据
            patterns = {
                'max_stress': r'最大应力:\s*([\d.]+e?[+-]?\d*)\s*Pa',
                'max_strain': r'在应变\s*([\d.]+)\s*%\s*时',
                'final_strain': r'最终应变:\s*([\d.]+)\s*%',
                'contact_time': r'接触时间:\s*([\d.]+)\s*秒',
                'total_particles': r'原始粒子数:\s*(\d+)',
                'filtered_particles': r'过滤后粒子数:\s*(\d+)'
            }
            
            for key, pattern in patterns.items():
                match = re.search(pattern, log_content)
                if match:
                    results[key] = float(match.group(1)) if key != 'total_particles' and key != 'filtered_particles' else int(match.group(1))
            
            return results
            
        except Exception as e:
            print(f"提取结果时出错: {str(e)}")
            return None
    
    def run_single_test(self, cohesion, scale_factor):
        """运行单个测试"""
        
        print(f"\n{'='*70}")
        print(f"运行测试: Cohesion={cohesion}, Scale Factor={scale_factor}")
        print(f"{'='*70}")
        
        try:
            # 修改源文件
            modified_content, dir_name, input_csv_path = self.modify_cpp_file(cohesion, scale_factor)
            
            # 检查输入文件
            if not self.check_input_file(input_csv_path):
                return False
            
            # 写入修改后的内容
            with open(self.source_file, 'w', encoding='utf-8') as f:
                f.write(modified_content)

            coe_str = self.format_cohesion_string(cohesion)
            factor_str = self.format_scale_string(scale_factor)
            
            print(f"✓ 已修改源文件")
            print(f"  - Cohesion: {cohesion}")
            print(f"  - Scale Factor: {scale_factor}")
            print(f"  - 输入文件: {input_csv_path}")
            print(f"  - 输出目录: {dir_name}")
            
            # 编译
            print("\n编译中...")
            compile_start = time.time()
            
            # 切换到build目录并编译
            os.chdir(self.build_dir)
            result = subprocess.run("ninja", shell=True, capture_output=True, text=True)
            
            if result.returncode != 0:
                print(f"✗ 编译失败!")
                print(f"错误信息: {result.stderr}")
                return False
            
            compile_time = time.time() - compile_start
            print(f"✓ 编译成功 (耗时: {compile_time:.1f}秒)")
            
            # 运行程序
            print("\n运行仿真...")
            run_start = time.time()
            
            # 创建日志目录
            log_dir = self.project_root / "compression_logs"
            log_dir.mkdir(exist_ok=True)
            log_file = log_dir / f"{dir_name}.log"
            
            # 运行程序并保存输出
            with open(log_file, 'w') as f:
                # 修改为实际的可执行文件名
                result = subprocess.run("./bin/SUPdemo_BulkCompress",  # 请确保这个名称正确
                                     shell=True, 
                                     stdout=f, 
                                     stderr=subprocess.STDOUT)
            
            run_time = time.time() - run_start
            
            if result.returncode == 0:
                print(f"✓ 仿真完成 (运行时间: {run_time:.1f}秒)")
                print(f"  日志文件: {log_file}")
                
                # 提取结果
                results = self.extract_results(log_file, dir_name)
                if results:
                    results['cohesion'] = cohesion
                    results['scale_factor'] = scale_factor
                    results['run_time'] = run_time
                    self.results_summary.append(results)
                    
                    # 打印关键结果
                    print(f"\n关键结果:")
                    if results['max_stress']:
                        print(f"  - 最大应力: {results['max_stress']:.2e} Pa")
                    if results['final_strain']:
                        print(f"  - 最终应变: {results['final_strain']:.1f}%")
                    if results['filtered_particles']:
                        print(f"  - 压缩粒子数: {results['filtered_particles']}")
                
                # 复制应力应变曲线
                self.copy_stress_strain_curve(dir_name, coe_str, factor_str)
                
                return True
            else:
                print(f"✗ 仿真运行失败!")
                print(f"  查看日志: {log_file}")
                return False
                
        except Exception as e:
            print(f"✗ 发生错误: {str(e)}")
            import traceback
            traceback.print_exc()
            return False
    
    def copy_stress_strain_curve(self, dir_name, coe_str, factor_str):
        """复制应力应变曲线到结果文件夹"""
        # 源文件路径
        source_dir = self.build_dir / dir_name
        source_file = source_dir / "stress_strain_curve.csv"
        
        if not source_file.exists():
            print(f"警告：找不到应力应变曲线文件: {source_file}")
            return
        
        # 创建目标文件夹
        target_folder = self.project_root / "compression_results" / f"f{factor_str}c{coe_str}"
        target_folder.mkdir(parents=True, exist_ok=True)
        
        # 复制文件
        shutil.copy2(source_file, target_folder / "stress_strain_curve.csv")
        
        # 也复制最终状态文件
        final_state = source_dir / "compression_final_state.csv"
        if final_state.exists():
            shutil.copy2(final_state, target_folder / "compression_final_state.csv")
        
        print(f"✓ 已复制结果文件到: {target_folder}")
    
    def save_summary(self):
        """保存结果汇总"""
        if not self.results_summary:
            return
        
        # 创建DataFrame
        df = pd.DataFrame(self.results_summary)
        
        # 保存为CSV
        summary_file = self.project_root / "compression_results" / "summary.csv"
        summary_file.parent.mkdir(exist_ok=True)
        df.to_csv(summary_file, index=False)
        
        print(f"\n结果汇总已保存到: {summary_file}")
        
        # 打印汇总表格
        print("\n结果汇总表:")
        print(df.to_string(index=False))
    
    def run_all_tests(self):
        """运行所有测试组合"""

        # 对scale_factors进行降序排序
        self.scale_factors = sorted(self.scale_factors, reverse=True)
        
        # 检查源文件是否存在
        if not self.source_file.exists():
            print(f"错误：找不到源文件 {self.source_file}")
            print(f"请确保源文件路径正确！")
            return
        
        # 创建备份
        print(f"备份原始文件到: {self.backup_file}")
        shutil.copy(self.source_file, self.backup_file)
        
        # 统计信息
        total_tests = len(self.cohesion_values) * len(self.scale_factors)
        successful = 0
        failed = 0
        skipped = 0
        
        print(f"\n准备运行 {total_tests} 个测试")
        print(f"Cohesion 值: {self.cohesion_values}")
        print(f"Scale Factor 值: {self.scale_factors}")
        
        start_time = time.time()
        
        try:
            # 运行所有组合
            for i, scale in enumerate(self.scale_factors):
                for j, cohesion in enumerate(self.cohesion_values):
                    test_num = i * len(self.cohesion_values) + j + 1
                    print(f"\n[{test_num}/{total_tests}]", end="")
                    
                    # 检查输入文件是否存在
                    input_path = self.get_input_csv_path(scale, cohesion)
                    if not Path(input_path).exists():
                        print(f" 跳过 (输入文件不存在): Scale={scale}, Cohesion={cohesion}")
                        skipped += 1
                        continue
                    
                    if self.run_single_test(cohesion, scale):
                        successful += 1
                    else:
                        failed += 1
                    
                    # 从备份恢复，为下一次测试准备
                    shutil.copy(self.backup_file, self.source_file)
        
        finally:
            # 清理：恢复原始文件
            print("\n\n恢复原始文件...")
            shutil.copy(self.backup_file, self.source_file)
            os.remove(self.backup_file)
            
            # 切回初始目录
            os.chdir(self.initial_cwd)
        
        # 保存汇总结果
        self.save_summary()
        
        # 总结
        total_time = time.time() - start_time
        print(f"\n{'='*70}")
        print(f"所有测试完成!")
        print(f"总耗时: {total_time:.1f}秒 ({total_time/60:.1f}分钟)")
        print(f"成功: {successful}")
        print(f"失败: {failed}")
        print(f"跳过: {skipped}")
        print(f"{'='*70}")
        
        # 显示输出文件位置
        print(f"\n输出文件位置:")
        print(f"仿真日志: {self.project_root}/compression_logs/")
        print(f"仿真数据: {self.build_dir}/CompressionOutput_*")
        print(f"压缩结果: {self.project_root}/compression_results/")
        print(f"结果汇总: {self.project_root}/compression_results/summary.csv")

def main():
    """主函数"""
    import argparse
    
    parser = argparse.ArgumentParser(description='批量运行压缩测试仿真')
    parser.add_argument('--dry-run', action='store_true',
                        help='只显示将要运行的参数组合，不实际运行')
    parser.add_argument('--cohesion', nargs='+', type=float,
                        help='自定义Cohesion值列表（默认: 0 0.05 0.1 0.2）')
    parser.add_argument('--scale', nargs='+', type=float,
                        help='自定义Scale Factor值列表（默认: 2.0 4.0）')
    parser.add_argument('--source', type=str,
                        help='源文件路径（相对于项目根目录）')
    parser.add_argument('--executable', type=str,
                        help='可执行文件名（不含路径）')
    
    args = parser.parse_args()
    
    runner = CompressionTestRunner()
    
    # 如果指定了源文件路径
    if args.source:
        runner.source_file = runner.project_root / args.source
    
    # 如果指定了可执行文件名
    if args.executable:
        runner.executable = runner.build_dir / "bin" / args.executable
    
    # 如果指定了自定义参数
    if args.cohesion:
        runner.cohesion_values = args.cohesion
    if args.scale:
        runner.scale_factors = args.scale
    
    if args.dry_run:
        print("预览模式 - 将运行以下参数组合：")
        print(f"源文件: {runner.source_file}")
        print(f"可执行文件: {runner.executable}")
        print("\n参数组合:")
        
        for scale in runner.scale_factors:
            for cohesion in runner.cohesion_values:
                coe_str = runner.format_cohesion_string(cohesion)
                factor_str = runner.format_scale_string(scale)
                dir_name = f"CompressionOutput_f{factor_str}c{coe_str}"
                input_csv = runner.get_input_csv_path(scale, cohesion)
                exists = "✓" if Path(input_csv).exists() else "✗"
                print(f"  Scale={scale}, Cohesion={cohesion} -> {dir_name}")
                print(f"    输入文件: {input_csv} [{exists}]")
                
        print(f"\n总计: {len(runner.cohesion_values) * len(runner.scale_factors)} 个测试")
    else:
        runner.run_all_tests()

if __name__ == "__main__":
    main()