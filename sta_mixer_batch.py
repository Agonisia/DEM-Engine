#!/usr/bin/env python3
# run_mixer_batch.py
# 批量运行SUP混拌器测试程序，编译串行执行，仿真可选并行

import os
import re
import subprocess
import shutil
import time
from pathlib import Path
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

class MixerTestRunner:
    def __init__(self):
        # 获取脚本所在目录（即DEM-Engine目录）
        self.project_root = Path(__file__).parent.resolve()
        
        # 源文件和构建设置
        self.source_file = self.project_root / "src/demo/SUPdemo_Mixer_Part2.cpp"  
        self.build_dir = self.project_root / "build"
        self.executable = self.build_dir / "bin/SUPdemo_Mixer_Part2"
        
        # 参数组合
        self.cohesion_values = [0, 0.1]  # 粘附力值
        self.scale_factors = [2.0, 4.0]  # 缩放因子
        
        # 备份文件
        self.backup_file = str(self.source_file) + ".backup"
        
        # 保存初始工作目录
        self.initial_cwd = Path.cwd()
        
        # 结果收集
        self.results_summary = []
        
        # 并行设置
        self.max_workers = multiprocessing.cpu_count() // 2  # 使用一半的CPU核心

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
    
    def get_input_csv_path(self, scale_factor):
        """根据scale_factor生成输入CSV文件路径
        注意: 输入文件名称格式固定为 SUPMixer_f{scale}se000.csv
        """
        factor_str = self.format_scale_string(scale_factor)
        
        # 构建输入文件路径 - 固定使用se000
        input_path = f"/home/huyuze/DEM-Data/SUPMixer_f{factor_str}se000.csv"
        return input_path
    
    def modify_and_compile(self, cohesion, scale_factor):
        """修改源文件并编译（必须串行执行）"""
        
        print(f"\n=== 编译: Scale={scale_factor}, Cohesion={cohesion} ===")
        
        # 读取备份文件内容
        with open(self.backup_file, 'r', encoding='utf-8') as f:
            content = f.read()
        
        # 1. 修改scale_factor
        scale_pattern = r'(const float my_scale_factor = )([\d.]+)f'
        content = re.sub(scale_pattern, rf'\g<1>{scale_factor}f', content)
        
        # 2. 修改particle_cohesion
        cohesion_pattern = r'(const float particle_cohesion = )([\d.]+)f'
        content = re.sub(cohesion_pattern, rf'\g<1>{cohesion}f', content)
        
        # 注意：C++代码中会自动根据scale_factor生成正确的输入文件名
        # int scale_int = static_cast<int>(my_scale_factor);
        # const std::string input_file = "SUPMixer_f" + std::to_string(scale_int) + "se000.csv";
        
        # 写入修改后的内容
        with open(self.source_file, 'w', encoding='utf-8') as f:
            f.write(content)
        
        print(f"  ✓ 源文件已修改")
        print(f"    - Scale Factor: {scale_factor}")
        print(f"    - Cohesion: {cohesion}")
        
        # 编译
        print(f"  编译中...")
        compile_start = time.time()
        
        os.chdir(self.build_dir)
        result = subprocess.run("ninja", shell=True, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"  ✗ 编译失败!")
            print(f"  错误信息: {result.stderr}")
            os.chdir(self.initial_cwd)
            return False
        
        compile_time = time.time() - compile_start
        print(f"  ✓ 编译成功 (耗时: {compile_time:.1f}秒)")
        
        # 等待文件系统同步
        print(f"  等待15秒以确保编译完全完成...")
        time.sleep(15)
        
        os.chdir(self.initial_cwd)
        return True
    
    def run_simulation(self, cohesion, scale_factor):
        """运行仿真（可以并行）"""
        
        coe_str = self.format_cohesion_string(cohesion)
        factor_str = self.format_scale_string(scale_factor)
        dir_name = f"SUPMixerOutput_f{factor_str}se{coe_str}"
        
        print(f"  运行仿真: {dir_name}")
        run_start = time.time()
        
        # 创建日志目录
        log_dir = self.project_root / "mixer_logs"
        log_dir.mkdir(exist_ok=True)
        log_file = log_dir / f"{dir_name}.log"
        
        # 运行程序
        os.chdir(self.build_dir)
        with open(log_file, 'w') as f:
            result = subprocess.run("./bin/SUPdemo_Mixer_Part2",
                                 shell=True, 
                                 stdout=f, 
                                 stderr=subprocess.STDOUT)
        os.chdir(self.initial_cwd)
        
        run_time = time.time() - run_start
        
        if result.returncode == 0:
            print(f"    ✓ 仿真完成 (运行时间: {run_time:.1f}秒)")
            
            # 提取结果
            results = self.extract_results(log_file, dir_name)
            results.update({
                'status': 'success',
                'cohesion': cohesion,
                'scale_factor': scale_factor,
                'run_time': run_time
            })
            
            # 复制结果文件
            # self.copy_result_files(dir_name, coe_str, factor_str)
            
            return results
        else:
            print(f"    ✗ 仿真运行失败! 查看日志: {log_file}")
            return {
                'status': 'failed',
                'cohesion': cohesion,
                'scale_factor': scale_factor,
                'reason': '运行失败'
            }
    
    def extract_results(self, log_file, dir_name):
        """从日志文件中提取关键结果"""
        try:
            with open(log_file, 'r') as f:
                log_content = f.read()
            
            # 提取关键信息
            results = {
                'directory': dir_name,
                'avg_kinetic_energy': None,
                'max_kinetic_energy': None,
                'final_kinetic_energy': None,
                'avg_torque': None,
                'max_torque': None,
                'particle_count': None
            }
            
            # 使用正则表达式提取数据
            patterns = {
                'avg_kinetic_energy': r'平均动能:\s*([\d.]+)\s*J',
                'max_kinetic_energy': r'最大动能:\s*([\d.]+)\s*J',
                'final_kinetic_energy': r'最终动能:\s*([\d.]+)\s*J',
                'avg_torque': r'平均扭矩\(绝对值\):\s*([\d.]+)\s*Nm',
                'max_torque': r'最大扭矩\(绝对值\):\s*([\d.]+)\s*Nm',
                'particle_count': r'成功加载\s*(\d+)\s*个粒子'
            }
            
            for key, pattern in patterns.items():
                match = re.search(pattern, log_content)
                if match:
                    if key == 'particle_count':
                        results[key] = int(match.group(1))
                    else:
                        results[key] = float(match.group(1))
            
            # 检查是否达到稳态
            if '系统接近稳态' in log_content:
                results['steady_state'] = True
            elif '系统尚未达到稳态' in log_content:
                results['steady_state'] = False
            else:
                results['steady_state'] = None
            
            return results
            
        except Exception as e:
            print(f"  提取结果时出错: {str(e)}")
            return {}
    
    # def copy_result_files(self, dir_name, coe_str, factor_str):
    #     """复制关键结果文件到结果文件夹"""
    #     # 源文件路径
    #     source_dir = self.build_dir / dir_name
        
    #     if not source_dir.exists():
    #         print(f"  警告：找不到输出目录: {source_dir}")
    #         return
        
    #     # 创建目标文件夹
    #     folder_name = f"f{factor_str}se{coe_str}"
    #     target_folder = self.project_root / "mixer_results" / folder_name
    #     target_folder.mkdir(parents=True, exist_ok=True)
        
    #     # 复制关键文件
    #     files_to_copy = [
    #         "kinetic_energy_history.csv",
    #         "simulation_statistics.txt"
    #     ]
        
    #     for filename in files_to_copy:
    #         source_file = source_dir / filename
    #         if source_file.exists():
    #             shutil.copy2(source_file, target_folder / filename)
        
    #     print(f"    ✓ 已复制结果文件到: {target_folder}")
    
    def run_all_tests(self):
        """运行所有测试组合（编译串行，运行可并行）"""
        
        # 检查源文件是否存在
        if not self.source_file.exists():
            print(f"错误：找不到源文件 {self.source_file}")
            return
        
        # 创建备份
        print(f"备份原始文件到: {self.backup_file}")
        shutil.copy(self.source_file, self.backup_file)
        
        # 统计信息
        total_tests = len(self.cohesion_values) * len(self.scale_factors)
        
        print(f"\n准备运行 {total_tests} 个测试")
        print(f"Cohesion 值: {self.cohesion_values}")
        print(f"Scale Factor 值: {self.scale_factors}")
        print(f"编译方式: 串行（每次编译后等待15秒）")
        print(f"运行方式: scale=1串行，其他并行")
        
        start_time = time.time()
        
        try:
            # 按scale_factor分组处理
            for scale in self.scale_factors:
                print(f"\n{'='*70}")
                print(f"处理 Scale Factor = {scale} 的测试组")
                print(f"{'='*70}")
                
                # 检查输入文件
                input_file = self.get_input_csv_path(scale)
                if not Path(input_file).exists():
                    print(f"⚠️  警告：输入文件不存在: {input_file}")
                    print(f"   跳过所有 scale={scale} 的测试")
                    continue
                
                print(f"✓ 输入文件存在: {input_file}")
                
                if scale == 1.0:
                    # Scale=1时，编译和运行都串行
                    print("执行模式: 串行（编译+运行）")
                    for cohesion in self.cohesion_values:
                        # 编译
                        if self.modify_and_compile(cohesion, scale):
                            # 运行
                            result = self.run_simulation(cohesion, scale)
                            self.results_summary.append(result)
                        else:
                            self.results_summary.append({
                                'status': 'failed',
                                'cohesion': cohesion,
                                'scale_factor': scale,
                                'reason': '编译失败'
                            })
                else:
                    # Scale!=1时，先串行编译所有配置，然后并行运行
                    print("执行模式: 串行编译 → 并行运行")
                    
                    # 第一步：串行编译所有cohesion配置
                    compiled_configs = []
                    for cohesion in self.cohesion_values:
                        if self.modify_and_compile(cohesion, scale):
                            # 编译成功后，复制可执行文件到特定名称
                            exec_name = f"SUPdemo_Mixer_s{int(scale)}c{self.format_cohesion_string(cohesion)}"
                            exec_path = self.build_dir / "bin" / exec_name
                            shutil.copy2(self.executable, exec_path)
                            compiled_configs.append((cohesion, scale, exec_name))
                            print(f"    ✓ 已保存可执行文件: {exec_name}")
                        else:
                            self.results_summary.append({
                                'status': 'failed',
                                'cohesion': cohesion,
                                'scale_factor': scale,
                                'reason': '编译失败'
                            })
                    
                    # 第二步：并行运行所有成功编译的配置
                    if compiled_configs:
                        print(f"\n并行运行 {len(compiled_configs)} 个仿真...")
                        with ProcessPoolExecutor(max_workers=min(self.max_workers, len(compiled_configs))) as executor:
                            futures = []
                            for cohesion, scale, exec_name in compiled_configs:
                                future = executor.submit(self.run_simulation_with_exec, 
                                                       cohesion, scale, exec_name)
                                futures.append(future)
                            
                            # 收集结果
                            for future in as_completed(futures):
                                result = future.result()
                                self.results_summary.append(result)
                                
                                # 显示进度
                                completed = sum(1 for f in futures if f.done())
                                print(f"  进度: {completed}/{len(futures)} 完成")
        
        finally:
            # 清理：恢复原始文件
            print("\n\n恢复原始文件...")
            shutil.copy(self.backup_file, self.source_file)
            os.remove(self.backup_file)
            
            # 清理临时可执行文件
            for file in (self.build_dir / "bin").glob("SUPdemo_Mixer_s*c*"):
                file.unlink(missing_ok=True)
            
            # 切回初始目录
            os.chdir(self.initial_cwd)
        
        # 保存汇总结果
        self.save_summary()
        
        # 统计结果
        successful = sum(1 for r in self.results_summary if r.get('status') == 'success')
        failed = sum(1 for r in self.results_summary if r.get('status') == 'failed')
        skipped = sum(1 for r in self.results_summary if r.get('status') == 'skipped')
        
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
        print(f"仿真日志: {self.project_root}/mixer_logs/")
        print(f"仿真数据: {self.build_dir}/SUPMixerOutput_*")
        print(f"结果文件: {self.project_root}/mixer_results/")
        print(f"结果汇总: {self.project_root}/mixer_results/summary.csv")
    
    def run_simulation_with_exec(self, cohesion, scale_factor, exec_name):
        """使用特定的可执行文件运行仿真（用于并行执行）"""
        
        coe_str = self.format_cohesion_string(cohesion)
        factor_str = self.format_scale_string(scale_factor)
        dir_name = f"SUPMixerOutput_f{factor_str}se{coe_str}"
        
        print(f"  运行仿真: {dir_name} (使用 {exec_name})")
        run_start = time.time()
        
        # 创建日志目录
        log_dir = self.project_root / "mixer_logs"
        log_dir.mkdir(exist_ok=True)
        log_file = log_dir / f"{dir_name}.log"
        
        # 运行特定的可执行文件
        os.chdir(self.build_dir)
        with open(log_file, 'w') as f:
            result = subprocess.run(f"./bin/{exec_name}",
                                 shell=True, 
                                 stdout=f, 
                                 stderr=subprocess.STDOUT)
        os.chdir(self.initial_cwd)
        
        run_time = time.time() - run_start
        
        if result.returncode == 0:
            print(f"    ✓ 仿真完成: {dir_name} (运行时间: {run_time:.1f}秒)")
            
            # 提取结果
            results = self.extract_results(log_file, dir_name)
            results.update({
                'status': 'success',
                'cohesion': cohesion,
                'scale_factor': scale_factor,
                'run_time': run_time
            })
            
            # 复制结果文件
            # self.copy_result_files(dir_name, coe_str, factor_str)
            
            return results
        else:
            print(f"    ✗ 仿真运行失败: {dir_name}")
            return {
                'status': 'failed',
                'cohesion': cohesion,
                'scale_factor': scale_factor,
                'reason': '运行失败'
            }
    
    def save_summary(self):
        """保存结果汇总"""
        if not self.results_summary:
            return
        
        # 只保存成功的结果
        successful_results = [r for r in self.results_summary if r.get('status') == 'success']
        
        if not successful_results:
            print("\n没有成功的测试结果可以保存")
            return
        
        # 创建DataFrame
        df = pd.DataFrame(successful_results)
        
        # 移除status列
        if 'status' in df.columns:
            df = df.drop('status', axis=1)
        
        # 保存为CSV
        summary_file = self.project_root / "mixer_results" / "summary.csv"
        summary_file.parent.mkdir(exist_ok=True)
        df.to_csv(summary_file, index=False)
        
        print(f"\n结果汇总已保存到: {summary_file}")
        
        # 打印汇总表格
        print("\n成功测试结果汇总:")
        # 选择要显示的关键列
        display_columns = ['scale_factor', 'cohesion', 'avg_kinetic_energy', 
                          'avg_torque', 'steady_state', 'run_time']
        display_columns = [col for col in display_columns if col in df.columns]
        print(df[display_columns].to_string(index=False))

def main():
    """主函数"""
    import argparse
    
    parser = argparse.ArgumentParser(description='批量运行SUP混拌器仿真')
    parser.add_argument('--dry-run', action='store_true',
                        help='只显示将要运行的参数组合，不实际运行')
    parser.add_argument('--cohesion', nargs='+', type=float,
                        help='自定义Cohesion值列表（默认: 0 0.05 0.1 0.2）')
    parser.add_argument('--scale', nargs='+', type=float,
                        help='自定义Scale Factor值列表（默认: 1.0 2.0 4.0）')
    parser.add_argument('--workers', type=int,
                        help='并行运行时使用的最大进程数')
    parser.add_argument('--source', type=str,
                        help='源文件路径（相对于项目根目录）')
    parser.add_argument('--executable', type=str,
                        help='可执行文件名（不含路径）')
    
    args = parser.parse_args()
    
    runner = MixerTestRunner()
    
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
    if args.workers:
        runner.max_workers = args.workers
    
    if args.dry_run:
        print("预览模式 - 将运行以下参数组合：")
        print(f"源文件: {runner.source_file}")
        print(f"可执行文件: {runner.executable}")
        print(f"最大并行进程数: {runner.max_workers}")
        print("\n参数组合:")
        
        for scale in runner.scale_factors:
            input_csv = runner.get_input_csv_path(scale)
            exists = "✓" if Path(input_csv).exists() else "✗"
            
            print(f"\nScale Factor = {scale}")
            print(f"  输入文件: {input_csv} [{exists}]")
            print(f"  执行模式: {'串行' if scale == 1.0 else '并行'}")
            
            for cohesion in runner.cohesion_values:
                coe_str = runner.format_cohesion_string(cohesion)
                factor_str = runner.format_scale_string(scale)
                dir_name = f"SUPMixerOutput_f{factor_str}se{coe_str}"
                print(f"    Cohesion={cohesion} -> {dir_name}")
                
        print(f"\n总计: {len(runner.cohesion_values) * len(runner.scale_factors)} 个测试")
    else:
        runner.run_all_tests()

if __name__ == "__main__":
    main()