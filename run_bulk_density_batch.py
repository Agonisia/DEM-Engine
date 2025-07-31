#!/usr/bin/env python3
# run_bulk_density_batch.py

import os
import re
import subprocess
import shutil
import time
from pathlib import Path

class BulkDensityRunner:
    def __init__(self):
        # 获取脚本所在目录（即DEM-Engine目录）
        self.project_root = Path(__file__).parent.resolve()
        self.source_file = self.project_root / "src/demo/SUPdemo_BulkDensity.cpp"
        self.build_dir = self.project_root / "build"
        self.executable = self.build_dir / "bin/SUPdemo_BulkDensity"
        
        # 参数组合
        self.cohesion_values = [0, 0.05, 0.1, 0.2]
        self.scale_factors = [4.0]
        
        # 创建备份
        self.backup_file = str(self.source_file) + ".backup"
        
        # 保存初始工作目录
        self.initial_cwd = Path.cwd()

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
        
    def modify_cpp_file(self, cohesion, scale_factor):
        """精确修改颗粒材料的Cohesion值和其他参数"""
        
        # 读取文件内容
        with open(self.source_file, 'r', encoding='utf-8') as f:
            content = f.read()
        
        # 1. 精确修改颗粒材料的Cohesion值
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
                # 修改这一行 - 使用 \\g<1> 语法避免歧义
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
        
        # 2. 修改scale_factor - 使用原始字符串避免问题
        scale_pattern = r'(const float my_scale_factor = )([\d.]+f)'
        content = re.sub(scale_pattern, rf'\g<1>{scale_factor}f', content)
        
        # 3. 修改输出目录
        # 生成目录名
        coe_str = self.format_cohesion_string(cohesion)
        factor_str = self.format_scale_string(scale_factor)
        dir_name = f"SUP_BulkDensity_coe{coe_str}_factor{factor_str}"
        
        # 替换输出目录 - 同样使用 \g<> 语法
        out_dir_pattern = r'(out_dir \+= "/)(SUP_BulkDensity_[^"]+)(")'
        content = re.sub(out_dir_pattern, rf'\g<1>{dir_name}\g<3>', content)
        
        return content, dir_name
    
    def run_single_test(self, cohesion, scale_factor):
        """运行单个测试"""
        
        print(f"\n{'='*60}")
        print(f"运行测试: Cohesion={cohesion}, Scale Factor={scale_factor}")
        print(f"{'='*60}")
        
        try:
            # 修改源文件
            modified_content, dir_name = self.modify_cpp_file(cohesion, scale_factor)
            
            # 写入修改后的内容
            with open(self.source_file, 'w', encoding='utf-8') as f:
                f.write(modified_content)

            coe_str = self.format_cohesion_string(cohesion)
            factor_str = self.format_scale_string(scale_factor)
            
            print(f"✓ 已修改源文件")
            print(f"  - Cohesion: {cohesion}")
            print(f"  - Scale Factor: {scale_factor}")
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
            
            # 创建日志目录（在项目根目录下）
            log_dir = self.project_root / "simulation_logs"
            log_dir.mkdir(exist_ok=True)
            log_file = log_dir / f"{dir_name}.log"
            
            # 运行程序并保存输出
            with open(log_file, 'w') as f:
                result = subprocess.run("./bin/SUPdemo_BulkDensity", 
                                     shell=True, 
                                     stdout=f, 
                                     stderr=subprocess.STDOUT)
            
            run_time = time.time() - run_start
            
            if result.returncode == 0:
                print(f"✓ 仿真完成 (运行时间: {run_time:.1f}秒)")
                print(f"  日志文件: {log_file}")

                # 创建结果文件夹并复制最后一帧
                self.copy_last_frame(dir_name, coe_str, factor_str)

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
    
    def copy_last_frame(self, dir_name, coe_str, factor_str):
        """复制最后一帧到专门的文件夹"""
        # 输出目录路径
        output_dir = Path("/root/autodl-tmp") / dir_name
        
        if not output_dir.exists():
            print(f"警告：输出目录不存在: {output_dir}")
            return
            
        # 查找所有输出文件
        pattern = f"SUPdemo_f{factor_str}c{coe_str}_*.csv"
        files = list(output_dir.glob(pattern))
        
        if not files:
            print(f"警告：未找到匹配的输出文件: {pattern}")
            return
            
        # 按文件名排序，找到最后一帧
        files.sort()
        last_frame = files[-1]
        
        # 创建目标文件夹
        target_folder = self.project_root / "results" / f"f{factor_str}c{coe_str}"
        target_folder.mkdir(parents=True, exist_ok=True)
        
        # 复制文件
        shutil.copy2(last_frame, target_folder / last_frame.name)
        print(f"✓ 已复制最后一帧到: {target_folder}")
    
    def run_all_tests(self):
        """运行所有测试组合"""
        
        # 检查源文件是否存在
        if not self.source_file.exists():
            print(f"错误：找不到源文件 {self.source_file}")
            print(f"当前目录: {Path.cwd()}")
            print(f"项目根目录: {self.project_root}")
            return
        
        # 创建备份
        print(f"备份原始文件到: {self.backup_file}")
        shutil.copy(self.source_file, self.backup_file)
        
        # 统计信息
        total_tests = len(self.cohesion_values) * len(self.scale_factors)
        successful = 0
        failed = 0
        
        print(f"\n准备运行 {total_tests} 个测试")
        print(f"Cohesion 值: {self.cohesion_values}")
        print(f"Scale Factor 值: {self.scale_factors}")
        
        start_time = time.time()
        
        try:
            # 运行所有组合
            for i, cohesion in enumerate(self.cohesion_values):
                for j, scale in enumerate(self.scale_factors):
                    test_num = i * len(self.scale_factors) + j + 1
                    print(f"\n[{test_num}/{total_tests}]", end="")
                    
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
        
        # 总结
        total_time = time.time() - start_time
        print(f"\n{'='*60}")
        print(f"所有测试完成!")
        print(f"总耗时: {total_time:.1f}秒 ({total_time/60:.1f}分钟)")
        print(f"成功: {successful}")
        print(f"失败: {failed}")
        print(f"{'='*60}")
        
        # 显示输出文件位置
        print(f"\n输出文件位置:")
        print(f"仿真日志: {self.project_root}/simulation_logs/")
        print(f"仿真数据: {self.project_root}/build/SUP_BulkDensity_*")
        print(f"最后帧结果: {self.project_root}/results/f*c*/")

def main():
    """主函数"""
    import argparse
    
    parser = argparse.ArgumentParser(description='批量运行BulkDensity仿真')
    parser.add_argument('--dry-run', action='store_true',
                        help='只显示将要运行的参数组合，不实际运行')
    parser.add_argument('--cohesion', nargs='+', type=float,
                        help='自定义Cohesion值列表（默认: 0 0.05 0.1 0.2 0.4）')
    parser.add_argument('--scale', nargs='+', type=float,
                        help='自定义Scale Factor值列表（默认: 2.0 4.0）')
    
    args = parser.parse_args()
    
    runner = BulkDensityRunner()
    
    # 如果指定了自定义参数
    if args.cohesion:
        runner.cohesion_values = args.cohesion
    if args.scale:
        runner.scale_factors = args.scale
    
    if args.dry_run:
        print("预览模式 - 将运行以下参数组合：")
        for cohesion in runner.cohesion_values:
            for scale in runner.scale_factors:
                if cohesion == 0:
                    coe_str = "0"
                else:
                    coe_str = f"{cohesion}".replace("0.", "")
                factor_str = str(int(scale))
                dir_name = f"SUP_BulkDensity_coe{coe_str}_factor{factor_str}"
                print(f"  Cohesion={cohesion}, Scale={scale} -> {dir_name}")
        print(f"\n总计: {len(runner.cohesion_values) * len(runner.scale_factors)} 个测试")
    else:
        runner.run_all_tests()

if __name__ == "__main__":
    main()