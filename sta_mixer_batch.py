from sta_mixer_force_torque import analyze_mixer_contacts
import os

# 批处理多个仿真目录
base_dir = '/home/peize/research/DEM-Engine/build'
simulation_dirs = ['SUPMixerOutput_Part2']

all_results = {}
for sim_dir in simulation_dirs:
    full_path = os.path.join(base_dir, sim_dir)
    output = f'results_{sim_dir}.csv'
    
    print(f"\nProcessing {sim_dir}...")
    df = analyze_mixer_contacts(
        data_directory=full_path,
        output_file=output,
        generate_plot=False,  # 批处理时可能不需要每个都绘图
        verbose=False  # 减少输出
    )
    all_results[sim_dir] = df

# 比较不同仿真的结果
for sim_name, df in all_results.items():
    avg_torque = df['mesh_torque_magnitude'].mean()
    print(f"{sim_name}: Average torque = {avg_torque:.6e} N·m")