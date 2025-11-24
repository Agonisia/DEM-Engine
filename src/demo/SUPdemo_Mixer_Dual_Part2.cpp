//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//
//	SPDX-License-Identifier: BSD-3-Clause

// =============================================================================
// SUP模型搅拌器仿真：Part2 - 读取大小粒子配置并运行搅拌仿真
// =============================================================================

#include <core/ApiVersion.h>
#include <core/utils/ThreadManager.h>
#include <DEM/API.h>
#include <DEM/HostSideHelpers.hpp>
#include <DEM/utils/Samplers.hpp>

#include <cstdio>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <iomanip>
#include <sstream>

using namespace deme;
using namespace std::filesystem;

int main() {
    // =========================================================================
    // 1. 仿真设置与SUP模型参数
    // =========================================================================
    
    // === SUP模型核心参数 ===
    const float my_scale_factor = 4.0f;                    // SUP缩放因子
    const float scale_force_index = 2.0f;                  // 力的缩放指数
    const float base_particle_mass = 0.0458f;              // 基准系统总质量：0.0458kg
    
    // === 大小粒子基准直径 ===
    const float base_small_diameter = 0.0003f;             // 小粒子基准直径：0.3mm
    const float base_large_diameter = 0.0005f;             // 大粒子基准直径：0.5mm（1.67倍）
    
    // === 粒子密度 ===
    const float particle_density = 1000.0f;                // 统一密度：1000 kg/m³
    const float particle_cohesion = 0.8f;                  // 颗粒间粘聚力
    
    // === 根据缩放因子自动计算的参数 ===
    const float small_diameter = base_small_diameter * my_scale_factor;
    const float small_radius = small_diameter / 2.0f;
    const float large_diameter = base_large_diameter * my_scale_factor;
    const float large_radius = large_diameter / 2.0f;
    
    // 小粒子质量（缩放后）
    const float particle_mass_small = particle_density * (4.0/3.0) * 3.14159265359 * 
                                      pow(small_radius, 3);
    // 大粒子质量（缩放后）
    const float particle_mass_large = particle_density * (4.0/3.0) * 3.14159265359 * 
                                      pow(large_radius, 3);
    
    // === 搅拌器参数 ===
    const float mixer_speed_rpm = 300.0f;                  // 搅拌器转速：300 RPM
    const float mixer_angular_velocity = mixer_speed_rpm * 2.0f * 3.14159265359f / 60.0f; // rad/s
    const float simulation_time = 3.0f;                    // 仿真时间：3秒
    
    // === 时间参数 ===
     const float step_size = (my_scale_factor == 1.0f) ? 5e-7f : 
                            (my_scale_factor == 2.0f) ? 1e-6f : 
                            2e-6f;                        // 时间步长
    const float time_per_frame = 
        (my_scale_factor == 1.0f) ? 5e-3f : 1e-3f;        // 每帧仿真时间
    
    // === 仿真域参数 ===
    const float domain_size = 0.1f;                        // 100mm × 100mm 水平尺寸
    const float wall_radius = 0.042f;                      // 圆柱墙体半径42mm
    
    // === Family ID定义（与Part1一致） ===
    const unsigned int FAMILY_MIXER = 1;
    const unsigned int FAMILY_PARTICLES_SMALL = 2;         // 小粒子
    const unsigned int FAMILY_PARTICLES_LARGE = 3;         // 大粒子
    
    // === 输入文件路径 ===
    const std::string input_dir = "/home/huyuze/DEM-Engine/DEM-Data/dualSize/";
    int scale_int = static_cast<int>(my_scale_factor);
    const std::string input_file = "SUPMixer_dualsize_f" + std::to_string(scale_int) + "se000.csv";

    // 输出SUP模型信息
    std::cout << "\n========== SUPdemo_Mixer_Part2 (大小粒子双层系统) ==========" << std::endl;
    std::cout << "缩放因子: " << my_scale_factor << std::endl;
    std::cout << "小粒子直径: " << small_diameter * 1000 << " mm" << std::endl;
    std::cout << "大粒子直径: " << large_diameter * 1000 << " mm" << std::endl;
    std::cout << "统一密度: " << particle_density << " kg/m³" << std::endl;
    std::cout << "小粒子质量: " << particle_mass_small * 1000 << " g" << std::endl;
    std::cout << "大粒子质量: " << particle_mass_large * 1000 << " g" << std::endl;
    std::cout << "搅拌器转速: " << mixer_speed_rpm << " RPM" << std::endl;
    std::cout << "仿真时间: " << simulation_time << " 秒" << std::endl;
    std::cout << "输入文件: " << input_dir + input_file << std::endl;
    std::cout << "===================================================================\n" << std::endl;

    // 创建求解器实例
    DEMSolver DEMSim;

    DEMSim.SetIntegrator(TIME_INTEGRATOR::CENTERED_DIFFERENCE);
    
    // 设置输出格式和内容
    DEMSim.SetVerbosity(INFO);
    DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV);
    DEMSim.SetOutputContent({"ABSV", "VEL", "ANG_VEL", "FAMILY"}); // 添加FAMILY以区分粒子类型
    DEMSim.SetContactOutputFormat(OUTPUT_FORMAT::CSV);
    DEMSim.SetContactOutputContent({"CNT_TYPE", "FORCE", "POINT"});

    // =========================================================================
    // 2. 材料属性设置
    // =========================================================================
    
    // 定义搅拌器材料
    auto mat_type_mixer = DEMSim.LoadMaterial({
        {"E", 1e8},         
        {"nu", 0.3},        
        {"CoR", 0.1},       
        {"mu", 0.3},        
        {"Crr", 0.0},       
        {"Cohesion", 0},
        {"scale_factor_l", my_scale_factor}
    });
    
    // 定义颗粒材料
    auto mat_type_particles = DEMSim.LoadMaterial({
        {"E", 1e8},         
        {"nu", 0.3},        
        {"CoR", 0.1},        
        {"mu", 0.3},
        {"Crr", 0.0},
        {"Cohesion", particle_cohesion},
        {"scale_factor_l", my_scale_factor}
    });
    
    // 设置材料相互作用属性
    DEMSim.SetMaterialPropertyPair("mu", mat_type_mixer, mat_type_particles, 0.3);
    DEMSim.SetMaterialPropertyPair("CoR", mat_type_mixer, mat_type_particles, 0.1);
    DEMSim.SetMaterialPropertyPair("Crr", mat_type_mixer, mat_type_particles, 0.0);
    DEMSim.SetMaterialPropertyPair("Cohesion", mat_type_mixer, mat_type_particles, 0.0);

    // 加载SUP接触力模型
    auto model_SUP = DEMSim.ReadContactForceModel("ForceModelSUP.cu");
    model_SUP->SetMustHaveMatProp({"E", "nu", "CoR", "mu", "Crr", "Cohesion", "scale_factor_l"});
    model_SUP->SetMustPairwiseMatProp({"CoR", "mu", "Crr", "Cohesion"});
    model_SUP->SetPerContactWildcards({"delta_time", "delta_tan_x", "delta_tan_y", "delta_tan_z"});

    // =========================================================================
    // 3. 仿真域设置
    // =========================================================================

    // 设置仿真域大小
    DEMSim.InstructBoxDomainDimension(
        {-domain_size/2, domain_size/2},            // X方向
        {-domain_size/2, domain_size/2},            // Y方向
        {0, 0.2f}                                   // Z方向：0-0.2m
    );

    // 设置边界条件
    DEMSim.InstructBoxDomainBoundingBC("top_open", mat_type_mixer);
    
    // 设置物理参数
    DEMSim.SetInitTimeStep(step_size);
    DEMSim.SetGravitationalAcceleration(make_float3(0, 0, -9.81));
    
    // 错误处理和安全限制
    DEMSim.SetErrorOutVelocity(20000.0);
    DEMSim.SetErrorOutAvgContacts(150);  
    
    // 接触检测和性能设置
    DEMSim.SetCDUpdateFreq(40);
    DEMSim.SetExpandSafetyAdder(2.0);

    // =========================================================================  
    // 4. 读取Part1输出的粒子配置  
    // =========================================================================  
    
    std::cout << "\n=== 读取粒子配置 ===" << std::endl;  
    
    // 读取粒子位置和方向  
    auto particle_xyz = DEMSim.ReadClumpXyzFromCsv(input_dir + input_file);  
    auto particle_quat = DEMSim.ReadClumpQuatFromCsv(input_dir + input_file);  
    
    // 准备加载粒子的向量  
    std::vector<float3> in_xyz;  
    std::vector<float4> in_quat;  
    std::vector<std::shared_ptr<DEMClumpTemplate>> in_types;  
    std::vector<unsigned int> in_families;  
    
    // 创建两种粒子模板（小粒子和大粒子）
    // 小粒子模板 - 使用SphereType
    auto small_sphere_template = DEMSim.LoadSphereType(  
        particle_mass_small,    // 小粒子质量  
        small_radius,  
        mat_type_particles  
    );
    
    // 大粒子模板 - 使用SphereType
    auto large_sphere_template = DEMSim.LoadSphereType(  
        particle_mass_large,    // 大粒子质量  
        large_radius,  
        mat_type_particles  
    );
    
    // 统计变量  
    int small_count = 0;  
    int large_count = 0;  
    
    // 处理读取的数据  
    for (const auto& [type_name, positions] : particle_xyz) {  
        auto& orientations = particle_quat[type_name];  
        
        // 根据type_name判断粒子类型
        // 假设Part1输出时，小粒子标记为"0000"，大粒子标记为"0001"
        // 需要根据实际输出调整
        bool is_small = (type_name == "0000");  // 小粒子
        bool is_large = (type_name == "0001");  // 大粒子
        
        for (size_t i = 0; i < positions.size(); i++) {  
            in_xyz.push_back(positions[i]);  
            in_quat.push_back(orientations[i]);  
            
            if (is_small) {  
                in_types.push_back(small_sphere_template);  
                in_families.push_back(FAMILY_PARTICLES_SMALL);  
                small_count++;  
            } else if (is_large) {  
                in_types.push_back(large_sphere_template);  
                in_families.push_back(FAMILY_PARTICLES_LARGE);  
                large_count++;  
            } else {  
                // 如果type_name不是预期的值，可能需要根据其他标准判断
                std::cerr << "警告: 未知的粒子类型 " << type_name << std::endl;  
            }  
        }  
        
        std::cout << "读取粒子类型 " << type_name << ": "   
                << positions.size() << " 个粒子" << std::endl;  
    }  
    
    std::cout << "  - 小粒子（" << small_diameter * 1000 << "mm）: " << small_count << " 个" << std::endl;  
    std::cout << "  - 大粒子（" << large_diameter * 1000 << "mm）: " << large_count << " 个" << std::endl;  
    std::cout << "  - 总计: " << in_xyz.size() << " 个粒子" << std::endl;  
    
    // 计算质量统计  
    float total_small_mass = small_count * particle_mass_small;  
    float total_large_mass = large_count * particle_mass_large;  
    float total_mass = total_small_mass + total_large_mass;  
    
    std::cout << "质量统计:" << std::endl;  
    std::cout << "  - 小粒子总质量: " << total_small_mass << " kg" << std::endl;  
    std::cout << "  - 大粒子总质量: " << total_large_mass << " kg" << std::endl;  
    std::cout << "  - 系统总质量: " << total_mass << " kg" << std::endl;  
    if (total_large_mass > 0) {  
        std::cout << "  - 质量比（小:大）: " << total_small_mass/total_large_mass << ":1" << std::endl;  
    }
    std::cout << "  - 粒子数比（小:大）: " << (float)small_count/large_count << ":1" << std::endl;
    std::cout << "  - 尺寸比（大/小）: " << large_diameter/small_diameter << std::endl;
    
    // 加载粒子到系统  
    DEMClumpBatch particle_batch(in_xyz.size());  
    particle_batch.SetTypes(in_types);  
    particle_batch.SetPos(in_xyz);  
    particle_batch.SetOriQ(in_quat);  
    particle_batch.SetFamilies(in_families);  // 设置Family信息  
    
    DEMSim.AddClumps(particle_batch);  
    
    std::cout << "成功加载 " << in_xyz.size() << " 个大小粒子到系统" << std::endl;

    // =========================================================================
    // 5. 创建搅拌器几何体
    // =========================================================================
    
    // 添加固定几何体 - 圆柱形边界墙
    auto walls = DEMSim.AddExternalObject();
    walls->AddCylinder(make_float3(0), make_float3(0, 0, 1), 
                      wall_radius, mat_type_mixer, 0);
    
    // 添加可动几何体 - 搅拌器网格
    auto mixer = DEMSim.AddWavefrontMeshObject(
        (GET_DATA_PATH() / "mesh/impeller.obj").string(), 
        mat_type_mixer
    );
    mixer->SetFamily(FAMILY_MIXER);
    
    std::cout << "搅拌器三角形网格数: " << mixer->GetNumTriangles() << std::endl;
    
    // 设置搅拌器的预定义角速度（绕z轴旋转）
    std::string omega_z_str = std::to_string(mixer_angular_velocity);
    DEMSim.SetFamilyPrescribedAngVel(FAMILY_MIXER, "0", "0", omega_z_str, true);  
    
    // 禁用搅拌器的自接触
    DEMSim.DisableContactBetweenFamilies(FAMILY_MIXER, FAMILY_MIXER);
    
    // 创建搅拌器的Tracker
    auto mixer_tracker = DEMSim.Track(mixer);

    // =========================================================================
    // 6. 创建输出目录和Inspector
    // =========================================================================
    
    // 创建输出目录（添加size_diff标识）
    std::ostringstream oss;
    oss << "SUPMixerOutput_SizeDiff_f" << (int)my_scale_factor 
        << "se" << std::setw(3) << std::setfill('0') << (int)(particle_cohesion * 100);
    std::string dir_name = oss.str();
    std::filesystem::path out_dir = std::filesystem::current_path() / dir_name;
    create_directory(out_dir);
    
    // 创建Inspector对象
    auto KE_inspector = DEMSim.CreateInspector("clump_kinetic_energy");
    auto max_v_inspector = DEMSim.CreateInspector("clump_max_absv");
    auto max_z_inspector = DEMSim.CreateInspector("clump_max_z");
    auto min_z_inspector = DEMSim.CreateInspector("clump_min_z");
    
    // 创建数据输出文件
    std::ofstream energy_file(out_dir / "kinetic_energy_history.csv");
    energy_file << "Time(s),Frame,TotalKE(J),MaxVelocity(m/s),"
                << "MaxZ(m),MinZ(m),AvgContacts,MixerTorque(Nm),"
                << "SmallParticles,LargeParticles" << std::endl;
    
    // 用于存储时序数据
    std::vector<float> time_history;
    std::vector<float> KE_history;
    std::vector<float> max_v_history;
    std::vector<float> mixer_torque_history;

    // =========================================================================
    // 7. 系统初始化
    // =========================================================================
    
    std::cout << "\n=== 初始化仿真系统 ===" << std::endl;
    DEMSim.Initialize();

    // =========================================================================
    // 8. 主仿真循环
    // =========================================================================
    
    std::cout << "\n=== 开始大小粒子混合仿真 ===" << std::endl;
    std::cout << "搅拌器角速度: " << mixer_angular_velocity << " rad/s" << std::endl;
    std::cout << "监测参数: 动能, 速度, 位置, 接触数, 搅拌器扭矩, 混合度" << std::endl;
    
    unsigned int frame_steps = (unsigned int)(time_per_frame / step_size);
    unsigned int curr_step = 0;
    unsigned int frame_count = 0;
    float current_time = 0.0f;
    
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    
    while (current_time < simulation_time) {
        // 执行一步仿真
        DEMSim.DoDynamics(step_size);
        
        // 每帧输出
        if (curr_step % frame_steps == 0) {
            std::cout << "\n--- 帧 " << frame_count << " (t = " << current_time << " s) ---" << std::endl;
            
            // 获取Inspector值
            float KE = KE_inspector->GetValue();
            float max_v = max_v_inspector->GetValue();
            float max_z = max_z_inspector->GetValue();
            float min_z = min_z_inspector->GetValue();
            float avg_contacts = DEMSim.GetAvgSphContacts();
            
            // 获取搅拌器扭矩
            float mixer_torque = 0.0f;
            std::vector<float3> forces, points;
            mixer_tracker->GetContactForcesForAll(points, forces);
            
            // 计算总扭矩（绕z轴）
            for (size_t i = 0; i < points.size(); i++) {
                float3 r = points[i]; // 相对于原点的位置向量
                float3 F = forces[i];
                // 扭矩的z分量: T_z = r_x * F_y - r_y * F_x
                mixer_torque += r.x * F.y - r.y * F.x;
            }
            
            // 存储历史数据
            time_history.push_back(current_time);
            KE_history.push_back(KE);
            max_v_history.push_back(max_v);
            mixer_torque_history.push_back(mixer_torque);
            
            // 写入CSV文件（实时写入）
            energy_file << std::fixed << std::setprecision(6)
                       << current_time << "," 
                       << frame_count << "," 
                       << KE << "," 
                       << max_v << "," 
                       << max_z << ","
                       << min_z << ","
                       << avg_contacts << ","
                       << mixer_torque << ","
                       << small_count << ","
                       << large_count << std::endl;
            energy_file.flush(); // 确保数据被写入
            
            // 在控制台显示
            std::cout << std::fixed << std::setprecision(4);
            std::cout << "系统动能: " << KE << " J" << std::endl;
            std::cout << "最大速度: " << max_v << " m/s" << std::endl;
            std::cout << "颗粒高度范围: [" << min_z << ", " << max_z << "] m" << std::endl;
            std::cout << "平均接触数: " << avg_contacts << std::endl;
            std::cout << "搅拌器扭矩: " << mixer_torque << " Nm" << std::endl;
            std::cout << "搅拌器接触点数: " << points.size() << std::endl;
            
            // 显示线程协作统计
            DEMSim.ShowThreadCollaborationStats();
            
            // 输出粒子数据文件
            char filename[200];
            sprintf(filename, "mixer_output_%04d.csv", frame_count);
            DEMSim.WriteSphereFile(out_dir / filename);

            // 输出接触力文件
            char contact_filename[200];
            sprintf(contact_filename, "mixer_contact_%04d.csv", frame_count);
            DEMSim.WriteContactFile(out_dir / contact_filename);
            
            // 输出mesh文件
            char meshfilename[200];
            sprintf(meshfilename, "mixer_mesh_%04d.vtk", frame_count);
            DEMSim.WriteMeshFile(out_dir / meshfilename);
            
            frame_count++;
        }
        
        curr_step++;
        current_time += step_size;
    }
    
    // 关闭数据文件
    energy_file.close();
    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_sec = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    
    // =========================================================================
    // 9. 后处理和统计
    // =========================================================================
   
    std::cout << "\n=== 大小粒子混合仿真完成 ===" << std::endl;
    std::cout << "总运行时间: " << time_sec.count() << " 秒" << std::endl;
    std::cout << "输出目录: " << out_dir << std::endl;
    std::cout << "数据文件: kinetic_energy_history.csv" << std::endl;
    
    // 计算并显示能量统计
    if (!KE_history.empty()) {
        float avg_KE = 0.0f;
        float max_KE = KE_history[0];
        float min_KE = KE_history[0];
        
        float avg_torque = 0.0f;
        float max_torque = std::abs(mixer_torque_history[0]);
        
        for (size_t i = 0; i < KE_history.size(); i++) {
            avg_KE += KE_history[i];
            max_KE = std::max(max_KE, KE_history[i]);
            min_KE = std::min(min_KE, KE_history[i]);
            
            avg_torque += std::abs(mixer_torque_history[i]);
            max_torque = std::max(max_torque, std::abs(mixer_torque_history[i]));
        }
        avg_KE /= KE_history.size();
        avg_torque /= mixer_torque_history.size();
        
        std::cout << "\n--- 动能统计 ---" << std::endl;
        std::cout << "平均动能: " << avg_KE << " J" << std::endl;
        std::cout << "最大动能: " << max_KE << " J" << std::endl;
        std::cout << "最小动能: " << min_KE << " J" << std::endl;
        std::cout << "最终动能: " << KE_history.back() << " J" << std::endl;
        
        std::cout << "\n--- 搅拌器扭矩统计 ---" << std::endl;
        std::cout << "平均扭矩(绝对值): " << avg_torque << " Nm" << std::endl;
        std::cout << "最大扭矩(绝对值): " << max_torque << " Nm" << std::endl;
        
        // 写入统计文件
        std::ofstream stats_file(out_dir / "simulation_statistics.txt");
        stats_file << "SUP Mixer Size Difference Simulation Statistics\n";
        stats_file << "================================================\n";
        stats_file << "Simulation Parameters:\n";
        stats_file << "  Simulation Time: " << simulation_time << " s\n";
        stats_file << "  Total Particles: " << in_xyz.size() << "\n";
        stats_file << "    - Small Particles (" << small_diameter * 1000 << " mm): " << small_count << "\n";
        stats_file << "    - Large Particles (" << large_diameter * 1000 << " mm): " << large_count << "\n";
        stats_file << "  Size Ratio (Large/Small): " << large_diameter/small_diameter << "\n";
        stats_file << "  Scale Factor: " << my_scale_factor << "\n";
        stats_file << "  Unified Density: " << particle_density << " kg/m³\n";
        stats_file << "  Small Particle Mass: " << particle_mass_small << " kg\n";
        stats_file << "  Large Particle Mass: " << particle_mass_large << " kg\n";
        stats_file << "  Total System Mass: " << total_mass << " kg\n";
        stats_file << "  Mass Ratio (Small:Large): " << total_small_mass/total_large_mass << ":1\n";
        stats_file << "  Number Ratio (Small:Large): " << (float)small_count/large_count << ":1\n";
        stats_file << "  Mixer Speed: " << mixer_speed_rpm << " RPM\n";
        stats_file << "  Time Step: " << step_size << " s\n";
        stats_file << "\nKinetic Energy Statistics:\n";
        stats_file << "  Average: " << avg_KE << " J\n";
        stats_file << "  Maximum: " << max_KE << " J\n";
        stats_file << "  Minimum: " << min_KE << " J\n";
        stats_file << "  Final: " << KE_history.back() << " J\n";
        stats_file << "\nMixer Torque Statistics:\n";
        stats_file << "  Average (abs): " << avg_torque << " Nm\n";
        stats_file << "  Maximum (abs): " << max_torque << " Nm\n";
        stats_file << "\nComputation Performance:\n";
        stats_file << "  Total Runtime: " << time_sec.count() << " s\n";
        stats_file << "  Real-time Factor: " << simulation_time / time_sec.count() << "\n";
        stats_file << "\nConfiguration Note:\n";
        stats_file << "  Small particles inserted first (bottom layer)\n";
        stats_file << "  Large particles inserted later (top layer)\n";
        stats_file.close();
        
        std::cout << "\n统计文件已保存至: simulation_statistics.txt" << std::endl;
    }
    
    // 显示详细统计信息
    std::cout << "\n--- 时间统计 ---" << std::endl;
    DEMSim.ShowTimingStats();
    DEMSim.ClearTimingStats();
    
    std::cout << "\n--- 内存统计 ---" << std::endl;
    DEMSim.ShowMemStats();
    
    std::cout << "\n--- 最终线程协作统计 ---" << std::endl;
    DEMSim.ShowThreadCollaborationStats();
    
    std::cout << "========================================" << std::endl;
    std::cout << "SUPdemo_Mixer_Part2 (Size Difference) 退出..." << std::endl;
    return 0;
}
