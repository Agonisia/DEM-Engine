//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//
//	SPDX-License-Identifier: BSD-3-Clause

// =============================================================================
// SUP模型搅拌器仿真：Part2 - 读取粒子配置并运行搅拌仿真
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

using namespace deme;
using namespace std::filesystem;

int main() {
    // =========================================================================
    // 1. 仿真设置与SUP模型参数
    // =========================================================================
    
    // === SUP模型核心参数（与Part1保持一致） ===
    const float my_scale_factor = 4.0f;                    // SUP缩放因子
    const float scale_force_index = 2.0f;                  // 力的缩放指数
    const int base_particle_count = 1600000;               // 基准粒子数量
    const float base_particle_diameter = 0.0005f;          // 基准粒子直径：0.5mm
    const float particle_density = 1000.0f;                // 颗粒密度：1000 kg/m³
    const float particle_cohesion = 0.0f;                  // 颗粒间粘聚力
    
    // === 根据缩放因子自动计算的参数 ===
    const float particle_diameter = base_particle_diameter * my_scale_factor;
    const float particle_radius = particle_diameter / 2.0f;
    const float particle_mass = particle_density * (4.0/3.0) * 3.14159265359 * 
                               pow(particle_radius, 3);
    
    // === 搅拌器参数 ===
    const float mixer_speed_rpm = 60.0f;                   // 搅拌器转速：60 RPM
    const float mixer_angular_velocity = mixer_speed_rpm * 2.0f * 3.14159265359f / 60.0f; // rad/s
    const float simulation_time = 10.0f;                   // 仿真时间：10秒
    
    // === 时间参数 ===
    const float step_size = 1e-6f;                         // 时间步长
    const float time_per_frame = 1e-2f;                    // 每帧仿真时间
    
    // === 仿真域参数 ===
    const float domain_size = 0.1f;                        // 100mm × 100mm 水平尺寸
    const float wall_radius = 0.042f;                      // 圆柱墙体半径42mm
    
    // === Family ID定义 ===
    const unsigned int FAMILY_MIXER = 1;
    const unsigned int FAMILY_PARTICLES = 2;
    
    // === 输入文件路径 ===
    const std::string input_dir = "/home/peize/research/DEME_data/250911/";
    const std::string input_file = "SUPMixer_f4c000.csv";  // 根据实际文件名调整
    
    // 输出SUP模型信息
    std::cout << "\n========== SUPdemo_Mixer_Part2 ==========" << std::endl;
    std::cout << "缩放因子: " << my_scale_factor << std::endl;
    std::cout << "粒子直径: " << particle_diameter * 1000 << " mm" << std::endl;
    std::cout << "搅拌器转速: " << mixer_speed_rpm << " RPM" << std::endl;
    std::cout << "仿真时间: " << simulation_time << " 秒" << std::endl;
    std::cout << "输入文件: " << input_dir + input_file << std::endl;
    std::cout << "========================================\n" << std::endl;

    // 创建求解器实例
    DEMSolver DEMSim;
    
    // 设置输出格式和内容
    DEMSim.SetVerbosity(INFO);
    DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV);
    DEMSim.SetOutputContent({"VEL", "ANG_VEL"}); 
    DEMSim.SetContactOutputFormat(OUTPUT_FORMAT::CSV);
    DEMSim.SetContactOutputContent({"CNT_TYPE", "FORCE", "TORQUE"});

    // =========================================================================
    // 2. 材料属性设置
    // =========================================================================
    
    // 定义搅拌器材料
    auto mat_type_mixer = DEMSim.LoadMaterial({
        {"E", 1e7},         
        {"nu", 0.3},        
        {"CoR", 0.1},       
        {"mu", 0.0},        
        {"Crr", 0.0},       
        {"Cohesion", 0},
        {"scale_factor_l", my_scale_factor}
    });
    
    // 定义颗粒材料
    auto mat_type_particles = DEMSim.LoadMaterial({
        {"E", 1e7},         
        {"nu", 0.3},        
        {"CoR", 0.1},        
        {"mu", 0.3},
        {"Crr", 0.0},
        {"Cohesion", particle_cohesion},
        {"scale_factor_l", my_scale_factor}
    });
    
    // 设置材料相互作用属性
    DEMSim.SetMaterialPropertyPair("mu", mat_type_mixer, mat_type_particles, 0.1);
    DEMSim.SetMaterialPropertyPair("CoR", mat_type_mixer, mat_type_particles, 0.0);
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
    DEMSim.SetErrorOutAvgContacts(50);
    
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
    
    // 创建粒子模板
    auto clump_template = DEMSim.LoadSphereType(
        particle_mass,
        particle_radius,
        mat_type_particles
    );
    
    // 处理读取的数据（假设只有一种粒子类型）
    for (const auto& [type_name, positions] : particle_xyz) {
        auto& orientations = particle_quat[type_name];
        
        for (size_t i = 0; i < positions.size(); i++) {
            in_xyz.push_back(positions[i]);
            in_quat.push_back(orientations[i]);
            in_types.push_back(clump_template);
        }
        
        std::cout << "读取粒子类型 " << type_name << ": " << positions.size() << " 个粒子" << std::endl;
    }
    
    // 加载粒子到系统
    DEMClumpBatch particle_batch(in_xyz.size());
    particle_batch.SetTypes(in_types);
    particle_batch.SetPos(in_xyz);
    particle_batch.SetOriQ(in_quat);
    
    particle_batch.SetFamily(FAMILY_PARTICLES);
    DEMSim.AddClumps(particle_batch);
    
    std::cout << "成功加载 " << in_xyz.size() << " 个粒子到系统" << std::endl;

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

    // =========================================================================
    // 6. 创建输出目录
    // =========================================================================
    
    // 创建输出目录
    path out_dir = current_path() / "SUPMixerOutput_Part2";
    create_directory(out_dir);

    // =========================================================================
    // 7. 系统初始化
    // =========================================================================
    
    std::cout << "\n=== 初始化仿真系统 ===" << std::endl;
    DEMSim.Initialize();

    // =========================================================================
    // 8. 主仿真循环
    // =========================================================================
    
    std::cout << "\n=== 开始搅拌仿真 ===" << std::endl;
    std::cout << "搅拌器角速度: " << mixer_angular_velocity << " rad/s" << std::endl;
    
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
            std::cout << "\n--- 帧 " << frame_count << " ---" << std::endl;
            
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
    
    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_sec = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    
    // =========================================================================
    // 9. 后处理和统计
    // =========================================================================
    
    // 最终状态输出
    DEMSim.WriteClumpFile(out_dir / "final_state.csv");
    DEMSim.WriteContactFile(out_dir / "final_contacts.csv");
    
    std::cout << "\n=== 仿真完成 ===" << std::endl;
    std::cout << "总运行时间: " << time_sec.count() << " 秒" << std::endl;
    std::cout << "输出目录: " << out_dir << std::endl;
    
    // 显示统计信息
    std::cout << "\n--- 时间统计 ---" << std::endl;
    DEMSim.ShowTimingStats();
    DEMSim.ClearTimingStats();
    
    std::cout << "\n--- 内存统计 ---" << std::endl;
    DEMSim.ShowMemStats();
    
    std::cout << "\n--- 最终线程协作统计 ---" << std::endl;
    DEMSim.ShowThreadCollaborationStats();
    
    std::cout << "========================================" << std::endl;
    std::cout << "SUPdemo_Mixer_Part2 退出..." << std::endl;
    return 0;
}
