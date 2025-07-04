//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//
//	SPDX-License-Identifier: BSD-3-Clause

// =============================================================================
// two particle overlap test using SUP contact model
// =============================================================================

#include <core/ApiVersion.h>
#include <core/utils/ThreadManager.h>
#include <DEM/API.h>
#include <DEM/HostSideHelpers.hpp>
#include <DEM/utils/Samplers.hpp>

#include <cstdio>
#include <chrono>
#include <filesystem>

using namespace deme;
using namespace std::filesystem;

int main() {
    // =========================================================================
    // 1. 仿真设置
    // =========================================================================
    
    // 全局参数
    const float my_scale_factor = 2.0f;                      // SUP缩放因子（2）
    const unsigned int FAMILY_PUSHER = 1;
    const unsigned int FAMILY_PARTICLES = 2;
    // 定义时间步长
    const float step_size = 10e-9; 
    // 定义每帧仿真时间
    const float time_per_frame = 10e-8;
    // 创建求解器实例
    DEMSolver DEMSim;
    
    // 设置输出格式和内容
    DEMSim.SetVerbosity(INFO);
    DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV);
    DEMSim.SetOutputContent({"VEL", "ANG_VEL", "ABS_ACC"}); 
    DEMSim.SetContactOutputFormat(OUTPUT_FORMAT::CSV);
    DEMSim.SetContactOutputContent({"POINT", "FORCE", "TORQUE"});
    
    // 求解器基本配置
    DEMSim.SetErrorOutAvgContacts(100);
    
    // 设置随机种子
    srand(52);

    // =========================================================================
    // 2. 材料属性
    // =========================================================================

    // 定义底部材料
    auto mat_type_walls = DEMSim.LoadMaterial({
        {"E", 1e7},         // 杨氏模量
        {"nu", 0.3},        // 泊松比
        {"CoR", 0.5},       // 恢复系数
        {"mu", 0.5},          // 滑动摩擦系数
        {"Crr", 0.03},      // 滚动摩擦系数
        {"Cohesion", 0.01}        // 粘聚力
    });

    // 定义颗粒材料
    auto mat_type_particles = DEMSim.LoadMaterial({
        {"E", 1e7},         
        {"nu", 0.3},        
        {"CoR", 0.5},        
        {"mu", 0.5},
        {"Crr", 0.03},
        {"Cohesion", 0.01}
    });

    
    // 设置材料相互作用属性 - 颗粒与墙壁
    DEMSim.SetMaterialPropertyPair("CoR", mat_type_walls, mat_type_particles, 0.5); 
    DEMSim.SetMaterialPropertyPair("mu", mat_type_walls, mat_type_particles, 0.5);
    DEMSim.SetMaterialPropertyPair("Crr", mat_type_walls, mat_type_particles, 0.03);
    DEMSim.SetMaterialPropertyPair("Cohesion", mat_type_walls, mat_type_particles, 0.01);
 
    // 加载SUP接触力模型
    auto model_SUP = DEMSim.ReadContactForceModel("ForceModelSUP.cu");
    model_SUP->SetMustHaveMatProp({"E", "nu", "CoR", "mu", "Crr", "Cohesion"});
    model_SUP->SetMustPairwiseMatProp({"CoR", "mu", "Crr", "Cohesion"});
    model_SUP->SetPerContactWildcards({"delta_time", "delta_tan_x", "delta_tan_y", "delta_tan_z", "scale_factor_l"});


    // =========================================================================
    // 3. 仿真域设置
    // =========================================================================

    // 颗粒基本参数
    float particle_diameter = 0.001f * my_scale_factor;                 // 1mm直径
    float particle_radius = particle_diameter / 2.0f;  // 半径
    float particle_density = 1000.0;                   // 颗粒密度：1000 kg/m³ 
    float particle_mass = particle_density * (4.0/3.0) * 3.14159265359 * pow(particle_radius, 3); // 颗粒质量

    // 定义仿真域
    float domain_size = 0.035;       // 35mm × 35mm 水平尺寸
    float plate_bottom = 0.f;       // Z坐标：底板位置

    // 设置仿真域大小
    DEMSim.InstructBoxDomainDimension(
        {-domain_size/2, domain_size/2},            // X方向
        {-domain_size/2, domain_size/2},            // Y方向
        {(plate_bottom), (plate_bottom + 0.1)}      // Z方向
    );

    // 设置边界条件 - 顶部开放
    DEMSim.InstructBoxDomainBoundingBC("top_open", mat_type_walls);

    // 设置物理参数
    DEMSim.SetInitTimeStep(step_size);
    DEMSim.SetGravitationalAcceleration(make_float3(0, 0, -9.81));

    // 设置最大速度（防止数值不稳定）
    DEMSim.SetMaxVelocity(25.0f);

    // 设置错误输出速度（用于处理初始沉降）
    DEMSim.SetErrorOutVelocity(50.0f);

    // =========================================================================
    // 4. 创建两个重叠的粒子
    // =========================================================================

    // 粒子模板
    std::shared_ptr<DEMClumpTemplate> clump_type;
    
    // 创建单球团簇模板
    std::vector<float> radii = {particle_radius};
    std::vector<float3> relPos = {make_float3(0, 0, 0)};
    
    float clump_mass = particle_mass;
    float3 MOI = make_float3(1, 1, 1) * 0.4 * clump_mass * pow(particle_radius, 2);
    
    clump_type = DEMSim.LoadClumpType(clump_mass, MOI, radii, relPos, mat_type_particles);

    // 设置两个粒子的位置 - 使它们重叠
    float overlap_distance = 0.7 * particle_diameter;  // 两个粒子中心距离为0.7倍直径（30%重叠）
    float overlap_amount = particle_diameter - overlap_distance;  // 实际重叠量
    
    float3 pos1 = make_float3(0, 0, 0.01f);  // 第一个粒子位置
    float3 pos2 = make_float3(overlap_distance, 0, 0.01f);  // 第二个粒子位置
    
    // 打印重叠信息
    std::cout << "\n=== 粒子重叠信息 ===" << std::endl;
    std::cout << "粒子直径: " << particle_diameter * 1000 << " mm" << std::endl;
    std::cout << "粒子半径: " << particle_radius * 1000 << " mm" << std::endl;
    std::cout << "两粒子中心距离: " << overlap_distance * 1000 << " mm" << std::endl;
    std::cout << "重叠量 (overlap): " << overlap_amount * 1000 << " mm" << std::endl;
    std::cout << "重叠百分比: " << (overlap_amount / particle_diameter) * 100 << "%" << std::endl;
    std::cout << "==================" << std::endl;

    // 创建DEMClumpBatch
    DEMClumpBatch particles(2);
    
    std::vector<std::shared_ptr<DEMClumpTemplate>> templates = {clump_type, clump_type};
    std::vector<float3> positions = {pos1, pos2};
    std::vector<float3> velocities = {make_float3(0, 0, 0), make_float3(0, 0, 0)};
    std::vector<unsigned int> families = {FAMILY_PARTICLES, FAMILY_PARTICLES};
    
    particles.SetTypes(templates);
    particles.SetPos(positions);
    particles.SetVel(velocities);
    particles.SetFamilies(families);
    
    DEMSim.SetFamilyFixed(FAMILY_PARTICLES);

    // =========================================================================
    // 5. SIMULATION INITIALIZATION 
    // =========================================================================
    
    // Initialize the simulation system
    DEMSim.Initialize();
    
    // Add particles
    auto clump_tracker = DEMSim.AddClumps(particles);
    DEMSim.UpdateClumps();
    
    // Set initial conditions for SUP model
    DEMSim.SetFamilyContactWildcardValueBoth(FAMILY_PARTICLES, "scale_factor_l", my_scale_factor);

    // =========================================================================
    // 6. OUTPUT SETUP 
    // =========================================================================
    
    // Create output directory
    path out_dir = current_path();
    out_dir += "/SUP_TwoParticleOverlap_Output";

    // Check if directory exists and clear it, otherwise create it
    if (exists(out_dir)) {
        std::cout << "Output directory exists, clearing contents..." << std::endl;
        for (const auto& entry : directory_iterator(out_dir)) {
            remove_all(entry.path());
        }
    } else {
        create_directory(out_dir);
        std::cout << "Created output directory: " << out_dir << std::endl;
    }

    // =========================================================================
    // 7. 仿真执行
    // =========================================================================
    
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    
    // 主仿真循环
    for (int frame = 0; frame < 2; frame++) {
        
        // 每10帧输出一次文件
        if (frame % 1 == 0) {
            // 生成输出文件名
            char outputfile[200], unifiedfile[200];
            sprintf(unifiedfile, "%s/SUPdemo_unified_%04d.csv", out_dir.c_str(), frame);                
            sprintf(outputfile, "%s/SUPdemo_output_%04d.vtk", out_dir.c_str(), frame);
            
            // 写入输出文件
            DEMSim.WriteSphereAndContactFile(std::string(outputfile));
            DEMSim.WriteUnifiedFile(std::string(unifiedfile));
            
            // 进度报告
            std::cout << "帧: " << frame << " - 输出文件已保存" << std::endl;
            
            // 在第一帧时详细输出粒子信息
            if (frame == 1) {
                std::cout << "\n=== 第0帧粒子状态 ===" << std::endl;
                std::cout << "粒子1位置: (" << pos1.x*1000 << ", " << pos1.y*1000 << ", " << pos1.z*1000 << ") mm" << std::endl;
                std::cout << "粒子2位置: (" << pos2.x*1000 << ", " << pos2.y*1000 << ", " << pos2.z*1000 << ") mm" << std::endl;
                std::cout << "预期接触力作用在X方向" << std::endl;
                std::cout << "====================" << std::endl;
            }
        }
        
        // 推进仿真
        DEMSim.DoDynamics(time_per_frame);
        DEMSim.ShowThreadCollaborationStats();
    }

    // =========================================================================
    // 8. 后处理
    // =========================================================================

    // 性能统计
    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_sec = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    std::cout << "\n仿真完成, 耗时: " << time_sec.count() << " 秒" << std::endl;
    
    DEMSim.ShowTimingStats();
    
    // 清理资源
    DEMSim.ClearTimingStats();

    std::cout << "----------------------------------------" << std::endl;
    DEMSim.ShowMemStats();
    std::cout << "----------------------------------------" << std::endl;

    std::cout << "SUP Two Particle Overlap Test 退出..." << std::endl;
    return 0;
}