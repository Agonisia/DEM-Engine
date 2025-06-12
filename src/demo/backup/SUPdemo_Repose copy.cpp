//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//
//	SPDX-License-Identifier: BSD-3-Clause

// =============================================================================
// 基于"An expression for the angle of repose of dry
// cohesive granular materials on Earth and in
// planetary environments" 设计的休止角实验
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
    const float my_scale_factor = 1.0f;                      // SUP缩放因子（暂时设为1）
    const unsigned int FAMILY_PLATE = 1;
    const unsigned int FAMILY_PARTICLES = 2;

    // 创建求解器实例
    DEMSolver DEMSim;
    
    // 设置输出格式和内容
    DEMSim.SetVerbosity(INFO);
    DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV);
    DEMSim.SetOutputContent({"VEL", "ANG_VEL"}); 
    DEMSim.SetContactOutputFormat(OUTPUT_FORMAT::CSV);
    DEMSim.SetContactOutputContent({"POINT", "FORCE", "TORQUE"});
    
    // 求解器基本配置
    DEMSim.SetErrorOutAvgContacts(1000);
    
    // 设置随机种子
    srand(52);

    // =========================================================================
    // 2. 材料属性
    // =========================================================================

    // 定义底部材料
    auto mat_type_walls = DEMSim.LoadMaterial({
        {"E", 1e9},        // 杨氏模量：63 GPa
        {"nu", 0.24},         // 泊松比：0.24
        {"CoR", 0.1},         // 恢复系数
        {"mu", 0.5},          // 滑动摩擦系数：0.5
        {"Crr", 0.05},        // 滚动摩擦系数：0.05
        {"Cohesion", 0}       // 与墙面无粘聚力
    });

    // 定义颗粒材料（玻璃珠）
    auto mat_type_particles = DEMSim.LoadMaterial({
        {"E", 1e9},        // 杨氏模量：63 GPa
        {"nu", 0.24},         // 泊松比：0.24
        {"CoR", 0.1},         // 恢复系数：0.1
        {"mu", 0.5},          // 滑动摩擦系数：0.5
        {"Crr", 0.05},        // 滚动摩擦系数：0.05
        {"Cohesion", 0.05}    // 表面能密度：0.05 J/m²
    });
    
    // 设置材料相互作用属性
    DEMSim.SetMaterialPropertyPair("CoR", mat_type_walls, mat_type_particles, 0.1);
    DEMSim.SetMaterialPropertyPair("Crr", mat_type_walls, mat_type_particles, 0.05);
    DEMSim.SetMaterialPropertyPair("mu", mat_type_walls, mat_type_particles, 0.5);
    DEMSim.SetMaterialPropertyPair("Cohesion", mat_type_walls, mat_type_particles, 0.0);

    // 加载SUP接触力模型
    auto model_SUP = DEMSim.ReadContactForceModel("ForceModelSUP.cu");
    model_SUP->SetMustHaveMatProp({"E", "nu", "CoR", "mu", "Crr", "Cohesion"});
    model_SUP->SetMustPairwiseMatProp({"CoR", "mu", "Crr", "Cohesion"});
    model_SUP->SetPerContactWildcards({"delta_time", "delta_tan_x", "delta_tan_y", "delta_tan_z", "scale_factor_l"});


    // =========================================================================
    // 3. 仿真域设置
    // =========================================================================

    // 颗粒基本参数
    float particle_diameter = 0.001f;                  // 1mm 直径颗粒
    float particle_radius = particle_diameter / 2.0f;  // 半径
    float particle_density = 2500.0;                   // 颗粒密度：2500 kg/m³
    float particle_mass = particle_density * (4.0/3.0) * 3.14159265359 * pow(particle_radius, 3); // 颗粒质量

    // 定义底板和仿真域参数
    float plate_size = 40 * particle_diameter;        // 40d = 40mm 底板尺寸（根据论文）
    float plate_thickness = particle_diameter;         // 底板厚度
    float domain_width = 150 * particle_diameter;     // 150mm 仿真域宽度（足够大以允许溢出）
    float domain_height = 150 * particle_diameter;    // 150mm 高度（容纳所有下落颗粒）
    float domain_bottom = -2 * particle_diameter;     // 底部稍低，避免穿透问题

    // 定义时间步长
    float step_size = 1e-6;  // 稍微增大时间步长，提高效率

    // 设置仿真域大小
    DEMSim.InstructBoxDomainDimension(
        {-domain_width/2, domain_width/2},     // X方向
        {-domain_width/2, domain_width/2},     // Y方向
        {domain_bottom, domain_height}         // Z方向
    );

    // 设置边界条件 - 全开放边界（允许颗粒从任意边界离开）
    DEMSim.InstructBoxDomainBoundingBC("top_open", mat_type_walls);

    // 设置物理参数
    DEMSim.SetInitTimeStep(step_size);
    DEMSim.SetGravitationalAcceleration(make_float3(0, 0, -9.81));

    // 设置最大速度（防止数值不稳定）
    DEMSim.SetMaxVelocity(10.0f);  // 降低到更合理的值

    // 设置错误输出速度（用于处理初始沉降）
    DEMSim.SetErrorOutVelocity(50.0f);  // 也降低这个值

    // =========================================================================
    // 4. GEOMETRY CREATION
    // =========================================================================

    // 创建底板（使用大球模拟平板）
    auto plate_template = DEMSim.LoadClumpType(
        1e10,  // 非常大的质量（相当于固定）
        make_float3(1e10, 1e10, 1e10),  // 非常大的惯性矩
        {plate_size/2},  // 半径
        {make_float3(0, 0, -plate_size/2)},  // 位置在底部
        mat_type_walls
    );
    auto plate = DEMSim.AddClumps({plate_template}, {make_float3(0, 0, 0)});
    plate->SetFamily(FAMILY_PLATE);
    DEMSim.SetFamilyFixed(FAMILY_PLATE);  // 固定底板

    // =========================================================================
    // 5. PARTICLE GENERATION
    // =========================================================================

    // --- 5.1 Particle/Clump Templates ---

    // Template generation parameters
    int num_template = 6;                              // Total number of random clump templates
    int min_sphere = 1;                                // Minimum number of spheres per clump
    int max_sphere = 1;                                // Maximum number of spheres per clump
    float min_rad = particle_radius * 0.8;             // 稍小于标准半径
    float max_rad = particle_radius * 0.8;             // 稍大于标准半径
    float min_relpos = -particle_radius * 0.5;         // 相对位置范围
    float max_relpos = particle_radius * 0.5;

    // Create array to store clump templates
    std::vector<std::shared_ptr<DEMClumpTemplate>> clump_types;

    // Generate random clump templates
    for (int i = 0; i < num_template; i++) {
        // 可以使用单球或多球团簇
        int num_sphere = min_sphere + rand() % (max_sphere - min_sphere + 1);
        
        // Define clump properties
        float clump_mass = particle_mass * num_sphere;  // 基于球数的质量
        float3 MOI = make_float3(1, 1, 1) * 0.4 * clump_mass * pow(particle_radius, 2);  // 简化的惯性矩
        
        std::vector<float> radii;
        std::vector<float3> relPos;
        
        // Generate sphere configurations (3D)
        float3 seed_pos = make_float3(0);
        for (int j = 0; j < num_sphere; j++) {
            // Random radius
            radii.push_back(((float)rand() / RAND_MAX) * (max_rad - min_rad) + min_rad);
            
            // Position relative to seed (full 3D)
            float3 tmp;
            if (j == 0) {
                tmp = make_float3(0, 0, 0);
            } else {
                tmp.x = ((float)rand() / RAND_MAX) * (max_relpos - min_relpos) + min_relpos;
                tmp.y = ((float)rand() / RAND_MAX) * (max_relpos - min_relpos) + min_relpos;
                tmp.z = ((float)rand() / RAND_MAX) * (max_relpos - min_relpos) + min_relpos;
            }
            tmp += seed_pos;
            relPos.push_back(tmp);
            
            // Update seed position
            seed_pos = relPos[rand() % (j + 1)];
        }
        
        // Create and store clump template
        auto clump_ptr = DEMSim.LoadClumpType(clump_mass, MOI, radii, relPos, mat_type_particles);
        clump_types.push_back(clump_ptr);
    }

    // --- 5.2 Particle Placement ---

    // Define particle generation parameters
    float insertion_radius = 10 * particle_diameter;   // 10d直径的圆形区域
    float initial_height = 20 * particle_diameter;      // 起始高度 20d
    int total_particles = 5000;                        // 总颗粒数
    int layers = 50;                                    // 分层数
    float layer_spacing = 2.0 * particle_diameter;      // 层间距

    std::vector<std::shared_ptr<DEMClumpTemplate>> input_pile_template_type;
    std::vector<float3> input_pile_xyz;

    // 分层生成颗粒（圆柱形分布，更接近实际投料）
    for (int layer = 0; layer < layers; layer++) {
        float z_height = initial_height + layer * layer_spacing;
        int particles_this_layer = total_particles / layers / 2;
        
        for (int p = 0; p < particles_this_layer; p++) {
            // 在圆形区域内随机分布
            float r = sqrt((float)rand() / RAND_MAX) * insertion_radius;
            float theta = 2 * M_PI * (float)rand() / RAND_MAX;
            
            float3 pos;
            pos.x = r * cos(theta);
            pos.y = r * sin(theta);
            pos.z = z_height;
            
            // 添加轻微的随机偏移，避免过于规则
            pos.x += ((float)rand() / RAND_MAX - 0.5) * particle_radius * 0.2;
            pos.y += ((float)rand() / RAND_MAX - 0.5) * particle_radius * 0.2;
            pos.z += ((float)rand() / RAND_MAX - 0.5) * particle_radius * 0.2;
            
            input_pile_xyz.push_back(pos);
            input_pile_template_type.push_back(clump_types[rand() % num_template]);
        }
    }

    // Add clumps to simulation
    auto the_pile = DEMSim.AddClumps(input_pile_template_type, input_pile_xyz);
    the_pile->SetFamily(FAMILY_PARTICLES);

    std::cout << "Generated " << input_pile_xyz.size() << " particles" << std::endl;

    // =========================================================================
    // 6. SIMULATION INITIALIZATION 
    // =========================================================================
    
    // Initialize the simulation system
    DEMSim.Initialize();
    
    // Set initial conditions for SUP model
    DEMSim.SetFamilyContactWildcardValueBoth(1, "scale_factor_l", my_scale_factor);

    // =========================================================================
    // 7. OUTPUT SETUP 
    // =========================================================================
    
    // Create output directory
    path out_dir = current_path();
    out_dir += "/SUP_Repose";
    create_directory(out_dir);

    path sub_dir = out_dir / "my_output";
    create_directory(sub_dir);

    // =========================================================================
    // 8. SIMULATION EXECUTION 
    // =========================================================================
    
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    
    // Main simulation loop
    for (int i = 0; i < 20; i++) {
        // 生成输出文件名
        char outputfile[200], unifiedfile[200];
        sprintf(unifiedfile, "%s/SUPdemo_unified_%04d.csv", out_dir.c_str(), i);                
        sprintf(outputfile, "%s/SUPdemo_output_%04d.vtk", sub_dir.c_str(), i);
        
        // 写入输出文件
        DEMSim.WriteSphereAndContactFile(std::string(outputfile));
        DEMSim.WriteUnifiedFile(std::string(unifiedfile));
        
        // Progress report
        std::cout << "Frame: " << i << std::endl;
        
        // Advance simulation by 0.1 seconds
        DEMSim.DoDynamics(1e-1); // 0.1 seconds time step
        DEMSim.ShowThreadCollaborationStats();
    }

    // =========================================================================
    // 9. POST-PROCESSING 
    // =========================================================================
    
    // Performance statistics
    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_sec = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    std::cout << time_sec.count() << " seconds (wall time) to finish the simulation" << std::endl;
    
    DEMSim.ShowTimingStats();
    
    // Clean up resources
    DEMSim.ClearTimingStats();
    
    std::cout << "----------------------------------------" << std::endl;
    DEMSim.ShowMemStats();
    std::cout << "----------------------------------------" << std::endl;
    
    std::cout << "DEMdemo_Repose2D_SUP exiting..." << std::endl;
    return 0;
}