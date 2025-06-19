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
    const float my_scale_factor = 1.0f;                      // SUP缩放因子（2）
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
        {"E", 1e7},         // 杨氏模量
        {"nu", 0.3},        // 泊松比
        {"CoR", 0.3},       // 恢复系数
        {"mu", 1},          // 滑动摩擦系数
        {"Crr", 0.06},      // 滚动摩擦系数
        {"SuE", 0.1}        // 粘聚力
    });

    // 定义颗粒材料
    auto mat_type_particles = DEMSim.LoadMaterial({
        {"E", 1e7},         
        {"nu", 0.3},        
        {"CoR", 0.7},        
        {"mu", 0.5},
        {"Crr", 0.03},
        {"SuE", 0.05}
    });
    
    // 设置材料相互作用属性
    DEMSim.SetMaterialPropertyPair("CoR", mat_type_walls, mat_type_particles, 0.3);
    DEMSim.SetMaterialPropertyPair("Crr", mat_type_walls, mat_type_particles, 0.5);
    DEMSim.SetMaterialPropertyPair("mu", mat_type_walls, mat_type_particles, 0.5);
    DEMSim.SetMaterialPropertyPair("SuE", mat_type_walls, mat_type_particles, 0.05);
 
    // DEMSim.UseFrictionalHertzianModel();

    // // 加载SUP接触力模型
    auto model_SUP = DEMSim.ReadContactForceModel("ForceModelSUP.cu");
    model_SUP->SetMustHaveMatProp({"E", "nu", "CoR", "mu", "Crr", "SuE"});
    model_SUP->SetMustPairwiseMatProp({"CoR", "mu", "Crr", "SuE"});
    model_SUP->SetPerContactWildcards({"delta_time", "delta_tan_x", "delta_tan_y", "delta_tan_z", "scale_factor_l"});


    // =========================================================================
    // 3. 仿真域设置
    // =========================================================================

    // 颗粒基本参数
    float particle_diameter = 0.0005f;  // 0.5mm直径
    float particle_radius = particle_diameter / 2.0f;  // 半径
    float particle_density = 1000.0;                   // 颗粒密度：1000 kg/m³ （塑料颗粒）
    float particle_mass = particle_density * (4.0/3.0) * 3.14159265359 * pow(particle_radius, 3); // 颗粒质量

    float plate_bottom = 0.f;          // Z坐标：底板位置

    // 定义底板和仿真域参数（根据论文）
    float domain_size = 50 * particle_diameter;       // 30d × 30d 水平尺寸

    // 定义时间步长
    float step_size = 5e-7; 

    // 设置仿真域大小（50 × 50 × 300）
    DEMSim.InstructBoxDomainDimension(
        {-domain_size/2, domain_size/2},            // X方向：-25d, 25d
        {-domain_size/2, domain_size/2},            // Y方向：-25d, 25d
        {(plate_bottom - 5) * particle_diameter, (plate_bottom + 395) * particle_diameter}      // Z方向：-5d, 295d
    );

    // 设置边界条件 - 顶部开放
    DEMSim.InstructBoxDomainBoundingBC("top_open", mat_type_walls);

    // 设置物理参数
    DEMSim.SetInitTimeStep(step_size);
    DEMSim.SetGravitationalAcceleration(make_float3(0, 0, -9.81));

    // 设置最大速度（防止数值不稳定）
    DEMSim.SetMaxVelocity(25.0f);  // 降低到更合理的值

    // 设置错误输出速度（用于处理初始沉降）
    DEMSim.SetErrorOutVelocity(50.0f);  // 也降低这个值

    // =========================================================================
    // 4. GEOMETRY CREATION
    // =========================================================================

    // 加载平面网格
    auto plate = DEMSim.AddWavefrontMeshObject(GetDEMEDataFile("mesh/plane_20by20.obj"), mat_type_walls);
    plate->Scale(0.005);

    // 设置为固定底板
    plate->SetFamily(FAMILY_PLATE);
    DEMSim.SetFamilyFixed(FAMILY_PLATE);

    // =========================================================================
    // 5. PARTICLE GENERATION
    // =========================================================================

        // --- 5.1 Particle/Clump Templates ---

        // Template generation parameters
        int num_template = 6;                              // Total number of random clump templates
        int min_sphere = 1;                                // Minimum number of spheres per clump
        int max_sphere = 1;                                // Maximum number of spheres per clump
        float min_rad = particle_radius;
        float max_rad = particle_radius;            
        float min_relpos = 0;
        float max_relpos = 0;

        // Create array to store clump templates
        std::vector<std::shared_ptr<DEMClumpTemplate>> clump_types;

        // Generate random clump templates
        for (int i = 0; i < num_template; i++) {
            // 可以使用单球或多球团簇
            int num_sphere = rand() % (max_sphere - min_sphere + 1) + 1;
            
            // Define clump properties
            float clump_mass = particle_mass * (float)num_sphere;  // 基于球数的质量
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

        // --- 5.2 Particle Placement with Dual Control Mode ---

        // ========== 控制模式选择 ==========
        enum class ControlMode {
            BY_COUNT,    // 按粒子数量控制
            BY_MASS      // 按总质量控制
        };
        
        // 选择控制模式（修改这里来切换模式）
        ControlMode control_mode = ControlMode::BY_COUNT;  // 可改为 ControlMode::BY_MASS
        
        // 控制参数
        int target_particle_count = 10000;      // 目标粒子数量（用于 BY_COUNT 模式）
        float target_total_mass = 0.001f;      // 目标总质量，单位：kg（用于 BY_MASS 模式）
        
        // 根据控制模式计算最大粒子数
        int max_particles;
        float expected_total_mass;
        
        switch (control_mode) {
            case ControlMode::BY_COUNT:
                max_particles = target_particle_count;
                expected_total_mass = target_particle_count * particle_mass;
                std::cout << "\n=== 使用粒子数量控制模式 ===" << std::endl;
                std::cout << "目标粒子数量: " << target_particle_count << std::endl;
                std::cout << "预期总质量: " << expected_total_mass * 1000 << " g" << std::endl;
                break;
                
            case ControlMode::BY_MASS:
                max_particles = (int)(target_total_mass / particle_mass);
                expected_total_mass = target_total_mass;
                std::cout << "\n=== 使用总质量控制模式 ===" << std::endl;
                std::cout << "目标总质量: " << target_total_mass * 1000 << " g" << std::endl;
                std::cout << "预期粒子数量: " << max_particles << std::endl;
                break;
        }
        
        std::cout << "单个粒子质量: " << particle_mass * 1000000 << " mg" << std::endl;
        std::cout << "单个粒子直径: " << particle_diameter * 1000 << " mm" << std::endl;

        // 插入区域参数
        float spacing = particle_diameter * 1.3f;                   // 粒子间距
        float insertion_size = 20.0f * particle_diameter;           // 16d × 16d 圆形插入区域
        float insertion_top = 190.0f * particle_diameter;           // 从仿真域顶部190d开始 
        float insertion_bottom = 10.0f * particle_diameter;         // 到底板上方10d为止
        float insertion_height = insertion_top - insertion_bottom;  // 插入区域高度

        // 估算区域容量
        float cylinder_volume = 3.14159f * pow(insertion_size/2, 2) * insertion_height;
        float packing_fraction = 0.6f;  // 随机堆积的典型体积分数
        float particle_volume = (4.0f/3.0f) * 3.14159f * pow(particle_radius, 3);
        int estimated_capacity = (int)(cylinder_volume * packing_fraction / particle_volume);
        
        std::cout << "\n插入区域信息:" << std::endl;
        std::cout << "  圆柱直径: " << insertion_size / particle_diameter << "d" << std::endl;
        std::cout << "  高度: " << insertion_height / particle_diameter << "d" << std::endl;
        std::cout << "  估计容量: " << estimated_capacity << " 个粒子" << std::endl;
        
        if (max_particles > estimated_capacity) {
            std::cout << "警告: 目标粒子数(" << max_particles << ")超过估计容量(" 
                    << estimated_capacity << ")，可能无法全部生成！" << std::endl;
        }

        // 设置泊松盘采样器
        PDSampler sampler(spacing);

        // 生成粒子的容器
        std::vector<std::shared_ptr<DEMClumpTemplate>> input_pile_template_type;
        std::vector<float3> input_pile_xyz;
        
        // 统计变量
        int particle_count = 0;
        float current_total_mass = 0.0f;
        int layer_count = 0;

        // 逐层在区域内生成粒子
        float layer_z = insertion_bottom;
        while (layer_z < insertion_top && particle_count < max_particles) {
            // 当前层的中心位置
            float3 layer_center = make_float3(0, 0, layer_z + spacing / 2);

            // 在圆柱形区域内采样
            auto layer_xyz = sampler.SampleCylinderZ(layer_center, insertion_size / 2, 0);
            
            // 如果没有生成粒子，说明空间已满，退出
            if (layer_xyz.empty()) {
                std::cout << "无法在 z=" << layer_z / particle_diameter << "d 处生成更多粒子，空间已满" << std::endl;
                break;
            }
            
            // 检查是否会超过最大粒子数
            if (particle_count + layer_xyz.size() > max_particles) {
                // 只取需要的粒子数
                int particles_needed = max_particles - particle_count;
                layer_xyz.resize(particles_needed);
                std::cout << "达到目标，最后一层只生成 " << particles_needed << " 个粒子" << std::endl;
            }
            
            // 为每个粒子分配类型
            for (size_t i = 0; i < layer_xyz.size(); i++) {
                input_pile_template_type.push_back(clump_types.at(i % num_template));
                current_total_mass += particle_mass;  // 更新当前总质量
            }
            
            // 添加到总列表
            input_pile_xyz.insert(input_pile_xyz.end(), layer_xyz.begin(), layer_xyz.end());
            
            // 更新统计
            particle_count += layer_xyz.size();
            layer_count++;
            
            // 进度显示（每10层显示一次）
            if (layer_count % 10 == 0) {
                float progress = (float)particle_count / max_particles * 100;
                std::cout << "进度: " << std::fixed << std::setprecision(1) << progress << "% "
                        << "- 层 " << layer_count << " (z=" << layer_z / particle_diameter << "d), "
                        << "已生成: " << particle_count << "/" << max_particles << " 粒子" << std::endl;
            }
            
            // 移动到下一层
            layer_z += spacing;
        }

        // Add clumps to simulation
        auto the_pile = DEMSim.AddClumps(input_pile_template_type, input_pile_xyz);
        the_pile->SetFamily(FAMILY_PARTICLES);

        // 输出最终统计信息
        std::cout << "\n=== 粒子生成完成 ===" << std::endl;
        std::cout << "生成层数: " << layer_count << std::endl;
        std::cout << "实际生成粒子数: " << particle_count << std::endl;
        std::cout << "实际总质量: " << current_total_mass * 1000 << " g" << std::endl;
        
        switch (control_mode) {
            case ControlMode::BY_COUNT:
                std::cout << "目标粒子数: " << target_particle_count << std::endl;
                std::cout << "完成率: " << (particle_count * 100.0 / target_particle_count) << "%" << std::endl;
                break;
                
            case ControlMode::BY_MASS:
                std::cout << "目标总质量: " << target_total_mass * 1000 << " g" << std::endl;
                std::cout << "完成率: " << (current_total_mass / target_total_mass * 100) << "%" << std::endl;
                break;
        }
        
        std::cout << "\n颗粒参数汇总：" << std::endl;
        std::cout << "  直径: " << particle_diameter * 1000 << " mm" << std::endl;
        std::cout << "  密度: " << particle_density << " kg/m³" << std::endl;
        std::cout << "  单个质量: " << particle_mass * 1000000 << " mg" << std::endl;
        std::cout << "  粒子间距: " << spacing / particle_diameter << "d" << std::endl;
        std::cout << "  堆积高度: " << (layer_count * spacing) / particle_diameter << "d" << std::endl;

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
    out_dir += "/SUP_Repose_Factor_1";
    // create_directory(out_dir);

    // Check if directory exists and clear it, otherwise create it
    if (exists(out_dir)) {
        // Directory exists, remove all files inside
        std::cout << "Output directory exists, clearing contents..." << std::endl;
        for (const auto& entry : directory_iterator(out_dir)) {
            remove_all(entry.path());
        }
    } else {
        // Directory doesn't exist, create it
        create_directory(out_dir);
        std::cout << "Created output directory: " << out_dir << std::endl;
    }


    // =========================================================================
    // 8. SIMULATION EXECUTION 
    // =========================================================================
    
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    
    // Main simulation loop
    for (int i = 0; i < 20; i++) {
        // 生成输出文件名
        char outputfile[200], unifiedfile[200], meshfile[200];
        sprintf(unifiedfile, "%s/SUPdemo_unified_%04d.csv", out_dir.c_str(), i);                
        sprintf(outputfile, "%s/SUPdemo_output_%04d.vtk", out_dir.c_str(), i);
        sprintf(meshfile, "%s/SUPdemo_mesh_%04d.vtk", out_dir.c_str(), i);
        
        // 写入输出文件
        DEMSim.WriteSphereAndContactFile(std::string(outputfile));
        DEMSim.WriteUnifiedFile(std::string(unifiedfile));
        DEMSim.WriteMeshFile(std::string(meshfile));
        
        // Progress report
        std::cout << "Frame: " << i << std::endl;
        
        // Advance simulation by 0.05 seconds
        DEMSim.DoDynamics(5e-2); // 0.05 seconds time step
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
    
    std::cout << "SUPdemo_Repose exiting..." << std::endl;
    return 0;
}