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
    const float my_scale_factor = 1.0f;                      // SUP缩放因子（1）
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
        {"Cohesion", 0.08}        // 粘聚力
    });

    // 定义颗粒材料
    auto mat_type_particles = DEMSim.LoadMaterial({
        {"E", 1e7},         
        {"nu", 0.3},        
        {"CoR", 0.7},        
        {"mu", 0.5},
        {"Crr", 0.03},
        {"Cohesion", 0.08}
    });
    
    // 设置材料相互作用属性
    DEMSim.SetMaterialPropertyPair("CoR", mat_type_walls, mat_type_particles, 0.3);
    DEMSim.SetMaterialPropertyPair("Crr", mat_type_walls, mat_type_particles, 0.5);
    DEMSim.SetMaterialPropertyPair("mu", mat_type_walls, mat_type_particles, 0.5);
    DEMSim.SetMaterialPropertyPair("Cohesion", mat_type_walls, mat_type_particles, 0.08);
 
    // DEMSim.UseFrictionalHertzianModel();

    // // 加载SUP接触力模型
    auto model_SUP = DEMSim.ReadContactForceModel("ForceModelSUP.cu");
    model_SUP->SetMustHaveMatProp({"E", "nu", "CoR", "mu", "Crr", "Cohesion"});
    model_SUP->SetMustPairwiseMatProp({"CoR", "mu", "Crr", "Cohesion"});
    model_SUP->SetPerContactWildcards({"delta_time", "delta_tan_x", "delta_tan_y", "delta_tan_z", "scale_factor_l"});


    // =========================================================================
    // 3. 仿真域设置
    // =========================================================================

    // 颗粒基本参数
    float particle_diameter = 0.0005f;  // 0.5mm直径
    float particle_radius = particle_diameter / 2.0f;  // 半径
    float particle_density = 1000.0;                   // 颗粒密度：1000 kg/m³ （塑料颗粒）
    float particle_mass = particle_density * (4.0/3.0) * 3.14159265359 * pow(particle_radius, 3); // 颗粒质量

    // 定义底板和仿真域
    float domain_size = 0.02;       // 20mm × 20mm 水平尺寸
    float plate_bottom = 0.f;       // Z坐标：底板位置

    // 定义时间步长
    float step_size = 5e-7; 

    // 设置仿真域大小（50 × 50 × 300）
    DEMSim.InstructBoxDomainDimension(
        {-domain_size/2, domain_size/2},            // X方向：-10, 10
        {-domain_size/2, domain_size/2},            // Y方向：-10, 10
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
    // 5. 粒子生成设置
    // =========================================================================

    // --- 5.1 粒子/团簇模板 ---

    // 模板生成参数
    int num_template = 6;                              // 随机团簇模板总数
    int min_sphere = 1;                                // 每个团簇最少球数
    int max_sphere = 1;                                // 每个团簇最多球数
    float min_rad = particle_radius;
    float max_rad = particle_radius;            
    float min_relpos = 0;
    float max_relpos = 0;

    // 创建数组存储团簇模板
    std::vector<std::shared_ptr<DEMClumpTemplate>> clump_types;

    // 生成随机团簇模板
    for (int i = 0; i < num_template; i++) {
        int num_sphere = rand() % (max_sphere - min_sphere + 1) + 1;
        
        // 定义团簇属性
        float clump_mass = particle_mass * (float)num_sphere;  // 基于球数的质量
        float3 MOI = make_float3(1, 1, 1) * 0.4 * clump_mass * pow(particle_radius, 2);  // 简化的惯性矩
        
        std::vector<float> radii;
        std::vector<float3> relPos;
        
        // 生成球配置 (3D)
        float3 seed_pos = make_float3(0);
        for (int j = 0; j < num_sphere; j++) {
            // 随机半径
            radii.push_back(((float)rand() / RAND_MAX) * (max_rad - min_rad) + min_rad);
            
            // 相对于种子的位置 (完整3D)
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
            
            // 更新种子位置
            seed_pos = relPos[rand() % (j + 1)];
        }
        
        // 创建并存储团簇模板
        auto clump_ptr = DEMSim.LoadClumpType(clump_mass, MOI, radii, relPos, mat_type_particles);
        clump_types.push_back(clump_ptr);
    }

    // --- 5.2 粒子放置与双重控制模式 ---

    // ========== 控制模式选择 ==========
    enum class ControlMode {
        BY_COUNT,    // 按粒子数量控制
        BY_MASS      // 按总质量控制
    };
    
    // 选择控制模式（修改这里来切换模式）
    ControlMode control_mode = ControlMode::BY_COUNT;  // 可改为 ControlMode::BY_MASS
    
    // 控制参数
    int target_particle_count = 8000;      // 目标粒子数量（用于 BY_COUNT 模式）
    float target_total_mass = 0.001f;       // 目标总质量，单位：kg（用于 BY_MASS 模式）
    
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
    float insertion_size = 20.0f * particle_diameter;           // 20d × 20d 圆形插入区域
    float insertion_top = 190.0f * particle_diameter;           // 从仿真域顶部190d开始 
    float insertion_bottom = 10.0f * particle_diameter;         // 到底板上方10d为止
    float insertion_height = insertion_top - insertion_bottom;  // 插入区域高度

    // 设置泊松盘采样器
    PDSampler sampler(spacing);

     // =========================================================================
    // 6. 分批插入粒子准备
    // =========================================================================
    
    // 计算每批插入的粒子数（分5批）
    const int num_batches = 5;
    int particles_per_batch = max_particles / num_batches;
    int remaining_particles = max_particles % num_batches;  // 剩余粒子加到最后一批
    
    // 存储每批的数据（使用vector的vector代替预分配的DEMClumpBatch）
    std::vector<std::vector<std::shared_ptr<DEMClumpTemplate>>> batch_templates(num_batches);
    std::vector<std::vector<float3>> batch_positions(num_batches);
    std::vector<std::vector<float3>> batch_velocities(num_batches);
    std::vector<std::vector<unsigned int>> batch_families(num_batches);
    
    // 预生成所有粒子位置和类型
    std::vector<std::shared_ptr<DEMClumpTemplate>> all_template_types;
    std::vector<float3> all_positions;
    
    // 统计变量
    int total_generated = 0;
    float current_total_mass = 0.0f;
    int layer_count = 0;

    // 逐层在区域内生成所有粒子位置
    float layer_z = insertion_bottom;
    while (layer_z < insertion_top && total_generated < max_particles) {
        // 当前层的中心位置
        float3 layer_center = make_float3(0, 0, layer_z + spacing / 2);

        // 在圆柱形区域内采样
        auto layer_xyz = sampler.SampleCylinderZ(layer_center, insertion_size / 2, 0);
        
        if (layer_xyz.empty()) {
            std::cout << "无法在 z=" << layer_z / particle_diameter << "d 处生成更多粒子，空间已满" << std::endl;
            break;
        }
        
        // 检查是否会超过最大粒子数
        if (total_generated + layer_xyz.size() > max_particles) {
            int particles_needed = max_particles - total_generated;
            layer_xyz.resize(particles_needed);
        }
        
        // 为每个粒子分配类型并存储
        for (size_t i = 0; i < layer_xyz.size(); i++) {
            all_template_types.push_back(clump_types.at(i % num_template));
            all_positions.push_back(layer_xyz[i]);
            current_total_mass += particle_mass;
        }
        
        total_generated += layer_xyz.size();
        layer_count++;
        layer_z += spacing;
    }

    // 将生成的粒子分配到各批次
    int particle_index = 0;
    for (int batch = 0; batch < num_batches; batch++) {
        int batch_size = (batch == num_batches - 1) ? 
                         particles_per_batch + remaining_particles : 
                         particles_per_batch;
        
        if (particle_index + batch_size > total_generated) {
            batch_size = total_generated - particle_index;
        }
        
        if (batch_size <= 0) continue;
        
        // 预分配空间
        batch_templates[batch].reserve(batch_size);
        batch_positions[batch].reserve(batch_size);
        batch_velocities[batch].reserve(batch_size);
        batch_families[batch].reserve(batch_size);
        
        // 准备当前批次的数据
        for (int i = 0; i < batch_size; i++) {
            if (particle_index >= total_generated) break;
            
            batch_templates[batch].push_back(all_template_types[particle_index]);
            batch_positions[batch].push_back(all_positions[particle_index]);
            batch_velocities[batch].push_back(make_float3(0, 0, -0.2));  // 初始向下速度
            batch_families[batch].push_back(FAMILY_PARTICLES);
            
            particle_index++;
        }
        
        std::cout << "批次 " << batch + 1 << " 准备完成，包含 " << batch_size << " 个粒子" << std::endl;
    }

    // 输出最终统计信息
    std::cout << "\n=== 粒子生成准备完成 ===" << std::endl;
    std::cout << "生成层数: " << layer_count << std::endl;
    std::cout << "实际准备粒子数: " << total_generated << std::endl;
    std::cout << "实际总质量: " << current_total_mass * 1000 << " g" << std::endl;
    
    switch (control_mode) {
        case ControlMode::BY_COUNT:
            std::cout << "目标粒子数: " << target_particle_count << std::endl;
            std::cout << "完成率: " << (total_generated * 100.0 / target_particle_count) << "%" << std::endl;
            break;
            
        case ControlMode::BY_MASS:
            std::cout << "目标总质量: " << target_total_mass * 1000 << " g" << std::endl;
            std::cout << "完成率: " << (current_total_mass / target_total_mass * 100) << "%" << std::endl;
            break;
    }
    std::cout << "================" << std::endl;

    // =========================================================================
    // 7. SIMULATION INITIALIZATION 
    // =========================================================================
    
    // Initialize the simulation system
    DEMSim.Initialize();
    
    // Set initial conditions for SUP model
    DEMSim.SetFamilyContactWildcardValueBoth(1, "scale_factor_l", my_scale_factor);

    // =========================================================================
    // 8. OUTPUT SETUP 
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
    // 9. 仿真执行（带分批插入）
    // =========================================================================
    
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    
    int current_batch = 0;
    
    // 主仿真循环
    for (int frame = 0; frame < 10; frame++) {
        // 在前5帧分批插入粒子
        if (frame < num_batches && current_batch < num_batches) {
            std::cout << "\n=== 第 " << frame + 1 << " 帧：插入批次 " << current_batch + 1 << " ===" << std::endl;
            
            // 获取当前批次的粒子数
            size_t batch_particle_count = batch_templates[current_batch].size();
            
            if (batch_particle_count > 0) {
                // 创建当前批次的DEMClumpBatch
                DEMClumpBatch current_batch_data(batch_particle_count);
                current_batch_data.SetTypes(batch_templates[current_batch]);
                current_batch_data.SetPos(batch_positions[current_batch]);
                current_batch_data.SetVel(batch_velocities[current_batch]);
                current_batch_data.SetFamilies(batch_families[current_batch]);
                
                // 添加当前批次
                auto clump_tracker = DEMSim.AddClumps(current_batch_data);
                std::cout << "成功插入 " << batch_particle_count << " 个粒子" << std::endl;
                
                // 更新仿真系统
                DEMSim.UpdateClumps();
                DEMSim.ClearCache();  // 清理缓存以释放内存
            }
            
            current_batch++;
        }
        
        // 生成输出文件名
        char outputfile[200], unifiedfile[200], meshfile[200];
        sprintf(unifiedfile, "%s/SUPdemo_unified_%04d.csv", out_dir.c_str(), frame);                
        sprintf(outputfile, "%s/SUPdemo_output_%04d.vtk", out_dir.c_str(), frame);
        sprintf(meshfile, "%s/SUPdemo_mesh_%04d.vtk", out_dir.c_str(), frame);
        
        // 写入输出文件
        DEMSim.WriteSphereAndContactFile(std::string(outputfile));
        DEMSim.WriteUnifiedFile(std::string(unifiedfile));
        DEMSim.WriteMeshFile(std::string(meshfile));
        
        // 进度报告
        std::cout << "帧: " << frame << " - 当前系统中粒子数: " << DEMSim.GetNumClumps() << std::endl;
        
        // 推进仿真0.05秒
        DEMSim.DoDynamics(5e-2);
        DEMSim.ShowThreadCollaborationStats();
    }

    // =========================================================================
    // 10. 后处理
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
    
    std::cout << "SUPdemo_Repose_BatchInsert 退出..." << std::endl;
    return 0;
}