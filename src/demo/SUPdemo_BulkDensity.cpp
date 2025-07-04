//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//
//	SPDX-License-Identifier: BSD-3-Clause

// =============================================================================
// bulk density demo using SUP contact model
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
    const float step_size = 5e-7; 
    // 定义每帧仿真时间
    const float time_per_frame = 5e-3;
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
        {"CoR", 0.5},       // 恢复系数
        {"mu", 0.5},          // 滑动摩擦系数
        {"Crr", 0},      // 滚动摩擦系数
        {"Cohesion", 0}        // 粘聚力
    });

    // 定义颗粒材料
    auto mat_type_particles = DEMSim.LoadMaterial({
        {"E", 1e7},         
        {"nu", 0.3},        
        {"CoR", 0.5},        
        {"mu", 0.5},
        {"Crr", 0},
        {"Cohesion", 0}
    });

    auto mat_type_smooth = DEMSim.LoadMaterial({
        {"E", 1e7},         
        {"nu", 0.3},        
        {"CoR", 0.1},        
        {"mu", 0},
        {"Crr", 0},
        {"Cohesion", 0}
    });
    
    // 设置材料相互作用属性 - 颗粒与墙壁
    DEMSim.SetMaterialPropertyPair("CoR", mat_type_walls, mat_type_particles, 0.5); 
    DEMSim.SetMaterialPropertyPair("mu", mat_type_walls, mat_type_particles, 0.5);
    DEMSim.SetMaterialPropertyPair("Crr", mat_type_walls, mat_type_particles, 0.03);
    DEMSim.SetMaterialPropertyPair("Cohesion", mat_type_walls, mat_type_particles, 0);

    // 设置材料相互作用属性 - 平滑圆筒
    DEMSim.SetMaterialPropertyPair("CoR", mat_type_smooth, mat_type_particles, 0.1);
    DEMSim.SetMaterialPropertyPair("mu", mat_type_smooth, mat_type_particles, 0.0);
    DEMSim.SetMaterialPropertyPair("Crr", mat_type_smooth, mat_type_particles, 0.0);
    DEMSim.SetMaterialPropertyPair("Cohesion", mat_type_smooth, mat_type_particles, 0.0);

    // 设置材料相互作用属性 - 平滑圆筒与墙壁
    DEMSim.SetMaterialPropertyPair("CoR", mat_type_smooth, mat_type_walls, 0.1);
    DEMSim.SetMaterialPropertyPair("mu", mat_type_smooth, mat_type_walls, 0.0);
    DEMSim.SetMaterialPropertyPair("Crr", mat_type_smooth, mat_type_walls, 0.0);
    DEMSim.SetMaterialPropertyPair("Cohesion", mat_type_smooth, mat_type_walls, 0.0);
 
    // 加载SUP接触力模型
    auto model_SUP = DEMSim.ReadContactForceModel("ForceModelSUP.cu");
    model_SUP->SetMustHaveMatProp({"E", "nu", "CoR", "mu", "Crr", "Cohesion"});
    model_SUP->SetMustPairwiseMatProp({"CoR", "mu", "Crr", "Cohesion"});
    model_SUP->SetPerContactWildcards({"delta_time", "delta_tan_x", "delta_tan_y", "delta_tan_z", "scale_factor_l"});


    // =========================================================================
    // 3. 仿真域设置
    // =========================================================================

    // 颗粒基本参数
    float particle_diameter = 0.0005f * my_scale_factor;                 // 0.5mm直径
    float particle_radius = particle_diameter / 2.0f;  // 半径
    float particle_density = 1000.0;                   // 颗粒密度：1000 kg/m³ （塑料颗粒）
    float particle_mass = particle_density * (4.0/3.0) * 3.14159265359 * pow(particle_radius, 3); // 颗粒质量

    // 定义仿真域
    float domain_size = 0.035;       // 30mm × 30mm 水平尺寸
    float plate_bottom = 0.f;       // Z坐标：底板位置

    // 设置仿真域大小（50 × 50 × 300）
    DEMSim.InstructBoxDomainDimension(
        {-domain_size/2, domain_size/2},            // X方向：-50mm, 50mm
        {-domain_size/2, domain_size/2},            // Y方向：-50mm, 50mm
        {(plate_bottom), (plate_bottom + 1)}      // Z方向：0, 100mm
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
    // 4. GEOMETRY CREATION
    // =========================================================================

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
    int target_particle_count = 32000;      // 目标粒子数量（用于 BY_COUNT 模式）
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
    float spacing = particle_diameter * 1.3f;      // 粒子间距
    float insertion_radius = 0.015;                 // 13mm 半径的插入区域
    float fixed_insertion_height = 0.00003;          // 每批粒子的起始高度
    // 设置泊松盘采样器
    PDSampler sampler(spacing);

    // =========================================================================
    // 6. 分批插入粒子准备
    // =========================================================================
    
    // 计算每批插入的粒子数（分16批）
    const int num_batches = 1;
    const int frames_between_batches = 9;
    int particles_per_batch = max_particles / num_batches;
    int remaining_particles = max_particles % num_batches;  // 剩余粒子加到最后一批
    float batch_interval_time = frames_between_batches * time_per_frame; // 每批次间隔时间

    // 存储每批的数据（使用vector的vector代替预分配的DEMClumpBatch）
    std::vector<std::vector<std::shared_ptr<DEMClumpTemplate>>> batch_templates(num_batches);
    std::vector<std::vector<float3>> batch_positions(num_batches);
    std::vector<std::vector<float3>> batch_velocities(num_batches);
    std::vector<std::vector<unsigned int>> batch_families(num_batches);
    
    // 预生成所有粒子位置和类型
    std::vector<std::shared_ptr<DEMClumpTemplate>> all_template_types;
    std::vector<float3> all_positions;

    // 计算所需的初速度
    // 使用运动学公式：s = v0*t + 0.5*g*t²
    // 我们希望粒子在batch_interval_time内至少下落一定距离
    float min_clearance = 10.0f * particle_diameter;  // 最小安全距离（5mm）
    float gravity = 9.81f;

    // v0 = (s - 0.5*g*t²) / t
    float required_initial_velocity = (min_clearance - 0.5f * gravity * batch_interval_time * batch_interval_time) / batch_interval_time;
    required_initial_velocity = -abs(required_initial_velocity);  // 确保向下

    std::cout << "计算的初速度: " << required_initial_velocity << " m/s" << std::endl;
    
    // 统计变量
    int total_generated = 0;
    float current_total_mass = 0.0f;
    int layer_count = 0;

    // ========== 固定起始高度分批生成 ==========


    std::cout << "开始为每批生成粒子位置..." << std::endl;

    for (int batch = 0; batch < num_batches; batch++) {
        int batch_size = (batch == num_batches - 1) ? 
                        particles_per_batch + remaining_particles : 
                        particles_per_batch;
        
        // 为当前批次生成粒子
        int batch_generated = 0;
        float current_z = fixed_insertion_height;  // 从固定高度开始
        int batch_layers = 0;
        
        std::cout << "\n批次 " << batch + 1 << " 开始生成..." << std::endl;
        
        while (batch_generated < batch_size) {
            float3 layer_center = make_float3(0, 0, current_z);
            
            // 在当前高度的圆形区域内采样
            auto layer_xyz = sampler.SampleCylinderZ(layer_center, insertion_radius, 0);
            
            if (layer_xyz.empty()) {
                std::cout << "  层 " << batch_layers << " 在 z=" << current_z/particle_diameter 
                        << "d 处无法生成粒子，尝试下一层" << std::endl;
                current_z += spacing;  // 向上移动一层
                continue;
            }
            
            // 添加当前层的粒子
            int layer_particles = 0;
            for (size_t i = 0; i < layer_xyz.size() && batch_generated < batch_size; i++) {
                all_template_types.push_back(clump_types.at((total_generated + batch_generated) % num_template));
                all_positions.push_back(layer_xyz[i]);
                batch_generated++;
                layer_particles++;
                current_total_mass += particle_mass;
            }
            
            std::cout << "  层 " << batch_layers << " (z=" << current_z/particle_diameter 
                    << "d): 生成 " << layer_particles << " 个粒子" << std::endl;
            
            batch_layers++;
            current_z += spacing;  // 向上移动到下一层
        }
        
        total_generated += batch_generated;
        
        std::cout << "批次 " << batch + 1 << " 完成: " << batch_generated 
                << " 个粒子，共 " << batch_layers << " 层" 
                << "，高度范围: " << fixed_insertion_height/particle_diameter << "d ~ " 
                << current_z/particle_diameter << "d" << std::endl;
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
            batch_velocities[batch].push_back(make_float3(0, 0, 0));  // 初始向下速度，先设为0
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
    DEMSim.SetFamilyContactWildcardValueBoth(FAMILY_PARTICLES, "scale_factor_l", my_scale_factor);

    // =========================================================================
    // 8. OUTPUT SETUP 
    // =========================================================================
    
    // Create output directory
    path out_dir = current_path();
    out_dir += "/SUP_BulkDensity_coe000hertz";
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
    // 9. 仿真执行
    // =========================================================================
    
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    
    int current_batch = 0;
    
    // 主仿真循环
    for (int frame = 0; frame < 45; frame++) {
        // 在前90帧分批插入粒子
        if (frame % frames_between_batches == 0 && current_batch < num_batches) {
            std::cout << "\n=== 第 " << frame << " 帧：插入批次 " << current_batch + 1 << " ===" << std::endl;
            
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
            }
            
            current_batch++;
        }
        
        // 每5帧输出一次文件
        if (frame % 1 == 0) {
            // 生成输出文件名
            char outputfile[200], unifiedfile[200];
            sprintf(unifiedfile, "%s/SUPdemo_unified_%04d.csv", out_dir.c_str(), frame);                
            sprintf(outputfile, "%s/SUPdemo_output_%04d.vtk", out_dir.c_str(), frame);
            
            // 写入输出文件
            DEMSim.WriteSphereAndContactFile(std::string(outputfile));
            DEMSim.WriteUnifiedFile(std::string(unifiedfile));
            
            // 进度报告
            std::cout << "帧: " << frame << " - 输出文件已保存 - 当前系统中粒子数: " << DEMSim.GetNumClumps() << std::endl;

        } else {
            // 非输出帧，只报告进度
            std::cout << "帧: " << frame << " - 当前系统中粒子数: " << DEMSim.GetNumClumps() << std::endl;
        }
        
        // 推进仿真0.005秒
        DEMSim.DoDynamics(time_per_frame); // 每帧推进0.005秒
        DEMSim.ShowThreadCollaborationStats();
    }

    // =========================================================================
    // 10. 创建检查器 (用于计算体积密度)
    // =========================================================================    
    // 首先创建最大Z坐标检查器，获取堆积高度
    auto max_z_finder = DEMSim.CreateInspector("clump_max_z");
    
    // 获取最高点高度
    float max_z = max_z_finder->GetValue();
    float pile_height = max_z - plate_bottom;  // 堆积总高度
    
    std::cout << "\n堆积高度信息：" << std::endl;
    std::cout << "  最高点Z坐标: " << max_z * 1000 << " mm" << std::endl;
    std::cout << "  堆积总高度: " << pile_height * 1000 << " mm" << std::endl;
    
    // 定义测量区域边界
    float measure_x_min = -domain_size/2 * 0.8f;  // 使用80%的区域避免边界效应
    float measure_x_max = domain_size/2 * 0.8f;
    float measure_y_min = -domain_size/2 * 0.8f;
    float measure_y_max = domain_size/2 * 0.8f;
    
    // Z方向：从堆积高度的20%到60%
    float measure_z_min = plate_bottom + pile_height * 0.2f;
    float measure_z_max = plate_bottom + pile_height * 0.6f;
    
    // 创建带条件的检查器 - 只计算测量区域内的粒子
    char condition[300];
    sprintf(condition, "return (X >= %f) && (X <= %f) && (Y >= %f) && (Y <= %f) && (Z >= %f) && (Z <= %f);",
            measure_x_min, measure_x_max, measure_y_min, measure_y_max, measure_z_min, measure_z_max);
    
    // 创建质量检查器（注意：我们不使用volume_inspector，因为它可能不工作）
    auto mass_inspector = DEMSim.CreateInspector("clump_mass", condition);
    
    // =========================================================================
    // 11. 计算体积密度
    // =========================================================================
    
    std::cout << "\n========== 体积密度计算 ==========" << std::endl;
    
    // 获取测量区域内的质量
    float matter_mass = mass_inspector->GetValue();
    
    // 计算测量区域体积
    float measure_region_volume = (measure_x_max - measure_x_min) * 
                                 (measure_y_max - measure_y_min) * 
                                 (measure_z_max - measure_z_min);
    
    // 使用质量和密度计算物质体积
    float matter_volume = matter_mass / particle_density;
    
    // 计算体积密度（体积分数）
    float bulk_density = matter_volume / measure_region_volume;
    
    // 计算质量密度
    float mass_density = matter_mass / measure_region_volume;
    
    // 计算空隙率和孔隙度
    float void_ratio = (measure_region_volume - matter_volume) / matter_volume;
    float porosity = (measure_region_volume - matter_volume) / measure_region_volume;
    
    // 输出详细结果
    std::cout << "\n测量区域参数：" << std::endl;
    std::cout << "  X方向: [" << measure_x_min * 1000 << ", " << measure_x_max * 1000 << "] mm" << std::endl;
    std::cout << "  Y方向: [" << measure_y_min * 1000 << ", " << measure_y_max * 1000 << "] mm" << std::endl;
    std::cout << "  Z方向: [" << measure_z_min * 1000 << ", " << measure_z_max * 1000 << "] mm" << std::endl;
    std::cout << "  Z方向范围: 堆积高度的20% ~ 60%" << std::endl;
    std::cout << "  测量区域体积: " << measure_region_volume * 1e6 << " mm³" << std::endl;
    
    std::cout << "\n测量结果：" << std::endl;
    std::cout << "  物质体积: " << matter_volume * 1e6 << " mm³" << std::endl;
    std::cout << "  物质质量: " << matter_mass * 1000 << " g" << std::endl;
    
    std::cout << "\n密度指标：" << std::endl;
    std::cout << "  体积密度（体积分数）: " << bulk_density << " (" << bulk_density * 100 << "%)" << std::endl;
    std::cout << "  质量密度: " << mass_density << " kg/m³" << std::endl;
    std::cout << "  空隙率: " << void_ratio << std::endl;
    std::cout << "  孔隙度: " << porosity * 100 << " %" << std::endl;
    
    std::cout << "\n验证计算：" << std::endl;
    std::cout << "  颗粒材料密度: " << particle_density << " kg/m³" << std::endl;
    std::cout << "  相对密度: " << (mass_density / particle_density) * 100 << " %" << std::endl;
    
    // 将结果写入文件（中文）
    std::ofstream density_file(out_dir / "bulk_density_results.txt");
    density_file << "体积密度分析结果\n";
    density_file << "================\n\n";
    
    density_file << "仿真参数：\n";
    density_file << "  颗粒直径: " << particle_diameter * 1000 << " mm\n";
    density_file << "  颗粒密度: " << particle_density << " kg/m³\n";
    density_file << "  颗粒总数: " << DEMSim.GetNumClumps() << " 个\n";
    density_file << "  缩放因子: " << my_scale_factor << "\n\n";
    
    density_file << "堆积信息：\n";
    density_file << "  最高点坐标: " << max_z * 1000 << " mm\n";
    density_file << "  堆积总高度: " << pile_height * 1000 << " mm\n\n";
    
    density_file << "测量区域：\n";
    density_file << "  X方向: [" << measure_x_min * 1000 << ", " << measure_x_max * 1000 << "] mm\n";
    density_file << "  Y方向: [" << measure_y_min * 1000 << ", " << measure_y_max * 1000 << "] mm\n";
    density_file << "  Z方向: [" << measure_z_min * 1000 << ", " << measure_z_max * 1000 << "] mm\n";
    density_file << "  Z方向说明: 堆积高度的20% ~ 60%区间\n";
    density_file << "  区域体积: " << measure_region_volume * 1e6 << " mm³\n\n";
    
    density_file << "测量值：\n";
    density_file << "  物质体积: " << matter_volume * 1e6 << " mm³\n";
    density_file << "  物质质量: " << matter_mass * 1000 << " g\n\n";
    
    density_file << "计算结果：\n";
    density_file << "  体积密度（体积分数）: " << bulk_density << " (" << bulk_density * 100 << "%)\n";
    density_file << "  质量密度: " << mass_density << " kg/m³\n";
    density_file << "  空隙率: " << void_ratio << "\n";
    density_file << "  孔隙度: " << porosity * 100 << " %\n";
    density_file << "  相对密度: " << (mass_density / particle_density) * 100 << " %\n\n";
    
    // 添加结果解释
    density_file << "结果解释：\n";
    density_file << "  - 体积密度表示固体颗粒占测量区域的体积比例\n";
    density_file << "  - 孔隙度表示空隙占测量区域的体积比例\n";
    density_file << "  - 相对密度表示实际密度与材料密度的比值\n";
    density_file << "  - 测量区域选择堆积中部(20%-60%)以获得代表性结果\n";
    
    density_file.close();
    
    std::cout << "\n结果已保存到: " << (out_dir / "bulk_density_results.txt") << std::endl;

    // 额外输出：根据体积密度判断堆积状态
    std::cout << "\n堆积状态分析：" << std::endl;
    if (bulk_density < 0.56) {
        std::cout << "  体积密度 < 0.56，属于非常松散的堆积" << std::endl;
    } else if (bulk_density < 0.60) {
        std::cout << "  体积密度在 0.56-0.60 之间，属于随机松散堆积" << std::endl;
    } else if (bulk_density < 0.64) {
        std::cout << "  体积密度在 0.60-0.64 之间，属于随机密实堆积" << std::endl;
    } else {
        std::cout << "  体积密度 > 0.64，接近密堆积状态" << std::endl;
    }

    // =========================================================================
    // 12. 后处理
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

    std::cout << "SUPdemo_Repose 退出..." << std::endl;
    return 0;
}
