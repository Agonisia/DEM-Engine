//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//
//	SPDX-License-Identifier: BSD-3-Clause

// =============================================================================
// SUP模型搅拌器仿真：双密度倒置配置粒子系统（明确分层版本）
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
    // 1. 仿真设置与SUP模型参数
    // =========================================================================
    
    // === SUP模型核心参数（集中管理） ===
    const float my_scale_factor = 1.0f;                    // SUP缩放因子
    const float scale_force_index = 2.0f;                  // 力的缩放指数
    const int base_particle_count = 1600000;               // 基准粒子数量设定（缩放因子为1时）
    const float base_particle_mass = 0.0458f;              // 基准系统总质量：0.0458kg（固定值）
    const float base_particle_diameter = 0.0005f;          // 基准粒子直径：0.5mm
    
    // === 两种粒子密度 ===
    const float particle_density_light = 1000.0f;          // 轻颗粒密度：1000 kg/m³
    const float particle_density_heavy = 2000.0f;          // 重颗粒密度：2000 kg/m³
    const float particle_cohesion = 0.0f;                  // 颗粒间粘聚力

    // === 时间参数 ===
    const float step_size = 1e-6f;                         // 时间步长
    const float time_per_frame = 1e-2f;                    // 每帧仿真时间

    // === 分批插入参数 ===
    const int num_batches = 16;                            // 分批插入的批次数
    const int frames_between_batches = 5;                  // 每批次之间的帧数
    const int frames_after_insertion = 50;                 // 插入完成后继续运行的帧数
    
    // === 根据缩放因子自动计算的参数 ===
    const float particle_diameter = base_particle_diameter * my_scale_factor;
    const float particle_radius = particle_diameter / 2.0f;
    
    // 轻粒子质量（密度1000）
    const float particle_mass_light = particle_density_light * (4.0/3.0) * 3.14159265359 * 
                               pow(particle_radius, 3);
    // 重粒子质量（密度2000）
    const float particle_mass_heavy = particle_density_heavy * (4.0/3.0) * 3.14159265359 * 
                               pow(particle_radius, 3);
    
    const int scale_factor_cubed = (int)pow(my_scale_factor, 3);
    
    // === 基于基准质量的粒子数计算 ===
    // 使用base_particle_mass作为系统总质量
    const float target_total_mass = base_particle_mass;  // 固定系统总质量：0.0458kg
    
    // 基于1:1质量比例计算粒子数
    const float target_light_mass = target_total_mass / 2.0f;
    const float target_heavy_mass = target_total_mass / 2.0f;
    
    // 计算所需的粒子数量
    const int target_light_particles = (int)(target_light_mass / particle_mass_light);
    const int target_heavy_particles = (int)(target_heavy_mass / particle_mass_heavy);
    const int actual_particle_count = target_light_particles + target_heavy_particles;
    
    // === 明确分层批次分配 ===
    // 根据粒子数比例计算批次分配
    const float light_ratio = (float)target_light_particles / actual_particle_count;
    const float heavy_ratio = (float)target_heavy_particles / actual_particle_count;
    
    // 计算轻粒子和重粒子应占用的批次数
    const int light_batches = (int)ceil(num_batches * light_ratio);  // 轻粒子批次（前几批）
    const int heavy_batches = num_batches - light_batches;           // 重粒子批次（后几批）
    
    // 每批次的粒子数
    const int light_particles_per_batch = target_light_particles / light_batches;
    const int light_particles_remainder = target_light_particles % light_batches;
    const int heavy_particles_per_batch = target_heavy_particles / heavy_batches;
    const int heavy_particles_remainder = target_heavy_particles % heavy_batches;

    // === 仿真域参数 ===
    const float domain_size = 0.1f;                        // 100mm × 100mm 水平尺寸
    const float wall_radius = 0.042f;                      // 圆柱墙体半径42mm
    const float plate_bottom = 0.0f;                       // Z坐标：底板位置
    const float insertion_radius = 0.038f;                 // 插入区域半径38mm
    const float insertion_height = 0.042f;                 // 起始高度 42mm
    const float spacing = particle_diameter * 1.5f;        // 粒子间距
    
    // === Family ID定义 ===
    const unsigned int FAMILY_MIXER = 1;
    const unsigned int FAMILY_PARTICLES_LIGHT = 2;         // 轻粒子
    const unsigned int FAMILY_PARTICLES_HEAVY = 3;         // 重粒子
    
    // 输出SUP模型信息
    std::cout << "\n========== SUP模型参数（倒置双密度系统 - 明确分层版） ==========" << std::endl;
    std::cout << "缩放因子 l: " << my_scale_factor << std::endl;
    std::cout << "基准系统总质量: " << base_particle_mass << " kg（固定值）" << std::endl;
    std::cout << "粒子直径: " << particle_diameter * 1000 << " mm" << std::endl;
    std::cout << "轻粒子单个质量: " << particle_mass_light * 1000 << " g" << std::endl;
    std::cout << "重粒子单个质量: " << particle_mass_heavy * 1000 << " g" << std::endl;
    std::cout << "质量比（轻:重）: 1:1" << std::endl;
    std::cout << "预计轻粒子数: " << target_light_particles << std::endl;
    std::cout << "预计重粒子数: " << target_heavy_particles << std::endl;
    std::cout << "预计总粒子数: " << actual_particle_count << std::endl;
    std::cout << "粒子数比（轻:重）: " << (float)target_light_particles/target_heavy_particles << ":1" << std::endl;
    std::cout << "\n=== 明确分层批次分配 ===" << std::endl;
    std::cout << "总批次数: " << num_batches << std::endl;
    std::cout << "轻粒子批次数（前" << light_batches << "批）: " << light_batches << std::endl;
    std::cout << "重粒子批次数（后" << heavy_batches << "批）: " << heavy_batches << std::endl;
    std::cout << "轻粒子每批: " << light_particles_per_batch << " 个（余数: " << light_particles_remainder << "）" << std::endl;
    std::cout << "重粒子每批: " << heavy_particles_per_batch << " 个（余数: " << heavy_particles_remainder << "）" << std::endl;
    std::cout << "说明: 轻粒子先插入形成下层，重粒子后插入形成上层（不稳定配置）" << std::endl;
    std::cout << "================================\n" << std::endl;

    // 创建求解器实例
    DEMSolver DEMSim;
    
    // 设置输出格式和内容
    DEMSim.SetVerbosity(INFO);
    DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV);
    DEMSim.SetOutputContent({"VEL", "ANG_VEL", "ACC"}); 
    DEMSim.SetNoForceRecord();

    // =========================================================================
    // 2. 材料属性设置
    // =========================================================================
    
    // 定义搅拌器材料
    auto mat_type_mixer = DEMSim.LoadMaterial({
        {"E", 1e7},         // 杨氏模量
        {"nu", 0.3},        // 泊松比
        {"CoR", 0.1},       // 恢复系数
        {"mu", 0.0},        // 滑动摩擦系数
        {"Crr", 0.0},       // 滚动摩擦系数
        {"Cohesion", 0},    // 粘聚力
        {"scale_factor_l", my_scale_factor}
    });
    
    // 定义统一的颗粒材料（通过质量差异实现密度差异）
    auto mat_type_particles = DEMSim.LoadMaterial({
        {"E", 1e7},
        {"nu", 0.3},
        {"CoR", 0.1},
        {"mu", 0.3},
        {"Crr", 0.0},
        {"Cohesion", particle_cohesion},
        {"scale_factor_l", my_scale_factor}
    });
    
    // 设置材料相互作用属性 - 搅拌器与颗粒
    DEMSim.SetMaterialPropertyPair("mu", mat_type_mixer, mat_type_particles, 0.1);
    DEMSim.SetMaterialPropertyPair("CoR", mat_type_mixer, mat_type_particles, 0.0);
    DEMSim.SetMaterialPropertyPair("Crr", mat_type_mixer, mat_type_particles, 0.0);
    DEMSim.SetMaterialPropertyPair("Cohesion", mat_type_mixer, mat_type_particles, 0.0);

    // === 加载SUP接触力模型 ===
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
        {plate_bottom, plate_bottom + 1.0f}         // Z方向：0-1m
    );

    // 设置边界条件
    DEMSim.InstructBoxDomainBoundingBC("top_open", mat_type_mixer);
    
    // 设置物理参数
    DEMSim.SetInitTimeStep(step_size);
    DEMSim.SetGravitationalAcceleration(make_float3(0, 0, -9.81));

    // 错误处理和安全限制
    DEMSim.SetErrorOutVelocity(20000.0);
    DEMSim.SetErrorOutAvgContacts(200);

    // 接触检测和性能设置
    DEMSim.SetCDUpdateFreq(40);
    DEMSim.SetExpandSafetyAdder(0.5);  // 考虑搅拌器叶片的高速运动
    
    // 高级接触检测设置
    DEMSim.SetCDNumStepsMaxDriftMultipleOfAvg(1.2);
    DEMSim.SetCDNumStepsMaxDriftAheadOfAvg(6);
    DEMSim.SetSortContactPairs(true);
    
    // 自适应bin尺寸设置
    DEMSim.SetAdaptiveBinSizeUpperProactivity(0.5);
    DEMSim.SetAdaptiveBinSizeLowerProactivity(0.15);

    // =========================================================================
    // 4. 几何体创建
    // =========================================================================
    
    // === 4.1 添加固定几何体 - 圆柱形边界墙 ===
    auto walls = DEMSim.AddExternalObject();
    walls->AddCylinder(make_float3(0), make_float3(0, 0, 1), 
                      wall_radius, mat_type_mixer, 0);
    
    // === 4.2 添加可动几何体 - 搅拌器网格 ===
    auto mixer = DEMSim.AddWavefrontMeshObject(
        (GET_DATA_PATH() / "mesh/impeller.obj").string(), 
        mat_type_mixer
    );
    std::cout << "搅拌器三角形网格数: " << mixer->GetNumTriangles() << std::endl;

    // 设置搅拌器的预定义角速度（绕z轴旋转，初始不旋转）
    DEMSim.SetFamilyPrescribedAngVel(FAMILY_MIXER, "0", "0", "0");

    // =========================================================================
    // 5. 粒子生成设置（明确分层的倒置双密度系统）
    // =========================================================================

    // --- 5.1 粒子/团簇模板 ---
    
    // 创建轻粒子模板（通过质量实现密度1000）
    std::vector<std::shared_ptr<DEMClumpTemplate>> clump_types_light;
    float3 light_MOI = make_float3(1, 1, 1) * 0.4 * particle_mass_light * pow(particle_radius, 2);
    auto light_clump_ptr = DEMSim.LoadClumpType(
        particle_mass_light,    // 较小的质量（密度1000）
        light_MOI,
        std::vector<float>(1, particle_radius),
        std::vector<float3>(1, make_float3(0, 0, 0)),
        mat_type_particles      // 使用统一的材料类型
    );
    clump_types_light.push_back(light_clump_ptr);
    
    // 创建重粒子模板（通过质量实现密度2000）
    std::vector<std::shared_ptr<DEMClumpTemplate>> clump_types_heavy;
    float3 heavy_MOI = make_float3(1, 1, 1) * 0.4 * particle_mass_heavy * pow(particle_radius, 2);
    auto heavy_clump_ptr = DEMSim.LoadClumpType(
        particle_mass_heavy,    // 较大的质量（密度2000）
        heavy_MOI,
        std::vector<float>(1, particle_radius),
        std::vector<float3>(1, make_float3(0, 0, 0)),
        mat_type_particles      // 使用统一的材料类型
    );
    clump_types_heavy.push_back(heavy_clump_ptr);

    // --- 5.2 明确分层的粒子放置 ---

    // 设置泊松盘采样器
    PDSampler sampler(spacing);

    // =========================================================================
    // 6. 分批插入粒子准备（明确分层）
    // =========================================================================

    // 存储每批的数据
    std::vector<std::vector<std::shared_ptr<DEMClumpTemplate>>> batch_templates(num_batches);
    std::vector<std::vector<float3>> batch_positions(num_batches);
    std::vector<std::vector<float3>> batch_velocities(num_batches);
    std::vector<std::vector<unsigned int>> batch_families(num_batches);
    
    // 计算所需的初速度
    float required_initial_velocity = -abs(0.2f);  // 固定初速度为-0.2m/s

    std::cout << "计算的初速度: " << required_initial_velocity << " m/s" << std::endl;
    
    // 统计变量
    int total_generated = 0;
    int total_light_particles = 0;
    int total_heavy_particles = 0;
    float current_total_mass = 0.0f;

    std::cout << "\n开始为每批生成粒子位置（明确分层）..." << std::endl;

    // ========== 第一阶段：生成轻粒子（前light_batches批） ==========
    std::cout << "\n=== 第一阶段：生成轻粒子（密度1000，形成下层） ===" << std::endl;
    
    for (int batch = 0; batch < light_batches; batch++) {
        // 计算当前批次应包含的轻粒子数
        int batch_target = light_particles_per_batch;
        if (batch == light_batches - 1) {
            batch_target += light_particles_remainder;  // 最后一批包含余数
        }
        
        // 为当前批次生成粒子
        int batch_generated = 0;
        float current_z = insertion_height;  // 从固定高度开始
        
        while (batch_generated < batch_target) {
            float3 layer_center = make_float3(0, 0, current_z);

            // 在当前高度的圆柱区域内采样
            auto layer_xyz = sampler.SampleCylinderZ(layer_center, insertion_radius, 0);
            
            if (layer_xyz.empty()) {
                current_z += spacing;  // 向上移动一层
                continue;
            }
            
            // 添加当前层的轻粒子
            for (size_t i = 0; i < layer_xyz.size() && batch_generated < batch_target; i++) {
                batch_templates[batch].push_back(clump_types_light[0]);
                batch_positions[batch].push_back(layer_xyz[i]);
                batch_velocities[batch].push_back(make_float3(0, 0, required_initial_velocity));
                batch_families[batch].push_back(FAMILY_PARTICLES_LIGHT);
                
                batch_generated++;
                total_light_particles++;
                current_total_mass += particle_mass_light;
            }
            
            current_z += spacing;  // 向上移动到下一层
        }
        
        total_generated += batch_generated;
        
        std::cout << "轻粒子批次 " << batch + 1 << "/" << light_batches 
                  << " 生成完成: " << batch_generated << " 个粒子" << std::endl;
    }
    
    std::cout << "第一阶段完成: 共生成 " << total_light_particles << " 个轻粒子" << std::endl;

    // ========== 第二阶段：生成重粒子（后heavy_batches批） ==========
    std::cout << "\n=== 第二阶段：生成重粒子（密度2000，形成上层） ===" << std::endl;
    
    for (int batch = light_batches; batch < num_batches; batch++) {
        // 计算当前批次应包含的重粒子数
        int batch_index = batch - light_batches;
        int batch_target = heavy_particles_per_batch;
        if (batch_index == heavy_batches - 1) {
            batch_target += heavy_particles_remainder;  // 最后一批包含余数
        }
        
        // 为当前批次生成粒子
        int batch_generated = 0;
        float current_z = insertion_height;  // 从固定高度开始
        
        while (batch_generated < batch_target) {
            float3 layer_center = make_float3(0, 0, current_z);

            // 在当前高度的圆柱区域内采样
            auto layer_xyz = sampler.SampleCylinderZ(layer_center, insertion_radius, 0);
            
            if (layer_xyz.empty()) {
                current_z += spacing;  // 向上移动一层
                continue;
            }
            
            // 添加当前层的重粒子
            for (size_t i = 0; i < layer_xyz.size() && batch_generated < batch_target; i++) {
                batch_templates[batch].push_back(clump_types_heavy[0]);
                batch_positions[batch].push_back(layer_xyz[i]);
                batch_velocities[batch].push_back(make_float3(0, 0, required_initial_velocity));
                batch_families[batch].push_back(FAMILY_PARTICLES_HEAVY);
                
                batch_generated++;
                total_heavy_particles++;
                current_total_mass += particle_mass_heavy;
            }
            
            current_z += spacing;  // 向上移动到下一层
        }
        
        total_generated += batch_generated;
        
        std::cout << "重粒子批次 " << batch_index + 1 << "/" << heavy_batches 
                  << " 生成完成: " << batch_generated << " 个粒子" << std::endl;
    }
    
    std::cout << "第二阶段完成: 共生成 " << total_heavy_particles << " 个重粒子" << std::endl;

    // 输出最终统计信息
    std::cout << "\n=== 粒子生成准备完成（明确分层倒置配置） ===" << std::endl;
    std::cout << "实际准备粒子数: " << total_generated << std::endl;
    std::cout << "  - 轻粒子数（密度1000，前" << light_batches << "批）: " << total_light_particles << std::endl;
    std::cout << "  - 重粒子数（密度2000，后" << heavy_batches << "批）: " << total_heavy_particles << std::endl;
    std::cout << "轻粒子总质量: " << total_light_particles * particle_mass_light << " kg" << std::endl;
    std::cout << "重粒子总质量: " << total_heavy_particles * particle_mass_heavy << " kg" << std::endl;
    std::cout << "实际总质量: " << current_total_mass << " kg" << std::endl;
    std::cout << "实际质量比（轻:重）: " << (total_light_particles * particle_mass_light) / (total_heavy_particles * particle_mass_heavy) << ":1" << std::endl;
    std::cout << "实际粒子数比（轻:重）: " << (float)total_light_particles / total_heavy_particles << ":1" << std::endl;
    
    float mass_error = abs((total_light_particles * particle_mass_light - total_heavy_particles * particle_mass_heavy)) 
                      / (current_total_mass) * 100.0f;
    std::cout << "质量平衡误差: " << mass_error << "%" << std::endl;
    
    std::cout << "\n重要说明: " << std::endl;
    std::cout << "1. 轻粒子（密度1000）将在前" << light_batches << "批中插入，形成下层" << std::endl;
    std::cout << "2. 重粒子（密度2000）将在后" << heavy_batches << "批中插入，形成上层" << std::endl;
    std::cout << "3. 这是一个不稳定配置（重在上，轻在下），用于研究混合过程" << std::endl;
    std::cout << "================" << std::endl;

    // =========================================================================
    // 7. SIMULATION INITIALIZATION 
    // =========================================================================
    
    // Initialize the simulation system
    DEMSim.Initialize();
    
    // =========================================================================
    // 8. OUTPUT SETUP 
    // =========================================================================
    
    std::ostringstream dir_name;
    // dir_name << "SUPMixerInverted_ClearLayers_f" << my_scale_factor 
    //          << "_light1000x" << light_batches << "batch"
    //          << "_heavy2000x" << heavy_batches << "batch";

    dir_name << "SUPMixerOutput_ClearLayers_f" << my_scale_factor 
             << "se" << std::setw(3) << std::setfill('0') << static_cast<int>(particle_cohesion * 100);
    
    path out_dir = current_path() / dir_name.str();
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
    int frames_total = num_batches * frames_between_batches + frames_after_insertion;

    // 主仿真循环
    for (int frame = 0; frame <= frames_total; frame++) {
        // 在前N帧分批插入粒子
        if (frame % frames_between_batches == 0 && current_batch < num_batches) {
            // 跳过空批次（额外的安全检查）
            if (batch_templates[current_batch].empty()) {
                current_batch++;
                continue;
            }
            
            // 获取当前最高点
            float current_max_z = 0.0f;
            if (current_batch > 0) {
                auto max_z_inspector = DEMSim.CreateInspector("clump_max_z");
                current_max_z = max_z_inspector->GetValue();
            }
            
            // 动态调整插入高度：在当前最高点上方留出安全距离
            float safety_margin = 5.0f * particle_diameter;  // 5倍粒径的安全距离
            float dynamic_insertion_height = std::max(insertion_height, 
                                                    current_max_z + safety_margin);
            
            // 更新所有粒子的Z坐标
            for (size_t i = 0; i < batch_positions[current_batch].size(); i++) {
                float z_offset = dynamic_insertion_height - insertion_height;
                batch_positions[current_batch][i].z += z_offset;
            }
            
            // 判断当前批次类型
            std::string batch_type = (current_batch < light_batches) ? "轻粒子" : "重粒子";
            int batch_phase_index = (current_batch < light_batches) ? 
                                    current_batch + 1 : 
                                    current_batch - light_batches + 1;
            int total_phase_batches = (current_batch < light_batches) ? 
                                      light_batches : 
                                      heavy_batches;
            
            std::cout << "\n=== 第 " << frame << " 帧：插入" << batch_type 
                     << "批次 " << batch_phase_index << "/" << total_phase_batches << " ===" << std::endl;
            
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
                std::cout << "成功插入 " << batch_particle_count << " 个" << batch_type << std::endl;
                std::cout << "动态插入高度: " << dynamic_insertion_height * 1000 << " mm" << std::endl;
                
                // 更新仿真系统
                DEMSim.UpdateClumps();
            }
            
            current_batch++;
        }

        // 每5帧输出一次文件
        if (frame % 5 == 0) {
            // 生成输出文件名
            char filename[200];
            sprintf(filename, "clear_layer_output_%04d.csv", frame);
            DEMSim.WriteClumpFile(out_dir / filename);
        
            // 进度报告
            std::cout << "帧: " << frame << " - 输出文件已保存 - 当前系统中粒子数: " << DEMSim.GetNumClumps() << std::endl;

        } else {
            // 非输出帧，只报告进度
            std::cout << "帧: " << frame << " - 当前系统中粒子数: " << DEMSim.GetNumClumps() << std::endl;
        }
        
        // 推进仿真
        DEMSim.DoDynamics(time_per_frame); // 每帧推进指定时间
        DEMSim.ShowThreadCollaborationStats();
    }

    // =========================================================================
    // 10. 后处理
    // =========================================================================

    // 性能统计
    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_sec = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    std::cout << "\n仿真完成, 耗时: " << time_sec.count() << " 秒" << std::endl;
    
    // 输出最终统计
    std::cout << "\n=== 最终统计（明确分层倒置配置） ===" << std::endl;
    std::cout << "最终粒子数: " << DEMSim.GetNumClumps() << std::endl;
    std::cout << "轻粒子数（前" << light_batches << "批插入）: " << total_light_particles << std::endl;
    std::cout << "重粒子数（后" << heavy_batches << "批插入）: " << total_heavy_particles << std::endl;
    std::cout << "最终系统总质量: " << current_total_mass << " kg" << std::endl;
    std::cout << "配置说明: 明确分层，轻粒子先插入形成下层，重粒子后插入形成上层" << std::endl;
    
    DEMSim.ShowTimingStats();
    
    // 清理资源
    DEMSim.ClearTimingStats();

    std::cout << "----------------------------------------" << std::endl;
    DEMSim.ShowMemStats();
    std::cout << "----------------------------------------" << std::endl;

    std::cout << "SUPdemo_Mixer_Part1 (ClearLayer) 退出..." << std::endl;
    return 0;
}
