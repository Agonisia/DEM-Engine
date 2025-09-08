//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//
//	SPDX-License-Identifier: BSD-3-Clause

// =============================================================================
// SUP模型搅拌器仿真：准备工作
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
    const float my_scale_factor = 4.0f;                    // SUP缩放因子
    const float scale_force_index = 2.0f;                  // 力的缩放指数
    const int base_particle_count = 1600000;               // 基准粒子数量（缩放因子为1时）
    const float base_particle_diameter = 0.0005f;          // 基准粒子直径：0.5mm
    const float base_particle_weight = 0.0458f;            // 基准粒子重量：0.0458kg
    const float particle_density = 1000.0f;                // 颗粒密度：1000 kg/m³
    const float particle_cohesion = 0.0f;               // 颗粒间粘聚力

    // === 时间参数 ===
    const float step_size = 1e-6f;                         // 时间步长
    const float time_per_frame = 1e-2f;                    // 每帧仿真时间

    // === 分批插入参数 ===
    const int num_batches = 16;                             // 分批插入的批次数
    const int frames_between_batches = 5;                  // 每批次之间的帧数
    const int frames_after_insertion = 50;                // 插入完成后继续运行的帧数
    
    // === 根据缩放因子自动计算的参数 ===
    const float particle_diameter = base_particle_diameter * my_scale_factor;
    const float particle_radius = particle_diameter / 2.0f;
    const float particle_mass = particle_density * (4.0/3.0) * 3.14159265359 * 
                               pow(particle_radius, 3);
    const int scale_factor_cubed = (int)pow(my_scale_factor, 3);
    const int actual_particle_count = base_particle_count / scale_factor_cubed;
    const float total_mass = base_particle_count * particle_density * (4.0/3.0) * 
                            3.14159265359 * pow(base_particle_diameter/2.0f, 3);
    const int particles_per_batch = actual_particle_count / num_batches;
    const int remaining_particles = actual_particle_count % num_batches;
    const float batch_interval_time = frames_between_batches * time_per_frame;

    // === 仿真域参数 ===
    const float domain_size = 0.1f;                      // 100mm × 100mm 水平尺寸
    const float wall_radius = 0.04f;                      // 圆柱墙体半径40mm
    const float plate_bottom = 0.0f;                       // Z坐标：底板位置
    const float insertion_radius = 0.038f;                 // 插入区域半径38mm
    const float insertion_height = 0.042f;                 // 起始高度 42mm
    const float spacing = particle_diameter * 1.5f;        // 粒子间距
    
    // === Family ID定义 ===
    const unsigned int FAMILY_MIXER = 1;
    const unsigned int FAMILY_PARTICLES = 2;
    
    // 输出SUP模型信息
    std::cout << "\n========== SUP模型参数 ==========" << std::endl;
    std::cout << "缩放因子 l: " << my_scale_factor << std::endl;
    std::cout << "基准粒子数 (l=1): " << base_particle_count << std::endl;
    std::cout << "实际粒子数: " << actual_particle_count << std::endl;
    std::cout << "粒子直径: " << particle_diameter * 1000 << " mm" << std::endl;
    std::cout << "单个粒子质量: " << particle_mass * 1000 << " g" << std::endl;
    std::cout << "系统总质量: " << total_mass << " kg" << std::endl;
    std::cout << "分批插入: " << num_batches << " 批" << std::endl;
    std::cout << "每批粒子数: " << particles_per_batch << std::endl;
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
    
    // 定义慢搅拌器材料
    auto mat_type_mixer = DEMSim.LoadMaterial({
        {"E", 1e7},         // 杨氏模量
        {"nu", 0.3},        // 泊松比
        {"CoR", 0.1},       // 恢复系数
        {"mu", 0.0},        // 滑动摩擦系数
        {"Crr", 0.0},       // 滚动摩擦系数
        {"Cohesion", 0},    // 粘聚力
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
    
    // 设置材料相互作用属性 - 慢搅拌器与颗粒
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
    DEMSim.SetErrorOutAvgContacts(50);

    // 接触检测和性能设置
    DEMSim.SetCDUpdateFreq(40);
    DEMSim.SetExpandSafetyAdder(2.0);  // 考虑慢搅拌器叶片的高速运动
    
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
    
    // === 4.2 添加可动几何体 - 慢搅拌器网格 ===
    auto mixer = DEMSim.AddWavefrontMeshObject(
        (GET_DATA_PATH() / "mesh/impeller.obj").string(), 
        mat_type_mixer
    );
    std::cout << "搅拌器三角形网格数: " << mixer->GetNumTriangles() << std::endl;

    // 设置搅拌器的预定义角速度（绕z轴旋转，初始不旋转）
    DEMSim.SetFamilyPrescribedAngVel(FAMILY_MIXER, "0", "0", "0");

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
    
    // 选择控制模式
    ControlMode control_mode = ControlMode::BY_COUNT;
    
    // 控制参数（自动使用SUP计算的值）
    int target_particle_count = actual_particle_count;
    float target_total_mass = total_mass;
    
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

    // 设置泊松盘采样器
    PDSampler sampler(spacing);

    // =========================================================================
    // 6. 分批插入粒子准备
    // =========================================================================

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
    float min_clearance = 5.0f * particle_diameter;  // 最小安全距离（5mm）
    float gravity = 9.81f;

    // v0 = (s - 0.5*g*t²) / t
    float required_initial_velocity = (min_clearance - 0.5f * gravity * batch_interval_time * batch_interval_time) / batch_interval_time;
    required_initial_velocity = -abs(0.2f);  // 固定初速度为-0.2m/s

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
        float current_z = insertion_height;  // 从固定高度开始
        int batch_layers = 0;
        
        // std::cout << "\n批次 " << batch + 1 << " 开始生成..." << std::endl;
        
        while (batch_generated < batch_size) {
            float3 layer_center = make_float3(0, 0, current_z);

            // 在当前高度的圆柱区域内采样
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
            
            // std::cout << "  层 " << batch_layers << " (z=" << current_z/particle_diameter 
            //         << "d): 生成 " << layer_particles << " 个粒子" << std::endl;
            
            batch_layers++;
            current_z += spacing;  // 向上移动到下一层
        }
        
        total_generated += batch_generated;
        
        // std::cout << "批次 " << batch + 1 << " 完成: " << batch_generated 
        //         << " 个粒子，共 " << batch_layers << " 层" 
        //         << "，高度范围: " << insertion_height/particle_diameter << "d ~ " 
        //         << current_z/particle_diameter << "d" << std::endl;
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
            batch_velocities[batch].push_back(make_float3(0, 0, required_initial_velocity));  // 初始向下速度
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
    
    // =========================================================================
    // 8. OUTPUT SETUP 
    // =========================================================================
    
    std::ostringstream dir_name;
    dir_name << "SUPMixerOutput_f" << my_scale_factor 
         << "c" << std::setw(3) << std::setfill('0') << static_cast<int>(particle_cohesion * 100);
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
        // 在前90帧分批插入粒子
        if (frame % frames_between_batches == 0 && current_batch < num_batches) {
            // 获取当前最高点
            float current_max_z = 0.0f;
            if (current_batch > 0) {
                auto max_z_inspector = DEMSim.CreateInspector("clump_max_z");
                current_max_z = max_z_inspector->GetValue();
            }
            
            // 动态调整插入高度：在当前最高点上方留出安全距离
            float safety_margin = 5.0f * particle_diameter;  // 20倍粒径的安全距离
            float dynamic_insertion_height = std::max(insertion_height, 
                                                    current_max_z + safety_margin);
            
            // 更新所有粒子的Z坐标
            for (size_t i = 0; i < batch_positions[current_batch].size(); i++) {
                float z_offset = dynamic_insertion_height - insertion_height;
                batch_positions[current_batch][i].z += z_offset;
            }
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

        // 每10帧输出一次文件
        if (frame % 5 == 0) {
            // 生成输出文件名
            char filename[200];
            sprintf(filename, "compression_output_%04d.csv", frame);
            DEMSim.WriteSphereFile(out_dir / filename);


            // 同时输出mesh文件
            char meshfilename[200];
            sprintf(meshfilename, "compression_plate_%04d.vtk", frame);
            DEMSim.WriteMeshFile(out_dir / meshfilename);
        
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

    std::cout << "SUPdemo_Mixer 退出..." << std::endl;
    return 0;
}
