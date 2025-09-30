//  Copyright (c) 2021, University of Wisconsin - Madison
//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//
//	SPDX-License-Identifier: BSD-3-Clause

// =============================================================================
// 压缩测试程序：读取SUP模型输出的粒子堆，过滤高位粒子，并进行压缩测试
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
#include <sstream>
#include <vector>
#include <iomanip>

using namespace deme;
using namespace std::filesystem;

// =============================================================================
// 自定义CSV读取函数 - 只读取X,Y,Z列
// =============================================================================
std::vector<float3> ReadParticlePositions(const std::string& filename) {
    std::vector<float3> positions;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "错误：无法打开文件 " << filename << std::endl;
        return positions;
    }
    
    std::string line;
    
    // 跳过标题行
    std::getline(file, line);
    
    // 读取数据行
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string token;
        float3 pos;
        int col = 0;
        
        // 只读取前三列（X,Y,Z）
        while (std::getline(ss, token, ',') && col < 3) {
            try {
                float value = std::stof(token);
                switch(col) {
                    case 0: pos.x = value; break;
                    case 1: pos.y = value; break;
                    case 2: pos.z = value; break;
                }
            } catch (const std::exception& e) {
                std::cerr << "解析错误在第 " << positions.size() + 1 << " 行: " << e.what() << std::endl;
                continue;
            }
            col++;
        }
        
        if (col == 3) {  // 确保读取了完整的3个坐标
            positions.push_back(pos);
        }
    }
    
    file.close();
    return positions;
}

// === 辅助函数：格式化参数 ===
std::string formatCohesion(float cohesion) {
    char buffer[4];
    sprintf(buffer, "%03d", (int)(cohesion * 100));
    return std::string(buffer);
}

std::string generateOutputDirName(float scale, float cohesion, int force_idx) {
    std::ostringstream ss;
    ss << "CompressionOutput_f" << (int)scale 
    << "c" << formatCohesion(cohesion)
    << "_fi" << force_idx;
    return ss.str();
}

int main() {
    // =========================================================================
    // 1. 仿真设置
    // =========================================================================
    
    // === SUP模型核心参数（集中管理） ===
    const float my_scale_factor = 2.0f;                    // SUP缩放因子
    const int scale_force_index = 3;                       // 力的缩放指数
    const float cohesion_value = 0.20f;                     // 粘附力
    const int base_particle_count = 1600000;               // 基准粒子数量（缩放因子为1时）
    const float base_particle_diameter = 0.0005f;          // 基准粒子直径：0.5mm
    const float particle_density = 1000.0f;                // 颗粒密度：1000 kg/m³
    
    // === 根据缩放因子自动计算的参数 ===
    const float particle_diameter = base_particle_diameter * my_scale_factor;
    const float particle_radius = particle_diameter / 2.0f;
    const float particle_mass = particle_density * (4.0/3.0) * 3.14159265359 * pow(particle_radius, 3);
    const int scale_factor_cubed = (int)pow(my_scale_factor, 3);
    const int actual_particle_count = base_particle_count / scale_factor_cubed;
    const float total_mass = base_particle_count * particle_density * (4.0/3.0) * 3.14159265359 * pow(base_particle_diameter/2.0f, 3);
    
    // === 压缩参数 ===
    const float filter_height = 0.1f;                      // 过滤高度：0.1m
    const float target_mass_g = 100.368f;                  // 目标质量：100.368g（根据因子1的运行结果设置）
    const float plate_gap = 0.0005f;                       // 压缩板与粒子初始间隙：0.5mm
    const float target_height = 0.095f;                     // 目标压缩高度：0.095m
    const float compression_velocity = 0.02f;              // 压缩速度：0.01 m/s
    const float plate_size = 0.040f;                       // 压缩板尺寸：40mm × 40mm
    const float plate_area = plate_size * plate_size;      // 压缩板面积
    
    // === 时间参数 ===
    const float step_size = 5e-7f;                         // 时间步长
    const float output_interval = 0.01f;                   // 输出间隔：0.01秒
    
    // === 输入文件路径 ===
    const std::string input_csv = "/root/DEM-Engine/results/f1c000/SUPdemo_f1c000_0150.csv";

    // 输出SUP模型信息
    std::cout << "\n========== SUP模型参数 ==========" << std::endl;
    std::cout << "缩放因子: " << my_scale_factor << std::endl;
    std::cout << "基准粒子数 (l=1): " << base_particle_count << std::endl;
    std::cout << "实际粒子数上限: " << actual_particle_count << std::endl;
    std::cout << "缩放比例: 1:" << scale_factor_cubed << std::endl;
    std::cout << "粒子直径: " << particle_diameter * 1000 << " mm" << std::endl;
    std::cout << "单个粒子质量: " << particle_mass * 1000000 << " mg" << std::endl;
    std::cout << "系统总质量: " << total_mass * 1000 << " g" << std::endl;
    std::cout << "================================\n" << std::endl;

    // 创建求解器实例
    DEMSolver DEMSim;
    
    // 设置输出格式和内容
    DEMSim.SetVerbosity(INFO);
    DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV);
    DEMSim.SetOutputContent({"VEL", "ANG_VEL"}); 
    DEMSim.SetContactOutputFormat(OUTPUT_FORMAT::CSV);
    DEMSim.SetContactOutputContent({"POINT", "FORCE", "TORQUE"});
    DEMSim.SetMeshOutputFormat(MESH_FORMAT::VTK);
    
    // 求解器基本配置
    DEMSim.SetErrorOutAvgContacts(1000);

    // 使用辅助函数：
    path out_dir = current_path();
    out_dir /= generateOutputDirName(my_scale_factor, cohesion_value, scale_force_index);
    create_directory(out_dir);

    // =========================================================================
    // 2. 材料属性
    // =========================================================================

    // 定义压缩板材料
    auto mat_type_plate = DEMSim.LoadMaterial({
        {"E", 1e7},         // 杨氏模量
        {"nu", 0.3},        // 泊松比
        {"CoR", 0.1},       // 恢复系数
        {"mu", 0.3},          // 滑动摩擦系数
        {"Crr", 0},         // 滚动摩擦系数
        {"Cohesion", 0},    // 粘聚力
        {"scale_factor_l", my_scale_factor}
    });

    // 定义颗粒材料
    auto mat_type_particles = DEMSim.LoadMaterial({
        {"E", 1e7},         
        {"nu", 0.3},        
        {"CoR", 0.1},        
        {"mu", 0.3},
        {"Crr", 0},
        {"Cohesion", cohesion_value},
        {"scale_factor_l", my_scale_factor}
    });
    
    // 设置材料相互作用属性 - 板与颗粒
    DEMSim.SetMaterialPropertyPair("CoR", mat_type_plate, mat_type_particles, 0.1); 
    DEMSim.SetMaterialPropertyPair("mu", mat_type_plate, mat_type_particles, 0.0);
    DEMSim.SetMaterialPropertyPair("Crr", mat_type_plate, mat_type_particles, 0.0);
    DEMSim.SetMaterialPropertyPair("Cohesion", mat_type_plate, mat_type_particles, 0.0);

    // 加载SUP接触力模型
    auto model_SUP = DEMSim.ReadContactForceModel("ForceModelSUP.cu");
    model_SUP->SetMustHaveMatProp({"E", "nu", "CoR", "mu", "Crr", "Cohesion", "scale_factor_l"});
    model_SUP->SetMustPairwiseMatProp({"CoR", "mu", "Crr", "Cohesion"});
    model_SUP->SetPerContactWildcards({"delta_time", "delta_tan_x", "delta_tan_y", "delta_tan_z"});

    // =========================================================================
    // 3. 仿真域设置
    // =========================================================================

    // 设置仿真域大小
    const float bottom_z = 0.0f;
    DEMSim.InstructBoxDomainDimension(
        {-plate_size/2, plate_size/2},          // X方向
        {-plate_size/2, plate_size/2},          // Y方向
        {0, 0.2f}                               // Z方向：0-0.2m
    );

    // 设置边界条件
    DEMSim.InstructBoxDomainBoundingBC("top_open", mat_type_particles);
    DEMSim.AddBCPlane(make_float3(0, 0, bottom_z), make_float3(0, 0, 1), mat_type_plate);

    // 设置物理参数
    DEMSim.SetInitTimeStep(step_size);
    DEMSim.SetGravitationalAcceleration(make_float3(0, 0, -9.81));
    // DEMSim.SetMaxVelocity(25.0f);
    DEMSim.SetErrorOutVelocity(10000.0f);

    // =========================================================================
    // 4. 读取和过滤粒子数据
    // =========================================================================

    std::cout << "\n=== 读取粒子数据 ===" << std::endl;
    std::cout << "输入文件: " << input_csv << std::endl;

    // 使用自定义函数读取粒子位置
    auto all_positions = ReadParticlePositions(input_csv);

    if (all_positions.empty()) {
        std::cerr << "错误：未能从文件读取任何粒子数据！" << std::endl;
        return 1;
    }

    std::cout << "成功读取 " << all_positions.size() << " 个粒子位置" << std::endl;

    // 创建粒子模板
    auto clump_template = DEMSim.LoadSphereType(
        particle_mass,
        particle_radius,
        mat_type_particles
    );

    // 过滤粒子
    std::vector<float3> filtered_positions;
    std::vector<std::shared_ptr<DEMClumpTemplate>> filtered_types;

    // 统计变量
    int total_particles = all_positions.size();
    int retained_particles = 0;
    float max_z_after_filter = -1e6f;
    float min_z_after_filter = 1e6f;
    float retained_mass = 0.0f;

    // 决定过滤策略
    if (my_scale_factor == 1.0f) {
        // 因子1：使用高度过滤
        std::cout << "使用高度过滤: " << filter_height << " m" << std::endl;
        
        for (const auto& pos : all_positions) {
            if (pos.z <= filter_height && retained_particles < actual_particle_count) {
                filtered_positions.push_back(pos);
                filtered_types.push_back(clump_template);
                retained_particles++;
                retained_mass += particle_mass;
                
                // 更新高度统计
                max_z_after_filter = std::max(max_z_after_filter, pos.z);
                min_z_after_filter = std::min(min_z_after_filter, pos.z);
            }
        }
        
        // 输出因子1的总质量，供其他因子参考
        std::cout << "\n*** 因子1过滤后总质量: " << retained_mass * 1000 << " g ***" << std::endl;
        std::cout << "*** 建议将此值用于其他缩放因子的质量控制 ***\n" << std::endl;
        
    } else {
        // 非因子1：使用质量控制

        const float target_mass = target_mass_g / 1000.0f;  // 转换为kg用于内部计算
    
        std::cout << "缩放因子=" << my_scale_factor 
              << "，使用质量控制，目标质量: " << target_mass_g << " g" << std::endl;
        
        // 将粒子按z坐标排序（从低到高）
        std::vector<std::pair<float, float3>> sorted_particles;
        for (const auto& pos : all_positions) {
            sorted_particles.push_back({pos.z, pos});
        }
        std::sort(sorted_particles.begin(), sorted_particles.end());
        
        // 从底部开始累加，直到达到目标质量或粒子数量限制
        for (const auto& [z, pos] : sorted_particles) {
            if (retained_mass >= target_mass || retained_particles >= actual_particle_count) {
                break;
            }
            
            filtered_positions.push_back(pos);
            filtered_types.push_back(clump_template);
            retained_particles++;
            retained_mass += particle_mass;
            
            // 更新高度统计
            max_z_after_filter = std::max(max_z_after_filter, pos.z);
            min_z_after_filter = std::min(min_z_after_filter, pos.z);
        }
    }

    std::cout << "\n=== 过滤统计 ===" << std::endl;
    std::cout << "原始粒子数: " << total_particles << std::endl;
    std::cout << "过滤后粒子数: " << retained_particles << std::endl;
    std::cout << "移除粒子数: " << total_particles - retained_particles << std::endl;
    std::cout << "保留率: " << (retained_particles * 100.0 / total_particles) << "%" << std::endl;
    std::cout << "过滤后总质量: " << retained_mass * 1000 << " g" << std::endl;
    std::cout << "过滤后最高点: " << max_z_after_filter << " m" << std::endl;
    std::cout << "过滤后最低点: " << min_z_after_filter << " m" << std::endl;
    std::cout << "================" << std::endl;


    // =========================================================================
    // 5. 加载粒子到系统
    // =========================================================================

    if (retained_particles == 0) {
        std::cerr << "错误：过滤后没有粒子！" << std::endl;
        return 1;
    }

    // 创建粒子批次
    DEMClumpBatch particle_batch(retained_particles);
    particle_batch.SetTypes(filtered_types);
    particle_batch.SetPos(filtered_positions);
    particle_batch.SetFamily(0);  // 粒子族为0
    
    // 添加到系统
    auto particles_tracker = DEMSim.AddClumps(particle_batch);
    std::cout << "成功加载 " << retained_particles << " 个粒子到系统" << std::endl;

    // =========================================================================
    // 6. 创建压缩板（使用mesh作为可视化，但实际用BC plane）
    // =========================================================================

    // 从过滤后的粒子位置直接计算最高点
    float initial_particle_height = max_z_after_filter;
    
    std::cout << "\n=== 初始状态 ===" << std::endl;
    std::cout << "粒子最高点（从数据计算）: " << initial_particle_height << " m" << std::endl;
    
    // 计算压缩板初始位置
    float initial_plate_height = initial_particle_height + plate_gap;
    
    std::cout << "压缩板初始位置: " << initial_plate_height << " m" << std::endl;
    std::cout << "初始间隙: " << plate_gap * 1000 << " mm" << std::endl;
    
    // 加载压缩板BC
    auto compression_plate = DEMSim.AddBCPlane(  
        make_float3(0, 0, initial_plate_height),  // 平面位置  
        make_float3(0, 0, -1),                    // 平面法向量（向下）  
        mat_type_plate                            // 材料  
    );  

    // compression_plate->Scale(make_float3(1.1f, 1.1f, 1.0f));  

    // 设置压缩板族（用于控制接触）
    compression_plate->SetFamily(1);
    compression_plate->SetMass(1e-10);  // 设置压缩板质量为1e-5kg

    // 创建跟踪器以监控压缩板状态
    auto plate_tracker = DEMSim.Track(compression_plate); 

    // =========================================================================
    // 7. 创建检查器
    // =========================================================================

    // 最高点检查器
    auto max_z_finder = DEMSim.CreateInspector("clump_max_z");
    
    // 压缩区域内的质量检查器
    auto mass_finder = DEMSim.CreateInspector("clump_mass", 
        "return (abs(X) <= 0.02) && (abs(Y) <= 0.02) && (Z >= 0);");
    
    // 速度检查器（监控系统稳定性）
    auto max_v_finder = DEMSim.CreateInspector("clump_max_absv");

    // =========================================================================
    // 8. 计算初始粒子高度（在初始化前）
    // =========================================================================
    
    // 计算压缩时间
    float compression_time = (initial_plate_height - target_height) / compression_velocity;
    
    // 设置时间依赖的压缩板速度表达式（必须在Initialize之前）
    // 在压缩时间内速度为-0.01 m/s，之后速度为0
    std::string vel_expr = "(t < " + std::to_string(compression_time) + ") ? " + 
                          std::to_string(-compression_velocity) + " : 0";
    DEMSim.SetFamilyPrescribedLinVel(1, "0", "0", vel_expr);
    
    // =========================================================================
    // 9. 系统初始化
    // =========================================================================
    
    std::cout << "\n=== 初始化仿真系统 ===" << std::endl;
    DEMSim.Initialize();
    
    // 设置SUP模型参数

    // DEMSim.SetContactWildcardValue("scale_factor_l", my_scale_factor);
    // DEMSim.SetContactWildcardValue("scale_force_index", scale_force_index);

    // 设置压缩板初始位置
    plate_tracker->SetPos(make_float3(0, 0, initial_plate_height));
    
    // 验证粒子最高点（使用检查器）
    float verified_max_z = max_z_finder->GetValue();
    float initial_max_v = max_v_finder->GetValue();
    std::cout << "粒子最高点（检查器验证）: " << verified_max_z << " m" << std::endl;
    std::cout << "粒子最大速度: " << initial_max_v << " m/s" << std::endl;

    // =========================================================================
    // 10. 压缩测试数据记录准备
    // =========================================================================
    
    // 数据记录容器
    std::vector<float> time_data;
    std::vector<float> position_data;
    std::vector<float> force_data;
    std::vector<float> stress_data;
    std::vector<float> strain_data;
    std::vector<float> velocity_data;
    
    // 打开输出文件
    std::ofstream stress_strain_file(out_dir / "stress_strain_curve.csv");
    stress_strain_file << "TimeFromContact,PlatePosition,Force,Stress,Strain,MaxVelocity\n";
    
    // 接触检测变量
    bool contact_started = false;
    float contact_height = 0.0f;
    float contact_time = 0.0f;
    const float contact_force_threshold = 0.01f;  // 接触力阈值：0.01N

    // =========================================================================
    // 11. 压缩过程
    // =========================================================================
    
    std::cout << "\n=== 开始压缩测试 ===" << std::endl;
    std::cout << "目标高度: " << target_height << " m" << std::endl;
    std::cout << "压缩距离: " << initial_plate_height - target_height << " m" << std::endl;
    std::cout << "压缩速度: " << compression_velocity << " m/s" << std::endl;
    std::cout << "预计时间: " << compression_time << " s" << std::endl;
    
    // 压缩循环变量
    unsigned int out_steps = (unsigned int)(output_interval / step_size);
    unsigned int curr_step = 0;
    unsigned int frame_count = 0;
    float current_time = 0.0f;
    
    // 运行时间稍微长一点，确保完全压缩到位
    float total_run_time = compression_time + 0.05f;  // 多运行0.05秒
    
    // 主压缩循环（基于时间控制）
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    
    while (current_time < total_run_time) {
        // 执行一步仿真
        DEMSim.DoDynamics(step_size);
        // DEMSim.SetContactWildcardValue("scale_factor_l", my_scale_factor);
        // DEMSim.SetContactWildcardValue("scale_force_index", scale_force_index);
        
        // 数据记录和输出
        if (curr_step % out_steps == 0) {
            // 获取当前状态
            float3 plate_pos = plate_tracker->Pos();
            float current_max_z = max_z_finder->GetValue();
            float current_max_v = max_v_finder->GetValue();
            float current_mass = mass_finder->GetValue();
            
            // 获取压缩板受力
            float3 plate_force = plate_tracker->ContactAcc() * plate_tracker->Mass();
            // printf("acc: %f, %f, %f\n", 
            //     plate_tracker->ContactAcc().x, plate_tracker->ContactAcc().y, plate_tracker->ContactAcc().z);
            float compression_force = abs(plate_force.z);
            
            // 检测接触
            if (!contact_started && compression_force > contact_force_threshold) {
                contact_started = true;
                contact_height = plate_pos.z;
                contact_time = current_time;
                std::cout << "\n*** 检测到接触！***" << std::endl;
                std::cout << "接触时间: " << contact_time << " s" << std::endl;
                std::cout << "接触高度: " << contact_height << " m" << std::endl;
                std::cout << "接触力: " << compression_force << " N" << std::endl;
                std::cout << "******************\n" << std::endl;
            }
            
            // 计算应力和应变（从接触开始）
            float stress = compression_force / plate_area;
            float strain = 0.0f;
            if (contact_started) {
                strain = (contact_height - plate_pos.z) / initial_particle_height;
            }
            
            // // 只在接触后记录应力应变数据
            if (contact_started) {
                // 记录数据
                time_data.push_back(current_time - contact_time);  // 相对于接触时刻的时间
                position_data.push_back(plate_pos.z);
                force_data.push_back(compression_force);
                stress_data.push_back(stress);
                strain_data.push_back(strain);
                velocity_data.push_back(current_max_v);
                
                // 写入文件
                stress_strain_file << current_time - contact_time << ","
                                  << plate_pos.z << ","
                                  << compression_force << ","
                                  << stress << ","
                                  << strain << ","
                                  << current_max_v << "\n";
            }
            
            // 输出仿真文件
            if (frame_count % 5 == 0) {
                char filename[200];
                sprintf(filename, "compression_output_%04d.csv", frame_count);
                DEMSim.WriteSphereFile(out_dir / filename);

                // 写入接触力文件
                char contact_filename[200];
                sprintf(contact_filename, "compression_contact_%04d.csv", frame_count);
                DEMSim.WriteContactFile(out_dir / contact_filename);

                // 同时输出mesh文件
                char meshfilename[200];
                sprintf(meshfilename, "compression_plate_%04d.vtk", frame_count);
                DEMSim.WriteMeshFile(out_dir / meshfilename);
            }
            
            // 进度报告
            if (!contact_started) {
                std::cout << "时间: " << current_time << " s"
                         << " | 压缩板位置: " << plate_pos.z << " m"
                         << " | 接触力: " << compression_force << " N"
                         << " | 等待接触..." << std::endl;
            } else {
                std::cout << "时间: " << current_time << " s"
                         << " | 压缩板位置: " << plate_pos.z << " m"
                         << " | 应变: " << strain * 100 << " %"
                         << " | 应力: " << stress << " Pa"
                         << " | 最大速度: " << current_max_v << " m/s" << std::endl;
            }
            
            frame_count++;
        }
        
        curr_step++;
        current_time += step_size;
    }
    
    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_sec = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    
    // 让系统最终稳定
    DEMSim.DoDynamicsThenSync(0.05f);
    
    // =========================================================================
    // 12. 最终数据输出
    // =========================================================================
    
    stress_strain_file.close();
    
    // 输出最终状态
    char final_filename[200];
    sprintf(final_filename, "compression_final_state.csv");
    DEMSim.WriteUnifiedFile(out_dir / final_filename);
    
    // 输出最终mesh状态
    char final_mesh_filename[200];
    sprintf(final_mesh_filename, "compression_plate_final.vtk");
    DEMSim.WriteMeshFile(out_dir / final_mesh_filename);
    
    // 输出统计信息
    std::cout << "\n=== 压缩测试完成 ===" << std::endl;
    std::cout << "总运行时间: " << time_sec.count() << " 秒" << std::endl;
    
    if (contact_started) {
        std::cout << "接触时间: " << contact_time << " 秒" << std::endl;
        std::cout << "实际压缩时间: " << (total_run_time - contact_time) << " 秒" << std::endl;
        std::cout << "最终应变: " << (contact_height - target_height) / initial_particle_height * 100 << " %" << std::endl;
        
        // 找出最大应力
        if (!stress_data.empty()) {
            float max_stress = 0.0f;
            size_t max_stress_index = 0;
            for (size_t i = 0; i < stress_data.size(); i++) {
                if (stress_data[i] > max_stress) {
                    max_stress = stress_data[i];
                    max_stress_index = i;
                }
            }
            std::cout << "最大应力: " << max_stress << " Pa (在应变 " 
                      << strain_data[max_stress_index] * 100 << "% 时)" << std::endl;
            
            // 绘制简单的应力应变曲线（控制台输出）
            std::cout << "\n=== 应力应变曲线预览 ===" << std::endl;
            std::cout << "应变(%) | 应力(Pa)" << std::endl;
            std::cout << "--------|----------" << std::endl;
            
            size_t step = std::max(size_t(1), strain_data.size() / 10);
            for (size_t i = 0; i < strain_data.size(); i += step) {
                std::cout << std::fixed << std::setprecision(1) 
                          << strain_data[i] * 100 << "%   | " 
                          << std::scientific << stress_data[i] << std::endl;
            }
        }
    } else {
        std::cout << "警告：未检测到接触！" << std::endl;
    }

    // =========================================================================
    // 13. 清理和统计
    // =========================================================================
    
    DEMSim.ShowTimingStats();
    DEMSim.ClearTimingStats();
    
    std::cout << "----------------------------------------" << std::endl;
    DEMSim.ShowMemStats();
    std::cout << "----------------------------------------" << std::endl;
    
    std::cout << "压缩测试程序退出..." << std::endl;
    return 0;
}