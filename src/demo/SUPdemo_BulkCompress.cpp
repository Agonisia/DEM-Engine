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

int main() {
    // =========================================================================
    // 1. 仿真设置
    // =========================================================================
    
    // === 核心参数（集中管理） ===
    const float my_scale_factor = 2.0f;                    // SUP缩放因子（需与原程序一致）
    const float particle_diameter = 0.001f;                // 粒子直径：1mm (0.5mm * 2)
    const float particle_radius = particle_diameter / 2.0f;
    const float particle_density = 1000.0f;                // 颗粒密度：1000 kg/m³
    const float particle_mass = particle_density * (4.0/3.0) * 3.14159265359 * pow(particle_radius, 3);
    
    // === 压缩参数 ===
    const float filter_height = 0.1f;                      // 过滤高度：0.1m
    const float plate_gap = 0.005f;                        // 压缩板与粒子初始间隙：5mm
    const float target_height = 0.09f;                     // 目标压缩高度：0.09m
    const float compression_velocity = 0.01f;              // 压缩速度：0.01 m/s
    const float plate_size = 0.040f;                       // 压缩板尺寸：40mm × 40mm
    const float plate_area = plate_size * plate_size;      // 压缩板面积
    
    // === 时间参数 ===
    const float step_size = 1e-6f;                         // 时间步长
    const float output_interval = 0.01f;                   // 输出间隔：0.01秒
    
    // === 输入文件路径 ===
    const std::string input_csv = "/home/peize/research/DEME_data/250722/results/f2c005/SUPdemo_f2c005_0150.csv";
    
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

    // 创建输出目录
    path out_dir = current_path();
    out_dir += "/CompressionTest_output";
    create_directory(out_dir);

    // =========================================================================
    // 2. 材料属性
    // =========================================================================

    // 定义压缩板材料
    auto mat_type_plate = DEMSim.LoadMaterial({
        {"E", 1e7},         // 杨氏模量
        {"nu", 0.3},        // 泊松比
        {"CoR", 0.1},       // 恢复系数
        {"mu", 0.3},        // 滑动摩擦系数
        {"Crr", 0},         // 滚动摩擦系数
        {"Cohesion", 0}     // 粘聚力
    });

    // 定义颗粒材料
    auto mat_type_particles = DEMSim.LoadMaterial({
        {"E", 1e7},         
        {"nu", 0.3},        
        {"CoR", 0.1},        
        {"mu", 0.3},
        {"Crr", 0},
        {"Cohesion", 0.05}
    });
    
    // 设置材料相互作用属性 - 板与颗粒
    DEMSim.SetMaterialPropertyPair("CoR", mat_type_plate, mat_type_particles, 0.1); 
    DEMSim.SetMaterialPropertyPair("mu", mat_type_plate, mat_type_particles, 0.3);
    DEMSim.SetMaterialPropertyPair("Crr", mat_type_plate, mat_type_particles, 0.0);
    DEMSim.SetMaterialPropertyPair("Cohesion", mat_type_plate, mat_type_particles, 0.0);

    // 加载SUP接触力模型
    auto model_SUP = DEMSim.ReadContactForceModel("ForceModelSUP.cu");
    model_SUP->SetMustHaveMatProp({"E", "nu", "CoR", "mu", "Crr", "Cohesion"});
    model_SUP->SetMustPairwiseMatProp({"CoR", "mu", "Crr", "Cohesion"});
    model_SUP->SetPerContactWildcards({"delta_time", "delta_tan_x", "delta_tan_y", "delta_tan_z", "scale_factor_l"});

    // =========================================================================
    // 3. 仿真域设置
    // =========================================================================

    // 设置仿真域大小
    DEMSim.InstructBoxDomainDimension(
        {-plate_size/2, plate_size/2},          // X方向
        {-plate_size/2, plate_size/2},          // Y方向
        {0, 0.2f}                               // Z方向：0-0.2m
    );

    // 设置边界条件
    DEMSim.InstructBoxDomainBoundingBC("top_open", mat_type_particles);

    // 设置物理参数
    DEMSim.SetInitTimeStep(step_size);
    DEMSim.SetGravitationalAcceleration(make_float3(0, 0, -9.81));
    DEMSim.SetMaxVelocity(25.0f);
    DEMSim.SetErrorOutVelocity(50.0f);

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
    
    // 过滤粒子（保留z <= filter_height的粒子）
    std::vector<float3> filtered_positions;
    std::vector<std::shared_ptr<DEMClumpTemplate>> filtered_types;
    
    // 创建粒子模板
    auto clump_template = DEMSim.LoadSphereType(
        particle_mass,
        particle_radius,
        mat_type_particles
    );
    
    // 统计变量
    int total_particles = all_positions.size();
    int retained_particles = 0;
    float max_z_after_filter = -1e6f;
    float min_z_after_filter = 1e6f;
    
    // 执行过滤
    for (const auto& pos : all_positions) {
        if (pos.z <= filter_height) {
            filtered_positions.push_back(pos);
            filtered_types.push_back(clump_template);
            retained_particles++;
            
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
    // 6. 创建压缩板（使用mesh）
    // =========================================================================

    // 加载压缩板mesh
    auto compression_plate = DEMSim.AddWavefrontMeshObject(
        (GET_DATA_PATH() / "mesh/plate40mm.obj").string(), 
        mat_type_plate
    );
    
    std::cout << "压缩板三角形数量: " << compression_plate->GetNumTriangles() << std::endl;
    
    // 设置压缩板族（用于控制接触）
    compression_plate->SetFamily(1);
    
    // 设置压缩板为固定族（我们将手动控制其位置）
    DEMSim.SetFamilyFixed(1);
    
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
    // 8. 系统初始化
    // =========================================================================
    
    std::cout << "\n=== 初始化仿真系统 ===" << std::endl;
    DEMSim.Initialize();
    
    // 设置SUP模型参数
    DEMSim.SetFamilyContactWildcardValueBoth(0, "scale_factor_l", my_scale_factor);

    // =========================================================================
    // 9. 设置压缩板初始位置
    // =========================================================================
    
    // 获取当前粒子堆最高点
    float initial_particle_height = max_z_finder->GetValue();
    float initial_max_v = max_v_finder->GetValue();
    
    std::cout << "\n=== 初始状态 ===" << std::endl;
    std::cout << "粒子最高点: " << initial_particle_height << " m" << std::endl;
    std::cout << "粒子最大速度: " << initial_max_v << " m/s" << std::endl;
    
    // 自动设置压缩板初始位置（在粒子最高点上方留出间隙）
    float initial_plate_height = initial_particle_height + plate_gap;
    plate_tracker->SetPos(make_float3(0, 0, initial_plate_height));
    
    std::cout << "压缩板初始位置: " << initial_plate_height << " m" << std::endl;
    std::cout << "初始间隙: " << plate_gap * 1000 << " mm" << std::endl;

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
    stress_strain_file << "Time,PlatePosition,Force,Stress,Strain,MaxVelocity\n";

    // =========================================================================
    // 11. 压缩过程
    // =========================================================================
    
    std::cout << "\n=== 开始压缩测试 ===" << std::endl;
    std::cout << "目标高度: " << target_height << " m" << std::endl;
    std::cout << "压缩距离: " << initial_plate_height - target_height << " m" << std::endl;
    std::cout << "压缩速度: " << compression_velocity << " m/s" << std::endl;
    std::cout << "预计时间: " << (initial_plate_height - target_height) / compression_velocity << " s" << std::endl;
    
    // 压缩循环变量
    unsigned int out_steps = (unsigned int)(output_interval / step_size);
    unsigned int curr_step = 0;
    unsigned int frame_count = 0;
    float current_time = 0.0f;
    float current_plate_z = initial_plate_height;
    
    // 主压缩循环
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    
    while (current_plate_z > target_height) {
        // 更新压缩板位置（手动控制）
        current_plate_z -= compression_velocity * step_size;
        plate_tracker->SetPos(make_float3(0, 0, current_plate_z));
        
        // 执行一步仿真
        DEMSim.DoDynamics(step_size);
        
        // 数据记录和输出
        if (curr_step % out_steps == 0) {
            // 获取当前状态
            float3 plate_pos = plate_tracker->Pos();
            float current_max_z = max_z_finder->GetValue();
            float current_max_v = max_v_finder->GetValue();
            float current_mass = mass_finder->GetValue();
            
            // 获取压缩板受力
            float3 plate_force = plate_tracker->ContactAcc() * plate_tracker->Mass();
            float compression_force = abs(plate_force.z);
            
            // 计算应力和应变
            float stress = compression_force / plate_area;
            float strain = (initial_particle_height - current_plate_z) / initial_particle_height;
            
            // 记录数据
            time_data.push_back(current_time);
            position_data.push_back(plate_pos.z);
            force_data.push_back(compression_force);
            stress_data.push_back(stress);
            strain_data.push_back(strain);
            velocity_data.push_back(current_max_v);
            
            // 写入文件
            stress_strain_file << current_time << ","
                              << plate_pos.z << ","
                              << compression_force << ","
                              << stress << ","
                              << strain << ","
                              << current_max_v << "\n";
            
            // 输出仿真文件
            if (frame_count % 5 == 0) {
                char filename[200];
                sprintf(filename, "compression_output_%04d.csv", frame_count);
                DEMSim.WriteUnifiedFile(out_dir / filename);
                
                // 同时输出mesh文件
                char meshfilename[200];
                sprintf(meshfilename, "compression_plate_%04d.vtk", frame_count);
                DEMSim.WriteMeshFile(out_dir / meshfilename);
            }
            
            // 进度报告
            std::cout << "时间: " << current_time << " s"
                     << " | 压缩板位置: " << plate_pos.z << " m"
                     << " | 应变: " << strain * 100 << " %"
                     << " | 应力: " << stress << " Pa"
                     << " | 最大速度: " << current_max_v << " m/s" << std::endl;
            
            frame_count++;
        }
        
        curr_step++;
        current_time += step_size;
    }
    
    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_sec = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    
    // =========================================================================
    // 12. 卸载过程
    // =========================================================================
    
    std::cout << "\n=== 开始卸载 ===" << std::endl;
    
    float unload_time = 0.5f;
    unsigned int unload_steps = (unsigned int)(unload_time / step_size);
    
    for (unsigned int i = 0; i < unload_steps; i++) {
        // 向上移动压缩板
        current_plate_z += compression_velocity * step_size;
        plate_tracker->SetPos(make_float3(0, 0, current_plate_z));
        
        // 执行仿真
        DEMSim.DoDynamics(step_size);
        
        // 定期输出状态
        if (i % out_steps == 0) {
            float3 plate_pos = plate_tracker->Pos();
            std::cout << "卸载中... 压缩板位置: " << plate_pos.z << " m" << std::endl;
            
            // 如果压缩板已经远离粒子，可以停止
            if (plate_pos.z > initial_particle_height + 0.02f) {
                std::cout << "压缩板已充分远离粒子，停止卸载" << std::endl;
                break;
            }
        }
    }
    
    // 让系统最终稳定
    DEMSim.DoDynamicsThenSync(0.1f);
    
    // =========================================================================
    // 13. 最终数据输出
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
    std::cout << "压缩时间: " << time_sec.count() << " 秒" << std::endl;
    std::cout << "最终应变: " << (initial_particle_height - target_height) / initial_particle_height * 100 << " %" << std::endl;
    
    // 找出最大应力
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
    
    // =========================================================================
    // 14. 清理和统计
    // =========================================================================
    
    DEMSim.ShowTimingStats();
    DEMSim.ClearTimingStats();
    
    std::cout << "----------------------------------------" << std::endl;
    DEMSim.ShowMemStats();
    std::cout << "----------------------------------------" << std::endl;
    
    std::cout << "压缩测试程序退出..." << std::endl;
    return 0;
}