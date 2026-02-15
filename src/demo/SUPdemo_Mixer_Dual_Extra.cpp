// =============================================================================
// SUP模型搅拌器仿真：断点续算版本
// 从mixer_output_XXXX.csv读取状态继续仿真
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
#include <map>
#include <algorithm>
#include <iomanip>

using namespace deme;
using namespace std::filesystem;

// =============================================================================
// 断点续算数据结构
// =============================================================================
struct ParticleState {
    float3 pos;           // 位置
    float3 vel;           // 线速度
    float3 angVel;        // 角速度
    float radius;         // 半径
    unsigned int family;  // Family ID
};

// =============================================================================
// 从mixer_output CSV文件读取粒子状态
// CSV格式: X,Y,Z,r,absv,v_x,v_y,v_z,w_x,w_y,w_z,family
// =============================================================================
std::vector<ParticleState> ReadCheckpointCSV(const std::string& filepath) {
    std::vector<ParticleState> particles;
    std::ifstream file(filepath);
    
    if (!file.is_open()) {
        std::cerr << "错误: 无法打开断点文件 " << filepath << std::endl;
        return particles;
    }
    
    std::string line;
    // 跳过标题行
    std::getline(file, line);
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        std::stringstream ss(line);
        std::string token;
        std::vector<float> values;
        
        while (std::getline(ss, token, ',')) {
            try {
                values.push_back(std::stof(token));
            } catch (...) {
                std::cerr << "警告: 解析数值失败: " << token << std::endl;
                values.push_back(0.0f);
            }
        }
        
        // 确保有足够的列: X,Y,Z,r,absv,v_x,v_y,v_z,w_x,w_y,w_z,family (12列)
        if (values.size() >= 12) {
            ParticleState p;
            p.pos = make_float3(values[0], values[1], values[2]);
            p.radius = values[3];
            // values[4] = absv (不需要，可从v_x,v_y,v_z计算)
            p.vel = make_float3(values[5], values[6], values[7]);
            p.angVel = make_float3(values[8], values[9], values[10]);
            p.family = static_cast<unsigned int>(values[11]);
            
            particles.push_back(p);
        }
    }
    
    file.close();
    std::cout << "从断点文件读取 " << particles.size() << " 个粒子" << std::endl;
    return particles;
}

// =============================================================================
// 查找最新的断点文件
// =============================================================================
std::pair<std::string, int> FindLatestCheckpoint(const std::string& output_dir, 
                                                   const std::string& prefix = "mixer_output_") {
    int max_frame = -1;
    std::string latest_file;
    
    for (const auto& entry : std::filesystem::directory_iterator(output_dir)) {
        std::string filename = entry.path().filename().string();
        
        // 检查是否匹配 mixer_output_XXXX.csv 格式
        if (filename.find(prefix) == 0 && filename.find(".csv") != std::string::npos) {
            // 提取帧号
            size_t start = prefix.length();
            size_t end = filename.find(".csv");
            std::string frame_str = filename.substr(start, end - start);
            
            try {
                int frame = std::stoi(frame_str);
                if (frame > max_frame) {
                    max_frame = frame;
                    latest_file = entry.path().string();
                }
            } catch (...) {
                continue;
            }
        }
    }
    
    return {latest_file, max_frame};
}

// =============================================================================
// 主函数
// =============================================================================
int main(int argc, char* argv[]) {
    // =========================================================================
    // 1. 断点续算参数
    // =========================================================================
    
    // 命令行参数或默认值
    std::string checkpoint_dir = "";
    int start_frame = 520;  // -1表示自动查找最新
    
    if (argc >= 2) {
        checkpoint_dir = argv[1];
    }
    if (argc >= 3) {
        start_frame = std::stoi(argv[2]);
    }
    
    // 如果未指定目录，使用默认
    if (checkpoint_dir.empty()) {
        checkpoint_dir = "/media/huyuze/Fanxiang/DEME_Data/backup251115/SUPMixerOutput_SizeDiff_f1se080";  // 根据实际修改
    }
    
    // =========================================================================
    // 2. SUP模型参数（与原仿真保持一致）
    // =====================================================s====================
    
    const float my_scale_factor = 1.0f;
    const float base_small_diameter = 0.0003f;
    const float base_large_diameter = 0.0005f;
    const float particle_density = 1000.0f;
    const float particle_cohesion = 0.8f;
    
    const float small_diameter = base_small_diameter * my_scale_factor;
    const float small_radius = small_diameter / 2.0f;
    const float large_diameter = base_large_diameter * my_scale_factor;
    const float large_radius = large_diameter / 2.0f;
    
    const float particle_mass_small = particle_density * (4.0/3.0) * 3.14159265359 * 
                                      pow(small_radius, 3);
    const float particle_mass_large = particle_density * (4.0/3.0) * 3.14159265359 * 
                                      pow(large_radius, 3);
    
    const float mixer_speed_rpm = 300.0f;
    const float mixer_angular_velocity = mixer_speed_rpm * 2.0f * 3.14159265359f / 60.0f;
    const float additional_time = 0.6f;  // 从断点继续运行的时长（根据需要修改）
    
    const float step_size = (my_scale_factor == 1.0f) ? 5e-7f : 
                            (my_scale_factor == 2.0f) ? 1e-6f : 2e-6f;
    const float time_per_frame = (my_scale_factor == 1.0f) ? 5e-3f : 1e-3f;
    
    const float domain_size = 0.1f;
    const float wall_radius = 0.042f;
    
    const unsigned int FAMILY_MIXER = 1;
    const unsigned int FAMILY_PARTICLES_SMALL = 2;
    const unsigned int FAMILY_PARTICLES_LARGE = 3;

    // =========================================================================
    // 3. 查找并加载断点
    // =========================================================================
    
    std::cout << "\n========== SUP Mixer 断点续算 ==========" << std::endl;
    
    std::string checkpoint_file;
    if (start_frame < 0) {
        auto [file, frame] = FindLatestCheckpoint(checkpoint_dir);
        checkpoint_file = file;
        start_frame = frame;
    } else {
        std::ostringstream oss;
        oss << checkpoint_dir << "/mixer_output_" 
            << std::setw(4) << std::setfill('0') << start_frame << ".csv";
        checkpoint_file = oss.str();
    }
    
    // =========================================================================
    // 注意：帧数选取建议
    // 搅拌器转速300 RPM，每帧5ms，每40帧完成一圈旋转回到初始角度
    // 建议选择 40 的整数倍帧号（40, 80, 120, 160, ...）作为断点
    // 这样搅拌器角度与初始状态一致，无需额外处理搅拌器旋转同步
    // =========================================================================
    
    if (checkpoint_file.empty() || start_frame < 0) {
        std::cerr << "错误: 未找到有效的断点文件" << std::endl;
        return 1;
    }
    
    std::cout << "断点文件: " << checkpoint_file << std::endl;
    std::cout << "起始帧号: " << start_frame << std::endl;
    
    // 计算起始时间
    float start_time = start_frame * time_per_frame;
    std::cout << "起始时间: " << start_time << " s" << std::endl;
    
    // 读取粒子状态
    std::vector<ParticleState> particles = ReadCheckpointCSV(checkpoint_file);
    if (particles.empty()) {
        std::cerr << "错误: 断点文件中没有粒子数据" << std::endl;
        return 1;
    }

    // =========================================================================
    // 4. 创建求解器并配置
    // =========================================================================
    
    DEMSolver DEMSim;
    DEMSim.SetIntegrator(TIME_INTEGRATOR::CENTERED_DIFFERENCE);
    DEMSim.SetVerbosity(INFO);
    DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV);
    DEMSim.SetOutputContent({"ABSV", "VEL", "ANG_VEL", "FAMILY"});
    DEMSim.SetContactOutputFormat(OUTPUT_FORMAT::CSV);
    DEMSim.SetContactOutputContent({"CNT_TYPE", "FORCE", "POINT"});

    // =========================================================================
    // 5. 材料属性（与原仿真一致）
    // =========================================================================
    
    auto mat_type_mixer = DEMSim.LoadMaterial({
        {"E", 1e8}, {"nu", 0.3}, {"CoR", 0.1}, {"mu", 0.3},
        {"Crr", 0.0}, {"Cohesion", 0}, {"scale_factor_l", my_scale_factor}
    });
    
    auto mat_type_particles = DEMSim.LoadMaterial({
        {"E", 1e8}, {"nu", 0.3}, {"CoR", 0.1}, {"mu", 0.3},
        {"Crr", 0.0}, {"Cohesion", particle_cohesion}, {"scale_factor_l", my_scale_factor}
    });
    
    DEMSim.SetMaterialPropertyPair("mu", mat_type_mixer, mat_type_particles, 0.3);
    DEMSim.SetMaterialPropertyPair("CoR", mat_type_mixer, mat_type_particles, 0.1);
    DEMSim.SetMaterialPropertyPair("Crr", mat_type_mixer, mat_type_particles, 0.0);
    DEMSim.SetMaterialPropertyPair("Cohesion", mat_type_mixer, mat_type_particles, 0.0);

    auto model_SUP = DEMSim.ReadContactForceModel("ForceModelSUP.cu");
    model_SUP->SetMustHaveMatProp({"E", "nu", "CoR", "mu", "Crr", "Cohesion", "scale_factor_l"});
    model_SUP->SetMustPairwiseMatProp({"CoR", "mu", "Crr", "Cohesion"});
    model_SUP->SetPerContactWildcards({"delta_time", "delta_tan_x", "delta_tan_y", "delta_tan_z"});

    // =========================================================================
    // 6. 仿真域设置
    // =========================================================================
    
    DEMSim.InstructBoxDomainDimension(
        {-domain_size/2, domain_size/2},
        {-domain_size/2, domain_size/2},
        {0, 0.2f}
    );
    DEMSim.InstructBoxDomainBoundingBC("top_open", mat_type_mixer);
    DEMSim.SetInitTimeStep(step_size);
    DEMSim.SetGravitationalAcceleration(make_float3(0, 0, -9.81));
    DEMSim.SetErrorOutVelocity(20000.0);
    DEMSim.SetErrorOutAvgContacts(150);
    DEMSim.SetCDUpdateFreq(40);
    DEMSim.SetExpandSafetyAdder(2.0);

    // =========================================================================
    // 7. 根据断点数据创建粒子
    // =========================================================================
    
    std::cout << "\n=== 从断点恢复粒子 ===" << std::endl;
    
    // 创建粒子模板
    auto small_sphere_template = DEMSim.LoadSphereType(
        particle_mass_small, small_radius, mat_type_particles);
    auto large_sphere_template = DEMSim.LoadSphereType(
        particle_mass_large, large_radius, mat_type_particles);
    
    // 准备粒子数据
    std::vector<float3> in_xyz;
    std::vector<float3> in_vel;
    std::vector<float3> in_angVel;
    std::vector<float4> in_quat;  // 使用单位四元数
    std::vector<std::shared_ptr<DEMClumpTemplate>> in_types;
    std::vector<unsigned int> in_families;
    
    int small_count = 0, large_count = 0;
    float radius_threshold = (small_radius + large_radius) / 2.0f;
    
    for (const auto& p : particles) {
        in_xyz.push_back(p.pos);
        in_vel.push_back(p.vel);
        in_angVel.push_back(p.angVel);
        in_quat.push_back(make_float4(1, 0, 0, 0));  // 单位四元数
        
        // 根据半径或family判断粒子类型
        bool is_small = (p.radius < radius_threshold) || (p.family == FAMILY_PARTICLES_SMALL);
        
        if (is_small) {
            in_types.push_back(small_sphere_template);
            in_families.push_back(FAMILY_PARTICLES_SMALL);
            small_count++;
        } else {
            in_types.push_back(large_sphere_template);
            in_families.push_back(FAMILY_PARTICLES_LARGE);
            large_count++;
        }
    }
    
    std::cout << "恢复粒子统计:" << std::endl;
    std::cout << "  - 小粒子: " << small_count << " 个" << std::endl;
    std::cout << "  - 大粒子: " << large_count << " 个" << std::endl;
    std::cout << "  - 总计: " << particles.size() << " 个" << std::endl;
    
    // 创建粒子批次并设置初始速度
    DEMClumpBatch particle_batch(in_xyz.size());
    particle_batch.SetTypes(in_types);
    particle_batch.SetPos(in_xyz);
    particle_batch.SetOriQ(in_quat);
    particle_batch.SetVel(in_vel);       // 设置初始线速度
    particle_batch.SetAngVel(in_angVel); // 设置初始角速度
    particle_batch.SetFamilies(in_families);
    
    DEMSim.AddClumps(particle_batch);
    
    std::cout << "成功恢复 " << in_xyz.size() << " 个粒子（含速度状态）" << std::endl;

    // =========================================================================
    // 8. 创建搅拌器（与原仿真一致）
    // =========================================================================
    
    auto walls = DEMSim.AddExternalObject();
    walls->AddCylinder(make_float3(0), make_float3(0, 0, 1), 
                      wall_radius, mat_type_mixer, 0);
    
    auto mixer = DEMSim.AddWavefrontMeshObject(
        (GET_DATA_PATH() / "mesh/impeller.obj").string(), mat_type_mixer);
    mixer->SetFamily(FAMILY_MIXER);
    
    std::cout << "搅拌器三角形网格数: " << mixer->GetNumTriangles() << std::endl;
    
    // 搅拌器从初始位置开始正常旋转（与原仿真一致）
    std::string omega_z_str = std::to_string(mixer_angular_velocity);
    DEMSim.SetFamilyPrescribedAngVel(FAMILY_MIXER, "0", "0", omega_z_str, true);
    DEMSim.DisableContactBetweenFamilies(FAMILY_MIXER, FAMILY_MIXER);
    
    auto mixer_tracker = DEMSim.Track(mixer);

    // =========================================================================
    // 9. 输出目录和Inspector
    // =========================================================================
    
    // 创建新的输出目录（与原脚本方式一致）
    std::ostringstream dir_oss;
    dir_oss << "SUPMixerOutput_SizeDiff_f" << (int)my_scale_factor 
            << "se" << std::setw(3) << std::setfill('0') << (int)(particle_cohesion * 100)
            << "_continue_from_" << start_frame;
    std::string dir_name = dir_oss.str();
    std::filesystem::path out_dir = std::filesystem::current_path() / dir_name;
    std::filesystem::create_directory(out_dir);
    
    std::cout << "输出目录: " << out_dir << std::endl;
    
    auto KE_inspector = DEMSim.CreateInspector("clump_kinetic_energy");
    auto max_v_inspector = DEMSim.CreateInspector("clump_max_absv");
    auto max_z_inspector = DEMSim.CreateInspector("clump_max_z");
    auto min_z_inspector = DEMSim.CreateInspector("clump_min_z");
    
    // 以追加模式打开能量历史文件
    std::ofstream energy_file(out_dir / "kinetic_energy_history.csv", std::ios::app);

    // =========================================================================
    // 10. 初始化并继续仿真
    // =========================================================================
    
    std::cout << "\n=== 初始化仿真系统 ===" << std::endl;
    DEMSim.Initialize();
    
    std::cout << "\n=== 从断点继续仿真 ===" << std::endl;
    std::cout << "起始时间: " << start_time << " s" << std::endl;
    float end_time = start_time + additional_time;
    std::cout << "目标时间: " << end_time << " s（继续运行 " << additional_time << " s）" << std::endl;
    
    unsigned int frame_steps = (unsigned int)(time_per_frame / step_size);
    unsigned int curr_step = 0;
    unsigned int frame_count = start_frame + 1;  // 从下一帧开始
    float current_time = start_time;
    
    std::chrono::high_resolution_clock::time_point start = 
        std::chrono::high_resolution_clock::now();
    
    while (current_time < end_time) {
        DEMSim.DoDynamics(step_size);
        
        if (curr_step % frame_steps == 0) {
            std::cout << "\n--- 帧 " << frame_count 
                      << " (t = " << current_time << " s) ---" << std::endl;
            
            float KE = KE_inspector->GetValue();
            float max_v = max_v_inspector->GetValue();
            float max_z = max_z_inspector->GetValue();
            float min_z = min_z_inspector->GetValue();
            float avg_contacts = DEMSim.GetAvgSphContacts();
            
            float mixer_torque = 0.0f;
            std::vector<float3> forces, points;
            mixer_tracker->GetContactForcesForAll(points, forces);
            for (size_t i = 0; i < points.size(); i++) {
                mixer_torque += points[i].x * forces[i].y - points[i].y * forces[i].x;
            }
            
            // 写入CSV
            energy_file << std::fixed << std::setprecision(6)
                       << current_time << "," << frame_count << ","
                       << KE << "," << max_v << "," << max_z << ","
                       << min_z << "," << avg_contacts << ","
                       << mixer_torque << "," << small_count << ","
                       << large_count << std::endl;
            energy_file.flush();
            
            std::cout << "动能: " << KE << " J, 最大速度: " << max_v << " m/s" << std::endl;
            
            // 输出文件
            char filename[200];
            sprintf(filename, "mixer_output_%04d.csv", frame_count);
            DEMSim.WriteSphereFile(out_dir / filename);
            
            sprintf(filename, "mixer_contact_%04d.csv", frame_count);
            DEMSim.WriteContactFile(out_dir / filename);
            
            sprintf(filename, "mixer_mesh_%04d.vtk", frame_count);
            DEMSim.WriteMeshFile(out_dir / filename);
            
            frame_count++;
        }
        
        curr_step++;
        current_time += step_size;
    }
    
    energy_file.close();
    
    std::chrono::high_resolution_clock::time_point end = 
        std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_sec = 
        std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    
    std::cout << "\n=== 断点续算完成 ===" << std::endl;
    std::cout << "续算运行时间: " << time_sec.count() << " 秒" << std::endl;
    std::cout << "仿真时间范围: " << start_time << " -> " << current_time << " s" << std::endl;
    
    DEMSim.ShowTimingStats();
    
    return 0;
}