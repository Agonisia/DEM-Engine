// 使用plane_20by20.obj作为底板的休止角实验

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
    float my_scale_factor = 1.0f;                      // 缩放因子（暂时设为1）
    float particle_diameter = 0.001f;                  // 1mm 直径颗粒
    float particle_radius = particle_diameter / 2.0f;  // 半径
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
    // 2. 材料属性 - 基于论文表1
    // =========================================================================
    
    // 定义底板材料
    auto mat_type_plate = DEMSim.LoadMaterial({
        {"E", 6.3e10},        // 杨氏模量：63 GPa
        {"nu", 0.24},         // 泊松比：0.24
        {"CoR", 0.1},         // 恢复系数
        {"mu", 0.5},          // 滑动摩擦系数：0.5
        {"Crr", 0.05},        // 滚动摩擦系数：0.05
        {"Cohesion", 0}       // 与墙面无粘聚力
    });

    // 定义颗粒材料（玻璃珠）
    auto mat_type_particles = DEMSim.LoadMaterial({
        {"E", 6.3e10},        // 杨氏模量：63 GPa
        {"nu", 0.24},         // 泊松比：0.24
        {"CoR", 0.1},         // 恢复系数：0.1
        {"mu", 0.5},          // 滑动摩擦系数：0.5
        {"Crr", 0.05},        // 滚动摩擦系数：0.05
        {"Cohesion", 0.05}    // 表面能密度：0.05 J/m²
    });
    
    // 设置材料相互作用属性
    DEMSim.SetMaterialPropertyPair("CoR", mat_type_plate, mat_type_particles, 0.1);
    DEMSim.SetMaterialPropertyPair("Crr", mat_type_plate, mat_type_particles, 0.05);
    DEMSim.SetMaterialPropertyPair("mu", mat_type_plate, mat_type_particles, 0.5);
    DEMSim.SetMaterialPropertyPair("Cohesion", mat_type_plate, mat_type_particles, 0.0);

    // =========================================================================
    // 3. 仿真域设置
    // =========================================================================
    
    // 底板尺寸要求：56mm × 56mm
    float required_plate_size = 56.0f * particle_diameter;  // 56mm
    float current_plane_size = 20.0f;  // plane_20by20.obj 的尺寸
    
    // 计算需要的缩放因子
    float plate_scale_factor = required_plate_size / current_plane_size;
    std::cout << "底板缩放因子: " << plate_scale_factor << std::endl;
    
    // 域尺寸
    float domain_size = 100.0f * particle_diameter;
    float domain_height = 300.0f * particle_diameter;
    
    // 时间步长（考虑scale factor）
    float step_size = 5e-6 * my_scale_factor;
    
    // 设置仿真域大小
    DEMSim.InstructBoxDomainDimension(
        {-domain_size/2, domain_size/2}, 
        {-domain_size/2, domain_size/2}, 
        {-particle_diameter, domain_height}
    );
    
    // 设置边界条件 - 水平方向开放边界
    DEMSim.InstructBoxDomainBoundingBC("top_open", mat_type_plate);
    
    // 设置物理参数
    DEMSim.SetInitTimeStep(step_size);
    DEMSim.SetGravitationalAcceleration(make_float3(0, 0, -9.81));
    
    // 设置最大速度
    DEMSim.SetMaxVelocity(25.0f);

    // =========================================================================
    // 4. 使用plane_20by20.obj作为底板
    // =========================================================================
    
    // 加载平面网格
    auto plate = DEMSim.AddWavefrontMeshObject(GetDEMEDataFile("mesh/plane_20by20.obj"), mat_type_plate);
    
    // 缩放到所需尺寸
    plate->Scale(plate_scale_factor);
    
    // 设置为固定底板
    plate->SetFamily(FAMILY_PLATE);
    DEMSim.SetFamilyFixed(FAMILY_PLATE);

    // =========================================================================
    // 5. 颗粒生成 - 基于论文规格
    // =========================================================================

    // 创建球形颗粒模板
    float particle_density = 2500.0f;  // 玻璃密度 kg/m³
    float particle_mass = particle_density * (4.0f/3.0f) * 3.14159f * std::pow(particle_radius, 3);

    // 计算需要的颗粒数
    // 插入速率：3.35 g/s，持续10秒 = 33.5 g总质量
    float total_mass = 3.35e-3f * 10.0f;  // kg
    int total_particles = static_cast<int>(total_mass / particle_mass);
    
    std::cout << "总颗粒数：" << total_particles << std::endl;
    std::cout << "单颗粒质量：" << particle_mass << " kg" << std::endl;
    
    // 定义插入区域（14 × 14 × 260 mm³）
    float insertion_width = 14.0f * particle_diameter;
    float insertion_height = 260.0f * particle_diameter;
    float insertion_bottom = 50.0f * particle_diameter;  // 底板上方50mm
    
    // 在插入区域内随机生成颗粒
    std::vector<float3> particle_positions;
    
    // 使用泊松盘采样获得更均匀的分布
    PDSampler sampler(particle_diameter * 2.2f);  // 间距稍大于直径
    
    // 分批次生成以模拟连续插入
    int batches = 20;
    int particles_per_batch = total_particles / batches;
    float batch_height = insertion_height / batches;
    
    for (int batch = 0; batch < batches; batch++) {
        float z_center = insertion_bottom + insertion_height/2 + 
                        (batch - batches/2.0f) * batch_height / 2.0f;
        
        // 在该批次高度范围内采样
        float3 sample_center = make_float3(0, 0, z_center);
        float3 sample_size = make_float3(insertion_width/2, insertion_width/2, batch_height/2);
        
        auto batch_positions = sampler.SampleBox(sample_center, sample_size);
        
        // 只保留在圆柱形区域内的颗粒
        for (auto& pos : batch_positions) {
            if (pos.x * pos.x + pos.y * pos.y <= (insertion_width/2) * (insertion_width/2)) {
                particle_positions.push_back(pos);
                if (particle_positions.size() >= total_particles) break;
            }
        }
        
        if (particle_positions.size() >= total_particles) break;
    }
    
    // 如果颗粒数不够，在顶部追加
    // while (particle_positions.size() < total_particles) {
    //     float x = ((float)rand() / RAND_MAX - 0.5f) * insertion_width;
    //     float y = ((float)rand() / RAND_MAX - 0.5f) * insertion_width;
    //     float z = insertion_bottom + insertion_height - particle_diameter * 2.0f;
        
    //     if (x*x + y*y <= (insertion_width/2) * (insertion_width/2)) {
    //         particle_positions.push_back(make_float3(x, y, z));
    //     }
    // }
    
    std::cout << "实际生成颗粒数：" << particle_positions.size() << std::endl;
    
    // 创建单球体的clump模板（考虑scale factor）
    std::vector<float> radii = {particle_radius * my_scale_factor};
    std::vector<float3> relPos = {make_float3(0, 0, 0)};
    float3 MOI = make_float3(2.0f/5.0f) * particle_mass * particle_radius * particle_radius;

    auto sphere_template = DEMSim.LoadClumpType(
        particle_mass,
        MOI,
        radii,
        relPos,
        mat_type_particles
    );
    
    // 添加所有颗粒
    std::vector<std::shared_ptr<DEMClumpTemplate>> templates(particle_positions.size(), sphere_template);
    auto particles = DEMSim.AddClumps(templates, particle_positions);
    particles->SetFamily(FAMILY_PARTICLES);
    
    // 设置初始向下速度（基于论文：0.2 m/s）
    particles->SetVel(make_float3(0, 0, -0.2f));

    // =========================================================================
    // 6. 仿真初始化
    // =========================================================================
    
    DEMSim.Initialize();

    // =========================================================================
    // 7. 输出设置
    // =========================================================================
    
    path out_dir = current_path();
    out_dir += "/SUP_ReposeAngle_Factor1";
    create_directory(out_dir);

    path sub_dir = out_dir / "my_output";
    create_directory(sub_dir);

    // =========================================================================
    // 8. 仿真执行
    // =========================================================================

    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();

    // 主仿真循环
    float sim_time = 15.0f;  // 总仿真时间15秒
    float output_interval = 0.5f;  // 每0.5秒输出一次
    int total_frames = static_cast<int>(sim_time / output_interval);

    // 创建检查器监测堆高
    auto max_z_inspector = DEMSim.CreateInspector("clump_max_z");
    auto min_z_inspector = DEMSim.CreateInspector("clump_min_z");

    for (int frame = 0; frame < total_frames; frame++) {
        
        // 生成输出文件名
        char unifiedfile[200], outputfile[200], meshfile[200];
        sprintf(unifiedfile, "%s/particles_%04d.csv", out_dir.c_str(), frame);
        sprintf(outputfile, "%s/output_%04d.vtk", sub_dir.c_str(), frame);
        sprintf(meshfile, "%s/plate_%04d.vtk", out_dir.c_str(), frame);
        
        // 写入输出文件
        DEMSim.WriteUnifiedFile(std::string(unifiedfile));
        DEMSim.WriteSphereAndContactFile(std::string(outputfile));
        DEMSim.WriteMeshFile(std::string(meshfile));
        
        // 获取堆积高度信息
        float max_z = max_z_inspector->GetValue();
        float min_z = min_z_inspector->GetValue();
        
        // 进度报告
        float current_time = frame * output_interval;
        std::cout << "帧: " << frame 
                  << ", 时间: " << current_time << " 秒"
                  << ", 堆高: " << (max_z - min_z) / particle_diameter << "d" 
                  << std::endl;
        
        // 推进仿真
        DEMSim.DoDynamics(output_interval);
        DEMSim.ShowThreadCollaborationStats();
    }

    // =========================================================================
    // 9. 后处理
    // =========================================================================

    // 性能统计
    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_sec = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    std::cout << time_sec.count() << " 秒（实际时间）完成仿真" << std::endl;

    DEMSim.ShowTimingStats();

    std::cout << "----------------------------------------" << std::endl;
    DEMSim.ShowMemStats();
    std::cout << "----------------------------------------" << std::endl;

    std::cout << "休止角实验完成。" << std::endl;

    return 0;
}