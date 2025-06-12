//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//
//	SPDX-License-Identifier: BSD-3-Clause

// =============================================================================
// This demo features a 2D repose angle test using SUP contact model. Particles 
// flow through a mesh-represented funnel and form a pile that has an apparent 
// angle. The motion is constrained to the XZ plane.
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
    // 1. SIMULATION SETUP
    // =========================================================================
    
    // Global parameters
    float my_scale_factor = 4;          // Scale factor for SUP model
    const unsigned int FAMILY_FUNNEL = 1;
    const unsigned int FAMILY_PARTICLES = 2;

    // Create solver instance
    DEMSolver DEMSim;
    
    // Set output format and content
    DEMSim.SetVerbosity(INFO);
    DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV);
    DEMSim.SetOutputContent({"VEL", "ANG_VEL"}); 
    DEMSim.SetContactOutputFormat(OUTPUT_FORMAT::CSV);
    DEMSim.SetContactOutputContent({"POINT", "FORCE", "TORQUE"});
    
    // Solver basic configuration
    // If you don't need individual force information, then this option makes the solver run a bit faster.
    // DEMSim.SetNoForceRecord();
    DEMSim.SetErrorOutAvgContacts(1000);
    
    // Set random seed for reproducibility
    srand(52);

    // =========================================================================
    // 2. MATERIAL PROPERTIES
    // =========================================================================
    
    // Define material types and their properties for SUP model
    auto mat_type_walls = DEMSim.LoadMaterial({
        {"E", 1e8}, {"nu", 0.3}, {"CoR", 0.3}, {"mu", 1}, {"Crr", 0},
        {"Cohesion", 0}
    });

    auto mat_type_particles = DEMSim.LoadMaterial({
        {"E", 1e9}, {"nu", 0.3}, {"CoR", 0.7}, {"mu", 0.50}, {"Crr", 0.03},
        {"Cohesion", 0.1}
    });
    
    // Set material interaction properties
    DEMSim.SetMaterialPropertyPair("CoR", mat_type_walls, mat_type_particles, 0.3);
    DEMSim.SetMaterialPropertyPair("Crr", mat_type_walls, mat_type_particles, 0.5);
    DEMSim.SetMaterialPropertyPair("mu", mat_type_walls, mat_type_particles, 0.5);
    DEMSim.SetMaterialPropertyPair("Cohesion", mat_type_walls, mat_type_particles, 0.0);
    
    // Load SUP contact force model
    auto model_SUP = DEMSim.ReadContactForceModel("ForceModelSUP.cu");
    model_SUP->SetMustHaveMatProp({"E", "nu", "CoR", "mu", "Crr", "Cohesion"});
    model_SUP->SetMustPairwiseMatProp({"CoR", "mu", "Crr", "Cohesion"});
    model_SUP->SetPerContactWildcards({"delta_time", "delta_tan_x", "delta_tan_y", "delta_tan_z", "scale_factor_l"});

    // =========================================================================
    // 3. SIMULATION DOMAIN
    // =========================================================================
    
    // Domain dimension parameters
    float funnel_bottom = 0.f;          // Z-coordinate of funnel bottom
    float step_size = 5e-6 * my_scale_factor;  // Time step size
    
    // Set simulation domain size
    DEMSim.InstructBoxDomainDimension({-10, 10}, {-10, 10}, {funnel_bottom - 10.f, funnel_bottom + 20.f});
    
    // Set boundary conditions
    DEMSim.InstructBoxDomainBoundingBC("top_open", mat_type_walls);
    
    // Set physical parameters
    DEMSim.SetInitTimeStep(step_size);
    DEMSim.SetGravitationalAcceleration(make_float3(0, 0, -9.81));
    
    // Set maximum velocity (for solver's reference in async contact detection)
    DEMSim.SetMaxVelocity(25.);

    // Set error out velocity to handle initial settling
    DEMSim.SetErrorOutVelocity(100.);

    // =========================================================================
    // 4. GEOMETRY CREATION
    // =========================================================================
    
    // Add fixed geometries
    auto funnel = DEMSim.AddWavefrontMeshObject(GetDEMEDataFile("mesh/funnel.obj"), mat_type_walls);
    funnel->Scale(0.15);

    funnel->SetFamily(FAMILY_FUNNEL);
    DEMSim.SetFamilyFixed(FAMILY_FUNNEL);

    // =========================================================================
    // 5. PARTICLE GENERATION
    // =========================================================================
    
    // --- 5.1 Particle/Clump Templates ---
    
    // Template generation parameters
    int num_template = 6;               // Total number of random clump templates
    int min_sphere = 1;                 // Minimum number of spheres per clump
    int max_sphere = 5;                 // Maximum number of spheres per clump
    float min_rad = 0.02 * my_scale_factor;     // Minimum radius of component spheres
    float max_rad = 0.02 * my_scale_factor;     // Maximum radius of component spheres
    float min_relpos = -0.01 * my_scale_factor; // Minimum relative position of spheres
    float max_relpos = 0.01 * my_scale_factor;  // Maximum relative position of spheres
    
    // Create array to store clump templates
    std::vector<std::shared_ptr<DEMClumpTemplate>> clump_types;
    
    // Generate random clump templates
    for (int i = 0; i < num_template; i++) {
        // For this demo, using single sphere clumps for 2D constraint
        int num_sphere = 1;
        
        // Define clump properties (all in SI units)
        float mass = 0.1 * (float)num_sphere * std::pow(my_scale_factor, 3);
        float3 MOI = make_float3(2e-5 * (float)num_sphere, 1.5e-5 * (float)num_sphere, 1.8e-5 * (float)num_sphere) *
                     50. * std::pow(my_scale_factor, 5);
        std::vector<float> radii;
        std::vector<float3> relPos;
        
        // Generate sphere configurations
        float3 seed_pos = make_float3(0);
        for (int j = 0; j < num_sphere; j++) {
            // Random radius
            radii.push_back(((float)rand() / RAND_MAX) * (max_rad - min_rad) + min_rad);
            
            // Position relative to seed (constrained to XZ plane)
            float3 tmp;
            if (j == 0) {
                tmp = make_float3(0, 0, 0);
            } else {
                tmp.x = ((float)rand() / RAND_MAX) * (max_relpos - min_relpos) + min_relpos;
                tmp.y = 0.0;  // 2D constraint - no Y displacement
                tmp.z = ((float)rand() / RAND_MAX) * (max_relpos - min_relpos) + min_relpos;
            }
            tmp += seed_pos;
            relPos.push_back(tmp);
            
            // Update seed position
            int choose_from = rand() % (j + 1);
            seed_pos = relPos.at(choose_from);
        }
        
        // Create and store clump template
        auto clump_ptr = DEMSim.LoadClumpType(mass, MOI, radii, relPos, mat_type_particles);
        clump_types.push_back(clump_ptr);
    }
    
    // --- 5.2 Particle Placement ---
    
    // Define particle filling region
    float spacing = max_rad * 2.0;
    float fill_width = 4.f;
    float fill_height = 2.f * fill_width;
    float fill_bottom = funnel_bottom + fill_width + spacing + fill_height / 2.0;
    
    // Set up Poisson Disk Sampling
    PDSampler sampler(spacing);
    std::vector<std::shared_ptr<DEMClumpTemplate>> input_pile_template_type;
    std::vector<float3> input_pile_xyz;
    
    // Sample positions in box region (2D constraint)
    float3 sample_center = make_float3(0, 0, fill_bottom);
    float3 sample_size = make_float3(fill_width, 0, fill_height / 2);  // Y size is 0 for 2D
    auto layer_xyz = sampler.SampleBox(sample_center, sample_size);
    
    // Assign clump types to positions
    unsigned int num_clumps = layer_xyz.size();
    for (unsigned int i = 0; i < num_clumps; i++) {
        input_pile_template_type.push_back(clump_types.at(i % num_template));
    }
    input_pile_xyz.insert(input_pile_xyz.end(), layer_xyz.begin(), layer_xyz.end());
    
    // Add clumps to simulation
    auto the_pile = DEMSim.AddClumps(input_pile_template_type, input_pile_xyz);
    the_pile->SetFamily(FAMILY_PARTICLES);

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
    out_dir += "/SUP_ReposeAngle2D";
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
        char filename[200], meshfile[200], contactfile[200], outputfile[200], unifiedfile[200];
        sprintf(filename, "%s/SUPdemo_sphere_%04d.csv", out_dir.c_str(), i);
        sprintf(contactfile, "%s/SUPdemo_contacts_%04d.csv", out_dir.c_str(), i);
        sprintf(unifiedfile, "%s/SUPdemo_unified_%04d.csv", sub_dir.c_str(), i);                
        sprintf(meshfile, "%s/SUPdemo_funnel_%04d.vtk", out_dir.c_str(), i);
        sprintf(outputfile, "%s/SUPdemo_output_%04d.vtk", sub_dir.c_str(), i);
        
        // 写入输出文件
        DEMSim.WriteSphereFile(std::string(filename));
        DEMSim.WriteMeshFile(std::string(meshfile));
        DEMSim.WriteContactFile(std::string(contactfile));
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