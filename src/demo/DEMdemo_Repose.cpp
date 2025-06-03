//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//
//	SPDX-License-Identifier: BSD-3-Clause

// =============================================================================
// This demo features a repose angle test. Particles flow through a mesh-represented 
// funnel and form a pile that has an apparent angle.
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
    float scaling = 2;                  // Scale factor for geometry
    
    // Create solver instance
    DEMSolver DEMSim;
    
    // Set contact model
    DEMSim.UseFrictionalHertzianModel();
    
    // Set output format and content
    DEMSim.SetVerbosity(INFO);
    DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV);
    DEMSim.SetOutputContent({"XYZ", "QUAT", "VEL", "ANG_VEL"});
    DEMSim.SetContactOutputFormat(OUTPUT_FORMAT::CSV);
    DEMSim.SetContactOutputContent({"OWNER", "GEO_ID", "FORCE", "POINT", "TORQUE"});
    
    // Solver basic configuration
    // If you don't need individual force information, then this option makes the solver run a bit faster.
    // DEMSim.SetNoForceRecord();
    
    // Set random seed for reproducibility
    srand(42);

    // =========================================================================
    // 2. MATERIAL PROPERTIES 
    // =========================================================================
    
    // Define material types and their properties
    auto mat_type_walls = DEMSim.LoadMaterial({{"E", 1e8}, {"nu", 0.3}, {"CoR", 0.3}, {"mu", 1}});
    auto mat_type_particles = DEMSim.LoadMaterial({{"E", 1e9}, {"nu", 0.3}, {"CoR", 0.7}, {"mu", 1}});
    
    // Set material interaction properties
    // Without this line, CoR between wall material and granular material would be 0.5 (average of the two)
    DEMSim.SetMaterialPropertyPair("CoR", mat_type_walls, mat_type_particles, 0.3);

    // =========================================================================
    // 3. SIMULATION DOMAIN 
    // =========================================================================
    
    // Domain dimension parameters
    float funnel_bottom = 0.f;          // Z-coordinate of funnel bottom
    
    // Set simulation domain size
    DEMSim.InstructBoxDomainDimension({-10, 10}, {-10, 10}, {funnel_bottom - 10.f, funnel_bottom + 20.f});
    
    // Set boundary conditions
    DEMSim.InstructBoxDomainBoundingBC("top_open", mat_type_walls);
    
    // Set physical parameters
    DEMSim.SetInitTimeStep(5e-6);
    DEMSim.SetGravitationalAcceleration(make_float3(0, 0, -9.81));
    
    // Set maximum velocity (for solver's reference in async contact detection)
    DEMSim.SetMaxVelocity(25.);
    
    // Advanced bin size configuration (optional)
    // You usually don't have to worry about initial bin size. In very rare cases, init bin size is so bad that auto bin
    // size adaption is ineffective, and you should notice in that case kT runs extremely slow. In that case setting
    // init bin size may save the simulation.
    // DEMSim.SetInitBinSize(min_rad * 6);

    // =========================================================================
    // 4. GEOMETRY CREATION 
    // =========================================================================
    
    // Add fixed geometries
    auto funnel = DEMSim.AddWavefrontMeshObject(GetDEMEDataFile("mesh/funnel.obj"), mat_type_walls);
    funnel->Scale(0.15);
    
    /* Commented out ground creation
    // First create clump type 0 for representing the ground
    float ground_sp_r = 0.02;
    auto template_ground = DEMSim.LoadSphereType(0.5, ground_sp_r, mat_type_walls);
    */

    // =========================================================================
    // 5. PARTICLE GENERATION 
    // =========================================================================
    
    // --- 5.1 Particle/Clump Templates ---
    
    // Template generation parameters
    int num_template = 6;               // Total number of random clump templates
    int min_sphere = 1;                 // Minimum number of spheres per clump
    int max_sphere = 5;                 // Maximum number of spheres per clump
    float min_rad = 0.01 * scaling;     // Minimum radius of component spheres
    float max_rad = 0.02 * scaling;     // Maximum radius of component spheres
    float min_relpos = -0.01 * scaling; // Minimum relative position of spheres
    float max_relpos = 0.01 * scaling;  // Maximum relative position of spheres
    
    // Create array to store clump templates
    std::vector<std::shared_ptr<DEMClumpTemplate>> clump_types;
    
    // Generate random clump templates
    for (int i = 0; i < num_template; i++) {
        // Decide the number of spheres in this clump
        int num_sphere = rand() % (max_sphere - min_sphere + 1) + 1;
        
        // Define clump properties (all in SI units)
        float mass = 0.1 * (float)num_sphere * std::pow(scaling, 3);
        float3 MOI = make_float3(2e-5 * (float)num_sphere, 1.5e-5 * (float)num_sphere, 1.8e-5 * (float)num_sphere) *
                     50. * std::pow(scaling, 5);
        std::vector<float> radii;
        std::vector<float3> relPos;
        
        // Generate sphere configurations
        float3 seed_pos = make_float3(0);
        for (int j = 0; j < num_sphere; j++) {
            // Random radius
            radii.push_back(((float)rand() / RAND_MAX) * (max_rad - min_rad) + min_rad);
            
            // Position relative to seed
            float3 tmp;
            if (j == 0) {
                // First sphere at origin
                tmp = make_float3(0, 0, 0);
            } else {
                // Random position relative to seed
                tmp.x = ((float)rand() / RAND_MAX) * (max_relpos - min_relpos) + min_relpos;
                tmp.y = ((float)rand() / RAND_MAX) * (max_relpos - min_relpos) + min_relpos;
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
    float spacing = 0.08 * scaling;     // Spacing between particles
    float fill_width = 5.f;             // Width of the fill region
    float fill_height = 2.f * fill_width; // Height of the fill region
    float fill_bottom = funnel_bottom + fill_width + spacing; // Bottom of the fill region
    
    // Set up Poisson Disk Sampling
    PDSampler sampler(spacing);
    std::vector<std::shared_ptr<DEMClumpTemplate>> input_pile_template_type;
    std::vector<float3> input_pile_xyz;
    
    // Create particles in horizontal layers
    float layer_z = 0;
    while (layer_z < fill_height) {
        // Center of current layer
        float3 sample_center = make_float3(0, 0, fill_bottom + layer_z + spacing / 2);
        
        // Sample positions in cylindrical region
        auto layer_xyz = sampler.SampleCylinderZ(sample_center, fill_width, 0);
        unsigned int num_clumps = layer_xyz.size();
        
        // Assign clump types to positions
        for (unsigned int i = 0; i < num_clumps; i++) {
            input_pile_template_type.push_back(clump_types.at(i % num_template));
        }
        
        // Add layer positions to overall pile
        input_pile_xyz.insert(input_pile_xyz.end(), layer_xyz.begin(), layer_xyz.end());
        
        // Move to next layer
        layer_z += spacing;
    }
    
    // Add clumps to simulation
    auto the_pile = DEMSim.AddClumps(input_pile_template_type, input_pile_xyz);

    // =========================================================================
    // 6. SIMULATION INITIALIZATION 
    // =========================================================================
    
    // Initialize the simulation system
    DEMSim.Initialize();

    // =========================================================================
    // 7. OUTPUT SETUP 
    // =========================================================================
    
    // Create output directory
    path out_dir = current_path();
    out_dir += "/DemoOutput_Repose_SUP";
    create_directory(out_dir);

    // =========================================================================
    // 8. SIMULATION EXECUTION 
    // =========================================================================
    
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    
    // Main simulation loop
    for (int i = 0; i < 22; i++) {
        // Generate output filenames
        char filename[200], meshfile[200], contactfile[200];
        sprintf(filename, "%s/DEMdemo_output_%04d.csv", out_dir.c_str(), i);
        sprintf(meshfile, "%s/DEMdemo_funnel_%04d.vtk", out_dir.c_str(), i);
        sprintf(contactfile, "%s/DEMdemo_contacts_%04d.csv", out_dir.c_str(), i);
        
        // Write output files
        DEMSim.WriteSphereFile(std::string(filename));
        DEMSim.WriteMeshFile(std::string(meshfile));
        DEMSim.WriteContactFile(std::string(contactfile));
        
        // Progress report
        std::cout << "Frame: " << i << std::endl;
        
        // Advance simulation by 0.1 seconds
        DEMSim.DoDynamics(1e-1);
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
    DEMSim.ClearTimingStats();
    
    std::cout << "----------------------------------------" << std::endl;
    DEMSim.ShowMemStats();
    std::cout << "----------------------------------------" << std::endl;
    
    std::cout << "DEMdemo_Repose exiting..." << std::endl;
    return 0;
}