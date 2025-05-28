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
    
    // Create solver and configure basic settings
    DEMSolver DEMSim;
    DEMSim.UseFrictionalHertzianModel();
    DEMSim.SetVerbosity(INFO);
    DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV);
    DEMSim.SetOutputContent({"XYZ", "QUAT", "VEL", "ANG_VEL"});
    DEMSim.SetContactOutputFormat(OUTPUT_FORMAT::CSV);
    DEMSim.SetContactOutputContent({"OWNER", "GEO_ID", "FORCE", "POINT", "TORQUE"});

    // If you don't need individual force information, then this option makes the solver run a bit faster.
    // DEMSim.SetNoForceRecord();

    // Set random seed for reproducibility
    srand(42);

    // =========================================================================
    // 2. MATERIAL PROPERTIES
    // =========================================================================
    
    // Define materials (E, nu, CoR, mu...)
    auto mat_type_walls = DEMSim.LoadMaterial({{"E", 1e8}, {"nu", 0.3}, {"CoR", 0.3}, {"mu", 1}});
    auto mat_type_particles = DEMSim.LoadMaterial({{"E", 1e9}, {"nu", 0.3}, {"CoR", 0.7}, {"mu", 1}});
    
    // Set interface properties between walls and particles
    // Without this line, CoR between wall material and granular material would be 0.5 (average of the two)
    DEMSim.SetMaterialPropertyPair("CoR", mat_type_walls, mat_type_particles, 0.3);

    // =========================================================================
    // 3. SIMULATION DOMAIN SETUP
    // =========================================================================

    // =========================================================================
    // 3. GEOMETRY CREATION
    // =========================================================================
    
    // Add mesh-represented funnel (fixed by default)
    auto funnel = DEMSim.AddWavefrontMeshObject(GetDEMEDataFile("mesh/funnel.obj"), mat_type_walls);
    funnel->Scale(0.15);
    float funnel_bottom = 0.f;  // Z-coordinate of funnel bottom

    /* Commented out ground creation
    // First create clump type 0 for representing the ground
    float ground_sp_r = 0.02;
    auto template_ground = DEMSim.LoadSphereType(0.5, ground_sp_r, mat_type_walls);
    */

    // =========================================================================
    // 4. PARTICLE CREATION
    // =========================================================================
    
    // Define parameters for clump generation
    float scaling = 2;                  // Scale factor for the simulation
    int num_template = 6;               // Total number of random clump templates to generate
    int min_sphere = 1;                 // Minimum number of spheres per clump
    int max_sphere = 5;                 // Maximum number of spheres per clump
    float min_rad = 0.01 * scaling;     // Minimum radius of component spheres
    float max_rad = 0.02 * scaling;     // Maximum radius of component spheres
    float min_relpos = -0.01 * scaling; // Minimum relative position of spheres in clump
    float max_relpos = 0.01 * scaling;  // Maximum relative position of spheres in clump

    // Make an array to store these generated clump templates
    std::vector<std::shared_ptr<DEMClumpTemplate>> clump_types;

    // Generate random clump templates
    for (int i = 0; i < num_template; i++) {
        // First decide the number of spheres that live in this clump
        int num_sphere = rand() % (max_sphere - min_sphere + 1) + 1;

        // Allocate the clump template definition arrays (all in SI units)
        float mass = 0.1 * (float)num_sphere * std::pow(scaling, 3);
        float3 MOI = make_float3(2e-5 * (float)num_sphere, 1.5e-5 * (float)num_sphere, 1.8e-5 * (float)num_sphere) *
                     50. * std::pow(scaling, 5);
        std::vector<float> radii;
        std::vector<float3> relPos;

        // Randomly generate clump template configurations
        // The relPos of a sphere is always seeded from one of the already-generated spheres
        float3 seed_pos = make_float3(0);
        for (int j = 0; j < num_sphere; j++) {
            // Generate random radius within specified range
            radii.push_back(((float)rand() / RAND_MAX) * (max_rad - min_rad) + min_rad);
            
            // Generate position relative to seed position
            float3 tmp;
            if (j == 0) {
                // First sphere is at the origin of the clump
                tmp.x = 0;
                tmp.y = 0;
                tmp.z = 0;
            } else {
                // Subsequent spheres are positioned randomly relative to the seed
                tmp.x = ((float)rand() / RAND_MAX) * (max_relpos - min_relpos) + min_relpos;
                tmp.y = ((float)rand() / RAND_MAX) * (max_relpos - min_relpos) + min_relpos;
                tmp.z = ((float)rand() / RAND_MAX) * (max_relpos - min_relpos) + min_relpos;
            }
            tmp += seed_pos;
            relPos.push_back(tmp);

            // Select a new seed position from one of the previously generated spheres
            int choose_from = rand() % (j + 1);
            seed_pos = relPos.at(choose_from);
        }

        // Create clump template and add to types array
        // LoadClumpType returns a shared_ptr that points to this template so you may modify it
        // Material can be vector or a material shared ptr, which will be applied to all component spheres
        auto clump_ptr = DEMSim.LoadClumpType(mass, MOI, radii, relPos, mat_type_walls);
        clump_types.push_back(clump_ptr);
    }

    // Define parameters for particle fill region
    float spacing = 0.08 * scaling;     // Spacing between particles
    float fill_width = 5.f;             // Width of the fill region
    float fill_height = 2.f * fill_width; // Height of the fill region
    float fill_bottom = funnel_bottom + fill_width + spacing; // Bottom of the fill region
    
    // Initialize particle sampler with defined spacing
    PDSampler sampler(spacing);
    
    // Use a PDSampler-based clump generation process
    std::vector<std::shared_ptr<DEMClumpTemplate>> input_pile_template_type;
    std::vector<float3> input_pile_xyz;
    
    // Create particles in horizontal layers until reaching desired height
    float layer_z = 0;
    while (layer_z < fill_height) {
        // Center of the current layer
        float3 sample_center = make_float3(0, 0, fill_bottom + layer_z + spacing / 2);
        
        // Sample positions for clumps in a cylindrical region for this layer
        auto layer_xyz = sampler.SampleCylinderZ(sample_center, fill_width, 0);
        unsigned int num_clumps = layer_xyz.size();
        
        // Select from available clump types in a round-robin fashion
        for (unsigned int i = 0; i < num_clumps; i++) {
            input_pile_template_type.push_back(clump_types.at(i % num_template));
        }
        
        // Add this layer's positions to the overall pile
        input_pile_xyz.insert(input_pile_xyz.end(), layer_xyz.begin(), layer_xyz.end());
        
        // Move up to the next layer
        layer_z += spacing;
    }
    
    // Add all clumps to the simulation
    // Note: AddClumps can be called multiple times before initialization to add more clumps to the system
    auto the_pile = DEMSim.AddClumps(input_pile_template_type, input_pile_xyz);

    // =========================================================================
    // 5. SIMULATION DOMAIN SETUP
    // =========================================================================
    
    // Define simulation domain boundaries
    DEMSim.InstructBoxDomainDimension({-10, 10}, {-10, 10}, {funnel_bottom - 10.f, funnel_bottom + 20.f});
    DEMSim.InstructBoxDomainBoundingBC("top_open", mat_type_walls);
    
    // =========================================================================
    // 6. SIMULATION PARAMETERS CONFIGURATION
    // =========================================================================
    
    // Basic physics settings
    DEMSim.SetInitTimeStep(5e-6);
    DEMSim.SetGravitationalAcceleration(make_float3(0, 0, -9.81));
    
    // Max velocity setting for the solver's reference
    // The solver won't take into account velocities larger than this when doing async-ed contact detection
    DEMSim.SetMaxVelocity(25.);
    
    // Advanced bin size configuration (commented out)
    // You usually don't have to worry about initial bin size. In very rare cases, init bin size is so bad that auto bin
    // size adaption is ineffective, and you should notice in that case kT runs extremely slow. In that case setting
    // init bin size may save the simulation.
    // DEMSim.SetInitBinSize(min_rad * 6);
    
    // Initialize the simulation system
    DEMSim.Initialize();

    // =========================================================================
    // 7. OUTPUT CONFIGURATION
    // =========================================================================
    
    // Create output directory
    path out_dir = current_path();
    out_dir += "/DemoOutput_Repose_SUP";
    create_directory(out_dir);

    // =========================================================================
    // 8. SIMULATION EXECUTION LOOP
    // =========================================================================
    
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();

    // Run simulation for 22 -1 frames
    for (int i = 0; i < 22; i++) {
        // Generate output filenames
        char filename[200], meshfile[200], contactfile[200];;
        sprintf(filename, "%s/DEMdemo_output_%04d.csv", out_dir.c_str(), i);
        sprintf(meshfile, "%s/DEMdemo_funnel_%04d.vtk", out_dir.c_str(), i);
        sprintf(contactfile, "%s/DEMdemo_contacts_%04d.csv", out_dir.c_str(), i);
        
        // Write output files
        DEMSim.WriteSphereFile(std::string(filename));
        DEMSim.WriteMeshFile(std::string(meshfile));
        DEMSim.WriteContactFile(std::string(contactfile));
        
        // Report simulation progress
        std::cout << "Frame: " << i << std::endl;
        
        // Advance simulation by one frame (0.1 seconds)
        DEMSim.DoDynamics(1e-1);
        DEMSim.ShowThreadCollaborationStats();
    }

    // =========================================================================
    // 9. PERFORMANCE REPORTING
    // =========================================================================
    
    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_sec = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    std::cout << time_sec.count() << " seconds (wall time) to finish the simulation" << std::endl;

    DEMSim.ShowTimingStats();
    DEMSim.ClearTimingStats();

    std::cout << "DEMdemo_Repose exiting..." << std::endl;
    return 0;
}