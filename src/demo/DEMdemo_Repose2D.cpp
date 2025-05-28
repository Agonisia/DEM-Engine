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
    DEMSim.SetVerbosity(INFO);
    DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV);
    DEMSim.SetOutputContent({"XYZ", "QUAT", "VEL", "ANG_VEL"});
    DEMSim.SetContactOutputFormat(OUTPUT_FORMAT::CSV);
    DEMSim.SetContactOutputContent({"OWNER", "GEO_ID", "FORCE", "POINT", "TORQUE",});

    // If you don't need individual force information, then this option makes the solver run a bit faster.
    // DEMSim.SetNoForceRecord();

    // Set random seed for reproducibility
    srand(52);
    float my_scale_factor= 5; // define sup factor here, before time step
    DEMSim.SetErrorOutAvgContacts(200);

    // =========================================================================
    // 2. MATERIAL PROPERTIES
    // =========================================================================
    // Define materials (E, nu, CoR, mu, Crr...)
    auto mat_type_walls = DEMSim.LoadMaterial({{"E", 1e8}, {"nu", 0.3}, {"CoR", 0.3}, {"mu", 1}});
    auto mat_type_particles = DEMSim.LoadMaterial({{"E", 1e9}, {"nu", 0.3}, {"CoR", 0.7}, {"mu", 0.50}, {"Crr", 0.03}});
    
    // Set interface properties between walls and particles
    // Without this line, CoR between wall material and granular material would be 0.5 (average of the two)
    DEMSim.SetMaterialPropertyPair("CoR", mat_type_walls, mat_type_particles, 0.3);
    DEMSim.SetMaterialPropertyPair("Crr", mat_type_walls, mat_type_particles, 0.5);
    DEMSim.SetMaterialPropertyPair("mu", mat_type_walls, mat_type_particles, 0.5);

    // =========================================================================
    // 3. SIMULATION DOMAIN SETUP
    // =========================================================================
    
    // Define simulation dimensions and parameters
    float scaling = 5;                  // Scale factor for the simulation
    float step_size = 5e-6 * my_scale_factor;             // Time step size
    float funnel_bottom = 0.f;          // Z-coordinate of funnel bottom
    
    // Configure domain boundaries
    DEMSim.InstructBoxDomainDimension({-10, 10}, {-10, 10}, {funnel_bottom - 10.f, funnel_bottom + 20.f});
    DEMSim.InstructBoxDomainBoundingBC("top_open", mat_type_walls);
    
    // Basic physics settings
    DEMSim.SetInitTimeStep(step_size);
    DEMSim.SetGravitationalAcceleration(make_float3(0, 0, -9.81));
    
    // Max velocity setting (generally just for the solver's reference)
    // The solver won't take into account velocities larger than this when doing async contact detection
    DEMSim.SetMaxVelocity(25.);

    // =========================================================================
    // 4. GEOMETRY CREATION
    // =========================================================================
    
    // Add mesh-represented funnel (fixed by default)
    auto funnel = DEMSim.AddWavefrontMeshObject(GetDEMEDataFile("mesh/funnel.obj"), mat_type_walls);
    funnel->Scale(0.15);

    // =========================================================================
    // 5. PARTICLE CREATION
    // =========================================================================
    
    // Configure clump template generation parameters
    int num_template = 6;           // Total number of random clump templates to generate
    int min_sphere = 1;             // Minimum number of spheres per clump
    int max_sphere = 5;             // Maximum number of spheres per clump
    float min_rad = 0.01 * scaling; // Minimum radius of component spheres
    float max_rad = 0.02 * scaling; // Maximum radius of component spheres
    float min_relpos = -0.01 * scaling; // Minimum relative position of spheres in clump
    float max_relpos = 0.01 * scaling;  // Maximum relative position of spheres in clump

    // Set up SUP force model
    auto model_SUP = DEMSim.ReadContactForceModel("ForceModelSUP.cu");
    model_SUP->SetMustHaveMatProp({"E", "nu", "CoR", "mu", "Crr"});
    model_SUP->SetMustPairwiseMatProp({"CoR", "mu", "Crr"});
    model_SUP->SetPerContactWildcards({"delta_time", "delta_tan_x", "delta_tan_y", "delta_tan_z", "scale_factor_l"});


    // Define dimensions for particle filling
    float spacing = max_rad * 2.0;
    float fill_width = 5.f;
    float fill_height = 5.f * fill_width;
    float fill_bottom = funnel_bottom + fill_width + spacing + fill_height / 2.0;

    // Create array to store generated clump templates
    std::vector<std::shared_ptr<DEMClumpTemplate>> clump_types;

    // Generate random clump templates
    for (int i = 0; i < num_template; i++) {
        // Decide the number of spheres in this clump
        // int num_sphere = rand() % (max_sphere - min_sphere + 1) + 1;
        int num_sphere = 1;
        
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
                // Subsequent spheres are positioned relative to the seed
                tmp.x = ((float)rand() / RAND_MAX) * (max_relpos - min_relpos) + min_relpos;
                tmp.y = 0.0; // 2D constraint
                tmp.z = ((float)rand() / RAND_MAX) * (max_relpos - min_relpos) + min_relpos;
            }
            tmp += seed_pos;
            relPos.push_back(tmp);

            // Select a new seed position from one of the previously generated spheres
            int choose_from = rand() % (j + 1);
            seed_pos = relPos.at(choose_from);
        }

        // Create clump template and add to types array
        auto clump_ptr = DEMSim.LoadClumpType(mass, MOI, radii, relPos, mat_type_walls);
        clump_types.push_back(clump_ptr);
    }

    // Set up particle distribution using Poisson Disk Sampling
    PDSampler sampler(spacing);
    std::vector<std::shared_ptr<DEMClumpTemplate>> input_pile_template_type;
    std::vector<float3> input_pile_xyz;

    // Sample positions for clumps in a box region
    float3 sample_center = make_float3(0, 0, fill_bottom);
    float3 sample_size = make_float3(fill_width, 0, fill_height / 2);
    auto layer_xyz = sampler.SampleBox(sample_center, sample_size);
    
    // Assign clump types to sampled positions
    unsigned int num_clumps = layer_xyz.size();
    for (unsigned int i = 0; i < num_clumps; i++) {
        input_pile_template_type.push_back(clump_types.at(i % num_template));
    }
    input_pile_xyz.insert(input_pile_xyz.end(), layer_xyz.begin(), layer_xyz.end());

    // Add clumps to the simulation (can be called multiple times before initialization)
    auto the_pile = DEMSim.AddClumps(input_pile_template_type, input_pile_xyz);

    // =========================================================================
    // 6. SIMULATION PARAMETERS CONFIGURATION
    // =========================================================================
    
    // Initialize the simulation system
    DEMSim.Initialize();
    DEMSim.SetFamilyContactWildcardValueAll(1, "scale_factor_l", my_scale_factor);

    // =========================================================================
    // 7. OUTPUT CONFIGURATION
    // =========================================================================
    
    // Create output directory
    path out_dir = current_path();
    out_dir += "/DemoOutput_Repose2D_SUP";
    create_directory(out_dir);

    // =========================================================================
    // 8. SIMULATION EXECUTION LOOP
    // =========================================================================
    
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();

    // Run simulation for 22 - 1 frames
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
