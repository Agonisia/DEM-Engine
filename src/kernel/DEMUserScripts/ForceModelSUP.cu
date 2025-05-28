// Based on FullHertzianForceModel.cu, modified for Scaled-Up Particle(SUP) model
// Reference: "Inter-particle torque scaling in coarse grained DEM with rolling resistance 
// and particle size distributions" - Hu et al., Powder Technology 438 (2024)

// Acquire scale_factor_l
float l = scale_factor_l;
// Ensure l is valid, default to 1.0 if not properly set or <= 0.
if (l < 1e-5f) { // Using a small epsilon
    l = 1.0f;
}

float overlap_s = overlapDepth;
// If no contact, clear history and exit
if (overlap_s > 0) {
    // ========================================================================
    // SUP Step 1: Input Scaling (Convert Scaled-up variables to Original variables)
    // Based on Table 1 from the paper
    // ========================================================================
    
    // Overlap scaling: δ_O = δ_S / l
    float overlap_o = overlap_s / l;

    if (overlap_o > 0.f) {
        // Material properties from the original model structure
        float E_cnt, G_cnt, CoR_cnt, mu_cnt, Crr_cnt;
        {
            float E_A_orig = E[bodyAMatType]; 
            float nu_A_orig = nu[bodyAMatType];
            float E_B_orig = E[bodyBMatType];
            float nu_B_orig = nu[bodyBMatType];
            matProxy2ContactParam<float>(E_cnt, G_cnt, E_A_orig, nu_A_orig, E_B_orig, nu_B_orig);
            CoR_cnt = CoR[bodyAMatType][bodyBMatType];
            mu_cnt = mu[bodyAMatType][bodyBMatType];
            Crr_cnt = Crr[bodyAMatType][bodyBMatType];
        }

        // Radius scaling: R_O = R_S / l
        float R_s_A = ARadius;
        float R_s_B = BRadius;
        float R_o_A = R_s_A / l;
        float R_o_B = R_s_B / l;

        // Contact point scaling: position_O = position_S / l
        float3 locCPA_s = locCPA;
        float3 locCPB_s = locCPB;
        float3 locCPA_o = locCPA_s / l;
        float3 locCPB_o = locCPB_s / l;

        // Rotational velocity scaling: ω_O = ω_S / l (CORRECTED from ω_S * l)
        // According to Table 1: ω_O = l·ω_S means ω_O = ω_S / l in the original scale
        float3 ARotVel_s = ARotVel;
        float3 BRotVel_s = BRotVel;
        float3 ARotVel_o = ARotVel_s / l;  // CORRECTED
        float3 BRotVel_o = BRotVel_s / l;  // CORRECTED

        // Calculate rotational velocities at contact points in original scale
        float3 rotVelCPA_o_local = cross(ARotVel_o, locCPA_o);
        float3 rotVelCPB_o_local = cross(BRotVel_o, locCPB_o);
        applyOriQToVector3<float, deme::oriQ_t>(rotVelCPA_o_local.x, rotVelCPA_o_local.y, rotVelCPA_o_local.z, AOriQ.w, AOriQ.x, AOriQ.y, AOriQ.z);
        applyOriQToVector3<float, deme::oriQ_t>(rotVelCPB_o_local.x, rotVelCPB_o_local.y, rotVelCPB_o_local.z, BOriQ.w, BOriQ.x, BOriQ.y, BOriQ.z);

        // Mass scaling: m_O = m_S / l³
        float mass_s_A = AOwnerMass;
        float mass_s_B = BOwnerMass;
        float mass_o_A = mass_s_A / (l * l * l);
        float mass_o_B = mass_s_B / (l * l * l);

        if (mass_o_A <= 0.f) {
            mass_o_A = 1e-12f
        };
        if (mass_o_B <= 0.f) {
            mass_o_B = 1e-12f
        };

        // Time step scaling: Δt_O = Δt_S / l
        float ts_s = ts;
        float ts_o = ts_s / l;

        // Tangential displacement history scaling
        float3 delta_tan_s = make_float3(delta_tan_x, delta_tan_y, delta_tan_z);
        float3 delta_tan_o = delta_tan_s / l;

        float delta_time_s = delta_time;
        float delta_time_o = delta_time_s / l;

        // Variables for force calculation at original scale
        float mass_eff_o, sqrt_Rd_o, beta_o;
        float R_star_o;
        float3 vrel_tan_o;

        // Initialize force components at original scale
        float3 F_normal_o_vec = make_float3(0.f, 0.f, 0.f);
        float3 F_tangential_o_vec = make_float3(0.f, 0.f, 0.f);
        float3 torque_only_force_o = make_float3(0.f, 0.f, 0.f);

        // ========================================================================
        // SUP Step 2: Calculate forces at original scale
        // ========================================================================
        {
            // Calculate relative velocity at original scale
            const float3 velB2A_o = (ALinVel + rotVelCPA_o_local) - (BLinVel + rotVelCPB_o_local);
            const float projection_o = dot(velB2A_o, B2A);
            vrel_tan_o = velB2A_o - projection_o * B2A;

            // Update tangential displacement history
            {
                delta_tan_o += ts_o * vrel_tan_o;
                const float disp_proj_o = dot(delta_tan_o, B2A);
                delta_tan_o -= disp_proj_o * B2A;
                delta_time_o += ts_o;
            }

            // Calculate effective mass
            mass_eff_o = (mass_o_A * mass_o_B) / (mass_o_A + mass_o_B);
            if (mass_o_A <= 1e-12f && mass_o_B <= 1e-12f) {
                mass_eff_o = 1e-6f;
            } else if (mass_o_A <= 1e-12f) {
                mass_eff_o = mass_o_B;
            } else if (mass_o_B <= 1e-12f) {
                mass_eff_o = mass_o_A;
            }

            // Calculate effective radius
            if (R_o_A <= 1.0e-12f && R_o_B <= 1.0e-12f) { R_star_o = 1.0e-6f; }
            else if (R_o_A <= 1.0e-12f) { R_star_o = R_o_B; }
            else if (R_o_B <= 1.0e-12f) { R_star_o = R_o_A; }
            else { R_star_o = (R_o_A * R_o_B) / (R_o_A + R_o_B); }

            // Calculate normal force using Hertzian contact model
            sqrt_Rd_o = sqrtf(overlap_o * R_star_o);
            const float Sn_o = 2.f * E_cnt * sqrt_Rd_o;

            const float loge_o = (CoR_cnt < DEME_TINY_FLOAT) ? logf(DEME_TINY_FLOAT) : logf(CoR_cnt);
            beta_o = loge_o / sqrtf(loge_o * loge_o + deme::PI_SQUARED);

            const float k_n_o = deme::TWO_OVER_THREE * Sn_o;
            const float gamma_n_o = deme::TWO_TIMES_SQRT_FIVE_OVER_SIX * beta_o * sqrtf(Sn_o * mass_eff_o);

            F_normal_o_vec = (k_n_o * overlap_o + gamma_n_o * projection_o) * B2A;
        }

        // Calculate rolling resistance if enabled
        if (Crr_cnt > 0.0f) {
            bool should_add_rolling_resistance_o = true;
            {
                const float R_eff_o = R_star_o;
                const float kn_simple_o = deme::FOUR_OVER_THREE * E_cnt * sqrtf(R_eff_o);
                const float gn_simple_o = -2.f * sqrtf(deme::FIVE_OVER_THREE * mass_eff_o * E_cnt) * beta_o * powf(R_eff_o, 0.25f);
                const float d_coeff_o = gn_simple_o / (2.f * sqrtf(kn_simple_o * mass_eff_o));

                if (d_coeff_o < 1.0f) {
                    float t_collision_o = deme::PI * sqrtf(mass_eff_o / (kn_simple_o * (1.f - d_coeff_o * d_coeff_o)));
                    if (delta_time_o <= t_collision_o) {
                        should_add_rolling_resistance_o = false;
                    }
                }
            }
            if (should_add_rolling_resistance_o) {
                const float3 v_rot_o_global = rotVelCPB_o_local - rotVelCPA_o_local;
                const float v_rot_o_mag = length(v_rot_o_global);
                if (v_rot_o_mag > DEME_TINY_FLOAT) {
                    torque_only_force_o = (v_rot_o_global / v_rot_o_mag) * (Crr_cnt * length(F_normal_o_vec));
                }
            }
        }

        // Calculate tangential force if friction is enabled
        if (mu_cnt > 0.0f) {
            const float kt_o = 8.f * G_cnt * sqrt_Rd_o;
            const float gt_o = -deme::TWO_TIMES_SQRT_FIVE_OVER_SIX * beta_o * sqrtf(mass_eff_o * kt_o);
            float3 tangent_force_trial_o = -kt_o * delta_tan_o - gt_o * vrel_tan_o;

            const float ft_o_mag_trial = length(tangent_force_trial_o);
            if (ft_o_mag_trial > DEME_TINY_FLOAT) {
                const float ft_max_o = length(F_normal_o_vec) * mu_cnt;
                if (ft_o_mag_trial > ft_max_o) {
                    // Sliding friction
                    F_tangential_o_vec = (ft_max_o / ft_o_mag_trial) * tangent_force_trial_o;
                    if (fabs(kt_o) > DEME_TINY_FLOAT) {
                        delta_tan_o = (F_tangential_o_vec + gt_o * vrel_tan_o) / (-kt_o);
                    } else {
                        delta_tan_o = make_float3(0.f,0.f,0.f);
                    }
                } else {
                    // Static friction
                    F_tangential_o_vec = tangent_force_trial_o;
                }
            } else {
                F_tangential_o_vec = make_float3(0.f, 0.f, 0.f);
            }
        }

        // ========================================================================
        // SUP Step 3: Scale forces back to scaled-up system
        // According to Eq. (23): F_IS = l²·F_IO
        // ========================================================================
        float l_sq = l * l;
        float3 F_total_o_vec = F_normal_o_vec + F_tangential_o_vec + torque_only_force_o;
        force = F_total_o_vec * l_sq;

        // Note: If torques are needed separately for rotation calculation,
        // they should also be scaled by l² according to Eq. (25): M_IS = l²·M_IO
        // Example (if needed):
        // float3 torque_o = cross(locCPA_o, F_total_o_vec);
        // torque = torque_o * l_sq;

        // ========================================================================
        // SUP Step 4: Update history variables in scaled-up system
        // ========================================================================
        delta_tan_s = delta_tan_o * l;
        delta_tan_x = delta_tan_s.x;
        delta_tan_y = delta_tan_s.y;
        delta_tan_z = delta_tan_s.z;

        delta_time_s = delta_time_o * l;
        delta_time = delta_time_s;
    } else { 
        // No contact at original scale
        delta_time = 0;
        delta_tan_x = 0;
        delta_tan_y = 0;
        delta_tan_z = 0;
    }
} else {
    // No contact at scaled-up scale
    delta_time = 0;
    delta_tan_x = 0;
    delta_tan_y = 0;
    delta_tan_z = 0;
}