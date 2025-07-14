// 基于 FullHertzianForceModel.cu 修改，用于缩放颗粒（SUP）模型
// 参考文献：
// "Inter-particle torque scaling in coarse grained DEM with rolling resistance 
// and particle size distributions" - Hu et al., Powder Technology 438 (2024)
// "implified Johnson-Kendall-Roberts (SJKR) Contact Model - Implementation in PFC"
// - C.J. Coetzee, 2020

// VERSION_250603: 移除了内部时间步缩放
// VERSION_250604: 添加了JKR粘附力模型
// VERSION_250630：修改JKR模型为SJKR-F
// VERSION_250701：修改JKR模型为更简单的计算

// 获取缩放因子 l
float l = scale_factor_l;
// 确保 l 有效，如果未正确设置或 <= 0，则默认为 1.0
if (l < 1e-5f) { // 使用小的 epsilon 值
    l = 1.0f;
}

float overlap_s = overlapDepth;
// 限制重叠量不超过粒子半径的十分之二
float max_overlap_A = ARadius * 0.2f;
float max_overlap_B = BRadius * 0.2f;
float max_overlap = fminf(max_overlap_A, max_overlap_B);
overlap_s = fminf(overlap_s, max_overlap);

// 如果没有接触，清除历史记录并退出
if (overlap_s > 0) {
    // ========================================================================
    // SUP 步骤 1：输入缩放（将缩放后的变量转换为原始变量）
    // 基于论文中的表 1
    // ========================================================================
    
    // 重叠量缩放：δ_O = δ_S / l
    float overlap_o = overlap_s / l;

    if (overlap_o > 0.f) {
        // printf("overlap_s: %f, overlap_o: %f\n", overlap_s, overlap_o);
        // 从原始模型结构中获取材料属性
        float E_cnt, G_cnt, CoR_cnt, mu_cnt, Crr_cnt, gamma_surf;
        {
            float E_A_orig = E[bodyAMatType]; 
            float nu_A_orig = nu[bodyAMatType];
            float E_B_orig = E[bodyBMatType];
            float nu_B_orig = nu[bodyBMatType];
            matProxy2ContactParam<float>(E_cnt, G_cnt, E_A_orig, nu_A_orig, E_B_orig, nu_B_orig);
            CoR_cnt = CoR[bodyAMatType][bodyBMatType];
            mu_cnt = mu[bodyAMatType][bodyBMatType];
            Crr_cnt = Crr[bodyAMatType][bodyBMatType];
            gamma_surf = Cohesion[bodyAMatType][bodyBMatType]; // 表面能
            
        }

        // 半径缩放：R_O = R_S / l
        float R_s_A = ARadius;
        float R_s_B = BRadius;
        float R_o_A = R_s_A / l;
        float R_o_B = R_s_B / l;

        // 接触点缩放：position_O = position_S / l
        float3 locCPA_s = locCPA;
        float3 locCPB_s = locCPB;
        float3 locCPA_o = locCPA_s / l;
        float3 locCPB_o = locCPB_s / l;

        // 旋转速度缩放：ω_O = ω_S * l
        float3 ARotVel_s = ARotVel;
        float3 BRotVel_s = BRotVel;
        float3 ARotVel_o = ARotVel_s * l;
        float3 BRotVel_o = BRotVel_s * l;

        // 计算原始尺度下接触点的旋转速度
        float3 rotVelCPA_o_local = cross(ARotVel_o, locCPA_o);
        float3 rotVelCPB_o_local = cross(BRotVel_o, locCPB_o);
        applyOriQToVector3<float, deme::oriQ_t>(rotVelCPA_o_local.x, rotVelCPA_o_local.y, rotVelCPA_o_local.z, AOriQ.w, AOriQ.x, AOriQ.y, AOriQ.z);
        applyOriQToVector3<float, deme::oriQ_t>(rotVelCPB_o_local.x, rotVelCPB_o_local.y, rotVelCPB_o_local.z, BOriQ.w, BOriQ.x, BOriQ.y, BOriQ.z);

        // 质量缩放：m_O = m_S / l³
        float mass_s_A = AOwnerMass;
        float mass_s_B = BOwnerMass;
        float mass_o_A = mass_s_A / (l * l * l);
        float mass_o_B = mass_s_B / (l * l * l);

        if (mass_o_A <= 0.f) {
            mass_o_A = 1e-12f;
        }
        if (mass_o_B <= 0.f) {
            mass_o_B = 1e-12f;
        }

        // 切向位移历史缩放（不进行时间步缩放）
        float3 delta_tan_s = make_float3(delta_tan_x, delta_tan_y, delta_tan_z);
        float3 delta_tan_o = delta_tan_s / l;

        float delta_time_s = delta_time;
        float delta_time_o = delta_time_s / l;

        // 原始尺度下力计算的变量
        float mass_eff_o, sqrt_Rd_o, beta_o;
        float R_star_o;
        float3 vrel_tan_o;

        // 初始化原始尺度下的力分量
        float3 F_normal_o_vec = make_float3(0.f, 0.f, 0.f);
        float3 F_tangential_o_vec = make_float3(0.f, 0.f, 0.f);
        float3 F_hertz_o_vec = make_float3(0.f, 0.f, 0.f);
        float3 F_adhesion_o_vec = make_float3(0.f, 0.f, 0.f);
        float3 torque_only_force_o = make_float3(0.f, 0.f, 0.f);

        // ========================================================================
        // SUP 步骤 2：在原始尺度下计算力（使用 SJKR-F 模型）
        // ========================================================================
        {
            // 计算原始尺度下的相对速度
            const float3 velB2A_o = (ALinVel + rotVelCPA_o_local) - (BLinVel + rotVelCPB_o_local);
            const float projection_o = dot(velB2A_o, B2A);
            vrel_tan_o = velB2A_o - projection_o * B2A;

            // 更新切向位移历史（使用 ts 而不是 ts_o）
            {
                delta_tan_o += ts * vrel_tan_o;  // 直接使用 ts（不缩放）
                const float disp_proj_o = dot(delta_tan_o, B2A);
                delta_tan_o -= disp_proj_o * B2A;
                delta_time_o += ts;  // 直接使用 ts（不缩放）
            }

            // 计算有效质量
            mass_eff_o = (mass_o_A * mass_o_B) / (mass_o_A + mass_o_B);
            if (mass_o_A <= 1e-12f && mass_o_B <= 1e-12f) {
                mass_eff_o = 1e-6f;
            } else if (mass_o_A <= 1e-12f) {
                mass_eff_o = mass_o_B;
            } else if (mass_o_B <= 1e-12f) {
                mass_eff_o = mass_o_A;
            }

            // 计算有效半径
            if (R_o_A <= 1.0e-12f && R_o_B <= 1.0e-12f) {
                R_star_o = 1.0e-6f;
            } else if (R_o_A <= 1.0e-12f) {
                R_star_o = R_o_B;
            } else if (R_o_B <= 1.0e-12f) {
                R_star_o = R_o_A;
            } else {
                R_star_o = (R_o_A * R_o_B) / (R_o_A + R_o_B);
            }

            sqrt_Rd_o = sqrtf(overlap_o * R_star_o); // 计算 a = sqrt(R_e * δ)，接触半径
            
            // 计算阻尼系数
            const float loge_o = (CoR_cnt < DEME_TINY_FLOAT) ? logf(DEME_TINY_FLOAT) : logf(CoR_cnt);
            beta_o = loge_o / sqrtf(loge_o * loge_o + deme::PI_SQUARED); // beta_o = loge / sqrt(loge² + π²)

            // ========================================================================
            // JKR简化模型实现
            // 根据文档中的方程：F = (4/3)√(R_e)(E*)δ^(3/2) + (4)√(πγE*(a³))n
            // ========================================================================

            // printf("E_cnt: %f, G_cnt: %f, CoR_cnt: %f, mu_cnt: %f, Crr_cnt: %f, gamma_surf: %f\n", 
            //        E_cnt, G_cnt, CoR_cnt, mu_cnt, Crr_cnt, gamma_surf);
            // printf("overlap_o: %f, R_star_o: %f, sqrt_Rd_o(a): %f\n", 
            //        overlap_o, R_star_o, sqrt_Rd_o);
            
            // 计算 SJKR-F 法向力
            // 第一项：Hertz 弹性力
            // F_H = (4/3)√(R_e)(E*)δ^(3/2)
            float F_hertz = (4.0f/3.0f) * E_cnt * sqrtf(R_star_o) * powf(overlap_o, 1.5f);
            // printf("hertz force: %f\n", F_hertz);
            
            // 第二项：粘附力（总是存在，即使 gamma_surf = 0）
            float F_adhesion = 0.0f;
            if (gamma_surf > 0.0f) {
                // F_JKR = 4√(πγE*)a³n
                F_adhesion = 4.0f * sqrtf(deme::PI * gamma_surf * E_cnt * sqrt_Rd_o * sqrt_Rd_o * sqrt_Rd_o) ;
                // printf("jkr force: %f\n", F_adhesion);
            }
            
            // SJKR-F 总法向力 = Hertz力 + JKR粘附力
            float F_normal_mag = F_hertz + F_adhesion;
            // printf("normal force: %f\n", F_normal_mag);

            // 使用Hertz模型的阻尼计算方式
            const float Sn_o = 2.0f * E_cnt * sqrt_Rd_o;
            const float k_n_o = deme::TWO_OVER_THREE * Sn_o;
            const float gamma_n_o = deme::TWO_TIMES_SQRT_FIVE_OVER_SIX * beta_o * sqrtf(Sn_o * mass_eff_o);

            // 总法向力 = SJKR-F力 + 阻尼力
            F_normal_o_vec = (F_normal_mag + gamma_n_o * projection_o) * B2A;

            // SJKR-F 模型的接触分离条件
            // 根据文档，SJKR-F 在 δ = 0 时激活/失活接触
            if (overlap_o <= 0.0f) {
                // 接触断开，清除所有力和历史
                F_normal_o_vec = make_float3(0.0f, 0.0f, 0.0f);
                delta_tan_o = make_float3(0.0f, 0.0f, 0.0f);
                delta_time_o = 0.0f;
            }
        }

        // 如果启用，计算滚动阻力
        if (Crr_cnt > 0.0f && length(F_normal_o_vec) > DEME_TINY_FLOAT) {
            // Figure out if we should apply rolling resistance force
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
            
            // If should, then compute it (using Schwartz model)
            if (should_add_rolling_resistance_o) {
                const float3 v_rot_o = rotVelCPB_o_local - rotVelCPA_o_local;
                const float v_rot_o_mag = length(v_rot_o);

                if (v_rot_o_mag > DEME_TINY_FLOAT) {
                    torque_only_force_o = (v_rot_o / v_rot_o_mag) * (Crr_cnt * length(F_normal_o_vec));
                    // printf("torque force: %f, %f, %f\n", torque_only_force_o.x, torque_only_force_o.y, torque_only_force_o.z);
                }
            }
        }

        // 如果启用摩擦，计算切向力
        if (mu_cnt > 0.0f && length(F_normal_o_vec) > DEME_TINY_FLOAT) {
            const float kt_o = 8.f * G_cnt * sqrt_Rd_o;
            const float gt_o = -deme::TWO_TIMES_SQRT_FIVE_OVER_SIX * beta_o * sqrtf(mass_eff_o * kt_o);
            float3 tangent_force_o = -kt_o * delta_tan_o - gt_o * vrel_tan_o;
            const float ft_o = length(tangent_force_o);
            
            if (ft_o > DEME_TINY_FLOAT) {
                // Reverse-engineer to get tangential displacement
                const float ft_max_o = length(F_normal_o_vec) * mu_cnt;
                if (ft_o > ft_max_o) {
                    tangent_force_o = (ft_max_o / ft_o) * tangent_force_o;
                    delta_tan_o = (tangent_force_o + gt_o * vrel_tan_o) / (-kt_o);
                }
            } else {
                tangent_force_o = make_float3(0.f, 0.f, 0.f);
            }
            
            // Use F_tangential_o_vec to collect tangent_force
            F_tangential_o_vec = tangent_force_o;
            // printf("tangent force: %f, %f, %f\n", F_tangential_o_vec.x, F_tangential_o_vec.y, F_tangential_o_vec.z);
        }

        // ========================================================================
        // SUP 步骤 3：将力缩放回缩放后的系统
        // 根据方程 (23)：F_IS = l²·F_IO
        // ========================================================================
        float l_sq = l * l;
        float3 F_total_o_vec = F_normal_o_vec + F_tangential_o_vec + torque_only_force_o;

        // 注意：JKR 模型中的粘附力已经包含在 F_normal_o_vec 中，
        // 所以不需要额外添加凝聚力项

        force = F_total_o_vec * l_sq;

        // 注意：如果需要单独计算扭矩用于旋转计算，
        // 根据方程 (25)：M_IS = l²·M_IO，它们也应该按 l² 缩放

        // ========================================================================
        // SUP 步骤 4：更新缩放后系统中的历史变量
        // ========================================================================
        delta_tan_s = delta_tan_o * l;
        delta_tan_x = delta_tan_s.x;
        delta_tan_y = delta_tan_s.y;
        delta_tan_z = delta_tan_s.z;

        delta_time_s = delta_time_o * l;
        delta_time = delta_time_s;
    } else { 
        // 原始尺度下没有接触
        delta_time = 0;
        delta_tan_x = 0;
        delta_tan_y = 0;
        delta_tan_z = 0;
    }
} else {
    // 缩放尺度下没有接触
    delta_time = 0;
    delta_tan_x = 0;
    delta_tan_y = 0;
    delta_tan_z = 0;
}
