// 基于 FullHertzianForceModel.cu 修改，用于缩放颗粒（SUP）模型
// 参考文献："Inter-particle torque scaling in coarse grained DEM with rolling resistance 
// and particle size distributions" - Hu et al., Powder Technology 438 (2024)

// VERSION_250603: 根据导师反馈移除了内部时间步缩放，添加了系数力计算
// VERSION_250604: 添加了JKR粘附力模型

// 获取缩放因子 l
float l = scale_factor_l;
// 确保 l 有效，如果未正确设置或 <= 0，则默认为 1.0
if (l < 1e-5f) { // 使用小的 epsilon 值
    l = 1.0f;
}

float overlap_s = overlapDepth;
// 如果没有接触，清除历史记录并退出
if (overlap_s > 0) {
    // ========================================================================
    // SUP 步骤 1：输入缩放（将缩放后的变量转换为原始变量）
    // 基于论文中的表 1
    // ========================================================================
    
    // 重叠量缩放：δ_O = δ_S / l
    float overlap_o = overlap_s / l;

    if (overlap_o > 0.f) {
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
            gamma_surf = Cohesion[bodyAMatType][bodyBMatType]; // 将 Cohesion 系数解释为表面能
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
        float3 torque_only_force_o = make_float3(0.f, 0.f, 0.f);

        // ========================================================================
        // SUP 步骤 2：在原始尺度下计算力（使用 JKR 模型）
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
            if (R_o_A <= 1.0e-12f && R_o_B <= 1.0e-12f) { R_star_o = 1.0e-6f; }
            else if (R_o_A <= 1.0e-12f) { R_star_o = R_o_B; }
            else if (R_o_B <= 1.0e-12f) { R_star_o = R_o_A; }
            else { R_star_o = (R_o_A * R_o_B) / (R_o_A + R_o_B); }

            sqrt_Rd_o = sqrtf(overlap_o * R_star_o);
            
            // 计算阻尼系数
            const float loge_o = (CoR_cnt < DEME_TINY_FLOAT) ? logf(DEME_TINY_FLOAT) : logf(CoR_cnt);
            beta_o = loge_o / sqrtf(loge_o * loge_o + deme::PI_SQUARED);

            // JKR 模型实现
            if (gamma_surf > 0.0f && overlap_o > -0.1f * R_star_o) {
                // JKR 模型：求解接触半径 a
                // 方程：a^4 - 2*R_star*delta_n*a^2 - (4*PI*gamma*R_star^2/E_star)*a + R_star^2*delta_n^2 = 0
                
                // 平衡接触半径（零载荷时）
                float a0 = powf(9.0f * deme::PI * gamma_surf * R_star_o * R_star_o / E_cnt, 1.0f/3.0f);
                
                // 初始猜测值
                float a = sqrtf(R_star_o * fmaxf(overlap_o, 0.0f)) + 0.5f * a0;
                
                // 牛顿-拉夫逊迭代求解接触半径
                for (int i = 0; i < 10; i++) {
                    float a2 = a * a;
                    float a3 = a2 * a;
                    float a4 = a2 * a2;
                    
                    float f = a4 - 2.0f * R_star_o * overlap_o * a2 
                             - 4.0f * deme::PI * gamma_surf * R_star_o * R_star_o * a / E_cnt 
                             + R_star_o * R_star_o * overlap_o * overlap_o;
                    
                    float df = 4.0f * a3 - 4.0f * R_star_o * overlap_o * a 
                              - 4.0f * deme::PI * gamma_surf * R_star_o * R_star_o / E_cnt;
                    
                    if (fabsf(df) > DEME_TINY_FLOAT) {
                        float da = -f / df;
                        a = a + da;
                        
                        // 收敛检查
                        if (fabsf(da) < 1e-6f * a) break;
                    }
                    
                    // 确保 a 为正值且合理
                    if (a < 0.1f * a0) a = 0.1f * a0;
                    if (a > 10.0f * sqrtf(R_star_o * fmaxf(overlap_o, 0.0f))) {
                        a = 10.0f * sqrtf(R_star_o * fmaxf(overlap_o, 0.0f));
                    }
                }
                
                // JKR 法向力计算
                float a3 = a * a * a;
                float sqrt_a3 = sqrtf(a3);
                
                // 法向力 = 弹性项 - 粘附项
                float F_elastic = (4.0f * E_cnt * a3) / (3.0f * R_star_o);
                float F_adhesion = 4.0f * sqrtf(deme::PI * gamma_surf * E_cnt) * sqrt_a3;
                float F_normal_mag = F_elastic - F_adhesion;
                
                // 添加阻尼（基于接触半径）
                float gamma_n_o = 0.0f;
                if (a > DEME_TINY_FLOAT) {
                    float Sn_jkr = 2.0f * E_cnt * a / sqrtf(R_star_o);
                    gamma_n_o = deme::TWO_TIMES_SQRT_FIVE_OVER_SIX * beta_o * sqrtf(Sn_jkr * mass_eff_o);
                }
                
                // 总法向力
                F_normal_o_vec = (F_normal_mag + gamma_n_o * projection_o) * B2A;
                
                // 分离条件检查（简化模型：当法向重叠为负时断开接触）
                if (overlap_o < 0.0f) {
                    F_normal_o_vec = make_float3(0.0f, 0.0f, 0.0f);
                    delta_tan_o = make_float3(0.0f, 0.0f, 0.0f);
                    delta_time_o = 0.0f;
                }
                
            } else {
                // 无粘附的 Hertzian 模型（原始代码）
                const float Sn_o = 2.f * E_cnt * sqrt_Rd_o;
                const float k_n_o = deme::TWO_OVER_THREE * Sn_o;
                const float gamma_n_o = deme::TWO_TIMES_SQRT_FIVE_OVER_SIX * beta_o * sqrtf(Sn_o * mass_eff_o);
                F_normal_o_vec = (k_n_o * overlap_o + gamma_n_o * projection_o) * B2A;
            }
        }

        // 如果启用，计算滚动阻力
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

        // 如果启用摩擦，计算切向力
        if (mu_cnt > 0.0f && length(F_normal_o_vec) > DEME_TINY_FLOAT) {
            const float kt_o = 8.f * G_cnt * sqrt_Rd_o;
            const float gt_o = -deme::TWO_TIMES_SQRT_FIVE_OVER_SIX * beta_o * sqrtf(mass_eff_o * kt_o);
            float3 tangent_force_trial_o = -kt_o * delta_tan_o - gt_o * vrel_tan_o;

            const float ft_o_mag_trial = length(tangent_force_trial_o);
            if (ft_o_mag_trial > DEME_TINY_FLOAT) {
                // 使用 JKR 模型的有效法向力
                const float ft_max_o = length(F_normal_o_vec) * mu_cnt;
                if (ft_o_mag_trial > ft_max_o) {
                    // 滑动摩擦
                    F_tangential_o_vec = (ft_max_o / ft_o_mag_trial) * tangent_force_trial_o;
                    if (fabs(kt_o) > DEME_TINY_FLOAT) {
                        delta_tan_o = (F_tangential_o_vec + gt_o * vrel_tan_o) / (-kt_o);
                    } else {
                        delta_tan_o = make_float3(0.f,0.f,0.f);
                    }
                } else {
                    // 静摩擦
                    F_tangential_o_vec = tangent_force_trial_o;
                }
            } else {
                F_tangential_o_vec = make_float3(0.f, 0.f, 0.f);
            }
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