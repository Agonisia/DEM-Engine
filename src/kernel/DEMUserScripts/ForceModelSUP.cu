// 基于 FullHertzianForceModel.cu 修改，使用SUP模型 + CDT滚动阻力 + 论文切向力模型
// 参考文献：Hu et al., Powder Technology 438 (2024)

// VERSION_250603: 移除了内部时间步缩放
// VERSION_250604: 添加了JKR粘附力模型
// VERSION_250630：修改JKR模型为SJKR-F
// VERSION_250701：修改JKR模型为更简单的计算
// VERSION_250825：力的缩放指数可调节
// VERSION_250826: 修改了滚动阻力和切向力的计算公式

// 获取缩放因子
float l = scale_factor_l[bodyAMatType]; // SUP缩放因子
float index = 2; // 力的缩放指数，通常为2

if (l < 1e-5f) {
    l = 1.0f;
}

float overlap_s = overlapDepth;
// 限制重叠量不超过粒子半径的十分之二
float max_overlap_A = ARadius * 0.2f;
float max_overlap_B = BRadius * 0.2f;
float max_overlap = fminf(max_overlap_A, max_overlap_B);
overlap_s = fminf(overlap_s, max_overlap);

if (overlap_s > 0) {
    // ========================================================================
    // SUP 步骤 1：输入缩放（将缩放后的变量转换为原始变量）
    // ========================================================================
    
    // 重叠量缩放：δ_O = δ_S / l
    float overlap_o = overlap_s / l;

    if (overlap_o > 0.f) {
        // 获取材料属性
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
            gamma_surf = Cohesion[bodyAMatType][bodyBMatType];
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
        applyOriQToVector3<float, deme::oriQ_t>(rotVelCPA_o_local.x, rotVelCPA_o_local.y, rotVelCPA_o_local.z, 
                                                AOriQ.w, AOriQ.x, AOriQ.y, AOriQ.z);
        applyOriQToVector3<float, deme::oriQ_t>(rotVelCPB_o_local.x, rotVelCPB_o_local.y, rotVelCPB_o_local.z, 
                                                BOriQ.w, BOriQ.x, BOriQ.y, BOriQ.z);

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

        // 切向位移历史缩放
        float3 delta_tan_s = make_float3(delta_tan_x, delta_tan_y, delta_tan_z);
        float3 delta_tan_o = delta_tan_s / l;

        float delta_time_s = delta_time;
        float delta_time_o = delta_time_s / l;

        // 有效质量计算
        float mass_eff_o = (mass_o_A * mass_o_B) / (mass_o_A + mass_o_B);
        if (mass_o_A <= 1e-12f && mass_o_B <= 1e-12f) {
            mass_eff_o = 1e-6f;
        } else if (mass_o_A <= 1e-12f) {
            mass_eff_o = mass_o_B;
        } else if (mass_o_B <= 1e-12f) {
            mass_eff_o = mass_o_A;
        }
        
        // 有效半径计算
        float R_star_o;
        if (R_o_A <= 1.0e-12f && R_o_B <= 1.0e-12f) {
            R_star_o = 1.0e-6f;
        } else if (R_o_A <= 1.0e-12f) {
            R_star_o = R_o_B;
        } else if (R_o_B <= 1.0e-12f) {
            R_star_o = R_o_A;
        } else {
            R_star_o = (R_o_A * R_o_B) / (R_o_A + R_o_B);
        }

        float sqrt_Rd_o = sqrtf(overlap_o * R_star_o);
        const float loge_o = (CoR_cnt < DEME_TINY_FLOAT) ? logf(DEME_TINY_FLOAT) : logf(CoR_cnt);
        float beta_o = loge_o / sqrtf(loge_o * loge_o + deme::PI_SQUARED);

        // 初始化力
        float3 F_normal_o_vec = make_float3(0.f, 0.f, 0.f);
        float3 F_tangential_o_vec = make_float3(0.f, 0.f, 0.f);
        float3 torque_only_force_o = make_float3(0.f, 0.f, 0.f);

        // ========================================================================
        // SUP 步骤 2：在原始尺度下计算力（使用论文中的模型）
        // ========================================================================
        
        // 计算相对速度
        const float3 velB2A_o = (ALinVel + rotVelCPA_o_local) - (BLinVel + rotVelCPB_o_local);
        const float projection_o = dot(velB2A_o, B2A);
        float3 vrel_tan_o = velB2A_o - projection_o * B2A;

        // 时间步缩放
        float ts_o = ts / l;

        // ============ 法向力计算（使用SJKR-F模型）============
        float F_hertz = (4.0f/3.0f) * E_cnt * sqrtf(R_star_o) * powf(overlap_o, 1.5f);
        
        float F_adhesion = 0.0f;
        if (gamma_surf > 0.0f) {
            F_adhesion = 4.0f * sqrtf(deme::PI * gamma_surf * E_cnt * sqrt_Rd_o * sqrt_Rd_o * sqrt_Rd_o);
        }
        
        float F_normal_mag = F_hertz + F_adhesion;
        
        // 法向阻尼
        const float Sn_o = 2.0f * E_cnt * sqrt_Rd_o;
        const float k_n_o = deme::TWO_OVER_THREE * Sn_o;
        const float gamma_n_o = deme::TWO_TIMES_SQRT_FIVE_OVER_SIX * beta_o * sqrtf(Sn_o * mass_eff_o);
        
        F_normal_o_vec = (F_normal_mag + gamma_n_o * projection_o) * B2A;
        
        if (overlap_o <= 0.0f) {
            F_normal_o_vec = make_float3(0.0f, 0.0f, 0.0f);
            delta_tan_o = make_float3(0.0f, 0.0f, 0.0f);
            delta_time_o = 0.0f;
        }

        // ============ 切向力计算（使用论文中的修正Cundall-Strack模型）============
        if (mu_cnt > 0.0f && length(F_normal_o_vec) > DEME_TINY_FLOAT) {
            // 更新切向位移历史（论文方程11中的积分）
            delta_tan_o += ts_o * vrel_tan_o;
            
            // 移除法向分量
            const float disp_proj_o = dot(delta_tan_o, B2A);
            delta_tan_o -= disp_proj_o * B2A;
            
            // 切向刚度和阻尼（论文第3页）
            const float kt_o = 8.0f * G_cnt * sqrt_Rd_o;
            const float gt_o = -deme::TWO_TIMES_SQRT_FIVE_OVER_SIX * beta_o * sqrtf(mass_eff_o * kt_o);
            
            // 计算切向力（论文方程11）
            float3 tangent_force_o = -kt_o * delta_tan_o - gt_o * vrel_tan_o;
            const float ft_o = length(tangent_force_o);
            
            if (ft_o > DEME_TINY_FLOAT) {
                // 库仑摩擦限制
                const float ft_max_o = length(F_normal_o_vec) * mu_cnt;
                if (ft_o > ft_max_o) {
                    // 达到滑动状态，调整切向力和位移
                    tangent_force_o = (ft_max_o / ft_o) * tangent_force_o;
                    // 反向计算切向位移以满足摩擦限制
                    delta_tan_o = (tangent_force_o + gt_o * vrel_tan_o) / (-kt_o);
                }
            } else {
                tangent_force_o = make_float3(0.f, 0.f, 0.f);
            }
            
            F_tangential_o_vec = tangent_force_o;
            delta_time_o += ts_o;
        }

        // ============ CDT滚动阻力模型（论文方程21）============
        if (Crr_cnt > 0.0f && length(F_normal_o_vec) > DEME_TINY_FLOAT) {
            // 计算相对角速度
            const float3 omega_rel_o = ARotVel_o - BRotVel_o;
            const float omega_rel_mag_o = length(omega_rel_o);
            
            if (omega_rel_mag_o > DEME_TINY_FLOAT) {
                // CDT模型：直接计算等效力
                // torque_only_force = -(ω_rel/|ω_rel|) * μ_r * F_N
                // 负号通过系数实现
                float coeff = -Crr_cnt * length(F_normal_o_vec) / omega_rel_mag_o;
                
                torque_only_force_o = make_float3(
                    omega_rel_o.x * coeff,
                    omega_rel_o.y * coeff,
                    omega_rel_o.z * coeff
                );
            }
        }

        // ========================================================================
        // SUP 步骤 3：将力缩放回缩放后的系统
        // ========================================================================
        
        float l_force = powf(l, (float)index);
        float3 F_total_o_vec = F_normal_o_vec + F_tangential_o_vec + torque_only_force_o;
        
        force = F_total_o_vec * l_force;


        //     printf("Final force = (%f, %f, %f), l_force = %f (l=%f, index=%d)\n",
        //    force.x, force.y, force.z,
        //    l_force, l, index);
        

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
        delta_time = 0;
        delta_tan_x = 0;
        delta_tan_y = 0;
        delta_tan_z = 0;
    }
} else {
    delta_time = 0;
    delta_tan_x = 0;
    delta_tan_y = 0;
    delta_tan_z = 0;
}
