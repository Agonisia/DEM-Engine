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
// 限制重叠量不超过粒子半径的十分之一
float max_overlap_A = ARadius * 0.1f;
float max_overlap_B = BRadius * 0.1f;
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

        // 线速度不缩放：v_O = v_S
        float3 ALinVel_o = ALinVel;
        float3 BLinVel_o = BLinVel;

        // 角速度缩放：ω_O = ω_S * l
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
        const float3 velB2A_o = (ALinVel_o + rotVelCPA_o_local) - (BLinVel_o + rotVelCPB_o_local);
        const float projection_o = dot(velB2A_o, B2A);
        float3 vrel_tan_o = velB2A_o - projection_o * B2A;

        //===========================修改1：修正时间步缩放===========================
        //float ts_o = ts;  // 修正：不使用时间步缩放ts_o
        //===========================修改结束===========================

        // 计算JKR接触半径
        float a_jkr = sqrt_Rd_o;  // 默认值（纯Hertz接触半径）
        if (gamma_surf > 0.0f) {
            // ===========================Cardano方法计算JKR接触半径（接近LIGGGHTS版本）===========================
            float sqrt_Rd = sqrtf(R_star_o * overlap_o);
            
            // LIGGGHTS的精确公式参数
            float reff = R_star_o;
            float deltan = overlap_o;
            float se = gamma_surf;
            float Yeff = E_cnt;
            
            float c0 = reff * reff * deltan * deltan;
            float c1 = -4.0f * deme::PI * se * reff * reff / (Yeff + 1e-10f);
            float c2 = -2.0f * reff * deltan;
            
            float P = -c2*c2/12.0f - c0;
            float Q = -c2*c2*c2/108.0f + c0*c2/3.0f - c1*c1/8.0f;
            
            // 使用fmaxf避免负值（类似LIGGGHTS的abs）
            float discriminant = Q*Q/4.0f + P*P*P/27.0f;
            discriminant = fmaxf(discriminant, 0.0f);
            
            float Uterm2 = sqrtf(discriminant);
            
            // 立方根的无分支处理
            float arg = -Q/2.0f + Uterm2;
            float sign_arg = (arg >= 0.0f) ? 1.0f : -1.0f;
            float U = sign_arg * powf(fabsf(arg) + 1e-10f, 0.333333f);
            
            // 无分支版本的s计算
            float P_nonzero = (fabsf(P) > 1e-10f) ? 1.0f : 0.0f;
            float s_with_P = -5.0f*c2/6.0f + U - P/(3.0f*U + 1e-10f);
            float sign_Q = (Q >= 0.0f) ? 1.0f : -1.0f;
            float s_without_P = -5.0f*c2/6.0f - sign_Q * powf(fabsf(Q) + 1e-10f, 0.333333f);
            float s = P_nonzero * s_with_P + (1.0f - P_nonzero) * s_without_P;
            
            // 继续LIGGGHTS算法
            float w_arg = c2 + 2.0f*s;
            w_arg = fmaxf(w_arg, 0.0f);
            float w = sqrtf(w_arg);
            
            float lambda = c1/(2.0f*w + 1e-10f);
            
            float aterm1 = w*w - 4.0f*(c2 + s + lambda);
            aterm1 = fmaxf(aterm1, 0.0f);
            float aterm2 = sqrtf(aterm1);
            
            float a = 0.5f * (w + aterm2);
            
            a_jkr = fmaxf(a, sqrt_Rd);
            // ===========================Cardano计算结束===========================
        }
        
        // 基于JKR接触半径计算力
        float F_repulsion = 0.0f;
        float F_sJKR = 0.0f;
        float F_normal_mag = 0.0f;
        
        if (a_jkr > 0.0f) {
            float a3_jkr = a_jkr * a_jkr * a_jkr;
            
            // Repulsion力（基于JKR接触半径）
            F_repulsion = (4.0f/3.0f) * E_cnt * a3_jkr / R_star_o;
            
            // JKR粘附力（吸引，负值）
            if (gamma_surf > 0.0f) {
                F_sJKR = -4.0f * sqrtf(deme::PI * gamma_surf * E_cnt * a3_jkr);
            }
            
            // 总接触力
            F_normal_mag = F_repulsion + F_sJKR;
            
            // 防止非物理的过度吸引力（可选的截断处理）
            if (F_normal_mag < 0.0f && overlap_o < 1e-8f) {
                F_normal_mag = 0.0f;
            }
        }
        
        // 法向阻尼（基于JKR接触半径）
        const float Sn_o = 2.0f * E_cnt * a_jkr;
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
            delta_tan_o += (ts * vrel_tan_o)/l;
            
            // 移除法向分量
            const float disp_proj_o = dot(delta_tan_o, B2A);
            delta_tan_o -= disp_proj_o * B2A;
            
            // 切向刚度基于JKR接触半径
            const float kt_o = 8.0f * G_cnt * a_jkr;
            
            const float gt_o = -deme::TWO_TIMES_SQRT_FIVE_OVER_SIX * beta_o * sqrtf(mass_eff_o * kt_o);
            
            // 计算切向力（论文方程11）
            float3 tangent_force_o = -kt_o * delta_tan_o - gt_o * vrel_tan_o;
            const float ft_o = length(tangent_force_o);
            
            if (ft_o > DEME_TINY_FLOAT) {
                // ===== 添加JKR摩擦增强 =====
                // 计算JKR粘附对摩擦的贡献
                float Fc_o = 0.0f;
                if (gamma_surf > 0.0f) {
                    Fc_o = 3.0f * deme::PI * gamma_surf * R_star_o;
                }
                
                // 库仑摩擦限制（包含JKR增强）
                const float ft_max_o = (length(F_normal_o_vec) + 2.0f * Fc_o) * mu_cnt;
                // ===== 结束JKR修改 =====
                
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
            delta_time_o += ts;
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
