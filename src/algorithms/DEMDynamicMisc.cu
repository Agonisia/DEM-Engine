//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//
//	SPDX-License-Identifier: BSD-3-Clause

#include <algorithms/DEMStaticDeviceSubroutines.h>
#include <algorithms/DEMStaticDeviceUtilities.cuh>

#include <core/utils/GpuError.h>
#include <kernel/DEMHelperKernels.cuh>

namespace deme {

__global__ void getContactForcesConcerningOwners_impl(float3* d_points,
                                                      float3* d_forces,
                                                      float3* d_torques,
                                                      unsigned long long* d_numUsefulCnt,
                                                      bodyID_t* d_ownerIDs,
                                                      size_t IDListSize,
                                                      DEMSimParams* simParams,
                                                      DEMDataDT* granData,
                                                      size_t numCnt,
                                                      bool need_torque,
                                                      bool torque_in_local) {
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < numCnt) {
        // ====================================================================================
        // 1. 获取接触对信息
        // ====================================================================================
        bodyID_t geoA = granData->idGeometryA[i];
        bodyID_t ownerA = granData->ownerClumpBody[geoA];
        bodyID_t geoB = granData->idGeometryB[i];
        contact_t typeB = granData->contactType[i];
        bodyID_t ownerB = DEME_GET_GEO_OWNER_ID(geoB, typeB);

        // ====================================================================================
        // 2. 检查是否涉及关心的物体
        // ====================================================================================
        bool AorB;  // true for A, false for B
        if (cuda_binary_search<bodyID_t, ssize_t>(d_ownerIDs, ownerA, 0, IDListSize - 1)) {
            AorB = true;
        } else if (cuda_binary_search<bodyID_t, ssize_t>(d_ownerIDs, ownerB, 0, IDListSize - 1)) {
            AorB = false;
        } else {
            return;
        }

        // ====================================================================================
        // 3. 获取接触力和力矩相关的力
        // ====================================================================================
        float3 force = granData->contactForces[i];          // 主接触力 (作用在A上)
        float3 torque_force;                                // 产生纯力矩的等效力
        if (need_torque) {
            torque_force = granData->contactTorque_convToForce[i];
        }
        
        // 检查力的大小
        {
            float mag = (need_torque) ? length(force) + length(torque_force) : length(force);
            if (mag < DEME_TINY_FLOAT)
                return;
        }

        // 获取写入索引
        unsigned long long writeIndex = atomicAdd(d_numUsefulCnt, 1);
        
        // ====================================================================================
        // 4. 根据关心的物体，设置相关变量
        // ====================================================================================
        float3 cntPnt_local;  // 接触点在局部坐标系的位置
        bodyID_t ownerID;
        
        if (AorB) {
            // 关心物体A
            cntPnt_local = granData->contactPointGeometryA[i];
            ownerID = ownerA;
            // force已经是作用在A上的，不需要改变
        } else {
            // 关心物体B
            cntPnt_local = granData->contactPointGeometryB[i];
            ownerID = ownerB;
            // B受到的力与A相反
            force = -force;
            if (need_torque) {
                torque_force = -torque_force;
            }
        }

        // ====================================================================================
        // 5. 获取物体的姿态四元数（只声明一次）
        // ====================================================================================
        float4 oriQ;
        oriQ.w = granData->oriQw[ownerID];
        oriQ.x = granData->oriQx[ownerID];
        oriQ.y = granData->oriQy[ownerID];
        oriQ.z = granData->oriQz[ownerID];

        // ====================================================================================
        // 6. 计算力矩（如果需要）
        // ====================================================================================
        float3 total_torque;
        if (need_torque) {
            // 步骤 6.1: 将主接触力从全局转换到局部坐标系
            float3 force_local = force;
            applyOriQToVector3<float, deme::oriQ_t>(force_local.x, force_local.y, force_local.z, 
                                                    oriQ.w, -oriQ.x, -oriQ.y, -oriQ.z);
            
            // 步骤 6.2: 计算主接触力产生的力矩 (M = r × F)
            float3 main_torque_local = cross(cntPnt_local, force_local);
            
            // 步骤 6.3: 将等效力从全局转换到局部坐标系
            float3 torque_force_local = torque_force;
            applyOriQToVector3<float, deme::oriQ_t>(torque_force_local.x, torque_force_local.y, torque_force_local.z,
                                                    oriQ.w, -oriQ.x, -oriQ.y, -oriQ.z);
            
            // 步骤 6.4: 计算等效力产生的力矩
            float3 other_torque_local = cross(cntPnt_local, torque_force_local);
            
            // 步骤 6.5: 叠加得到总力矩（局部坐标系）
            total_torque = main_torque_local + other_torque_local;
            
            // 步骤 6.6: 根据用户需求决定输出坐标系
            if (!torque_in_local) {
                // 转换回全局坐标系
                applyOriQToVector3<float, deme::oriQ_t>(total_torque.x, total_torque.y, total_torque.z,
                                                        oriQ.w, oriQ.x, oriQ.y, oriQ.z);
            }
        }

        // ====================================================================================
        // 7. 计算接触点的全局坐标
        // ====================================================================================
        double3 CoM;
        voxelID_t voxel = granData->voxelID[ownerID];
        subVoxelPos_t subVoxX = granData->locX[ownerID];
        subVoxelPos_t subVoxY = granData->locY[ownerID];
        subVoxelPos_t subVoxZ = granData->locZ[ownerID];
        
        voxelIDToPosition<double, voxelID_t, subVoxelPos_t>(CoM.x, CoM.y, CoM.z, voxel, subVoxX, subVoxY, subVoxZ,
                                                            simParams->nvXp2, simParams->nvYp2, simParams->voxelSize,
                                                            simParams->l);
        CoM.x += simParams->LBFX;
        CoM.y += simParams->LBFY;
        CoM.z += simParams->LBFZ;
        
        // 将局部接触点转换到全局坐标系
        float3 cntPnt_global = cntPnt_local;
        applyFrameTransformLocalToGlobal<float3, double3, float4>(cntPnt_global, CoM, oriQ);
        
        // ====================================================================================
        // 8. 写入输出数组
        // ====================================================================================
        d_points[writeIndex] = cntPnt_global;
        d_forces[writeIndex] = force;
        if (need_torque) {
            d_torques[writeIndex] = total_torque;
        }
    }
}

void getContactForcesConcerningOwners(float3* d_points,
                                      float3* d_forces,
                                      float3* d_torques,
                                      size_t* d_numUsefulCnt,
                                      bodyID_t* d_ownerIDs,
                                      size_t IDListSize,
                                      DEMSimParams* simParams,
                                      DEMDataDT* granData,
                                      size_t numCnt,
                                      bool need_torque,
                                      bool torque_in_local,
                                      cudaStream_t& this_stream) {
    size_t blocks_needed = (numCnt + DEME_MAX_THREADS_PER_BLOCK - 1) / DEME_MAX_THREADS_PER_BLOCK;
    getContactForcesConcerningOwners_impl<<<blocks_needed, DEME_MAX_THREADS_PER_BLOCK, 0, this_stream>>>(
        d_points, d_forces, d_torques, reinterpret_cast<unsigned long long*>(d_numUsefulCnt), d_ownerIDs, IDListSize,
        simParams, granData, numCnt, need_torque, torque_in_local);
    DEME_GPU_CALL(cudaStreamSynchronize(this_stream));
}

}  // namespace deme
