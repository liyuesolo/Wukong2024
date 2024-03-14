#include "../include/EoLRodSim.h"
#include "../autodiff/EoLRodBendingAndTwisting.h"



T EoLRodSim::add3DBendingAndTwistingEnergy(bool bending, bool twisting)
{
    T energy = 0.0;
    for (auto& rod : Rods)
    {
        T energy_current = energy;
        rod->iterate3Nodes([&](int node_i, int node_j, int node_k, int second, bool is_crossing){
            if (is_crossing)
                return;
            TV xi, xj, xk, Xi, Xj, Xk;
            rod->x(node_i, xi); rod->x(node_j, xj); rod->x(node_k, xk);
            rod->X(node_i, Xi); rod->X(node_j, Xj); rod->X(node_k, Xk);

            int left_rod = (rod->closed && second == 0) ? rod->numSeg() - 1 : second - 1;

            T theta0 = rod->reference_angles[left_rod];
            T theta1 = rod->reference_angles[second];
            TV referenceNormal1 = rod->reference_frame_us[left_rod];
            TV referenceNormal2 = rod->reference_frame_us[second];

            TV referenceTangent1 = rod->prev_tangents[left_rod];
            TV referenceTangent2 = rod->prev_tangents[second];


            Matrix<T, 2, 2> B = bending ? rod->bending_coeffs : Matrix<T, 2, 2>::Zero();
            T kt = twisting ? rod->kt : 0.0;

            T reference_twist = rod->reference_twist[second];
            energy += computeRodBendingAndTwistEnergy(B, kt, 0.0, referenceNormal1, referenceTangent1,
                referenceNormal2, referenceTangent2, reference_twist, xk, xi, xj, Xk, Xi, Xj, theta0, theta1);
            });
        
    }
    return energy;
}

void EoLRodSim::add3DBendingAndTwistingForce(Eigen::Ref<VectorXT> residual)
{
    VectorXT residual_cp = residual;
    for (auto& rod : Rods)
    {
        rod->iterate3NodesWithOffsets([&](int node_i, int node_j, int node_k, 
            Offset offset_i, Offset offset_j, Offset offset_k, int second, bool is_crossing)
        {
            if (is_crossing)
                return;
            TV xi, xj, xk, Xi, Xj, Xk, dXi, dXj, dXk;
            rod->x(node_i, xi); rod->x(node_j, xj); rod->x(node_k, xk);
            rod->XdX(node_i, Xi, dXi); rod->XdX(node_j, Xj, dXj); rod->XdX(node_k, Xk, dXk);
            
            int left_rod = (rod->closed && second == 0) ? rod->numSeg() - 1 : second - 1;

            T theta0 = rod->reference_angles[left_rod];
            T theta1 = rod->reference_angles[second];
            TV referenceNormal1 = rod->reference_frame_us[left_rod];
            TV referenceNormal2 = rod->reference_frame_us[second];

            TV referenceTangent1 = rod->prev_tangents[left_rod];
            TV referenceTangent2 = rod->prev_tangents[second];

            Matrix<T, 2, 2> B = rod->bending_coeffs;
            T reference_twist = rod->reference_twist[second];

            Vector<T, 20> F;
            computeRodBendingAndTwistEnergyGradient(B, rod->kt, 0.0, referenceNormal1, referenceTangent1,
                referenceNormal2, referenceTangent2,  reference_twist, xk, xi, xj, Xk, Xi, Xj, theta0, theta1, F);
            
            F *= -1.0;

            residual.template segment<3>(offset_k[0]) += F.template segment<3>(0);
            residual.template segment<3>(offset_i[0]) += F.template segment<3>(3);
            residual.template segment<3>(offset_j[0]) += F.template segment<3>(3 + 3);

            residual(offset_k[3]) += F.template segment<3>(3*3).dot(dXk);
            residual(offset_i[3]) += F.template segment<3>(3*3 + 3).dot(dXi);
            residual(offset_j[3]) += F.template segment<3>(3*3 + 2*3).dot(dXj);

            int dof_theta0 = (rod->closed && second == 0) ? rod->theta_dof_start_offset + rod->numSeg() - 1 : rod->theta_dof_start_offset + second - 1;
            int dof_theta1 = rod->theta_dof_start_offset + second;
            residual[dof_theta0] += F[18];
            residual[dof_theta1] += F[19];
            // residual.template segment<2>(rod->theta_dof_start_offset + (second - 1)) += 
            //     F.template segment<2>(18);
        });
    }
    if(print_force_mag)
        std::cout << "bending and twist norm: " << (residual - residual_cp).norm() << std::endl;
}


void EoLRodSim::add3DBendingAndTwistingK(std::vector<Entry>& entry_K)
{
    for (auto& rod : Rods)
    {
        rod->iterate3NodesWithOffsets([&](int node_i, int node_j, int node_k, 
            Offset offset_i, Offset offset_j, Offset offset_k, int second, bool is_crossing)
        {
            if (is_crossing)
                return;
            TV xi, xj, xk, Xi, Xj, Xk, dXi, dXj, dXk;
            rod->x(node_i, xi); rod->x(node_j, xj); rod->x(node_k, xk);
            rod->XdX(node_i, Xi, dXi); rod->XdX(node_j, Xj, dXj); rod->XdX(node_k, Xk, dXk);

            std::vector<int> nodes = { node_k, node_i, node_j };
            std::vector<TV> dXdu = { dXk, dXi, dXj };

            std::vector<Offset> offsets = { offset_k, offset_i, offset_j };

            int dof_theta0 = (rod->closed && second == 0) ? rod->theta_dof_start_offset + rod->numSeg() - 1 : rod->theta_dof_start_offset + second - 1;
            int dof_theta1 = rod->theta_dof_start_offset + second;
            
            std::vector<int> theta_dofs = {dof_theta0, dof_theta1};

            int left_rod = (rod->closed && second == 0) ? rod->numSeg() - 1 : second - 1;

            T theta0 = rod->reference_angles[left_rod];
            T theta1 = rod->reference_angles[second];
            TV referenceNormal1 = rod->reference_frame_us[left_rod];
            TV referenceNormal2 = rod->reference_frame_us[second];

            TV referenceTangent1 = rod->prev_tangents[left_rod];
            TV referenceTangent2 = rod->prev_tangents[second];

            Matrix<T, 2, 2> B = rod->bending_coeffs;
            T reference_twist = rod->reference_twist[second];
            Matrix<T, 20, 20> J;
            
            computeRodBendingAndTwistEnergyHessian(B, rod->kt, 0.0, referenceNormal1, referenceTangent1,
                referenceNormal2, referenceTangent2,  reference_twist, xk, xi, xj, Xk, Xi, Xj, theta0, theta1, J);

            J *= -1.0;
            for(int k = 0; k < nodes.size(); k++)
                for(int l = 0; l < nodes.size(); l++)
                    for(int i = 0; i < 3; i++)
                        for (int j = 0; j < 3; j++)
                            {
                                entry_K.push_back(Eigen::Triplet<T>(offsets[k][i], offsets[l][j], -J(k*3 + i, l * 3 + j)));

                                entry_K.push_back(Eigen::Triplet<T>(offsets[k][i], offsets[l][3], -J(k*3 + i, 3 * 3 + l * 3 + j) * dXdu[l][j]));

                                entry_K.push_back(Eigen::Triplet<T>(offsets[k][3], offsets[l][j], -J(3 * 3 + k * 3 + i, l * 3 + j) * dXdu[k][i]));

                                
                                entry_K.push_back(Eigen::Triplet<T>(offsets[k][3], 
                                                                    offsets[l][3], 
                                                                    -J(3 * 3 + k*3 + i, 3 * 3 + l * 3 + j) * dXdu[l][j] * dXdu[k][i]));

                            }
            for(int k = 0; k < nodes.size(); k++)
                {
                    
                    for (int j = 0; j < 2; j++)
                    {
                        for(int i = 0; i < 3; i++)
                        {
                            entry_K.push_back(Eigen::Triplet<T>(offsets[k][i], theta_dofs[j], -J(k*3 + i, 18 + j)));
                            entry_K.push_back(Eigen::Triplet<T>(theta_dofs[j], offsets[k][i], -J(18 + j, k * 3 + i)));

                            entry_K.push_back(Eigen::Triplet<T>(offsets[k][3], theta_dofs[j], -J(3 * 3 + k * 3 + i, 18 + j) * dXdu[k][i]));
                            entry_K.push_back(Eigen::Triplet<T>(theta_dofs[j], offsets[k][3], -J(18 + j, 3 * 3 + k * 3 + i) * dXdu[k][i]));
                        }
                    }
                }
            for (int i = 0; i < 2; i++)
                for (int j = 0; j < 2; j++)
                    entry_K.push_back(Eigen::Triplet<T>(theta_dofs[i], theta_dofs[j], -J(18 + i, 18 + j)));
        });
        
    }

    for (auto& rod : Rods)
    {   
        continue;
        rod->iterate3NodesWithOffsets([&](int node_i, int node_j, int node_k, 
            Offset offset_i, Offset offset_j, Offset offset_k, int second, bool is_crossing)
        {
            if (is_crossing)
                return;
            TV xi, xj, xk, Xi, Xj, Xk, dXi, dXj, dXk, ddXi, ddXj, ddXk;
            rod->x(node_i, xi); rod->x(node_j, xj); rod->x(node_k, xk);
            rod->XdXddX(node_i, Xi, dXi, ddXi); rod->XdXddX(node_j, Xj, dXj, ddXj); rod->XdXddX(node_k, Xk, dXk, ddXk);

            T theta0 = rod->reference_angles[second - 1];
            T theta1 = rod->reference_angles[second];

            TV referenceNormal1 = rod->reference_frame_us[second - 1];
            TV referenceNormal2 = rod->reference_frame_us[second];

            TV referenceTangent1 = rod->prev_tangents[second - 1];
            TV referenceTangent2 = rod->prev_tangents[second];

            Matrix<T, 2, 2> B = rod->bending_coeffs;
            T reference_twist = rod->reference_twist[second];
            Vector<T, 20> F;
            computeRodBendingAndTwistEnergyGradient(B, rod->kt, 0.0, referenceNormal1, referenceTangent1, 
                referenceNormal2, referenceTangent2, reference_twist, xk, xi, xj, Xk, Xi, Xj, theta0, theta1, F);
            
            F *= -1.0;

            for(int d = 0; d < 3; d++)
            {
                entry_K.push_back(Eigen::Triplet<T>(offset_i[3], 
                                    offset_i[3], -F[1*3 + 3*3 + d] * ddXi[d]));
                entry_K.push_back(Eigen::Triplet<T>(offset_j[3], 
                                    offset_j[3], -F[2*3 + 3*3 + d] * ddXj[d]));
                entry_K.push_back(Eigen::Triplet<T>(offset_k[3], 
                                    offset_k[3], -F[0*3 + 3*3 + d] * ddXk[d]));
                std::cout << ddXk[d] << std::endl;
            }
        });
    }
}