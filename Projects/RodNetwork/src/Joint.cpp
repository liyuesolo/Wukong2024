// #include "../include/RodNetwork.h"
// #include "../autodiff/EoLRodEulerAngleBendingAndTwisting.h"
// #include "../autodiff/EoLRodBendingAndTwisting.h"


// void RodNetwork::addJointBendingAndTwistingK(std::vector<Entry>& entry_K)
// {
//     for (auto& crossing : rod_crossings)
//     {
//         int node_i = crossing->node_idx;
        
//         std::vector<int> rods_involved = crossing->rods_involved;
//         std::unordered_map<int, int> on_rod_idx = crossing->on_rod_idx;

//         std::vector<Vector<T, 2>> sliding_ranges = crossing->sliding_ranges;
//         Vector<T, 3> omega = crossing->omega;
//         Matrix<T, 3, 3> omega_acc = crossing->rotation_accumulated;

//         int cnt = 0;
//         for (int rod_idx : rods_involved)
//         {
//             Offset offset_i, offset_j, offset_k;
//             TV xi, xj, xk, Xi, Xj, Xk, dXi, dXj, dXk;
            
//             auto rod = Rods[rod_idx];
//             rod->getEntry(node_i, offset_i);
//             rod->x(node_i, xi); rod->XdX(node_i, Xi, dXi);

//             int node_j = rod->nodeIdx(on_rod_idx[rod_idx] + 1);
//             int node_k = rod->nodeIdx(on_rod_idx[rod_idx] - 1);
//             int edge1 = 0, edge0 = 0;
//             if (node_j != -1)
//             {
//                 rod->getEntry(node_j, offset_j);
//                 rod->x(node_j, xj);
//                 rod->XdX(node_j, Xj, dXj);
//                 edge1 = on_rod_idx[rod_idx];
//                 if(rod->closed && on_rod_idx[rod_idx] == 0)
//                     edge1 = 0;
//             }
//             if (node_k != -1)
//             {
//                 rod->getEntry(node_k, offset_k);
//                 rod->x(node_k, xk);
//                 rod->XdX(node_k, Xk, dXk);
//                 edge0 = on_rod_idx[rod_idx]-1;
//                 if(rod->closed && on_rod_idx[rod_idx] == 0)
//                     edge0 = rod->numSeg() - 1;
//             }

//             T theta_edge0 = rod->reference_angles[edge0];
//             T theta_edge1 = rod->reference_angles[edge1];

//             Matrix<T, 2, 2> B = rod->bending_coeffs ;
//             T kt = rod->kt;
//             T reference_twist_edge0 = rod->reference_twist[edge0];
//             T reference_twist_edge1 = rod->reference_twist[edge1];

//             TV referenceNormal1 = rod->reference_frame_us[edge0];
//             TV referenceNormal2 = rod->reference_frame_us[edge1];

//             TV referenceTangent1 = rod->prev_tangents[edge0];
//             TV referenceTangent2 = rod->prev_tangents[edge1];            
            
//             TV rest_tangent1 = rod->rest_tangents[edge0];
//             TV rest_tangent2 = rod->rest_tangents[edge1];
//             TV rest_normal1 = rod->rest_normals[edge0];
//             TV rest_normal2 = rod->rest_normals[edge1];

//             if (crossing->is_fixed)
//             {
//                 Matrix<T, 16, 16> J0, J1;
//                 J0.setZero(); J1.setZero();
//                 // xi is the rigid body
                
//                 int dof_theta0 = rod->theta_dof_start_offset + edge0;
//                 int dof_theta1 = rod->theta_dof_start_offset + edge1;

//                 Offset omega_dof;
//                 omega_dof[0] = crossing->dof_offset;
//                 omega_dof[1] = crossing->dof_offset + 1;
//                 omega_dof[2] = crossing->dof_offset + 2;

//                 std::vector<Offset> offsets1 = { offset_i, offset_j };
//                 std::vector<Offset> offsets2 = { offset_i, offset_k };
//                 std::vector<TV> dXdu1 = { dXi, dXj };
//                 std::vector<TV> dXdu2 = { dXi, dXk };
                
//                 int entry_cnt = 0;
//                 if (node_j != -1)
//                 {
//                     computeRodEulerAngleBendingAndTwistEnergyRBFirstHessian(B, kt, 0.0, 
//                         rest_tangent1, rest_normal1, 
//                         referenceTangent2, referenceNormal2, 
//                         reference_twist_edge1, xi, xj, Xi, Xj, omega_acc, omega, theta_edge1, J0);
                    
//                     J0 *= 0.5;
//                     for(int k = 0; k < offsets1.size(); k++)
//                     {
//                         //dx/dx
//                         //dx/dX dX/dx
//                         for(int l = 0; l < offsets1.size(); l++)
//                             for (int i = 0; i < 3; i++)
//                                 for (int j = 0; j < 3; j++)
//                                 {
//                                     entry_K.emplace_back(offsets1[k][i], offsets1[l][j], J0(k*3 + i, l * 3 + j));
//                                     entry_K.emplace_back(offsets1[k][i], offsets1[l][3], J0(k*3 + i, 2 * 3 + l * 3 + j) * dXdu1[l][j]);
//                                     // entry_K.emplace_back(offsets1[l][3], offsets1[k][i], J0(2 * 3 + l * 3 + j, k*3 + i) * dXdu1[l][j]));
//                                     entry_K.emplace_back(offsets1[k][3], offsets1[l][j], J0(2 * 3 + k * 3 + i, l*3 + j) * dXdu1[k][i]);
//                                     entry_K.emplace_back(offsets1[k][3], offsets1[l][3], J0(2 * 3 + k*3 + i, 2 * 3 + l * 3 + j) * dXdu1[l][j] * dXdu1[k][i]);
//                                     entry_cnt+=4;
//                                 }
//                         for (int i = 0; i < 3; i++)
//                         {
//                             for (int j = 0; j < 3; j++)
//                             {
//                                 entry_K.emplace_back(offsets1[k][i], omega_dof[j], J0(k*3 + i, 4 * 3 + j));
//                                 entry_K.emplace_back(omega_dof[j], offsets1[k][i], J0(4 * 3 + j, k*3 + i));

//                                 entry_K.emplace_back(offsets1[k][3], omega_dof[j], J0(2 *3 + k*3 + i, 4 * 3 + j) * dXdu1[k][i]);
//                                 entry_K.emplace_back(omega_dof[j], offsets1[k][3], J0(4 * 3 + j, 2 *3 + k*3 + i) * dXdu1[k][i]);
//                                 entry_cnt+=4;
//                             }

//                             entry_K.emplace_back(offsets1[k][i], dof_theta1, J0(k*3 + i, 5 * 3));
//                             entry_K.emplace_back(dof_theta1, offsets1[k][i], J0(5 * 3, k*3 + i));

//                             entry_K.emplace_back(offsets1[k][3], dof_theta1, J0(2 * 3 + k*3 + i, 5 * 3) * dXdu1[k][i]);
//                             entry_K.emplace_back(dof_theta1, offsets1[k][3], J0(5 * 3, 2 * 3 + k*3 + i) * dXdu1[k][i]);
//                             entry_cnt+=4;
//                         }

//                     }
//                     for (int i = 0; i < 3; i++)
//                     {
//                         entry_K.emplace_back(omega_dof[i], dof_theta1, J0(4 * 3 + i, 5 * 3));
//                         entry_K.emplace_back(dof_theta1, omega_dof[i], J0(5 * 3, 4 * 3 + i));
//                         entry_cnt+=2;
//                         for (int j = 0; j < 3; j++)
//                         {
//                             entry_K.emplace_back(omega_dof[i], omega_dof[j], J0(4 * 3 + i, 4 * 3 + j));
//                             entry_cnt+=1;
//                         }
//                     }
//                     entry_K.emplace_back(dof_theta1, dof_theta1, J0(5 * 3, 5 * 3));
//                     entry_cnt+=1;

//                     // std::cout << 16 * 16 << " " << entry_cnt << std::endl;
//                     // std::getchar();
//                 }
//                 if (node_k != -1)
//                 {
//                     computeRodEulerAngleBendingAndTwistEnergyRBSecondHessian(B, kt, 0.0, 
//                         referenceTangent1, referenceNormal1,
//                         rest_tangent2, rest_normal2, 
//                         reference_twist_edge0, xi, xk, Xi, Xk, omega_acc, omega, theta_edge0, J1);
                    
//                     J1 *= 0.5;
//                     for(int k = 0; k < offsets1.size(); k++)
//                     {
//                         //dx/dx
//                         //dx/dX dX/dx
//                         for(int l = 0; l < offsets1.size(); l++)
//                             for (int i = 0; i < 3; i++)
//                                 for (int j = 0; j < 3; j++)
//                                 {
//                                     entry_K.emplace_back(offsets2[k][i], offsets2[l][j], J1(k*3 + i, l * 3 + j));
//                                     entry_K.emplace_back(offsets2[k][i], offsets2[l][3], J1(k*3 + i, 2 * 3 + l * 3 + j) * dXdu2[l][j]);
//                                     entry_K.emplace_back(offsets2[l][3], offsets2[k][i], J1(2 * 3 + l * 3 + j, k*3 + i) * dXdu2[l][j]);
//                                     entry_K.emplace_back(offsets2[k][3], offsets2[l][3], J1(2 * 3 + k*3 + i, 2 * 3 + l * 3 + j) * dXdu2[l][j] * dXdu2[k][i]);
//                                 }
//                         for (int i = 0; i < 3; i++)
//                         {
//                             for (int j = 0; j < 3; j++)
//                             {
//                                 entry_K.emplace_back(offsets2[k][i], omega_dof[j], J1(k*3 + i, 4 * 3 + j));
//                                 entry_K.emplace_back(omega_dof[j], offsets2[k][i], J1(4 * 3 + j, k*3 + i));
//                                 entry_K.emplace_back(offsets2[k][3], omega_dof[j], J1(2 *3 + k*3 + i, 4 * 3 + j) * dXdu2[k][i]);
//                                 entry_K.emplace_back(omega_dof[j], offsets2[k][3], J1(4 * 3 + j, 2 *3 + k*3 + i) * dXdu2[k][i]);
//                             }

//                             entry_K.emplace_back(dof_theta0, offsets2[k][i], J1(5 * 3, k*3 + i));
//                             entry_K.emplace_back(offsets2[k][i], dof_theta0, J1(k*3 + i, 5 * 3));

//                             entry_K.emplace_back(offsets2[k][3], dof_theta0, J1(2 * 3 + k*3 + i, 5 * 3) * dXdu2[k][i]);
//                             entry_K.emplace_back(dof_theta0, offsets2[k][3], J1(5 * 3, 2 * 3 + k*3 + i) * dXdu2[k][i]);
//                         }

//                     }
//                     for (int i = 0; i < 3; i++)
//                     {
//                         entry_K.emplace_back(omega_dof[i], dof_theta0, J1(4 * 3 + i, 5 * 3));
//                         entry_K.emplace_back(dof_theta0, omega_dof[i], J1(5 * 3, 4 * 3 + i));

//                         for (int j = 0; j < 3; j++)
//                         {
//                             entry_K.emplace_back(omega_dof[i], omega_dof[j], J1(4 * 3 + i, 4 * 3 + j));
//                         }
//                     }
//                     entry_K.emplace_back(dof_theta0, dof_theta0, J1(5 * 3, 5 * 3));
//                 }
//                 // std::cout << "here here" << std::endl;
//                 // std::getchar();
//             }
//             else
//             {
//                 if (node_k != -1 && node_j != -1)
//                 {
//                     // std::cout << node_i << std::endl;
//                     std::vector<int> nodes = { node_k, node_i, node_j };
//                     std::vector<TV> dXdu = { dXk, dXi, dXj };
//                     std::vector<Offset> offsets = { offset_k, offset_i, offset_j };

                    
//                     std::vector<int> theta_dofs = {rod->theta_dof_start_offset + edge0, rod->theta_dof_start_offset + edge1};

//                     Matrix<T, 20, 20> J;
                
//                     computeRodBendingAndTwistEnergyHessian(B, rod->kt, 0.0, referenceNormal1, referenceTangent1,
//                         referenceNormal2, referenceTangent2,  reference_twist_edge1, xk, xi, xj, Xk, Xi, Xj, theta_edge0, theta_edge1, J);

//                     J *= -1.0;
//                     for(int k = 0; k < nodes.size(); k++)
//                         for(int l = 0; l < nodes.size(); l++)
//                             for(int i = 0; i < 3; i++)
//                                 for (int j = 0; j < 3; j++)
//                                     {
//                                         entry_K.emplace_back(offsets[k][i], offsets[l][j], -J(k*3 + i, l * 3 + j));

//                                         entry_K.emplace_back(offsets[k][i], offsets[l][3], -J(k*3 + i, 3 * 3 + l * 3 + j) * dXdu[l][j]);
                                        
//                                         entry_K.emplace_back(offsets[k][3], offsets[l][j], -J(3 * 3 + k * 3 + i, l * 3 + j) * dXdu[k][i]);

                                        
//                                         entry_K.emplace_back(offsets[k][3], 
//                                                                             offsets[l][3], 
//                                                                             -J(3 * 3 + k*3 + i, 3 * 3 + l * 3 + j) * dXdu[l][j] * dXdu[k][i]);

//                                     }
//                     for(int k = 0; k < nodes.size(); k++)
//                         for (int j = 0; j < 2; j++)
//                             for(int i = 0; i < 3; i++)
//                             {
//                                 entry_K.emplace_back(offsets[k][i], theta_dofs[j], -J(k*3 + i, 18 + j));
//                                 entry_K.emplace_back(theta_dofs[j], offsets[k][i], -J(18 + j, k * 3 + i));

//                                 entry_K.emplace_back(offsets[k][3], theta_dofs[j], -J(3 * 3 + k * 3 + i, 18 + j) * dXdu[k][i]);
//                                 entry_K.emplace_back(theta_dofs[j], offsets[k][3], -J(18 + j, 3 * 3 + k * 3 + i) * dXdu[k][i]);
//                             }
//                     for (int i = 0; i < 2; i++)
//                         for (int j = 0; j < 2; j++)
//                             entry_K.emplace_back(theta_dofs[i], theta_dofs[j], -J(18 + i, 18 + j));
//                 }
                
//             }
            

//             cnt ++;
//         }
//     }
// }


// void RodNetwork::addJointBendingAndTwistingForce(Eigen::Ref<VectorXT> residual)
// {
//     VectorXT residual_cp = residual;
//     for (auto& crossing : rod_crossings)
//     {
//         int node_i = crossing->node_idx;
        
//         std::vector<int> rods_involved = crossing->rods_involved;
//         std::unordered_map<int, int> on_rod_idx = crossing->on_rod_idx;

//         std::vector<Vector<T, 2>> sliding_ranges = crossing->sliding_ranges;
//         Vector<T, 3> omega = crossing->omega;
//         // Vector<T, 4> omega_acc = crossing->omega_acc;
//         Matrix<T, 3, 3> omega_acc = crossing->rotation_accumulated;
//         // Vector<T, 3> omega_acc = crossing->omega_acc;

//         int cnt = 0;
//         for (int rod_idx : rods_involved)
//         {
//             Offset offset_i, offset_j, offset_k;
//             TV xi, xj, xk, Xi, Xj, Xk, dXi, dXj, dXk;
            
//             auto rod = Rods[rod_idx];
//             rod->getEntry(node_i, offset_i);
//             rod->x(node_i, xi); rod->XdX(node_i, Xi, dXi);

//             int node_j = rod->nodeIdx(on_rod_idx[rod_idx] + 1);
//             int node_k = rod->nodeIdx(on_rod_idx[rod_idx] - 1);
//             int edge1 = 0, edge0 = 0;
//             if (node_j != -1)
//             {
//                 rod->getEntry(node_j, offset_j);
//                 rod->x(node_j, xj);
//                 rod->XdX(node_j, Xj, dXj);
//                 edge1 = on_rod_idx[rod_idx];
//                 if(rod->closed && on_rod_idx[rod_idx] == 0)
//                     edge1 = 0;
//             }
//             if (node_k != -1)
//             {
//                 rod->getEntry(node_k, offset_k);
//                 rod->x(node_k, xk);
//                 rod->XdX(node_k, Xk, dXk);
//                 edge0 = on_rod_idx[rod_idx]-1;
//                 if(rod->closed && on_rod_idx[rod_idx] == 0)
//                     edge0 = rod->numSeg() - 1;

//             }

//             // std::cout << " node i " << node_i << " node j " << node_j << " node k " << node_k << std::endl;
//             // std::cout << " on_rod_idx[rod_idx] " << on_rod_idx[rod_idx] << " " << edge0 << " " << edge1 << std::endl;

//             T theta_edge0 = rod->reference_angles[edge0];
//             T theta_edge1 = rod->reference_angles[edge1];

//             Matrix<T, 2, 2> B = rod->bending_coeffs ;
//             T kt = rod->kt;
//             T reference_twist_edge0 = rod->reference_twist[edge0];
//             T reference_twist_edge1 = rod->reference_twist[edge1];

//             TV referenceNormal1 = rod->reference_frame_us[edge0];
//             TV referenceNormal2 = rod->reference_frame_us[edge1];

//             TV referenceTangent1 = rod->prev_tangents[edge0];
//             TV referenceTangent2 = rod->prev_tangents[edge1];       

//             TV rest_tangent1 = rod->rest_tangents[edge0];
//             TV rest_tangent2 = rod->rest_tangents[edge1];
            
//             TV rest_normal1 = rod->rest_normals[edge0];
//             TV rest_normal2 = rod->rest_normals[edge1];
            
//             if (crossing->is_fixed)
//             {
//                 Vector<T, 16> F0, F1;
//                 F0.setZero(); F1.setZero();
//                 // xi is the rigid body
                
//                 if (node_j != -1)
//                 {
                    
//                     computeRodEulerAngleBendingAndTwistEnergyRBFirstGradient(B, kt, 0.0, 
//                         rest_tangent1, rest_normal1, 
//                         referenceTangent2, referenceNormal2, 
//                         reference_twist_edge1, xi, xj, Xi, Xj, omega_acc, omega, theta_edge1, F0);

//                     // std::cout << node_i << " " << node_j << " " << node_k << std::endl;
//                     // std::cout << xi.transpose() << " "<< xj.transpose() << " "<< xk.transpose() << std::endl;
//                     // std::cout << Xi.transpose() << " "<< Xj.transpose() << " "<< Xk.transpose() << std::endl;
//                     // std::cout << referenceTangent2.transpose() << " " << referenceTangent1.transpose() << " "<< referenceNormal1.transpose() << " "<< referenceNormal2.transpose() << std::endl;
//                     // std::cout << "F0" << std::endl;
//                     // std::cout << F0 << std::endl;
//                     // std::getchar();

//                     F0 *= -0.5;
//                     residual.template segment<3>(offset_i[0]) += F0.template segment<3>(0);
//                     residual.template segment<3>(offset_j[0]) += F0.template segment<3>(3);
//                     residual(offset_i[3]) += F0.template segment<3>(2*3).dot(dXi);
//                     residual(offset_j[3]) += F0.template segment<3>(3*3).dot(dXj);
//                     residual.template segment<3>(crossing->dof_offset) += F0.template segment<3>(4*3);
//                     residual[rod->theta_dof_start_offset + edge1] += F0[5*3];
//                 }
//                 if (node_k != -1)                    
//                 {
                    

//                     computeRodEulerAngleBendingAndTwistEnergyRBSecondGradient(B, kt, 0.0, 
//                         referenceTangent1, referenceNormal1,
//                         rest_tangent2, rest_normal2, 
//                         reference_twist_edge0, xi, xk, Xi, Xk, omega_acc, omega, theta_edge0, F1);

//                     // std::cout << node_i << " " << node_k << " " << node_k << std::endl;
//                     // std::cout << xi.transpose() << " "<< xj.transpose() << " "<< xk.transpose() << std::endl;
//                     // std::cout << Xi.transpose() << " "<< Xj.transpose() << " "<< Xk.transpose() << std::endl;
//                     // std::cout << referenceTangent2.transpose() << " " << referenceTangent1.transpose() << " "<< referenceNormal1.transpose() << " "<< referenceNormal2.transpose() << std::endl;
//                     // std::cout << "F1" << std::endl;
//                     // std::cout << F1 << std::endl;
//                     // std::getchar();
                    
//                     F1 *= -0.5;
//                     residual.template segment<3>(offset_i[0]) += F1.template segment<3>(0);
//                     residual.template segment<3>(offset_k[0]) += F1.template segment<3>(3);
//                     residual(offset_i[3]) += F1.template segment<3>(2*3).dot(dXi);
//                     residual(offset_k[3]) += F1.template segment<3>(3*3).dot(dXk);
//                     residual.template segment<3>(crossing->dof_offset) += F1.template segment<3>(4*3);
//                     residual[rod->theta_dof_start_offset + edge0] += F1[5*3];
//                 }
                
//             }
//             else
//             {
//                 if (node_k != -1 && node_j != -1)
//                 {
//                     Vector<T, 20> F;
//                     computeRodBendingAndTwistEnergyGradient(B, rod->kt, 0.0, referenceNormal1, referenceTangent1,
//                         referenceNormal2, referenceTangent2,  reference_twist_edge1, xk, xi, xj, Xk, Xi, Xj, theta_edge0, theta_edge1, F);
//                     // std::cout << reference_twist_edge1 << " " << theta_edge0 << " " << theta_edge1 << std::endl;
//                     // std::cout << node_i << " " << node_j << " " << node_k << std::endl;
//                     // std::cout << xi.transpose() << " "<< xj.transpose() << " "<< xk.transpose() << std::endl;
//                     // std::cout << Xi.transpose() << " "<< Xj.transpose() << " "<< Xk.transpose() << std::endl;
//                     // std::cout << referenceTangent2.transpose() << " " << referenceTangent1.transpose() << " "<< referenceNormal1.transpose() << " "<< referenceNormal2.transpose() << std::endl;
//                     // std::cout << "F" << std::endl;
//                     // std::cout << F << std::endl;
//                     // std::getchar();

//                     F *= -1.0;

//                     residual.template segment<3>(offset_k[0]) += F.template segment<3>(0);
//                     residual.template segment<3>(offset_i[0]) += F.template segment<3>(3);
//                     residual.template segment<3>(offset_j[0]) += F.template segment<3>(3 + 3);

//                     residual(offset_k[3]) += F.template segment<3>(3*3).dot(dXk);
//                     residual(offset_i[3]) += F.template segment<3>(3*3 + 3).dot(dXi);
//                     residual(offset_j[3]) += F.template segment<3>(3*3 + 2*3).dot(dXj);

                    
//                     residual[rod->theta_dof_start_offset + edge0] += F[18];
//                     residual[rod->theta_dof_start_offset + edge1] += F[19];
                    
//                 }
//             }

//             cnt ++;
//         }
//     }
//     if (print_force_mag)
//         std::cout << "joint force norm: " << (residual - residual_cp).norm() << std::endl;
// }


// T RodNetwork::addJointBendingAndTwistingEnergy(bool bending, bool twisting)
// {
//     T energy = 0.0;
//     for (auto& crossing : rod_crossings)
//     {
//         int node_i = crossing->node_idx;
        
//         std::vector<int> rods_involved = crossing->rods_involved;
//         std::unordered_map<int, int> on_rod_idx = crossing->on_rod_idx;

//         std::vector<Vector<T, 2>> sliding_ranges = crossing->sliding_ranges;
//         Vector<T, 3> omega = crossing->omega;
//         Matrix<T, 3, 3> omega_acc = crossing->rotation_accumulated;
    
        
//         int cnt = 0;
//         for (int rod_idx : rods_involved)
//         {
//             Offset offset_i, offset_j, offset_k;
//             TV xi = TV::Zero(), xj = TV::Zero(), xk = TV::Zero(), 
//                 Xi = TV::Zero(), Xj = TV::Zero(), Xk = TV::Zero();

//             auto rod = Rods[rod_idx];
//             rod->getEntry(node_i, offset_i);
//             rod->x(node_i, xi); rod->X(node_i, Xi);

//             int node_j = rod->nodeIdx(on_rod_idx[rod_idx] + 1);
//             int node_k = rod->nodeIdx(on_rod_idx[rod_idx] - 1);
//             int edge1 = 0, edge0 = 0;
            
//             if (node_j != -1)
//             {
//                 rod->getEntry(node_j, offset_j);
//                 rod->x(node_j, xj);
//                 rod->X(node_j, Xj);
//                 edge1 = on_rod_idx[rod_idx];
//                 if(rod->closed && on_rod_idx[rod_idx] == 0)
//                     edge1 = 0;
//             }
//             if (node_k != -1)
//             {
//                 rod->getEntry(node_k, offset_k);
//                 rod->x(node_k, xk);
//                 rod->X(node_k, Xk);
//                 edge0 = on_rod_idx[rod_idx]-1;
//                 if(rod->closed && on_rod_idx[rod_idx] == 0)
//                     edge0 = rod->numSeg() - 1;
//             }

//             T theta_edge0 = rod->reference_angles[edge0];
//             T theta_edge1 = rod->reference_angles[edge1];

//             Matrix<T, 2, 2> B = bending ? rod->bending_coeffs : Matrix<T, 2, 2>::Zero();
//             T kt = twisting ? rod->kt : 0.0;

//             T reference_twist_edge0 = rod->reference_twist[edge0];
//             T reference_twist_edge1 = rod->reference_twist[edge1];

//             TV referenceNormal1 = rod->reference_frame_us[edge0];
//             TV referenceNormal2 = rod->reference_frame_us[edge1];

//             TV referenceTangent1 = rod->prev_tangents[edge0];
//             TV referenceTangent2 = rod->prev_tangents[edge1];       

//             TV rest_tangent1 = rod->rest_tangents[edge0];
//             TV rest_tangent2 = rod->rest_tangents[edge1];
//             TV rest_normal1 = rod->rest_normals[edge0];
//             TV rest_normal2 = rod->rest_normals[edge1];

//             if (crossing->is_fixed)
//             {
//                 // xi is the rigid body
//                 if (node_j != -1)
//                 {
//                     T E = 0.5 * computeRodEulerAngleBendingAndTwistEnergyRBFirst(B, kt, 0.0, 
//                         rest_tangent1, rest_normal1, 
//                         referenceTangent2, referenceNormal2, 
//                         reference_twist_edge1, xi, xj, Xi, Xj, omega_acc, omega, theta_edge1);
//                     energy += E;
                    
//                     // std::cout << E << std::endl;
//                     // std::getchar();
//                 }
//                 if (node_k != -1)
//                 {
//                     T E = 0.5 * computeRodEulerAngleBendingAndTwistEnergyRBSecond(B, kt, 0.0, 
//                         referenceTangent1, referenceNormal1,
//                         rest_tangent2, rest_normal2, 
//                         reference_twist_edge0, xi, xk, Xi, Xk, omega_acc, omega, theta_edge0);
//                     energy += E;
//                 }

                
//             }
//             else
//             {
//                 if (node_j != -1 && node_k != -1)
//                 {
//                     T E = computeRodBendingAndTwistEnergy(B, rod->kt, 0.0, referenceNormal1, referenceTangent1,
//                         referenceNormal2, referenceTangent2, reference_twist_edge1, xk, xi, xj, Xk, Xi, Xj, theta_edge0, theta_edge1);
                    
//                     energy += E;
//                 }
//             }
            
//             cnt++;
//         }
//     }
//     return energy;
// }