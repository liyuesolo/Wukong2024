#include "../include/EoLRodSim.h"
#include "../autodiff/EoLRodStretchingEnergy.h"


void EoLRodSim::addStretchingK(std::vector<Entry>& entry_K)
{
    for (auto& rod : Rods)
    {
        rod->iterateSegmentsWithOffset([&](int node_i, int node_j, Offset offset_i, Offset offset_j, int rod_idx)
        {
            TV xi, xj, Xi, Xj, dXi, dXj;
            rod->x(node_i, xi); rod->x(node_j, xj);
            rod->XdX(node_i, Xi, dXi); rod->XdX(node_j, Xj, dXj);

            std::vector<TV> x(4);
            x[0] = xi; x[1] = xj; x[2] = Xi; x[3] = Xj;

            // T J[12][12];
            // memset(J, 0, sizeof(J));

            Matrix<T, 12, 12> J;
            computeEoLRodStretchingEnergyHessian(rod->ks, Xi, Xj, xi, xj, J);
            J *= -1.0;
            
            std::vector<int> nodes = { node_i, node_j };
            std::vector<TV> dXdu = { dXi, dXj };

            std::vector<Offset> offsets = { offset_i, offset_j };
            
            for(int k = 0; k < nodes.size(); k++)
                for(int l = 0; l < nodes.size(); l++)
                    for(int i = 0; i < 3; i++)
                        for (int j = 0; j < 3; j++)
                            {
                                // entry_K.push_back(Entry(offsets[k][i], offsets[l][j], -J[k*3 + i][l * 3 + j]));
                                
                                // entry_K.push_back(Entry(offsets[k][i], offsets[l][3], -J[k*3 + i][2 * 3 + l * 3 + j] * dXdu[l][j]));

                                // entry_K.push_back(Entry(offsets[k][3], offsets[l][j] , -J[2 * 3 + k * 3 + i][l * 3 + j] * dXdu[k][i]));

                                // // dX/du ^T d2e/dX2 dXdu                    
                                // entry_K.push_back(Entry(offsets[k][3], 
                                //                         offsets[l][3], 
                                //                         -J[2 * 3 + k*3 + i][2 * 3 + l * 3 + j] * dXdu[l][j] * dXdu[k][i]));

                                entry_K.push_back(Entry(offsets[k][i], offsets[l][j], -J(k*3 + i, l * 3 + j)));
                                
                                entry_K.push_back(Entry(offsets[k][i], offsets[l][3], -J(k*3 + i, 2 * 3 + l * 3 + j) * dXdu[l][j]));

                                entry_K.push_back(Entry(offsets[k][3], offsets[l][j] , -J(2 * 3 + k * 3 + i, l * 3 + j) * dXdu[k][i]));

                                // dX/du ^T d2e/dX2 dXdu                    
                                entry_K.push_back(Entry(offsets[k][3], 
                                                        offsets[l][3], 
                                                        -J(2 * 3 + k*3 + i, 2 * 3 + l * 3 + j) * dXdu[l][j] * dXdu[k][i]));
                            }
        });   
    }
    // for (auto& rod : Rods)
    // {
    //     rod->iterateSegmentsWithOffset([&](int node_i, int node_j, Offset offset_i, Offset offset_j, int rod_idx)
    //     {
    //         // std::cout << "node i " << node_i << " node j " << node_j << std::endl;
    //         TV xi, xj, Xi, Xj, dXi, dXj, ddXi, ddXj;
    //         rod->x(node_i, xi); rod->x(node_j, xj);
    //         rod->XdXddX(node_i, Xi, dXi, ddXi); rod->XdXddX(node_j, Xj, dXj, ddXj);

    //         std::vector<TV> x(4);
    //         x[0] = xi; x[1] = xj; x[2] = Xi; x[3] = Xj;

    //         Vector<T, 3 * 4> F;
    //         F.setZero();
    //         if constexpr (3 == 2)
    //         {
    //             // #include "Maple/YarnStretchingF.mcg"
    //         }
    //         else if constexpr (3 == 3)
    //         {
    //             // #include "Maple/RodStretching3DF.mcg"
    //             computeEoLRodStretchingEnergyGradient(rod->ks, Xi, Xj, xi, xj, F);
    //         }

    //         for(int d = 0; d < 3; d++)
    //         {
    //             entry_K.push_back(Eigen::Triplet<T>(offset_i[3], 
    //                                 offset_i[3], -F[0*3 + 2*3 + d] * ddXi[d]));
    //             entry_K.push_back(Eigen::Triplet<T>(offset_j[3], 
    //                                 offset_j[3], -F[1*3 + 2*3 + d] * ddXj[d]));
    //         }    
    //     });
    // }
}

void EoLRodSim::addStretchingForce(Eigen::Ref<VectorXT> residual)
{
    
    VectorXT residual_cp = residual;
    for (auto& rod : Rods)
    {
        rod->iterateSegmentsWithOffset([&](int node_i, int node_j, Offset offset_i, Offset offset_j, int rod_idx)
        {
            // std::cout << "node i " << node_i << " node j " << node_j << std::endl;
            // std::cout << "node i " << offset_i.transpose() << " node j " << offset_j.transpose() << std::endl;
            TV xi, xj, Xi, Xj, dXi, dXj;
            rod->x(node_i, xi); rod->x(node_j, xj);
            rod->XdX(node_i, Xi, dXi); rod->XdX(node_j, Xj, dXj);
            
            std::vector<TV> x(4);
            x[0] = xi; x[1] = xj; x[2] = Xi; x[3] = Xj;

            Vector<T, 3 * 4> F;
            F.setZero();
            computeEoLRodStretchingEnergyGradient(rod->ks, Xi, Xj, xi, xj, F);
            F *= -1.0;
            
            residual.template segment<3>(offset_i[0]) += F.template segment<3>(0);
            residual.template segment<3>(offset_j[0]) += F.template segment<3>(3);

            residual(offset_i[3]) += F.template segment<3>(2*3).dot(dXi);
            residual(offset_j[3]) += F.template segment<3>(2*3 + 3).dot(dXj);
        });
    }
    // std::cout << (residual - residual_cp) << std::endl;
    // std::getchar();
    if(print_force_mag)
        std::cout << "stretching norm: " << (residual - residual_cp).norm() << std::endl;
}


T EoLRodSim::addStretchingEnergy()
{
    T E = 0;
    for (auto& rod : Rods)
    {
        rod->iterateSegments([&](int node_i, int node_j, int rod_idx)
        {
            TV xi, xj, Xi, Xj;
            rod->x(node_i, xi); rod->x(node_j, xj);
            rod->X(node_i, Xi); rod->X(node_j, Xj);

            std::vector<TV> x(4);
            x[0] = xi; x[1] = xj; x[2] = Xi; x[3] = Xj;

            T V[1];
            E += stretchingEnergyLocal(rod->ks, Xi, Xj, xi, xj);
        });
    }
    return E;
}