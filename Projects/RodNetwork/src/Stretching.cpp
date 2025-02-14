#include "../include/RodNetwork.h"
#include "../autodiff/RodStretchingEnergy.h"


void RodNetwork::addStretchingHessian(std::vector<Entry>& entry_K)
{
    for (auto& rod : rods)
    {
        rod->iterateSegments([&](int node_i, int node_j, int rod_idx)
        {
            TV xi, xj, Xi, Xj, dXi, dXj;
            rod->x(node_i, xi); rod->x(node_j, xj);
            rod->X(node_i, Xi); rod->X(node_j, Xj);

            std::vector<TV> x(4);
            x[0] = xi; x[1] = xj; x[2] = Xi; x[3] = Xj;

    
            Matrix<T, 6, 6> d2edx2;
            computeRodStretchingEnergyHessian(rod->ks, Xi, Xj, xi, xj, d2edx2);

            addHessianEntry<3, 3>(entry_K, { node_i, node_j }, d2edx2, 0, 0);
            
        });   
    }
    
}

void RodNetwork::addStretchingForce(Eigen::Ref<VectorXT> residual)
{
    
    for (auto& rod : rods)
    {
        rod->iterateSegments([&](int node_i, int node_j, int rod_idx)
        {
            std::cout << "node i " << node_i << " node j " << node_j << std::endl;
            // std::cout << "node i " << offset_i.transpose() << " node j " << offset_j.transpose() << std::endl;
            TV xi, xj, Xi, Xj, dXi, dXj;
            rod->x(node_i, xi); rod->x(node_j, xj);
            rod->X(node_i, Xi); rod->X(node_j, Xj);
            
            std::cout << "xi " << xi.transpose() << " xj " << xj.transpose() << " Xi " << Xi.transpose() << " Xj " << Xj.transpose() << std::endl;
            std::cout << "|X| " << (Xi - Xj).norm() << " |x| " << (xi - xj).norm() << std::endl;


            std::vector<TV> x(4);
            x[0] = xi; x[1] = xj; x[2] = Xi; x[3] = Xj;

            Vector<T, 6> dedx;
            computeRodStretchingEnergyGradient(rod->ks, Xi, Xj, xi, xj, dedx);

            std::cout << dedx.transpose() << std::endl;
            if ((Xi - Xj).norm() < 1e-6)
                std::getchar();
            
            residual.segment<3>(node_i * 3) -= dedx.segment<3>(0);
            residual.segment<3>(node_j * 3) -= dedx.segment<3>(3);

        });
    }
}


T RodNetwork::addStretchingEnergy()
{
    T E = 0;
    for (auto& rod : rods)
    {
        rod->iterateSegments([&](int node_i, int node_j, int rod_idx)
        {
            TV xi, xj, Xi, Xj;
            rod->x(node_i, xi); rod->x(node_j, xj);
            rod->X(node_i, Xi); rod->X(node_j, Xj);

            std::vector<TV> x(4);
            x[0] = xi; x[1] = xj; x[2] = Xi; x[3] = Xj;

            T energy;
            computeRodStretchingEnergy(rod->ks, Xi, Xj, xi, xj, energy);
        });
    }
    return E;
}