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

            Matrix<T, 6, 6> d2edx2;
            computeRodStretchingEnergyHessian(rod->ks, Xi, Xj, xi, xj, d2edx2);

            addHessianEntry<3, 3>(entry_K, { node_i, node_j }, d2edx2, 0, 0);
            
        });   
    }
    
}

void RodNetwork::addStretchingForce(VectorXT& residual)
{
    
    for (auto& rod : rods)
    {
        rod->iterateSegments([&](int node_i, int node_j, int rod_idx)
        {
            // std::cout << "node i " << node_i << " node j " << node_j << std::endl;
            // std::cout << "node i " << offset_i.transpose() << " node j " << offset_j.transpose() << std::endl;
            TV xi, xj, Xi, Xj, dXi, dXj;
            rod->x(node_i, xi); rod->x(node_j, xj);
            rod->X(node_i, Xi); rod->X(node_j, Xj);

            Vector<T, 6> dedx;
            computeRodStretchingEnergyGradient(rod->ks, Xi, Xj, xi, xj, dedx);
            
            addForceEntry<3>(residual, { node_i, node_j }, -dedx, 0);

        });
        // break;
    }
}


T RodNetwork::addStretchingEnergy()
{
    T E = 0;
    for (auto& rod : rods)
    {
        rod->iterateSegments([&](int node_i, int node_j, int rod_idx)
        {
            // std::cout << "node i " << node_i << " node j " << node_j << std::endl;
            TV xi, xj, Xi, Xj;
            rod->x(node_i, xi); rod->x(node_j, xj);
            rod->X(node_i, Xi); rod->X(node_j, Xj);

            T energy;
            computeRodStretchingEnergy(rod->ks, Xi, Xj, xi, xj, energy);
            E += energy;
        });
        // break;
    }
    return E;
}