#include "../include/EoLRodSim.h"


void EoLRodSim::addRegularizingK(std::vector<Eigen::Triplet<T>>& entry_K)
{
    for (auto& crossing : rod_crossings)
    {
        if (crossing->is_fixed)
            continue;
        int node_idx = crossing->node_idx;
        std::vector<int> rods_involved = crossing->rods_involved;
        for (int rod_idx : rods_involved)
        {
            Offset offset;
            Rods[rod_idx]->getEntry(node_idx, offset);
            entry_K.push_back(Entry(offset[3], offset[3], ke));
        }
    }
}


void EoLRodSim::addRegularizingForce(Eigen::Ref<VectorXT> residual)
{
    
    VectorXT residual_cp = residual;
    for (auto& crossing : rod_crossings)
    {
        if (crossing->is_fixed)
            continue;
        
        int node_idx = crossing->node_idx;
        std::vector<int> rods_involved = crossing->rods_involved;
        for (int rod_idx : rods_involved)
        {
            Offset offset;
            Rods[rod_idx]->getEntry(node_idx, offset);
            T u, U;
            Rods[rod_idx]->u(node_idx, u);
            Rods[rod_idx]->U(node_idx, U);
            
            residual[offset[3]] += -ke * (u-U);
        }
    }
    

    if(print_force_mag)
        std::cout << "Eulerian penalty norm: " << (residual - residual_cp).norm() << std::endl;
}



T EoLRodSim::addRegularizingEnergy()
{
    T energy = 0.0;
    for (auto& crossing : rod_crossings)
    {
        if (crossing->is_fixed)
            continue;
        int node_idx = crossing->node_idx;
        std::vector<int> rods_involved = crossing->rods_involved;
        for (int rod_idx : rods_involved)
        {
            Offset offset;
            Rods[rod_idx]->getEntry(node_idx, offset);
            T u, U;
            Rods[rod_idx]->u(node_idx, u);
            Rods[rod_idx]->U(node_idx, U);
            energy += 0.5 * ke * std::pow(u - U, 2);
        }
    }


    return energy;
}
