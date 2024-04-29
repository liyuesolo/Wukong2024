#include "../include/EoLRodSim.h"

static double min_du = 0.01;


void EoLRodSim::addParallelContactK(std::vector<Entry>& entry_K)
{
    if (!add_contact_penalty)
        return;
    for (auto& crossing : rod_crossings)
    {
        if (crossing->is_fixed)
            continue;
        int node_idx = crossing->node_idx;
        std::vector<int> rods_involved = crossing->rods_involved;
        std::vector<Vector<T, 2>> sliding_ranges = crossing->sliding_ranges;

        int cnt = 0;
        for (int rod_idx : rods_involved)
        {
            Offset offset;
            Rods[rod_idx]->getEntry(node_idx, offset);
            T u, U, u_left, u_right;
            
            auto left_right = Rods[rod_idx]->neighbouringCrossingIndex(crossing->on_rod_idx[rod_idx]);
            Rods[rod_idx]->u(left_right.first, u_left);
            Rods[rod_idx]->u(left_right.second, u_right);
            Rods[rod_idx]->u(node_idx, u);

            Offset offset_left, offset_right;
            Rods[rod_idx]->getEntry(left_right.first, offset_left);
            Rods[rod_idx]->getEntry(left_right.second, offset_right);

            if (node_idx != 0 || !Rods[rod_idx]->closed)
            {
                if (u_right - u < min_du && offset_right[3] != offset[3])
                {
                    entry_K.emplace_back(offset[3], offset[3], k_yc);
                    entry_K.emplace_back(offset_right[3], offset_right[3], k_yc);
                    
                }
                if (u - u_left < min_du && offset_left[3] != offset[3])
                {
                    entry_K.emplace_back(offset[3], offset[3], k_yc);
                    entry_K.emplace_back(offset_left[3], offset_left[3], k_yc);   
                }
            }

            Rods[rod_idx]->U(node_idx, U);
            T delta_u = (u - U);
            Range range = sliding_ranges[cnt];
            // 0 is the positive side sliding range
            if(delta_u >= range[0])
            {
                entry_K.emplace_back(offset[3], offset[3], k_yc);
            }
            else if (delta_u <= -range[1])
            {
                entry_K.emplace_back(offset[3], offset[3], k_yc);
            }
            cnt++;
        }
    }
}


void EoLRodSim::parallelContactdfdp(Eigen::Ref<VectorXT> dfdp)
{
    dfdp.setZero();
    if (!add_contact_penalty)
        return;

    int offset_cnt = 0;
    for (auto& crossing : rod_crossings)
    {
        if (crossing->is_fixed)
        {
            offset_cnt += 2 * crossing->rods_involved.size();    
            continue;
        }
        int node_idx = crossing->node_idx;
        std::vector<int> rods_involved = crossing->rods_involved;
        std::vector<Vector<T, 2>> sliding_ranges = crossing->sliding_ranges;

        int cnt = 0;
        for (int rod_idx : rods_involved)
        {
            Offset offset;
            Rods[rod_idx]->getEntry(node_idx, offset);
            T u, U;
            Rods[rod_idx]->u(node_idx, u);
            Rods[rod_idx]->U(node_idx, U);
            T delta_u = (u - U);
            Range range = sliding_ranges[cnt];
            // std::cout << "du " << delta_u << " range " << range.transpose() << std::endl;
            if(delta_u >= range[0])
            {
                dfdp[offset[3] + (offset_cnt) * deformed_states.rows()] += k_yc;
                // std::cout << "larger " << std::endl;
            }
            
            else if (delta_u <= -range[1])
            {
                dfdp[offset[3] + (offset_cnt + 1) * deformed_states.rows()] += -k_yc;
                // std::cout << "smaller " << std::endl;             
                
            }
            cnt++;
            offset_cnt += 2;
        }
    }   
    // std::cout << dfdp.norm() << std::endl;
}


void EoLRodSim::addParallelContactForce(Eigen::Ref<VectorXT> residual)
{
    if (!add_contact_penalty)
        return;
    VectorXT residual_cp = residual;
    for (auto& crossing : rod_crossings)
    {
        if (crossing->is_fixed)
            continue;
        int node_idx = crossing->node_idx;
        std::vector<int> rods_involved = crossing->rods_involved;
        std::vector<Vector<T, 2>> sliding_ranges = crossing->sliding_ranges;
        std::unordered_map<int, int> on_rod_idx = crossing->on_rod_idx;

        int cnt = 0;
        for (int rod_idx : rods_involved)
        {
            Offset offset;
            Rods[rod_idx]->getEntry(node_idx, offset);
            T u, U, u_left, u_right;
            
            auto left_right = Rods[rod_idx]->neighbouringCrossingIndex(crossing->on_rod_idx[rod_idx]);
            
            Rods[rod_idx]->u(left_right.first, u_left);
            Rods[rod_idx]->u(left_right.second, u_right);
            Rods[rod_idx]->u(node_idx, u);

            Offset offset_left, offset_right;
            Rods[rod_idx]->getEntry(left_right.first, offset_left);
            Rods[rod_idx]->getEntry(left_right.second, offset_right);

            if (node_idx != 0 || !Rods[rod_idx]->closed)
            {
                if (u_right - u < min_du && offset_right[3] != offset[3])
                {
                    residual[offset[3]] += k_yc * (u_right - u - min_du);
                    residual[offset_right[3]] += -k_yc * (u_right - u - min_du);
                }
                if (u - u_left < min_du && offset_left[3] != offset[3])
                {
                    residual[offset[3]] += -k_yc * (u - u_left - min_du);
                    residual[offset_left[3]] += k_yc * (u - u_left - min_du);
                }
            }

            
            Rods[rod_idx]->U(node_idx, U);
            T delta_u = (u - U);
            Range range = sliding_ranges[cnt];
            // if(rod_idx == 1)
            //     std::cout << range.transpose() << " " << delta_u << std::endl;
            // 0 is the positive side sliding range
            if(delta_u >= range[0])
            {
                // std::cout<< delta_u  << " " << range[0] << std::endl;
                residual[offset[3]] += -k_yc * (delta_u - range[0]);
            }
            // 1 is the sliding range along the negative direction
            else if (delta_u <= -range[1])
            {
                // std::cout<< delta_u  << " " << range[1] << std::endl;
                residual[offset[3]] += -k_yc * (delta_u + range[1]);
            }
            cnt++;
        }
    }
    if (print_force_mag)
        std::cout << "contact force norm: " << (residual - residual_cp).norm() << std::endl;
}


T EoLRodSim::addParallelContactEnergy()
{
    if (!add_contact_penalty)
        return 0.0;
    T energy = 0.0;
    for (auto& crossing : rod_crossings)
    {
        if (crossing->is_fixed)
            continue;
        int node_idx = crossing->node_idx;
        std::vector<int> rods_involved = crossing->rods_involved;
        std::vector<Vector<T, 2>> sliding_ranges = crossing->sliding_ranges;

        int cnt = 0;
        for (int rod_idx : rods_involved)
        {
            Offset offset;
            Rods[rod_idx]->getEntry(node_idx, offset);
            T u, U, u_left, u_right;
            
            auto left_right = Rods[rod_idx]->neighbouringCrossingIndex(crossing->on_rod_idx[rod_idx]);
            Rods[rod_idx]->u(left_right.first, u_left);
            Rods[rod_idx]->u(left_right.second, u_right);
            Rods[rod_idx]->u(node_idx, u);

            Offset offset_left, offset_right;
            Rods[rod_idx]->getEntry(left_right.first, offset_left);
            Rods[rod_idx]->getEntry(left_right.second, offset_right);

            if (node_idx != 0 || !Rods[rod_idx]->closed)
            {
                if (u_right - u < min_du && offset_right[3] != offset[3])
                {
                    energy += 0.5 * k_yc * std::pow(u_right - u - min_du, 2);
                }
                if (u - u_left < min_du && offset_left[3] != offset[3])
                {
                    energy += 0.5 * k_yc * std::pow(u - u_left - min_du, 2);
                }
            }

            Rods[rod_idx]->U(node_idx, U);
            T delta_u = (u - U);
            Range range = sliding_ranges[cnt];
            // 0 is the positive side sliding range
            // std::cout << delta_u << " " << range.transpose() << std::endl;
            if(delta_u >= range[0])
            {
                // std::cout << delta_u << " " << u << " " << U << std::endl;
                energy += 0.5 * k_yc * std::pow(delta_u - range[0], 2);
            }
            else if (delta_u <= -range[1])
            {
                // std::cout << delta_u << " " << u << " " << U << std::endl;
                energy += 0.5 * k_yc * std::pow(delta_u + range[1], 2);
            }
            
            

            cnt++;
        }
    }
    return energy;
}