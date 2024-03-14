#include "../include/EoLRodSim.h"
#include "../autodiff/RotationPenalty.h"


void EoLRodSim::addPBCK(std::vector<Entry>& entry_K)
{
    
    iteratePBCStrainData([&](Offset offset_i, Offset offset_j, TV strain_dir, T Dij){

        TV xj = deformed_states.template segment<3>(offset_j[0]);
        TV xi = deformed_states.template segment<3>(offset_i[0]);

        T dij = (xj - xi).dot(strain_dir);
        TM Hessian = strain_dir * strain_dir.transpose();

        for(int i = 0; i < 3; i++)
        {
            for(int j = 0; j < 3; j++)
            {
                entry_K.push_back(Entry(offset_i[i], offset_i[j], k_strain * Hessian(i, j)));
                entry_K.push_back(Entry(offset_i[i], offset_j[j], -k_strain * Hessian(i, j)));
                entry_K.push_back(Entry(offset_j[i], offset_i[j], -k_strain * Hessian(i, j)));
                entry_K.push_back(Entry(offset_j[i], offset_j[j], k_strain * Hessian(i, j)));
            }
        }
    });
    

    if (add_rotation_penalty)
    {
        std::vector<TV2> data, dXdu, d2Xdu2;
        std::vector<Offset> offsets(4);
        buildRotationPenaltyData(data, offsets, dXdu, d2Xdu2);

        Matrix<T, 16, 16> J;
        computeRotationPenaltyEnergyHessian(kr, data[0], data[1], data[2], data[3],
                                            data[4], data[5], data[6], data[7], J);
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                for (int k = 0; k < 2; k++)
                    for (int l = 0; l < 2; l++)
                    {
                        entry_K.push_back(Entry(offsets[i][k], offsets[j][l], J(i*2 + k,j*2+l)));
                        entry_K.push_back(Eigen::Triplet<T>(offsets[k][i], offsets[l][3], J(k*2 + i, 4 * 2 + l * 2 + j) * dXdu[l][j]));

                        entry_K.push_back(Eigen::Triplet<T>(offsets[k][3], offsets[l][j], J(4 * 2 + k * 2 + i, l * 2 + j) * dXdu[k][i]));

                        
                        entry_K.push_back(Eigen::Triplet<T>(offsets[k][3], 
                                                            offsets[l][3], 
                                                            J(4 * 2 + k*2 + i, 4 * 2 + l * 2 + j) * dXdu[l][j] * dXdu[k][i]));
                    }
    }

    iteratePBCPairs([&](Offset offset_ref_i, Offset offset_ref_j, Offset offset_i, Offset offset_j)
    {
        
        TV xj = deformed_states.template segment<3>(offset_j[0]);
        TV xi = deformed_states.template segment<3>(offset_i[0]);

        TV xj_ref = deformed_states.template segment<3>(offset_ref_j[0]);
        TV xi_ref = deformed_states.template segment<3>(offset_ref_i[0]);

        std::vector<Offset> offsets = {offset_i, offset_j, offset_ref_i, offset_ref_j};
        std::vector<T> sign_J = {-1, 1, 1, -1};
        std::vector<T> sign_F = {1, -1, -1, 1};

        if ((offset_ref_j - offset_j).cwiseAbs().sum() < 1e-6 && (offset_ref_i - offset_i).cwiseAbs().sum() < 1e-6)
        {
            // entry_K.push_back(Entry(offset_i[3-1], offset_i[3-1], k_pbc));
            // entry_K.push_back(Entry(offset_j[3-1], offset_j[3-1], k_pbc));

            

            entry_K.push_back(Entry(offset_i[3-1], offset_i[3-1], k_pbc));
            entry_K.push_back(Entry(offset_j[3-1], offset_j[3-1], k_pbc));
            entry_K.push_back(Entry(offset_i[3-1], offset_j[3-1], -k_pbc));
            entry_K.push_back(Entry(offset_j[3-1], offset_i[3-1], -k_pbc));
            return;
        }

        for(int k = 0; k < 4; k++)
            for(int l = 0; l < 4; l++)
                for(int i = 0; i < 3; i++)
                    entry_K.push_back(Entry(offsets[k][i], offsets[l][i], -k_pbc *sign_F[k]*sign_J[l]));
    });
    
}



void EoLRodSim::addPBCForce(Eigen::Ref<VectorXT> residual)
{
    
    VectorXT residual_cp = residual;
    iteratePBCStrainData([&](Offset offset_i, Offset offset_j, TV strain_dir, T Dij){

        TV xj = deformed_states.template segment<3>(offset_j[0]);
        TV xi = deformed_states.template segment<3>(offset_i[0]);

        
        T dij = (xj - xi).dot(strain_dir);
        
        residual.template segment<3>(offset_i[0]) += k_strain * strain_dir * (dij - Dij);
        residual.template segment<3>(offset_j[0]) += -k_strain * strain_dir * (dij - Dij);
    });
    
    

    if (add_rotation_penalty)
    {
        std::vector<TV2> data, dXdu, d2Xdu2;
        std::vector<Offset> offsets(4);
        buildRotationPenaltyData(data, offsets, dXdu, d2Xdu2);
        Vector<T, 16> F;
        computeRotationPenaltyEnergyGradient(kr, data[0], data[1], data[2], data[3],
                                            data[4], data[5], data[6], data[7], F);
        for (int i = 0; i < 4; i++)
        {
            residual.template segment<2>(offsets[i][0]) += -F.template segment<2>(i*2);

            residual[offsets[i][3]] += -F.template segment<2>(4*2+ i*2).dot(dXdu[i]);
        }
    }

    // std::cout << (residual - residual_cp).norm() << std::endl;
    // std::getchar();
    
    iteratePBCPairs([&](Offset offset_ref_i, Offset offset_ref_j, Offset offset_i, Offset offset_j){
        
        TV xj = deformed_states.template segment<3>(offset_j[0]);
        TV xi = deformed_states.template segment<3>(offset_i[0]);

        TV xj_ref = deformed_states.template segment<3>(offset_ref_j[0]);
        TV xi_ref = deformed_states.template segment<3>(offset_ref_i[0]);


        if ((offset_ref_j - offset_j).cwiseAbs().sum() < 1e-6 && (offset_ref_i - offset_i).cwiseAbs().sum() < 1e-6)
        {
            // residual[offset_i[3 - 1]] += -k_pbc * xi[3-1];
            // residual[offset_j[3 - 1]] += -k_pbc * xj[3-1];
            T dx = xj[3-1] - xi[3-1];
            residual[offset_i[3 - 1]] += k_pbc * dx;
            residual[offset_j[3 - 1]] += -k_pbc * dx;
            return;
        }

        TV pair_dis_vec = xj - xi - (xj_ref - xi_ref);
        residual.template segment<3>(offset_i[0]) += k_pbc *pair_dis_vec;
        residual.template segment<3>(offset_j[0]) += -k_pbc *pair_dis_vec;
        residual.template segment<3>(offset_ref_i[0]) += -k_pbc *pair_dis_vec;
        residual.template segment<3>(offset_ref_j[0]) += k_pbc *pair_dis_vec;
    });

    // std::cout << (residual - residual_cp).norm() << std::endl;
    // std::getchar();

    if (print_force_mag)
        std::cout << "pbc force " << (residual - residual_cp).norm() << std::endl;
}


void EoLRodSim::buildRotationPenaltyData(
    std::vector<TV2>& data, std::vector<Offset>& offsets, 
    std::vector<TV2>& dXdu, std::vector<TV2>& d2Xdu2)
{
    auto ref0 = pbc_pairs_reference[0];
    auto ref1 = pbc_pairs_reference[1];

    data.resize(8); dXdu.resize(4); d2Xdu2.resize(4);

    data[0] = deformed_states.template segment<2>(ref0.first.first[0]);
    data[1] = deformed_states.template segment<2>(ref0.first.second[0]);
    data[2] = deformed_states.template segment<2>(ref1.first.first[0]);
    data[3] = deformed_states.template segment<2>(ref1.first.second[0]);

    offsets[0] = ref0.first.first;
    offsets[1] = ref0.first.second;
    offsets[2] = ref1.first.first;
    offsets[3] = ref1.first.second;
    
    TV pos, dpos, ddpos;
    T u = deformed_states[ref0.first.first[3]];
    Rods[ref0.second.first]->rest_state->getMaterialPos(u, pos, dpos, ddpos, true, true);
    data[4] = pos.template segment<2>(0);
    dXdu[0] = dpos.template segment<2>(0);
    d2Xdu2[0] = ddpos.template segment<2>(0);

    u = deformed_states[ref0.first.second[3]];
    Rods[ref0.second.second]->rest_state->getMaterialPos(u, pos, dpos, ddpos, true, true);
    data[5] = pos.template segment<2>(0);
    dXdu[1] = dpos.template segment<2>(0);
    d2Xdu2[1] = ddpos.template segment<2>(0);

    u = deformed_states[ref1.first.first[3]];
    Rods[ref1.second.first]->rest_state->getMaterialPos(u, pos, dpos, ddpos, true, true);
    data[6] = pos.template segment<2>(0);
    dXdu[2] = dpos.template segment<2>(0);
    d2Xdu2[2] = ddpos.template segment<2>(0);

    u = deformed_states[ref1.first.second[3]];
    Rods[ref1.second.second]->rest_state->getMaterialPos(u, pos, dpos, ddpos, true, true);
    data[7] = pos.template segment<2>(0);
    dXdu[3] = dpos.template segment<2>(0);
    d2Xdu2[3] = ddpos.template segment<2>(0);

}



T EoLRodSim::addPBCEnergy()
{
    T energy_pbc = 0.0;
    iteratePBCStrainData([&](Offset offset_i, Offset offset_j, TV strain_dir, T Dij){
        TV xj = deformed_states.template segment<3>(offset_j[0]);
        TV xi = deformed_states.template segment<3>(offset_i[0]);
        T dij = (xj - xi).dot(strain_dir);
        energy_pbc += 0.5 * k_strain * (dij - Dij) * (dij - Dij);
        
    });

    if (add_rotation_penalty)
    {
        std::vector<TV2> data, dXdu, d2Xdu2;
        std::vector<Offset> offsets(4);
        buildRotationPenaltyData(data, offsets, dXdu, d2Xdu2);
        T E_rotation = computeRotationPenaltyEnergy(kr, data[0], data[1], data[2], data[3],
                                                    data[4], data[5], data[6], data[7]); 
        energy_pbc += E_rotation;
    }
    

    iteratePBCPairs([&](Offset offset_ref_i, Offset offset_ref_j, Offset offset_i, Offset offset_j){
        
        
        TV xj = deformed_states.template segment<3>(offset_j[0]);
        TV xi = deformed_states.template segment<3>(offset_i[0]);

        
        TV xj_ref = deformed_states.template segment<3>(offset_ref_j[0]);
        TV xi_ref = deformed_states.template segment<3>(offset_ref_i[0]);

        
        if ((offset_ref_j - offset_j).cwiseAbs().sum() < 1e-6 && (offset_ref_i - offset_i).cwiseAbs().sum() < 1e-6)
        {            
            // energy_pbc += 0.5 * k_pbc * std::pow(xi[3-1], 2);
            // energy_pbc += 0.5 * k_pbc * std::pow(xj[3-1], 2);
            energy_pbc += 0.5 * k_pbc * std::pow(xj[3-1] - xi[3-1], 2);
            return;
        }
        
        TV pair_dis_vec = xj - xi - (xj_ref - xi_ref);
        energy_pbc += 0.5  *k_pbc * pair_dis_vec.dot(pair_dis_vec);
    });

    return energy_pbc;
}
