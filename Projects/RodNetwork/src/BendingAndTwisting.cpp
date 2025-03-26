#include "../autodiff/RodBendingAndTwisting.h"
#include "../include/RodNetwork.h"

T RodNetwork::addBendingAndTwistingEnergy()
{
    T energy = 0.0;
    for (auto& rod : rods)
    {
        T energy_current = energy;
        rod->iterate3Nodes(
            [&](int node_i, int node_j, int node_k, int node_idx_local)
            {
                TV xi, xj, xk, Xi, Xj, Xk;
                rod->x(node_i, xi);
                rod->x(node_j, xj);
                rod->x(node_k, xk);
                rod->X(node_i, Xi);
                rod->X(node_j, Xj);
                rod->X(node_k, Xk);

                int left_rod = node_idx_local - 1;

                T theta0 = rod->reference_angles[left_rod];
                T theta1 = rod->reference_angles[node_idx_local];
                TV referenceNormal1 = rod->reference_frame_us[left_rod];
                TV referenceNormal2 = rod->reference_frame_us[node_idx_local];

                TV referenceTangent1 = rod->prev_tangents[left_rod];
                TV referenceTangent2 = rod->prev_tangents[node_idx_local];

                Matrix<T, 2, 2> B = rod->bending_coeffs;
                T kt = rod->kt;

                T reference_twist = rod->reference_twist[node_idx_local];
                T bending_plus_twisting = 0.0;
                computeRodBendingAndTwistEnergy(
                    B, kt, 0.0, referenceNormal1, referenceTangent1,
                    referenceNormal2, referenceTangent2, reference_twist, xk,
                    xi, xj, Xk, Xi, Xj, theta0, theta1, bending_plus_twisting);
                energy += bending_plus_twisting;
            });
    }
    return energy;
}

void RodNetwork::addBendingAndTwistingForceEntries(VectorXT& residual)
{

    for (auto& rod : rods)
    {
        rod->iterate3Nodes(
            [&](int node_i, int node_j, int node_k, int node_idx_local)
            {
                TV xi, xj, xk, Xi, Xj, Xk;
                rod->x(node_i, xi);
                rod->x(node_j, xj);
                rod->x(node_k, xk);
                rod->X(node_i, Xi);
                rod->X(node_j, Xj);
                rod->X(node_k, Xk);

                int left_rod = node_idx_local - 1;

                T theta0 = rod->reference_angles[left_rod];
                T theta1 = rod->reference_angles[node_idx_local];
                TV referenceNormal1 = rod->reference_frame_us[left_rod];
                TV referenceNormal2 = rod->reference_frame_us[node_idx_local];

                TV referenceTangent1 = rod->prev_tangents[left_rod];
                TV referenceTangent2 = rod->prev_tangents[node_idx_local];

                Matrix<T, 2, 2> B = rod->bending_coeffs;
                T reference_twist = rod->reference_twist[node_idx_local];

                Vector<T, 11> dedx;
                computeRodBendingAndTwistEnergyGradient(
                    B, rod->kt, 0.0, referenceNormal1, referenceTangent1,
                    referenceNormal2, referenceTangent2, reference_twist, xk,
                    xi, xj, Xk, Xi, Xj, theta0, theta1, dedx);
                addForceEntry<3>(residual, {node_k, node_i, node_j},
                                 -dedx.segment<9>(0), 0);
                int dof_theta0 =
                    rod->theta_dof_start_offset + node_idx_local - 1;
                int dof_theta1 = rod->theta_dof_start_offset + node_idx_local;
                residual[dof_theta0] -= dedx[9];
                residual[dof_theta1] -= dedx[10];
            });
    }
}

void RodNetwork::addBendingAndTwistingHessianEntries(
    std::vector<Entry>& entry_K)
{
    for (auto& rod : rods)
    {
        rod->iterate3Nodes(
            [&](int node_i, int node_j, int node_k, int node_idx_local)
            {
                TV xi, xj, xk, Xi, Xj, Xk, dXi, dXj, dXk;
                rod->x(node_i, xi);
                rod->x(node_j, xj);
                rod->x(node_k, xk);
                rod->X(node_i, Xi);
                rod->X(node_j, Xj);
                rod->X(node_k, Xk);

                std::vector<int> nodes = {node_k, node_i, node_j};

                int dof_theta0 =
                    rod->theta_dof_start_offset + node_idx_local - 1;
                int dof_theta1 = rod->theta_dof_start_offset + node_idx_local;

                std::vector<int> theta_dofs = {dof_theta0, dof_theta1};

                int left_rod = node_idx_local - 1;

                T theta0 = rod->reference_angles[left_rod];
                T theta1 = rod->reference_angles[node_idx_local];
                TV referenceNormal1 = rod->reference_frame_us[left_rod];
                TV referenceNormal2 = rod->reference_frame_us[node_idx_local];

                TV referenceTangent1 = rod->prev_tangents[left_rod];
                TV referenceTangent2 = rod->prev_tangents[node_idx_local];

                Matrix<T, 2, 2> B = rod->bending_coeffs;
                T reference_twist = rod->reference_twist[node_idx_local];
                Matrix<T, 11, 11> d2edx2;
                computeRodBendingAndTwistEnergyHessian(
                    B, rod->kt, 0.0, referenceNormal1, referenceTangent1,
                    referenceNormal2, referenceTangent2, reference_twist, xk,
                    xi, xj, Xk, Xi, Xj, theta0, theta1, d2edx2);
                // projectBlockPD<11>(d2edx2);
                addHessianEntry<3, 3>(entry_K, nodes, d2edx2.block<9, 9>(0, 0),
                                      0, 0);
                addHessianEntry<1, 1>(entry_K, theta_dofs,
                                      d2edx2.block<2, 2>(9, 9), 0, 0);
                addJacobianEntry<3, 1>(entry_K, nodes, theta_dofs,
                                       d2edx2.block<9, 2>(0, 9), 0);
                addJacobianEntry<1, 3>(entry_K, theta_dofs, nodes,
                                       d2edx2.block<2, 9>(9, 0), 0);
            });
    }
}