#include <Eigen/SparseCore>
#include <fstream>
#include <chrono>
#include <iomanip>

#include <Spectra/SymEigsShiftSolver.h>
#include <Spectra/MatOp/SparseSymShiftSolve.h>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/SparseSymMatProd.h>

#include "../include/EoLRodSim.h"


void EoLRodSim::resetScene()
{ 
    deformed_states = rest_states; 
    for (auto& rod : Rods)
    {
        rod->reference_twist.setZero();
        rod->reference_angles.setZero();
    }
    for (auto& crossing : rod_crossings)
    {
        crossing->omega.setZero();
        crossing->rotation_accumulated.setIdentity();
    }

    for (auto& rod : Rods)
    {
        rod->setupBishopFrame();
    }
}


void EoLRodSim::computeBoundingBox(TV& bottom_left, TV& top_right) const
{
    bottom_left.setConstant(1e6);
    top_right.setConstant(-1e6);

    for (auto& rod : Rods)
    {
        int cnt = 0;
        for (int idx : rod->indices)
        {
            TV x;
            rod->x(idx, x);
            for (int d = 0; d < 3; d++)
            {
                top_right[d] = std::max(top_right[d], x[d]);
                bottom_left[d] = std::min(bottom_left[d], x[d]);
            }
        }
    }   
}

void EoLRodSim::fixRegionAvoidRod(std::function<bool(const TV&)> inside_region, int rod_idx)
{
    for (auto& rod : Rods)
    {
        if (rod->rod_id == rod_idx)
            continue;
        int cnt = 0;
        for (int idx : rod->indices)
        {
            TV x;
            rod->x(idx, x);
            Offset offset;
            
            if (inside_region(x))
            {
                rod->getEntry(idx, offset);
                for (int d = 0; d < 3; d++)
                {
                    dirichlet_dof[rod->reduced_map[offset[d]]] = 0.0;
                }
                if (cnt < rod->numSeg())
                    dirichlet_dof[rod->theta_reduced_dof_start_offset + cnt] = 0;
            }
            cnt++;
        }
    }
    for (auto& crossing : rod_crossings)
    {
        auto rod = Rods[crossing->rods_involved.front()];
        TV x;
        rod->x(crossing->node_idx, x);
        if (inside_region(x))
        {
            for (int d = 0; d < 3; d++)
                dirichlet_dof[crossing->reduced_dof_offset + d] = 0.0;
        }
    }
}

void EoLRodSim::fixRegion(std::function<bool(const TV&)> inside_region)
{
    for (auto& rod : Rods)
    {
        int cnt = 0;
        for (int idx : rod->indices)
        {
            TV x;
            rod->x(idx, x);
            Offset offset;
            
            if (inside_region(x))
            {
                rod->getEntry(idx, offset);
                for (int d = 0; d < 3; d++)
                {
                    dirichlet_dof[rod->reduced_map[offset[d]]] = 0.0;
                }
                if (cnt < rod->numSeg())
                    dirichlet_dof[rod->theta_reduced_dof_start_offset + cnt] = 0;
            }
            cnt++;
        }
    }
    for (auto& crossing : rod_crossings)
    {
        auto rod = Rods[crossing->rods_involved.front()];
        TV x;
        rod->x(crossing->node_idx, x);
        if (inside_region(x))
        {
            for (int d = 0; d < 3; d++)
                dirichlet_dof[crossing->reduced_dof_offset + d] = 0.0;
        }
    }
}


void EoLRodSim::fixRegionalDisplacement(
    std::function<bool(const TV&, TV&,Vector<bool, 3>&)> helper)
{
    for (auto& rod : Rods)
    {
        int cnt = 0;
        for (int idx : rod->indices)
        {
            TV x;
            rod->x(idx, x);
            Offset offset;
            TV dx;
            Vector<bool, 3> mask;
            if (helper(x, dx, mask))
            {
                rod->getEntry(idx, offset);
                for (int d = 0; d < 3; d++)
                    if (mask[d])
                        dirichlet_dof[rod->reduced_map[offset[d]]] = dx[d];
                
                if (cnt < rod->numSeg())
                    dirichlet_dof[rod->theta_reduced_dof_start_offset + cnt] = 0;
            }
            cnt++;
        }
    }
    for (auto& crossing : rod_crossings)
    {
        auto rod = Rods[crossing->rods_involved.front()];
        TV x;
        rod->x(crossing->node_idx, x);
        TV dx;
        Vector<bool, 3> mask;
        if (helper(x, dx, mask))
        {
            for (int d = 0; d < 3; d++)
                
                dirichlet_dof[crossing->reduced_dof_offset + d] = 0.0;
        }
    }
}


void EoLRodSim::fixCrossingLagrangian(int crossing_idx, 
    TV delta, Mask mask)
{
    auto crossing = rod_crossings[crossing_idx];
    auto rod = Rods[crossing->rods_involved.front()];
    Offset offset;
    rod->getEntry(crossing->node_idx, offset);
    for (int d = 0; d < 3; d++)
        if (mask[d])
            dirichlet_dof[rod->reduced_map[offset[d]]] = delta[d];
}


void EoLRodSim::fixCrossing()
{
    for(auto& crossing : rod_crossings)
    {
        int node_idx = crossing->node_idx;
        std::vector<int> rods_involved = crossing->rods_involved;

        int cnt = 0;    
        for (int rod_idx : rods_involved)
        {
            Offset offset;
            Rods[rod_idx]->getEntry(node_idx, offset);
            if (crossing->is_fixed)
                dirichlet_dof[Rods[rod_idx]->reduced_map[offset[3]]] = 0;
            else
            {
                for (int d = 0; d < 3; d++) 
                    dirichlet_dof[crossing->reduced_dof_offset + d] = 0;
                if ((crossing->sliding_ranges[cnt] - TV2::Zero()).norm() < 1e-6)
                    dirichlet_dof[Rods[rod_idx]->reduced_map[offset[3]]] = 0;
            }
            cnt++;
        }
        
    }
}


T EoLRodSim::computeTotalEnergy(Eigen::Ref<const VectorXT> dq, 
        bool verbose)
{
    VectorXT dq_projected = dq;
    
    if(!add_penalty && !run_diff_test)
        iterateDirichletDoF([&](int offset, T target)
        {
            dq_projected[offset] = target;
        });

   
    deformed_states = rest_states + W * dq_projected;

    for (auto& rod : Rods)
    {
        rod->reference_angles = deformed_states.template segment(rod->theta_dof_start_offset, 
            rod->indices.size() - 1);
    }
    for (auto& crossing : rod_crossings)
        crossing->omega = deformed_states.template segment<3>(crossing->dof_offset);

    T total_energy = 0;
    T E_stretching = 0, E_bending = 0, E_shearing = 0, E_twisting = 0, E_bending_twisting = 0,
        E_eul_reg = 0, E_pbc = 0, E_penalty = 0, E_contact = 0;
    
    if (add_stretching)
        E_stretching += addStretchingEnergy();
    if (add_bending && add_twisting)
    {
        E_bending_twisting = add3DBendingAndTwistingEnergy();
    }
    if (add_pbc_bending || add_pbc_twisting)
        E_bending_twisting += add3DPBCBendingAndTwistingEnergy(add_pbc_bending, add_pbc_twisting);
    if (add_rigid_joint)
        E_bending_twisting += addJointBendingAndTwistingEnergy();
    
    if (add_pbc)
        E_pbc += addPBCEnergy();
    if (add_eularian_reg)
        E_eul_reg += addRegularizingEnergy();
    if (add_contact_penalty)
        E_contact += addParallelContactEnergy();
    total_energy = E_stretching + E_bending + E_shearing + E_twisting + E_bending_twisting + E_eul_reg + E_pbc + E_penalty + E_contact;
    
    if (verbose)
        std::cout << "E_stretching " << E_stretching << " E_bending " 
        << E_bending << " E_twisting " << E_twisting << 
        " E_bending_twisting " << E_bending_twisting <<  
        " E_shearing " << E_shearing << " E_eul_reg " << E_eul_reg << 
        " E_pbc " << E_pbc << " E_penalty " << E_penalty << " E_contact " << E_contact << std::endl;
    return total_energy;
}



T EoLRodSim::computeResidual(Eigen::Ref<VectorXT> residual, Eigen::Ref<const VectorXT> dq)
{
    VectorXT dq_projected = dq;
    if(!add_penalty && !run_diff_test)
        iterateDirichletDoF([&](int offset, T target)
        {
            dq_projected[offset] = target;
        });

    deformed_states = rest_states + W * dq_projected;

    for (auto& rod : Rods)
    {
        rod->reference_angles = deformed_states.template segment(rod->theta_dof_start_offset, 
            rod->indices.size() - 1);
        if (!run_diff_test)
            rod->rotateReferenceFrameToLastNewtonStepAndComputeReferenceTwsit();
    }

    for (auto& crossing : rod_crossings)
    {
        crossing->omega = deformed_states.template segment<3>(crossing->dof_offset);
        // crossing->updateRotation(deformed_states.template segment<3>(crossing->dof_offset));
    }
        // crossing->updateRotation(deformed_states.template segment<3>(crossing->dof_offset));

    VectorXT full_residual(deformed_states.rows());
    full_residual.setZero();

    if (add_stretching)
        addStretchingForce(full_residual);
    if (add_bending && add_twisting)
    {
        add3DBendingAndTwistingForce(full_residual);
    }
    if (add_pbc_bending || add_pbc_twisting)
        add3DPBCBendingAndTwistingForce(full_residual, add_pbc_bending, add_pbc_twisting);
    if (add_rigid_joint)
        addJointBendingAndTwistingForce(full_residual);
    
    if (add_pbc)
        addPBCForce(full_residual);
    
    if (add_eularian_reg)
        addRegularizingForce(full_residual);
    if (add_contact_penalty)
        addParallelContactForce(full_residual);
    
    // std::cout << " full residual " << std::endl;
    // std::cout << full_residual << std::endl;
    // std::getchar();
    
    residual = W.transpose() * full_residual;
    
    if (!run_diff_test)
        iterateDirichletDoF([&](int offset, T target)
        {
            residual[offset] = 0;
        });
    
    return residual.norm();
}


void EoLRodSim::addStiffnessMatrix(std::vector<Eigen::Triplet<T>>& entry_K,
         Eigen::Ref<const VectorXT> dq)
{
    VectorXT dq_projected = dq;
    if(!add_penalty && !run_diff_test)
        iterateDirichletDoF([&](int offset, T target)
        {
            dq_projected[offset] = target;
        });

    deformed_states = rest_states + W * dq_projected;
    for (auto& rod : Rods)
    {
        rod->reference_angles = deformed_states.template segment(rod->theta_dof_start_offset, rod->indices.size() - 1);
    }
    for (auto& crossing : rod_crossings)
        crossing->omega = deformed_states.template segment<3>(crossing->dof_offset);

    if (add_stretching)
        addStretchingK(entry_K);
    if (add_bending && add_twisting)
    {
        add3DBendingAndTwistingK(entry_K);
    }
    if (add_pbc_bending || add_pbc_twisting)
        add3DPBCBendingAndTwistingK(entry_K, add_pbc_bending, add_pbc_twisting);
    if (add_rigid_joint)
        addJointBendingAndTwistingK(entry_K);
    
    if (add_pbc)
        addPBCK(entry_K);
    if (add_eularian_reg)
        addRegularizingK(entry_K);
    if (add_contact_penalty)
        addParallelContactK(entry_K);
}


void EoLRodSim::buildSystemDoFMatrix(
    Eigen::Ref<const VectorXT> dq, StiffnessMatrix& K)
{
    std::vector<Entry> entry_K;
    addStiffnessMatrix(entry_K, dq);
    StiffnessMatrix A(deformed_states.rows(), deformed_states.rows());

    A.setFromTriplets(entry_K.begin(), entry_K.end());
    // std::cout << A.rows() << " " << A.cols() << std::endl;
    // std::cout << A << std::endl;
    K = W.transpose() * A * W;
    
    if(!add_penalty && !run_diff_test)
        projectDirichletDoFSystemMatrix(K);

    // std::cout << K << std::endl;
    
    K.makeCompressed();
}


void EoLRodSim::projectDirichletDoFSystemMatrix(StiffnessMatrix& A)
{
    iterateDirichletDoF([&](int offset, T target)
    {
        A.row(offset) *= 0.0;
        A.col(offset) *= 0.0;
        A.coeffRef(offset, offset) = 1.0;
    });
}


void EoLRodSim::projectDirichletDoFMatrix(StiffnessMatrix& A, const std::unordered_map<int, T>& data)
{
    for (auto iter : data)
    {
        A.row(iter.first) *= 0.0;
        A.col(iter.first) *= 0.0;
        A.coeffRef(iter.first, iter.first) = 1.0;
    }

}

bool EoLRodSim::linearSolve(StiffnessMatrix& K, 
    Eigen::Ref<const VectorXT> residual, Eigen::Ref<VectorXT> ddq)
{
    
    StiffnessMatrix I(K.rows(), K.cols());
    I.setIdentity();

    StiffnessMatrix H = K;
    // Eigen::SimplicialLLT<StiffnessMatrix> solver;
    Eigen::SimplicialLDLT<StiffnessMatrix> solver;
    
    T mu = 10e-6;
    solver.analyzePattern(K);
    for (int i = 0; i < 50; i++)
    {
        solver.factorize(K);
        if (solver.info() == Eigen::NumericalIssue)
        {
            K = H + mu * I;        
            mu *= 10;
            continue;
        }
        ddq = solver.solve(residual);

        T dot_dx_g = ddq.normalized().dot(residual.normalized());

        VectorXT d_vector = solver.vectorD();
        int num_negative_eigen_values = 0;

        for (int i = 0; i < d_vector.size(); i++)
        {
            if (d_vector[i] < 0)
            {
                num_negative_eigen_values++;
                break;
            }
        
        }
        bool positive_definte = num_negative_eigen_values == 0;
        bool search_dir_correct_sign = dot_dx_g > 1e-6;
        bool solve_success = (K*ddq - residual).norm() < 1e-6 && solver.info() == Eigen::Success;

        if (positive_definte && search_dir_correct_sign && solve_success)
            break;
        else
        {
            K = H + mu * I;        
            mu *= 10;
        }
    }

    return true;
}




T EoLRodSim::lineSearchNewton(Eigen::Ref<VectorXT> dq, 
        Eigen::Ref<const VectorXT> residual, 
        int line_search_max)
{
    bool debug = true;
    bool verbose = false;
    
    VectorXT ddq(W.cols());
    ddq.setZero();

    StiffnessMatrix K;
    buildSystemDoFMatrix(dq, K);
    bool success = linearSolve(K, residual, ddq);
    if (!success)
        return 1e16;

    // T norm = ddq.cwiseAbs().maxCoeff();
    T norm = ddq.norm();
    // std::cout << norm << std::endl;
    // if (norm < 1e-6) return norm;
    T alpha = 1;
    T E0 = computeTotalEnergy(dq);
    // std::cout << "E0: " << E0 << std::endl;
    int cnt = 0;
    bool set_to_gradient = true;
    std::vector<T> Es;
    while(true)
    {
        VectorXT dq_ls = dq + alpha * ddq;
        T E1 = computeTotalEnergy(dq_ls, verbose);
        if (debug)
            Es.push_back(E1);
        // std::cout << "E1: " << E1 << std::endl;
        if (E1 - E0 < 0) {
            dq = dq_ls;
            // return 1e16;
            // testGradient(dq);
            // testHessian(dq);
            break;
        }
        alpha *= T(0.5);
        cnt += 1;
        if (cnt > 15)
        {
            
            // checkHessianPD(dq);
            // for(T a = 1.0; a > 1e-8; a*=0.5)
            // {
            //     T E_grad = computeTotalEnergy(dq + a * residual, verbose);
            //     if (E_grad < E0)
            //         std::cout << "take gradient" << std::endl;
            //         // dq_ls = dq + a * residual;
            //     // std::cout << "E_grad " << std::setprecision(15) << E_grad << std::endl;
            // }
            
            // std::cout << "E0 " << std::setprecision(15) << E0 << std::endl;
            // for (T E : Es)
            //     std::cout << std::setprecision(15) << E << " ";
            // std::cout << std::endl;
            // std::cout << residual.norm() << std::endl;
            // std::getchar();
            // testGradient(dq);
            // testHessian(dq);
            // return 1e16;
            // computeSmallestEigenVector(K, dq_ls);
            // std::cout << ddq.norm() << std::endl;
            // std::cout << ddq.normalized().dot(residual.normalized()) << std::endl;
            // std::cout << residual.norm() << std::endl;
            // std::getchar();
            // for (auto& crossing : rod_crossings)
            // // for (int i = 0; i < 10; i++)
            // {
            //     // auto crossing = rod_crossings[i];
            //     Offset off;
            //     Rods[crossing->rods_involved.front()]->getEntry(crossing->node_idx, off);
            //     T r = static_cast <T> (rand()) / static_cast <T> (RAND_MAX);
            //     int z_off = Rods[crossing->rods_involved.front()]->reduced_map[off[3-1]];
            //     dq_ls[z_off] += 0.001 * (r - 0.5) * unit;
            //     // break;
            //     // dq[z_off] += 0.001 * r * unit;
                
            // }
            dq = dq_ls;
            break;
        }
        if (cnt == line_search_max) 
            break;
            // return 1e16;
            
    }
    
    for (auto& crossing : rod_crossings)
    {
        Vector<T, 3> new_omega = dq.template segment<3>(crossing->reduced_dof_offset);
        crossing->updateRotation(new_omega);
        dq.template segment<3>(crossing->reduced_dof_offset).setZero();
    }
    
    // std::cout << "#ls: " << cnt << std::endl;
    return norm;   
}


bool EoLRodSim::staticSolve(Eigen::Ref<VectorXT> dq)
{
    int cnt = 0;
    T residual_norm = 1e10, dq_norm = 1e10;
    
    int max_newton_iter = 1000;
    // std::chrono::high_resolution_clock::time_point time = std::chrono::high_resolution_clock::now();
    while (true)
    {
        VectorXT residual(W.cols());
        residual.setZero();
        
        if (cnt == 0)
        {
            dq += perturb;
        }

        
        computeResidual(residual, dq);

        // auto f_time = std::chrono::duration_cast<std::chrono::milliseconds>( std::chrono::high_resolution_clock::now() - time ).count();

        residual_norm = residual.norm();
        if (verbose)
            std::cout << "residual_norm " << residual_norm << std::endl;
        
        // std::cout << residual_norm << std::endl;
        // std::getchar();

        if (residual_norm < newton_tol)
            break;
        
        dq_norm = lineSearchNewton(dq, residual, 50);

        // auto newton_time = std::chrono::duration_cast<std::chrono::milliseconds>( std::chrono::high_resolution_clock::now() - time ).count();
        // std::cout << f_time / 1000.0 << " " << newton_time / 1000.0 << std::endl;
        // std::getchar();
        if(cnt == max_newton_iter || dq_norm > 1e10)
            break;
        cnt++;
    }
    
    // if (verbose)
        std::cout << "# of newton solve: " << cnt << " exited with |g|: " << residual_norm << "|dq|: " << dq_norm  << std::endl;
    
    if (cnt == max_newton_iter || dq_norm > 1e10 || residual_norm > 1)
        return false;
    return true;
}


void EoLRodSim::computeSmallestEigenVector(const StiffnessMatrix& K, 
    Eigen::Ref<VectorXT> eigen_vector)
{
    int nmodes = 2;
    Spectra::SparseSymShiftSolve<T, Eigen::Upper> op(K);

    double shift = -1e-1;
    Spectra::SymEigsShiftSolver<T, Spectra::LARGEST_MAGN, Spectra::SparseSymShiftSolve<T, Eigen::Upper> > eigs(&op, nmodes, 20, shift);
    eigs.init();

    int nconv = eigs.compute();

    if (eigs.info() == Spectra::SUCCESSFUL)
    {
        Eigen::MatrixXd eigen_vectors = eigs.eigenvectors().real();
        // Eigen::VectorXd eigen_values = eigs.eigenvalues().real();
        eigen_vector = eigen_vectors.row(nmodes - 1);

    }
    else
    {
        eigen_vector = VectorXT::Zero(K.cols());
    }

    // Eigen::MatrixXd A_dense = K;

    // Eigen::EigenSolver<Eigen::MatrixXd> eigen_solver;
    
    // eigen_solver.compute(A_dense, /* computeEigenvectors = */ false);
    // auto eigen_values = eigen_solver.eigenvalues().real();
    // // auto eigen_vectors = eigen_solver.eigenvectors().real();
    // std::cout << eigen_values << std::endl;
}


void EoLRodSim::checkHessianPD(Eigen::Ref<const VectorXT> dq)
{
    bool compute_position_sub_hessian_evs = false;

    int nmodes = 10;

    StiffnessMatrix K;
    buildSystemDoFMatrix(dq, K);
    // K.conservativeResize(Rods[0]->theta_reduced_dof_start_offset, Rods[0]->theta_reduced_dof_start_offset);
    // K.conservativeResize(rod_crossings[0]->reduced_dof_offset, rod_crossings[0]->reduced_dof_offset);
    
    // std::cout << K << std::endl;
    // Eigen::SimplicialLLT<StiffnessMatrix> solver;
    // solver.compute(K);

    Spectra::SparseSymShiftSolve<T, Eigen::Upper> op(K);

    double shift = -1e-1;
    Spectra::SymEigsShiftSolver<T, Spectra::LARGEST_MAGN, Spectra::SparseSymShiftSolve<T, Eigen::Upper> > eigs(&op, nmodes, 2 * nmodes, shift);
    eigs.init();

    int nconv = eigs.compute();

    if (eigs.info() == Spectra::SUCCESSFUL)
    {
        Eigen::MatrixXd eigen_vectors = eigs.eigenvectors().real();
        Eigen::VectorXd eigen_values = eigs.eigenvalues().real();
        std::cout << eigen_values << std::endl;
        std::ofstream out("eigen_vectors.txt");
        out << eigen_vectors.rows() << " " << eigen_vectors.cols() << std::endl;
        for (int i = 0; i < eigen_vectors.rows(); i++)
        {
            for (int j = 0; j < eigen_vectors.cols(); j++)
                out << eigen_vectors(i, j) << " ";
            out << std::endl;
        }       
        out << std::endl;
        out.close();
    }
    else
    {
        std::cout << "Eigen decomposition failed" << std::endl;
    }

    if (compute_position_sub_hessian_evs)
    {

        Eigen::MatrixXd A_dense = K;
        int x_dof = Rods[0]->theta_reduced_dof_start_offset;
        // int x_dof = rod_crossings[0]->reduced_dof_offset;
        int total_dof = K.cols();
        Eigen::MatrixXd Kmm = A_dense.template block(0, 0, x_dof, x_dof);
        Eigen::MatrixXd Kms = A_dense.template block(0, x_dof, x_dof, total_dof - x_dof);
        Eigen::MatrixXd Kss = A_dense.template block(x_dof, x_dof, total_dof - x_dof, total_dof - x_dof);

        A_dense = Kmm - Kms * Kss.inverse() * Kms.transpose();

        
        
        // A_dense = Kss - Kms.transpose() * Kmm.inverse() * Kms;



        Eigen::EigenSolver<Eigen::MatrixXd> eigen_solver;
        // // // Eigen::VectorXcd eigen_vector;
        // // std::cout << A_dense << std::endl;
        // eigen_solver.compute(A_dense, /* computeEigenvectors = */ true);
        eigen_solver.compute(A_dense, /* computeEigenvectors = */ false);
        auto eigen_values = eigen_solver.eigenvalues();
        // auto eigen_vectors = eigen_solver.eigenvectors();

        Eigen::VectorXcd smallest_eigen_vector, second_smallest;
        T min_ev = 1e10;
        T second_min_ev = 1e10;
        for (int i = 0; i < A_dense.cols(); i++)
            if (eigen_values[i].real() < min_ev)
            {
                min_ev = eigen_values[i].real();
                // smallest_eigen_vector = eigen_vectors.col(i);        
            }
        
        // for (int i = 0; i < A_dense.cols(); i++)
        //     if (eigen_values[i].real() < second_min_ev && eigen_values[i].real() > min_ev)
        //     {
        //         second_min_ev = eigen_values[i].real();
        //         second_smallest = eigen_vectors.col(i);        
        //     }
        
        // Eigen::VectorXd ds = smallest_eigen_vector.real();
        // Eigen::VectorXd dm = -Kmm.colPivHouseholderQr().solve(Kms * ds);
        std::cout << min_ev << std::endl;
        
        std::vector<T> ev_all(A_dense.cols());
        for (int i = 0; i < A_dense.cols(); i++)
        {
            ev_all[i] = eigen_values[i].real();
        }
        std::sort(ev_all.begin(), ev_all.end());
        // for (int i = 0; i < 10; i++)
        for (int i = 0; i < ev_all.size(); i++)
            std::cout << ev_all[i] << " ";
        std::cout << std::endl;
        // std::cout << ev_all << std::endl;

        
        // std::cout << K.cols() << " " << dm.rows() << " " << ds.rows() << std::endl;
        // std::ofstream out("eigen_vector.txt");
        // // for (int i = 0; i < Kss.cols(); i++)
        // // {
        // //     out << double(smallest_eigen_vector[i].real()) << std::endl;
        // // }
        // for (int i = 0; i < dm.rows(); i++)
        // {
        //     out << dm[i] << std::endl;
        // }
        // for (int i = 0; i < ds.rows(); i++)
        // {
        //     out << ds[i] << std::endl;
        // }
        // out.close();
    }

}


bool EoLRodSim::forward(Eigen::Ref<VectorXT> dq)
{
    
    bool succeed = staticSolve(dq);
    if (!succeed)
        return false;
    
    VectorXT dq_projected = dq;

    iterateDirichletDoF([&](int offset, T target){
            dq_projected[offset] = target;
        });
    dq = dq_projected;
    deformed_states = rest_states + W * dq_projected;
    return true;
}


void EoLRodSim::staticSolveIncremental(int step)
{
    if (incremental_steps == 0)
    {
        std::cout << "[EoLRodSim.cpp] Zero incremetal step size" << std::endl;
        return;
    }
    resetScene();
    incremental_bc(*this, step);
    advanceOneStep(0);
    rest_states = deformed_states;
}


bool EoLRodSim::advanceOneStep(int step)
{
    // newton_tol = 1e-6;
    n_dof = W.cols();
    VectorXT dq(n_dof);
    dq.setZero();
    
    // checkHessianPD(dq);
    // return;

    staticSolve(dq);
    
    VectorXT dq_projected = dq;

    iterateDirichletDoF([&](int offset, T target){
            dq_projected[offset] = target;
        });
    deformed_states = rest_states + W * dq_projected;
    
    
    T e_total = 0.0;
    if (add_stretching)
        e_total += addStretchingEnergy();
    e_total += add3DBendingAndTwistingEnergy();
    if (add_pbc)
        e_total += add3DPBCBendingAndTwistingEnergy();
    e_total += addJointBendingAndTwistingEnergy();
    

    std::cout << "E total: " << e_total << std::endl;
    std::cout << "total stretching " << addStretchingEnergy() << std::endl;
    std::cout << "total bending " << add3DBendingAndTwistingEnergy(true, false) + 
                                    add3DPBCBendingAndTwistingEnergy(true, false) + 
                                    addJointBendingAndTwistingEnergy(true, false) << std::endl;

    std::cout << "total twist " << add3DBendingAndTwistingEnergy(false, true) + 
                                    add3DPBCBendingAndTwistingEnergy(false, true) + 
                                    addJointBendingAndTwistingEnergy(false, true) << std::endl;
    
    
    // checkHessianPD(dq);

    // VectorXT force(W.rows()); force.setZero();
    // addPBCForce(force);


    // TV fi = force.template segment<3>(pbc_pairs_reference[0].first.first[0]);
    // TV fj = force.template segment<3>(pbc_pairs_reference[1].first.first[0]);

    // // std::cout << fi.transpose() << " " << fj.transpose() << " " << std::endl;

    // for (auto& rod : Rods)
    // {
    //     rod->reference_angles = deformed_states.template segment(rod->theta_dof_start_offset, 
    //         rod->indices.size() - 1);
    //     VectorXT reference_twist = rod->reference_twist + rod->reference_angles;
    //     // std::cout << rod->reference_angles.transpose() << std::endl;
    //     // std::cout << reference_twist.transpose() << std::endl;
    // }
    return true;
}
