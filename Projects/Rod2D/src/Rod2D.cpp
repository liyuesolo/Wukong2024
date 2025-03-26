#include "../include/Rod2D.h"
#include <Eigen/CholmodSupport>

void Rod2D::initialize()
{
    int n_node = 10;           
    int angle = M_PI / 4;       
    T length = 1.0;                
    double segmentLength = length / (n_node - 1);

    undeformed.resize(n_node * 2);
    for (int i = 0; i < n_node; i++) 
    {
        double x = segmentLength * i * cos(angle);
        double y = segmentLength * i * sin(angle);
        undeformed[i * 2] = x;
        undeformed[i * 2 + 1] = y;
    }
    deformed = undeformed;
    xn = deformed;
    vn.resize(n_node * 2);
    vn.setZero();
    u.resize(n_node * 2);
    u.setZero();
    external_force.resize(n_node * 2);
    for (int i = 0; i < n_node; i++) 
    {
        external_force[i * 2] = 0;
        external_force[i * 2 + 1] = -9.8;
    }
    rod_indices.resize((n_node-1) * 2);
    for (int i = 0; i < n_node - 1; i++) 
    {
        rod_indices[i * 2] = i;
        rod_indices[i * 2 + 1] = i + 1; 
    }
}

T Rod2D::computeTotalEnergy()
{

    if (!run_diff_test)
        iterateDirichletDoF([&](int offset, T target) { u[offset] = target; });
    deformed = undeformed + u;
    T energy = 0.0;

    // if (dynamics)
    //    addInertialEnergy(energy);
    // energy -= u.dot(external_force);
    return energy;
}

T Rod2D::computeResidual(VectorXT& residual)
{
    if (!run_diff_test)
        iterateDirichletDoF([&](int offset, T target) { u[offset] = target; });
    deformed = undeformed + u;

    if (!run_diff_test)
        iterateDirichletDoF([&](int offset, T target)
                            { residual[offset] = 0; });
    return residual.norm();
}

bool Rod2D::linearSolve(StiffnessMatrix& K, const VectorXT& residual,
                        VectorXT& du)
{

    START_TIMING(LinearSolve)
    Eigen::CholmodSupernodalLLT<StiffnessMatrix, Eigen::Lower> solver;

    T alpha = 1e-6;
    if (!dynamics)
    {
        StiffnessMatrix H(K.rows(), K.cols());
        H.setIdentity();
        H.diagonal().array() = 1e-10;
        K += H;
    }
    solver.analyzePattern(K);

    int indefinite_count_reg_cnt = 0, invalid_search_dir_cnt = 0,
        invalid_residual_cnt = 0;
    int i = 0;
    T dot_dx_g = 0.0;
    for (; i < 50; i++)
    {
        solver.factorize(K);
        if (solver.info() == Eigen::NumericalIssue)
        {
            K.diagonal().array() += alpha;
            alpha *= 10;
            indefinite_count_reg_cnt++;
            continue;
        }

        du = solver.solve(residual);

        dot_dx_g = du.normalized().dot(residual.normalized());

        int num_negative_eigen_values = 0;
        int num_zero_eigen_value = 0;

        bool positive_definte = num_negative_eigen_values == 0;
        bool search_dir_correct_sign = dot_dx_g > 1e-6;
        if (!search_dir_correct_sign)
        {
            invalid_search_dir_cnt++;
        }

        bool solve_success = du.norm() < 1e3;

        if (!solve_success)
            invalid_residual_cnt++;

        if (positive_definte && search_dir_correct_sign && solve_success)
        {

            if (verbose)
            {
                std::cout << "	===== Linear Solve ===== " << std::endl;
                std::cout << "	nnz: " << K.nonZeros() << std::endl;

                std::cout << "	# regularization step " << i << " indefinite "
                          << indefinite_count_reg_cnt << " invalid search dir "
                          << invalid_search_dir_cnt << " invalid solve "
                          << invalid_residual_cnt << std::endl;
                std::cout << "	dot(search, -gradient) " << dot_dx_g
                          << std::endl;
                std::cout << "	======================== " << std::endl;
                FINISH_TIMING_PRINT(LinearSolve)
            }
            return true;
        }
        else
        {
            K.diagonal().array() += alpha;
            alpha *= 10;
        }
    }
    if (verbose)
    {
        std::cout << "	===== Linear Solve ===== " << std::endl;
        std::cout << "	nnz: " << K.nonZeros() << std::endl;
        // std::cout << "	 takes " << t.elapsed_sec() << "s" << std::endl;
        std::cout << "	# regularization step " << i << " indefinite "
                  << indefinite_count_reg_cnt << " invalid search dir "
                  << invalid_search_dir_cnt << " invalid solve "
                  << invalid_residual_cnt << std::endl;
        std::cout << "	dot(search, -gradient) " << dot_dx_g << std::endl;
        std::cout << "	======================== " << std::endl;
        FINISH_TIMING_PRINT(LinearSolve)
    }
    return false;
}

void Rod2D::buildSystemMatrix(StiffnessMatrix& K)
{
    int n_dof = deformed.rows();
    if (!run_diff_test)
        iterateDirichletDoF([&](int offset, T target) { u[offset] = target; });
    deformed = undeformed + u;
    std::vector<Entry> entries;

    K.resize(n_dof, n_dof);
    K.setFromTriplets(entries.begin(), entries.end());

    if (!run_diff_test)
        projectDirichletDoFMatrix(K, dirichlet_data);
    K.makeCompressed();
}

T Rod2D::lineSearchNewton(const VectorXT& residual)
{
    VectorXT du = residual;
    du.setZero();

    du = residual;
    StiffnessMatrix K(residual.rows(), residual.rows());
    buildSystemMatrix(K);
    bool success = linearSolve(K, residual, du);
    if (!success)
    {
        std::cout << "Linear Solve Failed" << std::endl;
        return 1e16;
    }

    T norm = du.norm();
    if (verbose)
        std::cout << "	|du | " << norm << std::endl;

    T E0 = computeTotalEnergy();
    T alpha = 1.0;
    int cnt = 0;
    VectorXT u_current = u;
    while (true)
    {
        u = u_current + alpha * du;
        T E1 = computeTotalEnergy();
        if (E1 - E0 < 0 || cnt > 10)
        {
            break;
        }
        alpha *= 0.5;
        cnt += 1;
    }
    return alpha * du.norm();
}

bool Rod2D::advanceOneStep(int step) 
{ 
    if (dynamics)
    {
        std::cout << "=================== Time STEP " << step * dt
                  << "s===================" << std::endl;
        bool finished = advanceOneTimeStep();
        updateDynamicStates();
        if (step * dt > simulation_duration)
            return true;
        return false;
    }
    else
    {
        std::cout << "===================STEP " << step
                  << "===================" << std::endl;
        VectorXT residual = external_force;
        T residual_norm = computeResidual(residual);
        
        std::cout << "[NEWTON] iter " << step << "/" << max_newton_iter
                  << ": residual_norm " << residual_norm
                  << " tol: " << newton_tol << std::endl;
        if (residual_norm < newton_tol || step == max_newton_iter)
        {
            return true;
        }

        T du_norm = 1e10;
        du_norm = lineSearchNewton(residual);
        return false;
    }
}

void Rod2D::updateDynamicStates()
{
    vn = (deformed - xn) / dt;
    xn = deformed;
}

bool Rod2D::advanceOneTimeStep()
{

    int iter = 0;
    while (true)
    {
        VectorXT residual = external_force;
        T residual_norm = computeResidual(residual);
        // residual_norms.push_back(residual_norm);
        if (verbose)
            std::cout << "[NEWTON] iter " << iter << "/" << max_newton_iter
                      << ": residual_norm " << residual_norm
                      << " tol: " << newton_tol << std::endl;
        if (residual_norm < newton_tol || iter == max_newton_iter)
        {
            std::cout << "[NEWTON] iter " << iter << "/" << max_newton_iter
                      << ": residual_norm " << residual_norm
                      << " tol: " << newton_tol << std::endl;
            break;
        }
        T du_norm = 1e10;
        du_norm = lineSearchNewton(residual);
        if (du_norm < 1e-10)
            break;
        iter++;
    }
    return true;
}
void Rod2D::projectDirichletDoFMatrix(StiffnessMatrix& A,
                                      const std::unordered_map<int, T>& data)
{
    for (auto iter : data)
    {
        A.row(iter.first) *= 0.0;
        A.col(iter.first) *= 0.0;
        A.coeffRef(iter.first, iter.first) = 1.0;
    }
}
