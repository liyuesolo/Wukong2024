#include <Eigen/CholmodSupport>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/massmatrix.h>
#include <igl/adjacency_list.h>
#include <Spectra/SymEigsShiftSolver.h>
#include <Spectra/MatOp/SparseSymShiftSolve.h>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <boost/math/policies/error_handling.hpp>

#include "../autodiff/CST3DShell.h"
#include "../include/IMLSContact.h"

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/heat_method_distance.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/surface_mesh.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;
using namespace boost::math::tools;

bool IMLSContact::advanceOneTimeStep()
{
    int iter = 0;
    while (true)
    {
        VectorXT residual = external_force;
        T residual_norm = computeResidual(residual);
        // residual_norms.push_back(residual_norm);
        if (verbose)
            std::cout << "[NEWTON] iter " << iter << "/" 
                << max_newton_iter << ": residual_norm " 
                << residual_norm << " tol: " << newton_tol << std::endl;
        if (residual_norm < newton_tol || iter == max_newton_iter)
        {
            std::cout << "[NEWTON] iter " << iter << "/" 
                << max_newton_iter << ": residual_norm " 
                << residual_norm << " tol: " << newton_tol << std::endl;
            break;
        }
        T du_norm = 1e10;
        du_norm = lineSearchNewton(residual);
        if (du_norm < 1e-10)
            break;
        iter ++;
        
    }

    return true;
}

bool IMLSContact::advanceOneStep(int step)
{
    if (dynamics)
    {
        Eigen::MatrixXd V;
        Eigen::MatrixXi F;
        std::cout << "=================== Time STEP " << step * dt << "s===================" << std::endl;
        bool finished = advanceOneTimeStep();
        vectorToIGLMatrix<T, 3>(deformed, V);
        vectorToIGLMatrix<int,3>(faces, F);
        igl::writeOBJ("Frame_"+std::to_string(step)+".obj",V,F);
        updateDynamicStates();
        if (step * dt > simulation_duration)
            return true;
        return false;
    }
    else
    {
        std::cout << "===================STEP " << step << "===================" << std::endl;
        VectorXT residual = external_force;
        T residual_norm = computeResidual(residual);
        residual_norms.push_back(residual_norm);
        std::cout << "[NEWTON] iter " << step << "/" 
            << max_newton_iter << ": residual_norm " 
            << residual_norm << " tol: " << newton_tol << std::endl;
        if (residual_norm < newton_tol || step == max_newton_iter)
        {
            return true;
        }

        T du_norm = 1e10;
        du_norm = lineSearchNewton(residual);
        return false;
    }
}

bool IMLSContact::linearSolve(StiffnessMatrix& K, const VectorXT& residual, VectorXT& du)
{
    START_TIMING(LinearSolve)
    Eigen::CholmodSupernodalLLT<StiffnessMatrix, Eigen::Lower> solver;
    
    T alpha = 1e-6;
    if (!dynamics)
    {
        StiffnessMatrix H(K.rows(), K.cols());
        H.setIdentity(); H.diagonal().array() = 1e-10;
        K += H;
    }
    solver.analyzePattern(K);
    // T time_analyze = t.elapsed_sec();
    // std::cout << "\t analyzePattern takes " << time_analyze << "s" << std::endl;

    // std::cout << K << std::endl;
    
    int indefinite_count_reg_cnt = 0, invalid_search_dir_cnt = 0, invalid_residual_cnt = 0;
    int i = 0;
    T dot_dx_g = 0.0;
    for (; i < 50; i++)
    {
        solver.factorize(K);
        // std::cout << "factorize" << std::endl;
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
        
        // bool solve_success = true;
        // bool solve_success = (K * du - residual).norm() / residual.norm() < 1e-6;
        bool solve_success = du.norm() < 1e3;
        
        if (!solve_success)
            invalid_residual_cnt++;
        // std::cout << "PD: " << positive_definte << " direction " 
        //     << search_dir_correct_sign << " solve " << solve_success << std::endl;

        if (positive_definte && search_dir_correct_sign && solve_success)
        {
            
            if (verbose)
            {
                std::cout << "\t===== Linear Solve ===== " << std::endl;
                std::cout << "\tnnz: " << K.nonZeros() << std::endl;
                // std::cout << "\t takes " << t.elapsed_sec() << "s" << std::endl;
                std::cout << "\t# regularization step " << i 
                    << " indefinite " << indefinite_count_reg_cnt 
                    << " invalid search dir " << invalid_search_dir_cnt
                    << " invalid solve " << invalid_residual_cnt << std::endl;
                std::cout << "\tdot(search, -gradient) " << dot_dx_g << std::endl;
                std::cout << "\t======================== " << std::endl;
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
        std::cout << "\t===== Linear Solve ===== " << std::endl;
        std::cout << "\tnnz: " << K.nonZeros() << std::endl;
        // std::cout << "\t takes " << t.elapsed_sec() << "s" << std::endl;
        std::cout << "\t# regularization step " << i 
            << " indefinite " << indefinite_count_reg_cnt 
            << " invalid search dir " << invalid_search_dir_cnt
            << " invalid solve " << invalid_residual_cnt << std::endl;
        std::cout << "\tdot(search, -gradient) " << dot_dx_g << std::endl;
        std::cout << "\t======================== " << std::endl;
        FINISH_TIMING_PRINT(LinearSolve)
    }
    return false;
}

T IMLSContact::computeTotalEnergy()
{
    deformed = undeformed + u;

    T energy = 0.0;
    addShellEnergy(energy);
    if (add_gravity)
        addShellGravitionEnergy(energy);
    if (dynamics)
        addInertialEnergy(energy);
    if (BARRIER_ENERGY)
        addFastRIMLSSCTestEnergy(energy);
    energy -= u.dot(external_force);
    return energy;
}

T IMLSContact::computeResidual(VectorXT& residual)
{
    deformed = undeformed + u;
    addShellForceEntry(residual);
    if (add_gravity)
        addShellGravitionForceEntry(residual);
    if (dynamics)
        addInertialForceEntry(residual);
    if (BARRIER_ENERGY)
        addFastRIMLSSCForceEntries(residual);

    if (!run_diff_test)
        iterateDirichletDoF([&](int offset, T target)
        {
            residual[offset] = 0;
        });
    return residual.norm();
}

void IMLSContact::computeLinearModes(MatrixXT& eigen_vectors, VectorXT& eigen_values)
{
    int nmodes = 10;
    StiffnessMatrix K;
    run_diff_test = true;
    buildSystemMatrix(K);
    run_diff_test = false;
    Spectra::SparseSymShiftSolve<T, Eigen::Lower> op(K);
    
    T shift = -1e-4;
    Spectra::SymEigsShiftSolver<Spectra::SparseSymShiftSolve<T, Eigen::Lower>> 
            eigs(op, nmodes, 2 * nmodes, shift);

    eigs.init();

    int nconv = eigs.compute(Spectra::SortRule::LargestMagn);

    if (eigs.info() == Spectra::CompInfo::Successful)
    {
        eigen_vectors = eigs.eigenvectors().real();
        eigen_values = eigs.eigenvalues().real();
    }

    MatrixXT tmp_vec = eigen_vectors;
    VectorXT tmp_val = eigen_values;
    for (int i = 0; i < nmodes; i++)
    {
        eigen_vectors.col(i) = tmp_vec.col(nmodes-i-1);
        eigen_values[i] = tmp_val[nmodes-i-1];
    }
    
}

void IMLSContact::buildSystemMatrix(StiffnessMatrix& K)
{
    deformed = undeformed + u;
    std::vector<Entry> entries;
    addShellHessianEntries(entries);
    if (add_gravity)
        addShellGravitionHessianEntry(entries);
    if (dynamics)
        addInertialHessianEntries(entries);
    if (BARRIER_ENERGY)
        addFastRIMLSSCHessianEntries(entries);
    int n_dof = deformed.rows();
    K.resize(n_dof, n_dof);
    K.setFromTriplets(entries.begin(), entries.end());
    if (!run_diff_test)
        projectDirichletDoFMatrix(K, dirichlet_data);
    
    K.makeCompressed();
}

void IMLSContact::projectDirichletDoFMatrix(StiffnessMatrix& A, const std::unordered_map<int, T>& data)
{
    for (auto iter : data)
    {
        A.row(iter.first) *= 0.0;
        A.col(iter.first) *= 0.0;
        A.coeffRef(iter.first, iter.first) = 1.0;
    }
}

T IMLSContact::lineSearchNewton(const VectorXT& residual)
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

    T alpha = 1.0;
    T norm = du.norm();
    if (verbose)
        std::cout << "\t|du | " << norm << std::endl;

    if(BARRIER_ENERGY)
    {
        //T ipc_step_size = computeCollisionFreeStepsize2D(_u, du);
        T ipc_step_size = computeCollisionFreeStepsizeUnsigned(u, du);
        //std::cout<<"IPC3D result: "<<ipc_step_size2<<std::endl;
        alpha = std::min(alpha, ipc_step_size);
        std::cout << "after barrier step size: " << alpha << std::endl;
    }
    
    T E0 = computeTotalEnergy();
    // std::cout << "obj: " << E0 << std::endl;
   
    int cnt = 0;
    VectorXT u_current = u;
    while (true)
    {
        u = u_current + alpha * du;
        T E1 = computeTotalEnergy();
        if (E1 - E0 < 0 || cnt > 10)
        {
            // if (cnt > 10)
                // std::cout << "cnt > 10" << " |du| " << norm << " |g| " << residual.norm() << std::endl;
            break;
        }
        alpha *= 0.5;
        cnt += 1;
    }

    return alpha * du.norm();
}

T IMLSContact::computeCollisionFreeStepsizeUnsigned(const VectorXT& _u, const VectorXT& du)
{

    // T ipc_step = 1;

    // std::vector<bool> vertex_candidates(num_nodes,1);
    
    // // std::cout<<accumulate(vertex_candidates.begin(), vertex_candidates.end(),0)<<std::endl;

    // Eigen::VectorXd max_step_size(num_nodes);
    // for(int i=0; i<num_nodes; ++i)
    // {
    //     max_step_size[i] = ipc_step;
    // }

    // double max_step = ipc_step;
    // std::vector<int> previous_sign(num_nodes,0);

    // double temp;
    // addFastRIMLSSCEnergy(temp);

    // for(int i=0; i<dist_info.size(); ++i)
    // {
    //     //std::cout<<"1111"<<dist_info[i].first<<" "<<dist_info[i].second<<std::endl;
    //     if(dist_info[i].first == 3216 || dist_info[i].first == 1551) std::cout<<"1111"<<dist_info[i].first<<" "<<dist_info[i].second<<std::endl;
    //     if(dist_info[i].first > 0 && dist_info[i].second > 0){
    //         previous_sign[dist_info[i].first] = 1;
    //     }
    //     else if(dist_info[i].first > 0 && dist_info[i].second < 0){
    //         previous_sign[dist_info[i].first] = -1;
    //     }
    // }

    // std::vector<int> potential_sign_pts;
    // std::vector<int> potential_signs;

    // VectorXa x = deformed+max_step*du;
    // VectorXa bk = deformed;
    // deformed = x;

    // addFastRIMLSSCEnergy(temp);
    // deformed = bk;


    // for(int i=0; i<dist_info.size(); ++i)
    // {
    //     int index = dist_info[i].first;
    //     double distance;
    //     evaluateUnsignedDistanceSqRIMLS(index,max_step,du,max_step,distance);
        
    //     if(dist_info[i].first > 0 && previous_sign[dist_info[i].first] == 0 && dist_info[i].second != 0)
    //     {
    //         potential_sign_pts.push_back(dist_info[i].first);
    //         if(dist_info[i].second > 0) potential_signs.push_back(1);
    //         else if(dist_info[i].second < 0) potential_signs.push_back(-1);
    //     }
    // }

    // // Check previous sign
    // int num_tests = potential_sign_pts.size();

    // tbb::parallel_for(0, num_tests , [&](int i)
    // //for(int i=0; i<num_tests; ++i)
    // {
    //     int index = potential_sign_pts[i];
    //     int sign_start = 0;
    //     int sign_end = potential_signs[i];

    //     double start = 0;
    //     double end = max_step;

    //     while((end-start)*du.norm()>1e-12)
    //     {
            
    //         double mid = (start+end)/2.;
    //         double distance;
    //         int sign_mid = 0;
    //         evaluateUnsignedDistanceSqRIMLS(index,mid,du,max_step,distance);
    //         if(distance > 0) sign_mid = 1;
    //         if(distance < 0) sign_mid = -1;

    //         //if(index == 1291) std::cout<<index<<" "<<start<<" "<<sign_start<<" "<<mid<<" "<<distance<<" "<<end<<" "<<sign_end<<std::endl;

    //         if(sign_mid == 0)
    //         {
    //             start = mid;
    //         }
    //         else if(sign_mid == sign_end)
    //         {
    //             end = mid;
    //         }
    //         else 
    //         {
    //             previous_sign[index] = sign_mid;
    //             break;
    //         }
    //         previous_sign[index] = sign_end;
    //     }
    // }
    // );

    // if(max_step < 1e-2)
    // {
    //     max_step = 0.1;
    //     for(int i=0; i<num_nodes; ++i)
    //     {
    //         max_step_size[i] = max_step;
    //     }
    // } 

    // int iter = 0;
    // double prev_step;

    // while(true)
    // {
    //     prev_step = max_step;
    //     VectorXa x = deformed+max_step*du;
    //     VectorXa bk = deformed;
    //     deformed = x;

    //     double temp;
    //     addFastRIMLSSCEnergy(temp);
        
    //     deformed = bk;

    //     bool is_step_ok = true;
    //     std::vector<int> collision_points;
        
    //     for(int i=0; i<dist_info.size(); ++i)
    //     {
    //         if(dist_info[i].first > 0){
    //             if(previous_sign[dist_info[i].first] == 0 || previous_sign[dist_info[i].first]*dist_info[i].second > 0) continue;
    //             std::cout<<dist_info[i].first<<" "<<dist_info[i].second<<" sign "<<previous_sign[dist_info[i].first]<<std::endl;
    //             collision_points.push_back(dist_info[i].first);
    //             is_step_ok = false;
    //         }
    //     }

    //     if(is_step_ok) break;
    //     int num_tests = collision_points.size();

    //     tbb::parallel_for(0, num_tests , [&](int i)
    //     {
    //         int k = collision_points[i];
    //         // double parameters[1] = {0.0};
    //         double max_step_local = std::min(max_step, max_step_size.minCoeff());

    //         auto func = [this,k,du,max_step_local,temp](double t)->std::tuple<double, double>
    //         {
    //             double temp;
    //             double fx = this->evaluateUnsignedDistanceSqRIMLS(k, t, du, max_step_local, temp);
    //             double dfx = this->evaluateUnsignedDistanceSqGradientRIMLS(k, t, du, max_step_local);
    //             return {fx,dfx};
    //         };

    //         double guess = max_step_local;
    //         double lower = 0.0;
    //         double upper = max_step_local;
    //         int digits = std::numeric_limits<double>::digits;
    //         boost::uintmax_t max_iter = 20;

    //         try {
    //             double root = boost::math::tools::newton_raphson_iterate(
    //                 func, guess, lower, upper, digits, max_iter);
    //             std::cout << "Root found: " << root << "\n";
    //             if (root > 0 && root < max_step_local) {
    //                 max_step_size(k) = 0.9 * root;
    //             } else {
    //                 max_step_size(k) = max_step_local;
    //             }
    //         }catch (const std::exception& e) {
    //             max_step_size(k) = max_step_local;
    //         }
    //         // ceres::GradientProblemSolver::Options options;
    //         // options.max_num_iterations = 100;
    //         // options.minimizer_progress_to_stdout = false;
    //         // options.gradient_tolerance = 1e-12;
    //         // options.function_tolerance = 1e-12;
    //         // options.parameter_tolerance = 1e-12;
    //         // options.line_search_direction_type = ceres::LBFGS;

    //         // ceres::GradientProblemSolver::Summary summary;

    //         // RosenbrockCost* cost_function = new RosenbrockCost();
    //         // cost_function->solver = this;
    //         // cost_function->index = k;
    //         // cost_function->du = du;
    //         // cost_function->ipc_stepsize = max_step_local;

    //         // ceres::GradientProblem problem(cost_function);

    //         // ceres::Solve(options, problem, parameters, &summary);

    //         // Update step size based on result
    //         // if (summary.final_cost < 1e-10 && parameters[0] > 0 && parameters[0] < max_step_local) {
    //         //     max_step_size(k) = 0.9 * parameters[0];
    //         // } else {
    //         //     max_step_size(k) = max_step_local;
    //         // }
    //     }
    //     );
        

    //     //max_step /= 2;

    //     Eigen::Index index;
    //     max_step = max_step_size.minCoeff(&index);
    //     // std::cout<<"Minimal Step Size Index: "<<index<<" max step "<<max_step<<std::endl;
    //     if(max_step == prev_step && iter>1)
    //     {
    //         max_step *= 0.5;
    //     }
    //     iter ++;
    // }

    // return max_step;

    double min_step = 1.0;
    int iter = 0;

    while(true)
    {
        VectorXa x = undeformed+_u+min_step*du;
        VectorXa bk = deformed;
        deformed = x;
        double temp;
        addFastRIMLSSCEnergy(temp);
        bool is_step_ok = true;
        for(int i=0; i<dist_info.size(); ++i)
        {
            if(dist_info[i].first > 0 && dist_info[i].second <= 0){
                std::cout<<dist_info[i].first<<" "<<dist_info[i].second<<std::endl;
                is_step_ok = false;
                break;
            }
        }
        if(is_step_ok || iter >= 30){break;} 
        if(iter >= 30){min_step = 0; break;} 
        iter++;
        min_step /= 2.0;
        deformed = bk;
    }
    
    return min_step;
}

double sign(double t)
{
    if (t<0) return -1.;
    else if (t>0) return 1.;
    else return 0.;
}

double IMLSContact::evaluateUnsignedDistanceSqRIMLS(int i, double t, const VectorXa& du, double ipc_stepsize, double& distance)
{
    //if(t >= 1 || t < 0) return 1;
    double t0 = t;
    t = min(t,ipc_stepsize);
    t = max(t,0.);
    VectorXa deformed_cur = deformed+t*du;
    // Tree local_accelerator;
    // local_accelerator.clear();

    // for(int i=0; i<num_nodes; ++i)
    // {
    //     if(is_surface_vertex[i] == 0) continue;
    //     Eigen::VectorXd point = deformed.segment<dim>(dim*i);
    //     local_accelerator.insert(std::make_tuple(Point_d(point[0], point[1], point[2]), i));
    // }

    std::vector<int> local_current_index;

    Vector3a xi = deformed_cur.segment<3>(3*i);

    //Fuzzy_sphere fs(Point_d(xi[0], xi[1], xi[2]), radius+radius*0.2, radius*0.2);

    // std::vector<std::tuple<Point_d, int>> points_query;
    // accelerator.search(std::back_inserter(points_query), fs);

    //for(int k=0; k<points_query.size(); ++k)
    for(int k=0; k<num_nodes; ++k)
    {
        Vector3a ck = deformed_cur.segment<3>(3*k);
        //int index = std::get<1>(points_query[k]);
        int index = k;
        if(is_surface_vertex[index] && (geodist_close_matrix.coeff(i,index) == 0) && (ck-xi).norm()<=1.2*radius)
        {
            local_current_index.push_back(index);
        }     
    }

    //if(local_current_index.size() == 0) return 1.+fabs(t-t0);
    if(local_current_index.size() == 0)
    {
        distance = 0;
        return 1.+fabs(t-t0);
    }


    AScalar fx = 0;
    AScalar gx = 0;
    VectorXa N(3); N.setZero();
    VectorXa GF(3); GF.setZero();
    VectorXa GW(3); GW.setZero();
    AScalar fx2 = 0;
    AScalar gx2 = 0;

    for(int k=0; k<local_current_index.size(); ++k)
    {
        Vector3a ck = deformed_cur.segment<3>(local_current_index[k]*3);
        int valence = nodes_adjacency[local_current_index[k]].size();
        if(valence > 11) std::cout<<"WTF!!!!"<<std::endl;
        VectorXa vs;

        vs.resize(33);
        for(int a=0; a<11; ++a)
        {
            vs.segment<3>(a*3) = ck;
        }


        for(int a=0; a<nodes_adjacency[local_current_index[k]].size(); ++a)
        {
            vs.segment<3>(a*3) = deformed_cur.segment<3>(nodes_adjacency[local_current_index[k]][a]*3);
        }

        if((ck-xi).norm()<=radius)
        {
            fx += fx_func11(xi, vs, ck, radius);
            gx += gx_func11(xi, ck, radius);
            GF(0) += sumGF1_func(xi, vs, ck, radius); GF(1) += sumGF2_func(xi, vs, ck, radius); GF(2) += sumGF3_func(xi, vs, ck, radius);
            GW(0) += sumGW1_func(xi, vs, ck, radius); GW(1) += sumGW2_func(xi, vs, ck, radius); GW(2) += sumGW3_func(xi, vs, ck, radius);
            N(0) += sumN1_func(xi, vs, ck, radius); N(1) += sumN2_func(xi, vs, ck, radius); N(2) += sumN3_func(xi, vs, ck, radius);
        }
    }


    if(abs(gx)<1e-6)
    {
        distance = 0;
        return 1.+fabs(t-t0);
    } 

    AScalar f = fx/gx;
    VectorXa grad_f = (GF-f*GW+N)/gx;

    for(int k=0; k<local_current_index.size(); ++k)
    {
        Vector3a ck = deformed_cur.segment<3>(local_current_index[k]*3);
        int valence = nodes_adjacency[local_current_index[k]].size();
        VectorXa vs;

        vs.resize(33);
        for(int a=0; a<11; ++a)
        {
            vs.segment<3>(a*3) = ck;
        }


        for(int a=0; a<nodes_adjacency[local_current_index[k]].size(); ++a)
        {
            vs.segment<3>(a*3) = deformed_cur.segment<3>(nodes_adjacency[local_current_index[k]][a]*3);
        }

        if((ck-xi).norm()<=radius)
        {
            fx2 += fx_func_2(xi, vs, ck, radius,f,grad_f,sigma_r,sigma_n);
            gx2 += gx_func_2(xi, vs, ck, radius,f,grad_f,sigma_r,sigma_n);
        }
    }

    if(abs(gx2)<=1e-6)
    {
        distance = 0;
        return 1.+fabs(t-t0);
    } 
    AScalar dist = fx2/(gx2*radius);
    distance = fx2/gx2;
    
    return fabs(dist)+fabs(t-t0);
}

double IMLSContact::evaluateUnsignedDistanceSqGradientRIMLS(int i, double t, const VectorXa& du, double ipc_stepsize)
{
    using namespace std;
    if(t >= ipc_stepsize || t <= 0.) return sign(t);
    t = min(t,ipc_stepsize);
    t = max(t,0.);
    VectorXa deformed_cur = deformed+t*du;
    // Tree local_accelerator;
    // local_accelerator.clear();

    // for(int i=0; i<num_nodes; ++i)
    // {
    //     if(is_surface_vertex[i] == 0) continue;
    //     Eigen::VectorXd point = deformed.segment<dim>(dim*i);
    //     local_accelerator.insert(std::make_tuple(Point_d(point[0], point[1], point[2]), i));
    // }

	std::vector<int> local_current_index;

    Vector3a xi = deformed_cur.segment<3>(3*i);

    //Fuzzy_sphere fs(Point_d(xi[0], xi[1], xi[2]), radius+radius*0.2, radius*0.2);

    // std::vector<std::tuple<Point_d, int>> points_query;
    // accelerator.search(std::back_inserter(points_query), fs);

    //for(int k=0; k<points_query.size(); ++k)
    for(int k=0; k<num_nodes; ++k)
    {
        Vector3a ck = deformed_cur.segment<3>(3*k);
        //int index = std::get<1>(points_query[k]);
        int index = k;
        if(is_surface_vertex[index] && (geodist_close_matrix.coeff(i,index) == 0) && (ck-xi).norm()<=radius)
        {
            local_current_index.push_back(index);
        }     
    }

    if(local_current_index.size() == 0) return sign(t);

    VectorXa gradient(num_nodes*dim);
    gradient.setZero();
    VectorXa sum_dfdx(dim*num_nodes);
    sum_dfdx.setZero();
    VectorXa sum_dgdx(dim*num_nodes);
    sum_dgdx.setZero();
    VectorXa sum_dGF1dx(dim*num_nodes);
    sum_dGF1dx.setZero();
    VectorXa sum_dGF2dx(dim*num_nodes);
    sum_dGF2dx.setZero();
    VectorXa sum_dGF3dx(dim*num_nodes);
    sum_dGF3dx.setZero();
    VectorXa sum_dGW1dx(dim*num_nodes);
    sum_dGW1dx.setZero();
    VectorXa sum_dGW2dx(dim*num_nodes);
    sum_dGW2dx.setZero();
    VectorXa sum_dGW3dx(dim*num_nodes);
    sum_dGW3dx.setZero();
    VectorXa sum_dN1dx(dim*num_nodes);
    sum_dN1dx.setZero();
    VectorXa sum_dN2dx(dim*num_nodes);
    sum_dN2dx.setZero();
    VectorXa sum_dN3dx(dim*num_nodes);
    sum_dN3dx.setZero();
    VectorXa sum_df2dx(dim*num_nodes);
    sum_df2dx.setZero();
    VectorXa sum_dg2dx(dim*num_nodes);
    sum_dg2dx.setZero();

    //std::cout<<i<<" out of "<<slave_nodes_3d[0].size()<<std::endl;
    AScalar fx = 0;
    AScalar gx = 0;
    VectorXa N(3); N.setZero();
    VectorXa GF(3); GF.setZero();
    VectorXa GW(3); GW.setZero();
    AScalar fx2 = 0;
    AScalar gx2 = 0;
    VectorXa ele_dfdx;
    VectorXa ele_dgdx;
    VectorXa ele_df2dx;
    VectorXa ele_dg2dx;
    VectorXa ele_dGF1dx,ele_dGF2dx,ele_dGF3dx;
    VectorXa ele_dGW1dx,ele_dGW2dx,ele_dGW3dx;
    VectorXa ele_dN1dx,ele_dN2dx,ele_dN3dx;
    //std::cout<<i<<" out of "<<slave_nodes_3d[0].size()<<std::endl;

    for(int k=0; k<local_current_index.size(); ++k)
    {
        Vector3a ck = deformed_cur.segment<3>(local_current_index[k]*3);
        //f(i == 1567) std::cout<<k<<" k out of "<<local_current_index.size()<<std::endl;
        int valence = nodes_adjacency[local_current_index[k]].size();
        Eigen::VectorXi valence_indices;
        VectorXa vs;

        ele_dfdx.resize(39);
        ele_dgdx.resize(39);
        vs.resize(33);
        valence_indices.resize(13);
        for(int a=0; a<11; ++a)
        {
            vs.segment<3>(a*3) = ck;
            valence_indices(a+2)=(local_current_index[k]);
        }


        valence_indices(0)=(i);
        valence_indices(1)=(local_current_index[k]);


        for(int a=0; a<nodes_adjacency[local_current_index[k]].size(); ++a)
        {
            vs.segment<3>(a*3) = deformed_cur.segment<3>(nodes_adjacency[local_current_index[k]][a]*3);
            valence_indices(a+2) = nodes_adjacency[local_current_index[k]][a];
        }

        if((ck-xi).norm()<=radius) {
            ele_dfdx = dfdx_func11(xi, vs, ck, radius);
            ele_dgdx = dgdx_func11(xi, ck, radius);
            ele_dGF1dx = dsumGF1dx_func(xi, vs, ck, radius); ele_dGF2dx = dsumGF2dx_func(xi, vs, ck, radius); ele_dGF3dx = dsumGF3dx_func(xi, vs, ck, radius);
            ele_dGW1dx = dsumGW1dx_func(xi, vs, ck, radius); ele_dGW2dx = dsumGW2dx_func(xi, vs, ck, radius); ele_dGW3dx = dsumGW3dx_func(xi, vs, ck, radius);
            ele_dN1dx = dsumN1dx_func(xi, vs, ck, radius); ele_dN2dx = dsumN2dx_func(xi, vs, ck, radius); ele_dN3dx = dsumN3dx_func(xi, vs, ck, radius);
            fx += fx_func11(xi, vs, ck, radius);
            gx += gx_func11(xi, ck, radius);
            GF(0) += sumGF1_func(xi, vs, ck, radius); GF(1) += sumGF2_func(xi, vs, ck, radius); GF(2) += sumGF3_func(xi, vs, ck, radius);
            GW(0) += sumGW1_func(xi, vs, ck, radius); GW(1) += sumGW2_func(xi, vs, ck, radius); GW(2) += sumGW3_func(xi, vs, ck, radius);
            N(0) += sumN1_func(xi, vs, ck, radius); N(1) += sumN2_func(xi, vs, ck, radius); N(2) += sumN3_func(xi, vs, ck, radius);
            IMLS_local_gradient_to_global_gradient(ele_dfdx,valence_indices,dim,sum_dfdx);
            IMLS_local_gradient_to_global_gradient(ele_dgdx,valence_indices,dim,sum_dgdx);
            IMLS_local_gradient_to_global_gradient(ele_dGF1dx,valence_indices,dim,sum_dGF1dx);
            IMLS_local_gradient_to_global_gradient(ele_dGF2dx,valence_indices,dim,sum_dGF2dx);
            IMLS_local_gradient_to_global_gradient(ele_dGF3dx,valence_indices,dim,sum_dGF3dx);
            IMLS_local_gradient_to_global_gradient(ele_dGW1dx,valence_indices,dim,sum_dGW1dx);
            IMLS_local_gradient_to_global_gradient(ele_dGW2dx,valence_indices,dim,sum_dGW2dx);
            IMLS_local_gradient_to_global_gradient(ele_dGW3dx,valence_indices,dim,sum_dGW3dx);
            IMLS_local_gradient_to_global_gradient(ele_dN1dx,valence_indices,dim,sum_dN1dx);
            IMLS_local_gradient_to_global_gradient(ele_dN2dx,valence_indices,dim,sum_dN2dx);
            IMLS_local_gradient_to_global_gradient(ele_dN3dx,valence_indices,dim,sum_dN3dx);

        }
    }
    if(abs(gx)<1e-6) return sign(t);

    AScalar f = fx/gx;
    VectorXa grad_f = (GF-f*GW+N)/gx;

    VectorXa dfdx = (gx*sum_dfdx-fx*sum_dgdx)/(pow(gx,2));
    // A = GF/gx
    VectorXa dA1dx = (gx*sum_dGF1dx-GF(0)*sum_dgdx)/(pow(gx,2)); VectorXa dA2dx = (gx*sum_dGF2dx-GF(1)*sum_dgdx)/(pow(gx,2)); VectorXa dA3dx = (gx*sum_dGF3dx-GF(2)*sum_dgdx)/(pow(gx,2));
    // B = GW/gx
    VectorXa dB1dx = (gx*sum_dGW1dx-GW(0)*sum_dgdx)/(pow(gx,2)); VectorXa dB2dx = (gx*sum_dGW2dx-GW(1)*sum_dgdx)/(pow(gx,2)); VectorXa dB3dx = (gx*sum_dGW3dx-GW(2)*sum_dgdx)/(pow(gx,2));
    // C = N/gx
    VectorXa dC1dx = (gx*sum_dN1dx-N(0)*sum_dgdx)/(pow(gx,2)); VectorXa dC2dx = (gx*sum_dN2dx-N(1)*sum_dgdx)/(pow(gx,2)); VectorXa dC3dx = (gx*sum_dN3dx-N(2)*sum_dgdx)/(pow(gx,2));

    VectorXa dgrad_f1dx = dA1dx - dfdx*(GW(0)/gx) - f*dB1dx + dC1dx;
    VectorXa dgrad_f2dx = dA2dx - dfdx*(GW(1)/gx) - f*dB2dx + dC2dx;
    VectorXa dgrad_f3dx = dA3dx - dfdx*(GW(2)/gx) - f*dB3dx + dC3dx;

    for(int k=0; k<local_current_index.size(); ++k)
    {
        Vector3a ck = deformed_cur.segment<3>(local_current_index[k]*3);
        int valence = nodes_adjacency[local_current_index[k]].size();
        Eigen::VectorXi valence_indices;
        VectorXa vs;

        ele_df2dx.resize(43); ele_dg2dx.resize(43);
        vs.resize(33);
        valence_indices.resize(13);
        for(int a=0; a<11; ++a)
        {
            vs.segment<3>(a*3) = ck;
            valence_indices(a+2)=(local_current_index[k]);
        }


        valence_indices(0)=(i);
        valence_indices(1)=(local_current_index[k]);


        for(int a=0; a<nodes_adjacency[local_current_index[k]].size(); ++a)
        {
            vs.segment<3>(a*3) = deformed_cur.segment<3>(nodes_adjacency[local_current_index[k]][a]*3);
            valence_indices(a+2) = nodes_adjacency[local_current_index[k]][a];
        }

        AScalar df2df = 0; AScalar dg2df = 0;
        Vector3a df2dgrad_f; Vector3a dg2dgrad_f;
        if((ck-xi).norm()<=radius)
        {
            // df2dx contains df2/df df2/dgrad_f1 ...
            VectorXa temp1 = dfdx_func_2(xi, vs, ck, radius,f,grad_f,sigma_r,sigma_n);
            VectorXa temp2 = dgdx_func_2(xi, vs, ck, radius,f,grad_f,sigma_r,sigma_n);
            df2df = temp1(39);  dg2df = temp2(39);
            df2dgrad_f = temp1.segment<3>(40); dg2dgrad_f = temp2.segment<3>(40);
            ele_df2dx = temp1.segment<39>(0); ele_dg2dx = temp2.segment<39>(0);
            fx2 += fx_func_2(xi, vs, ck, radius,f,grad_f,sigma_r,sigma_n);
            gx2 += gx_func_2(xi, vs, ck, radius,f,grad_f,sigma_r,sigma_n);
            IMLS_local_gradient_to_global_gradient(ele_df2dx,valence_indices,dim,sum_df2dx);
            IMLS_local_gradient_to_global_gradient(ele_dg2dx,valence_indices,dim,sum_dg2dx);
            sum_df2dx += (df2df * dfdx + df2dgrad_f(0) * dgrad_f1dx + df2dgrad_f(1) * dgrad_f2dx + df2dgrad_f(2) * dgrad_f3dx);
            sum_dg2dx += (dg2df * dfdx + dg2dgrad_f(0) * dgrad_f1dx + dg2dgrad_f(1) * dgrad_f2dx + dg2dgrad_f(2) * dgrad_f3dx);
        }
    }

    AScalar dist = fx2/gx2;

    if(abs(gx2) <= 1e-6) return sign(t);
               

    if(dist > 0)
        gradient = (gx2*sum_df2dx-fx2*sum_dg2dx)/(pow(gx2,2)*pow(radius,1));
    else
        gradient = -(gx2*sum_df2dx-fx2*sum_dg2dx)/(pow(gx2,2)*pow(radius,1));

    return gradient.dot(du);
}

void IMLSContact::initializeFromFile(const std::string& filename)
{
    MatrixXT V; MatrixXi F;
    std::vector<Eigen::Triplet<double>> tripletList;
    igl::readOBJ(filename, V, F);
    std::cout<<"Done"<<std::endl;

    TV min_corner = V.colwise().minCoeff();
    TV max_corner = V.colwise().maxCoeff();
    T bb_diag = (max_corner - min_corner).norm();

    V *= 1.0 / bb_diag;

    // V *= 0.5;

    search_radius = 0.2;
    radius = 0.01;

    auto rotationMatrixFromEulerAngle = [](T angle_z, T angle_y, T angle_x)
    {
        Eigen::Matrix3d R, yaw, pitch, roll;
        yaw.setZero(); pitch.setZero(); roll.setZero();
        yaw(0, 0) = cos(angle_z);	yaw(0, 1) = -sin(angle_z);
        yaw(1, 0) = sin(angle_z);	yaw(1, 1) = cos(angle_z);
        yaw(2, 2) = 1.0;
        //y rotation
        pitch(0, 0) = cos(angle_y); pitch(0, 2) = sin(angle_y);
        pitch(1, 1) = 1.0;
        pitch(2, 0) = -sin(angle_y); pitch(2, 2) = cos(angle_y);
        //x rotation
        roll(0, 0) = 1.0;
        roll(1, 1) = cos(angle_x); roll(1, 2) = -sin(angle_x);
        roll(2, 1) = sin(angle_x); roll(2, 2) = cos(angle_x);
        R = yaw * pitch * roll;
        return R;
    };

    // TM R = rotationMatrixFromEulerAngle(0, 0, M_PI_2);
    // for (int i = 0; i < V.rows(); i++)
    // {
    //     V.row(i) = (R * V.row(i).transpose()).transpose();
    // }
    

    iglMatrixFatten<T, 3>(V, undeformed);
    iglMatrixFatten<int, 3>(F, faces);

    min_corner = V.colwise().minCoeff();
    max_corner = V.colwise().maxCoeff();
    num_nodes = V.rows();

    // Build the dirichlet data structure
    deformed = undeformed;
    for (int i = 0; i < V.rows(); i++)
    {
        TV xi = deformed.segment<3>(i * 3);
        if ( fabs(xi[0] -  min_corner[0])<1e-3 || fabs(xi[0] -  max_corner[0])<1e-3)
        {
            dirichlet_data[i * 3] = 0.0;
            dirichlet_data[i * 3 + 1] = 0.0;
            dirichlet_data[i * 3 + 2] = 0.0;
        }
    }
    std::cout << dirichlet_data.size() << std::endl;

    // Initialize geodist matrix
    std::unique_ptr<SurfaceMesh> mesh;
    std::unique_ptr<VertexPositionGeometry> geometry;
    std::tie(mesh, geometry) = readSurfaceMesh(filename);

    HeatMethodDistanceSolver heatSolver(*geometry);

    for(int i=0; i<V.rows(); ++i)
    {
        nodes_adjacency.push_back(std::vector<int>());
        is_surface_vertex.push_back(0);
    }
        

    int num = 0;
    int node_start = 0;
    int face_start = 0;
    for(Vertex v:mesh->vertices())
    {
        VertexData<double> distToSource = heatSolver.computeDistance(v);
        is_surface_vertex[num+node_start] = 1;
        for(int j=0; j<distToSource.size(); ++j)
        {
            if(distToSource[j] <= search_radius)
            {
                //std::cout<<v.getIndex()+node_start<<" "<<j+node_start<<std::endl;
                //if(mesh_dim == 2)
                tripletList.push_back(Eigen::Triplet<double>(v.getIndex()+node_start,j+node_start,1));
                // else if(mesh_dim == 3)
                //     tripletList.push_back(Eigen::Triplet<double>(v.getIndex()+node_start,j+node_start,1));
            }
            // else if((V.row(J(num))-V.row(J(j))).norm() <= radius && N.row(j).dot(V.row(J(num))-V.row(J(j))) <= 0)
            // {
            //     //std::cout<<v.getIndex()+node_start<<" "<<j+node_start<<std::endl;
            //     //if(mesh_dim == 2)
            //     tripletList.push_back(Eigen::Triplet<double>(J(v.getIndex())+node_start,J(j)+node_start,1));
            //     // else if(mesh_dim == 3)
            //     //     tripletList.push_back(Eigen::Triplet<double>(v.getIndex()+node_start,j+node_start,1));
            // }
        }
        if(num == 0){
            std::cout<<distToSource.size()<<std::endl;
        } 
        num++;
    }
    std::cout<<num<<std::endl;

    // for(int i=0; i<V.rows(); ++i)
    // {
    //     if(fabs(V(i,0)-0)< 1e-5 || fabs(V(i,0)-2)< 1e-5)
    //     {
    //         is_surface_vertex[i+node_start] = 0;
    //     } 
    // }

    geodist_close_matrix.resize(V.rows(),F.rows());
    geodist_close_matrix.setFromTriplets(tripletList.begin(),tripletList.end());
        


    mesh_begin_vertices.push_back(node_start);
    std::vector<std::vector<int>> current_node_adjancy;
    igl::adjacency_list(F,current_node_adjancy,true);
    //igl::adjacency_list(F,current_node_adjancy,true);
    for(int j=0; j<current_node_adjancy.size(); ++j)
    {
        for(int k=0; k<current_node_adjancy[j].size(); ++k)
        {
            nodes_adjacency[j+node_start].push_back(current_node_adjancy[j][k] + node_start);
        }
        if(nodes_adjacency[j+node_start].size()>11) std::cout<<"Too Much!!!"<<std::endl;
    }

    double temp = 0;
    addFastRIMLSSCEnergy(temp);
    
    u = VectorXT::Zero(deformed.rows());
    external_force = VectorXT::Zero(deformed.rows());
    // external_force[5 * 3 + 2] = -10.0;
    buildHingeStructure();
    dynamics = true;
    add_gravity = true;
    use_consistent_mass_matrix = true;
    // E = 0.0;
    dt = 1e-1;
    simulation_duration = 1000000;
    
    hinge_stiffness.setConstant(10);
    if (dynamics)
    {
        initializeDynamicStates();
    }
    
}

void IMLSContact::buildHingeStructure()
{
    struct Hinge
	{
		Hinge()
		{
			for (int i = 0; i < 2; i++)
			{
				edge[i] = -1;
				flaps[i] = -1;
				tris[i] = -1;
			}
		}
		int edge[2];
		int flaps[2];
		int tris[2];
	};
	
	std::vector<Hinge> hinges_temp;
	
	hinges_temp.clear();
	std::map<std::pair<int, int>, int> edge2index;
	for (int i = 0; i < faces.size() / 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			int i1 = faces(3 * i + j);
			int i2 = faces(3 * i + (j + 1) % 3);
			int i1t = i1;
			int i2t = i2;
			bool swapped = false;
			if (i1t > i2t)
			{
				std::swap(i1t, i2t);
				swapped = true;
			}
			
			auto ei = std::make_pair(i1t, i2t);
			auto ite = edge2index.find(ei);
			if (ite == edge2index.end())
			{
				//insert new hinge
				edge2index[ei] = hinges_temp.size();
				hinges_temp.push_back(Hinge());
				Hinge& hinge = hinges_temp.back();
				hinge.edge[0] = i1t;
				hinge.edge[1] = i2t;
				int itmp = swapped ? 1 : 0;
				hinge.tris[itmp] = i;
				hinge.flaps[itmp] = faces(3 * i + (j + 2) % 3);
			}
			else
			{
				//hinge for this edge already exists, add missing information for this triangle
				Hinge& hinge = hinges_temp[ite->second];
				int itmp = swapped ? 1 : 0;
				hinge.tris[itmp] = i;
				hinge.flaps[itmp] = faces(3 * i + (j + 2) % 3);
			}
		}
	}
	//ordering for edges
	
	hinges.resize(hinges_temp.size(), Eigen::NoChange);
	int ii = 0;
	/*
      auto diff code takes
           x3
         /   \
        x2---x1
         \   /
           x0	

      hinge is 
           x2
         /   \
        x0---x1
         \   /
           x3	
    */
    for(Hinge & hinge : hinges_temp) {
		if ((hinge.tris[0] == -1) || (hinge.tris[1] == -1)) {
			continue; //skip boundary edges
		}
		hinges(ii, 2) = hinge.edge[0]; //x0
		hinges(ii, 1) = hinge.edge[1]; //x1
		hinges(ii, 3) = hinge.flaps[0]; //x2
		hinges(ii, 0) = hinge.flaps[1]; //x3
		++ii;
	}
	hinges.conservativeResize(ii, Eigen::NoChange);
    hinge_stiffness.resize(hinges.rows());
    hinge_stiffness.setOnes();
}

void IMLSContact::addShellInplaneEnergy(T& energy)
{
    iterateFaceSerial([&](int face_idx)
    {
        FaceVtx vertices = getFaceVtxDeformed(face_idx);
        FaceVtx undeformed_vertices = getFaceVtxUndeformed(face_idx);

        FaceIdx indices = faces.segment<3>(face_idx * 3);
        
        TV x0 = vertices.row(0); TV x1 = vertices.row(1); TV x2 = vertices.row(2);
        TV X0 = undeformed_vertices.row(0); 
        TV X1 = undeformed_vertices.row(1); 
        TV X2 = undeformed_vertices.row(2);

        T k_s = E * thickness / (1.0 - nu * nu);

        energy += compute3DCSTShellEnergy(nu, k_s, x0, x1, x2, X0, X1, X2);

        
    });
}

void IMLSContact::addShellBendingEnergy(T& energy)
{
    iterateHingeSerial([&](const HingeIdx& hinge_idx, int hinge_cnt){
        
        HingeVtx deformed_vertices = getHingeVtxDeformed(hinge_idx);
        HingeVtx undeformed_vertices = getHingeVtxUndeformed(hinge_idx);


        T k_bend = hinge_stiffness[hinge_cnt] * E * std::pow(thickness, 3) / (24 * (1.0 - std::pow(nu, 2)));

        TV x0 = deformed_vertices.row(0); TV x1 = deformed_vertices.row(1); TV x2 = deformed_vertices.row(2); TV x3 = deformed_vertices.row(3);
        TV X0 = undeformed_vertices.row(0); TV X1 = undeformed_vertices.row(1); TV X2 = undeformed_vertices.row(2); TV X3 = undeformed_vertices.row(3);

        energy += computeDSBendingEnergy(k_bend, x0, x1, x2, x3, X0, X1, X2, X3);

    });
}

void IMLSContact::addShellEnergy(T& energy)
{
    T in_plane_energy = 0.0;
    T bending_energy = 0.0;
    addShellInplaneEnergy(in_plane_energy);
    addShellBendingEnergy(bending_energy);    

    energy += in_plane_energy;
    energy += bending_energy;
}

void IMLSContact::addShellInplaneForceEntries(VectorXT& residual)
{
    iterateFaceSerial([&](int face_idx)
    {
        FaceVtx vertices = getFaceVtxDeformed(face_idx);
        FaceVtx undeformed_vertices = getFaceVtxUndeformed(face_idx);
        FaceIdx indices = faces.segment<3>(face_idx * 3);

        T k_s = E * thickness / (1.0 - nu * nu);
        TV x0 = vertices.row(0); TV x1 = vertices.row(1); TV x2 = vertices.row(2);
        TV X0 = undeformed_vertices.row(0); TV X1 = undeformed_vertices.row(1); TV X2 = undeformed_vertices.row(2);
        
        Vector9a dedx;
        compute3DCSTShellEnergyGradient(nu, k_s, x0, x1, x2, X0, X1, X2, dedx);
        
        addForceEntry<3>(residual, {indices[0], indices[1], indices[2]}, -dedx);
    });
}

void IMLSContact::addShellBendingForceEntries(VectorXT& residual)
{
    iterateHingeSerial([&](const HingeIdx& hinge_idx, int hinge_cnt){
                
        HingeVtx deformed_vertices = getHingeVtxDeformed(hinge_idx);
        HingeVtx undeformed_vertices = getHingeVtxUndeformed(hinge_idx);


        T k_bend = hinge_stiffness[hinge_cnt] * E * std::pow(thickness, 3) / (24 * (1.0 - std::pow(nu, 2)));

        Vector12a dedx;
        
        TV x0 = deformed_vertices.row(0); TV x1 = deformed_vertices.row(1); TV x2 = deformed_vertices.row(2); TV x3 = deformed_vertices.row(3);
        TV X0 = undeformed_vertices.row(0); TV X1 = undeformed_vertices.row(1); TV X2 = undeformed_vertices.row(2); TV X3 = undeformed_vertices.row(3);

        computeDSBendingEnergyGradient(k_bend, x0, x1, x2, x3, X0, X1, X2, X3, dedx);
        addForceEntry<3>(residual, 
            {hinge_idx[0], hinge_idx[1], hinge_idx[2], hinge_idx[3]}, -dedx);
    });
}
void IMLSContact::addShellForceEntry(VectorXT& residual)
{
    addShellInplaneForceEntries(residual);
    addShellBendingForceEntries(residual);
}

void IMLSContact::addShellInplaneHessianEntries(std::vector<Entry>& entries)
{
    iterateFaceSerial([&](int face_idx)
    {
        FaceVtx vertices = getFaceVtxDeformed(face_idx);
        FaceVtx undeformed_vertices = getFaceVtxUndeformed(face_idx);

        FaceIdx indices = faces.segment<3>(face_idx * 3);

        T k_s = E * thickness / (1.0 - nu * nu);

        TV x0 = vertices.row(0); TV x1 = vertices.row(1); TV x2 = vertices.row(2);
        TV X0 = undeformed_vertices.row(0); TV X1 = undeformed_vertices.row(1); TV X2 = undeformed_vertices.row(2);

        Matrix<T, 9, 9> hessian;
        compute3DCSTShellEnergyHessian(nu, k_s, x0, x1, x2, X0, X1, X2, hessian);


        addHessianEntry<3, 3>(entries, {indices[0], indices[1], indices[2]}, hessian);

    });
}

void IMLSContact::addShellBendingHessianEntries(std::vector<Entry>& entries)
{
    iterateHingeSerial([&](const HingeIdx& hinge_idx, int hinge_cnt){
        
        HingeVtx deformed_vertices = getHingeVtxDeformed(hinge_idx);
        HingeVtx undeformed_vertices = getHingeVtxUndeformed(hinge_idx);

        T k_bend = hinge_stiffness[hinge_cnt] * E * std::pow(thickness, 3) / (24 * (1.0 - std::pow(nu, 2)));
        
        Matrix<T, 12, 12> hess;
        TV x0 = deformed_vertices.row(0); TV x1 = deformed_vertices.row(1); TV x2 = deformed_vertices.row(2); TV x3 = deformed_vertices.row(3);
        TV X0 = undeformed_vertices.row(0); TV X1 = undeformed_vertices.row(1); TV X2 = undeformed_vertices.row(2); TV X3 = undeformed_vertices.row(3);

        computeDSBendingEnergyHessian(k_bend, x0, x1, x2, x3, X0, X1, X2, X3, hess);
        addHessianEntry<3, 3>(entries, 
                            {hinge_idx[0], hinge_idx[1], hinge_idx[2], hinge_idx[3]}, 
                            hess);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          

    });
}

void IMLSContact::addShellHessianEntries(std::vector<Entry>& entries)
{
    addShellInplaneHessianEntries(entries);
    addShellBendingHessianEntries(entries);
}

void IMLSContact::setHingeStiffness()
{
    int dir = 1;
    T eps = 0.02;
    TV min_corner, max_corner;
    computeBoundingBox(min_corner, max_corner);
    TV center = 0.5 * (min_corner + max_corner);
    iterateHingeSerial([&](const HingeIdx& hinge_idx, int hinge_cnt){
        
        HingeVtx undeformed_vertices = getHingeVtxUndeformed(hinge_idx);
        TV v0 = undeformed_vertices.row(1);
        TV v1 = undeformed_vertices.row(2);
        bool center_v0 = v0[dir] < center[dir] + eps && v0[dir] > center[dir] - eps;
        bool center_v1 = v1[dir] < center[dir] + eps && v1[dir] > center[dir] - eps;
        if (center_v0 && center_v1)
        {
            hinge_stiffness[hinge_cnt] = E;
            std::cout << "center" << std::endl;
        }
    });
}

void IMLSContact::computeBoundingBox(TV& min_corner, TV& max_corner)
{
    min_corner.setConstant(1e6);
    max_corner.setConstant(-1e6);
    int num_nodes = deformed.rows() / 3;
    for (int i = 0; i < num_nodes; i++)
    {
        for (int d = 0; d < 3; d++)
        {
            max_corner[d] = std::max(max_corner[d], undeformed[i * 3 + d]);
            min_corner[d] = std::min(min_corner[d], undeformed[i * 3 + d]);
        }
    }
}

void IMLSContact::addShellGravitionEnergy(T& energy)
{
    iterateFaceSerial([&](int face_idx)
    {
        FaceVtx vertices = getFaceVtxDeformed(face_idx);
        FaceVtx undeformed_vertices = getFaceVtxUndeformed(face_idx);

        FaceIdx indices = faces.segment<3>(face_idx * 3);

        TV x0 = vertices.row(0); TV x1 = vertices.row(1); TV x2 = vertices.row(2);
        TV X0 = undeformed_vertices.row(0); 
        TV X1 = undeformed_vertices.row(1); 
        TV X2 = undeformed_vertices.row(2);

        energy += compute3DCSTGravitationalEnergy(density, thickness, 
            gravity, x0, x1, x2, X0, X1, X2);

        
    });
}

void IMLSContact::addShellGravitionForceEntry(VectorXT& residual)
{
    iterateFaceSerial([&](int face_idx)
    {
        FaceVtx vertices = getFaceVtxDeformed(face_idx);
        FaceVtx undeformed_vertices = getFaceVtxUndeformed(face_idx);

        FaceIdx indices = faces.segment<3>(face_idx * 3);

        TV x0 = vertices.row(0); TV x1 = vertices.row(1); TV x2 = vertices.row(2);
        TV X0 = undeformed_vertices.row(0); 
        TV X1 = undeformed_vertices.row(1); 
        TV X2 = undeformed_vertices.row(2);

        Vector9a dedx;
        compute3DCSTGravitationalEnergyGradient(density, thickness, gravity, x0, x1, x2, X0, X1, X2, dedx);

        addForceEntry<3>(residual, {indices[0], indices[1], indices[2]}, -dedx);
    });
}

void IMLSContact::addShellGravitionHessianEntry(std::vector<Entry>& entries)
{
    iterateFaceSerial([&](int face_idx)
    {
        FaceVtx vertices = getFaceVtxDeformed(face_idx);
        FaceVtx undeformed_vertices = getFaceVtxUndeformed(face_idx);

        FaceIdx indices = faces.segment<3>(face_idx * 3);

        TV x0 = vertices.row(0); TV x1 = vertices.row(1); TV x2 = vertices.row(2);
        TV X0 = undeformed_vertices.row(0); 
        TV X1 = undeformed_vertices.row(1); 
        TV X2 = undeformed_vertices.row(2);

        
        Matrix<T, 9, 9> hessian;
        compute3DCSTGravitationalEnergyHessian(density, thickness, gravity, x0, x1, x2, X0, X1, X2, hessian);
    });
}

// ============================= DERIVATIVE TESTS =============================

void IMLSContact::checkTotalGradient(bool perturb)
{
    run_diff_test = true;

    int n_dof = deformed.rows();

    if (perturb)
    {
        VectorXT du(n_dof);
        du.setRandom();
        du *= 1.0 / du.norm();
        du *= 0.001;
        u += du;
    }

    std::cout << "======================== CHECK GRADIENT ========================" << std::endl;
    std::cout << "****** Only mismatching entries are printed ******" << std::endl;
    
    T epsilon = 1e-6;
    VectorXT gradient(n_dof);
    gradient.setZero();

    computeResidual(gradient);

    // std::cout << gradient.transpose() << std::endl;
    
    VectorXT gradient_FD(n_dof);
    gradient_FD.setZero();

    int cnt = 0;
    for(int dof_i = 0; dof_i < n_dof; dof_i++)
    {
        u(dof_i) += epsilon;
        // std::cout << W * dq << std::endl;
        T E0 = computeTotalEnergy();
        
        u(dof_i) -= 2.0 * epsilon;
        T E1 = computeTotalEnergy();
        u(dof_i) += epsilon;
        // std::cout << "E1 " << E1 << " E0 " << E0 << std::endl;
        gradient_FD(dof_i) = (E1 - E0) / (2.0 *epsilon);
        if( gradient_FD(dof_i) == 0 && gradient(dof_i) == 0)
            continue;
        if (std::abs( gradient_FD(dof_i) - gradient(dof_i)) < 1e-3 * std::abs(gradient(dof_i)))
            continue;
        std::cout << " dof " << dof_i << " " << gradient_FD(dof_i) << " " << gradient(dof_i) << std::endl;
        std::getchar();
        cnt++;   
    }
    std::cout << "Gradient Diff Test Finished" << std::endl;
    run_diff_test = false;
}

void IMLSContact::checkTotalGradientScale(bool perturb)
{
    
    run_diff_test = true;
    
    std::cout << "===================== Check Gradient Scale =====================" << std::endl;
    std::cout << "********************You Should Be Seeing 4s********************" << std::endl;
    
    int n_dof = deformed.rows();

    if (perturb)
    {
        VectorXT du(n_dof);
        du.setRandom();
        du *= 1.0 / du.norm();
        du *= 0.001;
        u += du;
    }

    VectorXT gradient(n_dof);
    gradient.setZero();
    computeResidual(gradient);
    gradient *= -1;
    T E0 = computeTotalEnergy();
    VectorXT dx(n_dof);
    dx.setRandom();
    dx *= 1.0 / dx.norm();
    dx *= 0.01;
    T previous = 0.0;
    VectorXT u_backup = u;
    for (int i = 0; i < 10; i++)
    {
        u = u_backup + dx;
        T E1 = computeTotalEnergy();
        T dE = E1 - E0;
        
        dE -= gradient.dot(dx);
        // std::cout << "dE " << dE << std::endl;
        if (i > 0)
        {
            std::cout << (previous/dE) << std::endl;
        }
        previous = dE;
        dx *= 0.5;
    }
    run_diff_test = false;
}


void IMLSContact::checkTotalHessian(bool perturb)
{
    std::cout << "======================== CHECK Hessian ========================" << std::endl;
    std::cout << "****** Only mismatching entries are printed ******" << std::endl;
    run_diff_test = true;
    T epsilon = 1e-5;
    int n_dof = deformed.rows();

    if (perturb)
    {
        VectorXT du(n_dof);
        du.setRandom();
        du *= 1.0 / du.norm();
        du *= 0.001;
        u += du;
    }
    StiffnessMatrix A(n_dof, n_dof);
    buildSystemMatrix(A);
    
    for(int dof_i = 0; dof_i < n_dof; dof_i++)
    {
        // std::cout << dof_i << std::endl;
        u(dof_i) += epsilon;
        VectorXT g0(n_dof), g1(n_dof);
        g0.setZero(); g1.setZero();
        
        computeResidual(g0); 
        
        u(dof_i) -= 2.0 * epsilon;
        
        computeResidual(g1); 
        
        u(dof_i) += epsilon;
        VectorXT row_FD = (g1 - g0) / (2.0 * epsilon);

        for(int i = 0; i < n_dof; i++)
        {
            
            if(A.coeff(i, dof_i) == 0 && row_FD(i) == 0)
                continue;
            if (std::abs( A.coeff(i, dof_i) - row_FD(i)) < 1e-3 * std::abs(row_FD(i)))
                continue;
            // std::cout << "node i: "  << std::floor(dof_i / T(dof)) << " dof " << dof_i%dof 
            //     << " node j: " << std::floor(i / T(dof)) << " dof " << i%dof 
            //     << " FD: " <<  row_FD(i) << " symbolic: " << A.coeff(i, dof_i) << std::endl;
            std::cout << "H(" << i << ", " << dof_i << ") " << " FD: " <<  row_FD(i) << " symbolic: " << A.coeff(i, dof_i) << std::endl;
            std::getchar();
        }
    }
    std::cout << "Hessian Diff Test Finished" << std::endl;
    run_diff_test = false;
}


void IMLSContact::checkTotalHessianScale(bool perturb)
{
    
    std::cout << "===================== check Hessian Scale =====================" << std::endl;
    std::cout << "********************You Should Be Seeing 4s********************" << std::endl;
    std::cout << "************** Unless your function is quadratic **************" << std::endl;
    run_diff_test = true;
    int n_dof = deformed.rows();

    if (perturb)
    {
        VectorXT du(n_dof);
        du.setRandom();
        du *= 1.0 / du.norm();
        du *= 0.001;
        u += du;
    }
    
    StiffnessMatrix A(n_dof, n_dof);
    
    buildSystemMatrix(A);

    VectorXT f0(n_dof);
    f0.setZero();
    computeResidual(f0);
    f0 *= -1;
    
    VectorXT dx(n_dof);
    dx.setRandom();
    dx *= 1.0 / dx.norm();
    for(int i = 0; i < n_dof; i++) dx[i] += 0.5;
    dx *= 0.001;
    T previous = 0.0;
    VectorXT u_backup = u;
    for (int i = 0; i < 10; i++)
    {
        VectorXT f1(n_dof);
        f1.setZero();
        u = u_backup + dx;
        computeResidual(f1);
        f1 *= -1;
        T df_norm = (f0 + (A * dx) - f1).norm();
        // std::cout << "df_norm " << df_norm << std::endl;
        if (i > 0)
        {
            std::cout << (previous/df_norm) << std::endl;
        }
        previous = df_norm;
        dx *= 0.5;
    }
    run_diff_test = false;
}
// ============================= DERIVATIVE TESTS END =============================


// ============================= Dynamics =============================
void IMLSContact::addInertialEnergy(T& energy)
{
    T kinetic_energy = 0.0;
    if (use_consistent_mass_matrix)
    {
        iterateFaceSerial([&](int face_idx)
        {
            FaceVtx vertices = getFaceVtxDeformed(face_idx);
            FaceVtx undeformed_vertices = getFaceVtxUndeformed(face_idx);
            FaceIdx indices = faces.segment<3>(face_idx * 3);
            Matrix<T, 9, 9> mass_mat;
            computeConsistentMassMatrix(undeformed_vertices, mass_mat);
            Vector9a x_n_plus_1_vec, xn_vec, vn_vec;
            for (int i = 0; i < 3; i++)
            {
                x_n_plus_1_vec.segment<3>(i * 3) = vertices.row(i);
                xn_vec.segment<3>(i * 3) = xn.segment<3>(indices[i] * 3);
                vn_vec.segment<3>(i * 3) = vn.segment<3>(indices[i] * 3);
            }

            T xTMx = x_n_plus_1_vec.transpose() * mass_mat * x_n_plus_1_vec;
            T xTMxn_vn_dt = 2.0 * x_n_plus_1_vec.transpose() * mass_mat * (xn_vec + vn_vec * dt);
            kinetic_energy += (xTMx - xTMxn_vn_dt) / (2.0 * std::pow(dt, 2));
            
        });
    }
    else
    {
        for (int i = 0; i < deformed.rows() / 3; i++)
        {
            TV x_n_plus_1 = deformed.segment<3>(i * 3);
            kinetic_energy += (density * mass_diagonal[i] * (x_n_plus_1.dot(x_n_plus_1)
                                                    - 2.0 * x_n_plus_1.dot(xn.segment<3>(i * 3) + vn.segment<3>(i * 3) * dt)
                                                    )) / (2.0 * std::pow(dt, 2));
        }
    }
    
    energy += kinetic_energy;
}

void IMLSContact::addInertialForceEntry(VectorXT& residual)
{
    if (use_consistent_mass_matrix)
    {
        iterateFaceSerial([&](int face_idx)
        {
            FaceVtx vertices = getFaceVtxDeformed(face_idx);
            FaceVtx undeformed_vertices = getFaceVtxUndeformed(face_idx);
            FaceIdx indices = faces.segment<3>(face_idx * 3);
            Matrix<T, 9, 9> mass_mat;
            computeConsistentMassMatrix(undeformed_vertices, mass_mat);
            Vector9a x_n_plus_1_vec, xn_vec, vn_vec;
            for (int i = 0; i < 3; i++)
            {
                x_n_plus_1_vec.segment<3>(i * 3) = vertices.row(i);
                xn_vec.segment<3>(i * 3) = xn.segment<3>(indices[i] * 3);
                vn_vec.segment<3>(i * 3) = vn.segment<3>(indices[i] * 3);
            }
            Vector9a dedx = mass_mat * (2.0 * x_n_plus_1_vec - 2.0 * (xn_vec + vn_vec * dt)) / (2.0 * std::pow(dt, 2));
            addForceEntry<3>(residual, {indices[0], indices[1], indices[2]}, -dedx);
        });
    }
    else
    {
        for (int i = 0; i < deformed.rows() / 3; i++)
        {
            TV x_n_plus_1 = deformed.segment<3>(i * 3);
            residual.segment<3>(i * 3) -= (density * mass_diagonal[i] * (2.0 * x_n_plus_1
                                                    - 2.0 * (xn.segment<3>(i * 3) + vn.segment<3>(i * 3) * dt)
                                                    )) / (2.0 * std::pow(dt, 2));
        }
    }
}

void IMLSContact::addInertialHessianEntries(std::vector<Entry>& entries)
{
    if (use_consistent_mass_matrix)
    {
        iterateFaceSerial([&](int face_idx)
        {
            FaceVtx vertices = getFaceVtxDeformed(face_idx);
            FaceVtx undeformed_vertices = getFaceVtxUndeformed(face_idx);
            Matrix<T, 9, 9> mass_mat;
            computeConsistentMassMatrix(undeformed_vertices, mass_mat);
            FaceIdx indices = faces.segment<3>(face_idx * 3);
            addHessianEntry<3, 3>(entries, {indices[0], indices[1], indices[2]}, mass_mat / std::pow(dt, 2));
        });
    }
    else
    {
        for (int i = 0; i < deformed.rows() / 3; i++)
        {
            TM hess = density * TM::Identity() * mass_diagonal[i] / std::pow(dt, 2);
            addHessianEntry<3, 3>(entries, {i}, hess);
        }
    }
}

void IMLSContact::updateDynamicStates()
{
    vn = (deformed - xn) / dt;
    xn = deformed;
}

void IMLSContact::initializeDynamicStates()
{
    if (!use_consistent_mass_matrix)
        computeMassMatrix();
    xn = undeformed;
    vn = VectorXT::Zero(undeformed.rows());
}

void IMLSContact::computeMassMatrix()
{
    MatrixXT V; MatrixXi F;
    vectorToIGLMatrix<T, 3>(undeformed, V);
    vectorToIGLMatrix<int, 3>(faces, F);
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, M);
    mass_diagonal.resize(deformed.rows() / 3);
    mass_diagonal.setZero();
    mass_diagonal = M.diagonal();   
}

void IMLSContact::computeConsistentMassMatrix(const FaceVtx& vtx_pos, Matrix<T, 9, 9>& mass_mat)
{
    double m[81];
    T t1 = vtx_pos(0, 0) - vtx_pos(2, 0);
    T t2 = vtx_pos(0, 1) - vtx_pos(2, 1);
    T t3 = vtx_pos(0, 2) - vtx_pos(2, 2);
    T t4 = vtx_pos(1, 0) - vtx_pos(2, 0);
    T t5 = vtx_pos(1, 1) - vtx_pos(2, 1);
    T t6 = vtx_pos(1, 2) - vtx_pos(2, 2);
    T t7 = t1 * t4 + t2 * t5 + t3 * t6;
    t1 = (pow(t1, 0.2e1) + pow(t2, 0.2e1) + pow(t3, 0.2e1)) * (pow(t4, 0.2e1) + pow(t5, 0.2e1) + pow(t6, 0.2e1)) - pow(t7, 0.2e1);
    t1 = sqrt(t1);
    t2 = t1 / 0.12e2;
    t1 = t1 / 0.24e2;
    m[0] = t2;
    m[1] = 0;
    m[2] = 0;
    m[3] = t1;
    m[4] = 0;
    m[5] = 0;
    m[6] = t1;
    m[7] = 0;
    m[8] = 0;
    m[9] = 0;
    m[10] = t2;
    m[11] = 0;
    m[12] = 0;
    m[13] = t1;
    m[14] = 0;
    m[15] = 0;
    m[16] = t1;
    m[17] = 0;
    m[18] = 0;
    m[19] = 0;
    m[20] = t2;
    m[21] = 0;
    m[22] = 0;
    m[23] = t1;
    m[24] = 0;
    m[25] = 0;
    m[26] = t1;
    m[27] = t1;
    m[28] = 0;
    m[29] = 0;
    m[30] = t2;
    m[31] = 0;
    m[32] = 0;
    m[33] = t1;
    m[34] = 0;
    m[35] = 0;
    m[36] = 0;
    m[37] = t1;
    m[38] = 0;
    m[39] = 0;
    m[40] = t2;
    m[41] = 0;
    m[42] = 0;
    m[43] = t1;
    m[44] = 0;
    m[45] = 0;
    m[46] = 0;
    m[47] = t1;
    m[48] = 0;
    m[49] = 0;
    m[50] = t2;
    m[51] = 0;
    m[52] = 0;
    m[53] = t1;
    m[54] = t1;
    m[55] = 0;
    m[56] = 0;
    m[57] = t1;
    m[58] = 0;
    m[59] = 0;
    m[60] = t2;
    m[61] = 0;
    m[62] = 0;
    m[63] = 0;
    m[64] = t1;
    m[65] = 0;
    m[66] = 0;
    m[67] = t1;
    m[68] = 0;
    m[69] = 0;
    m[70] = t2;
    m[71] = 0;
    m[72] = 0;
    m[73] = 0;
    m[74] = t1;
    m[75] = 0;
    m[76] = 0;
    m[77] = t1;
    m[78] = 0;
    m[79] = 0;
    m[80] = t2;

    for (int i = 0; i < 9; i++)
    {
        for (int j = 0; j < 9; j++)
        {
            mass_mat(i, j) = m[i * 9 + j] * density;
        }
        
    }
    
}
// ============================= Dynamics End =============================