#include <Eigen/CholmodSupport>
#include <igl/readOBJ.h>
#include <igl/readMSH.h>
#include <igl/edges.h>
#include <igl/boundary_facets.h>
#include <Spectra/SymEigsShiftSolver.h>
#include <Spectra/MatOp/SparseSymShiftSolve.h>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/SparseSymMatProd.h>

#include "../autodiff/FEMEnergy.h"
#include "../include/FEM3D.h"

bool FEM3D::advanceOneStep(int step)
{
    std::cout << "===================STEP " << step << "===================" << std::endl;
    VectorXT residual;
    residual.resize(deformed.rows());
    residual.setZero();
    T residual_norm = computeResidual(residual);
    residual_norms.push_back(residual_norm);
    std::cout << "[NEWTON] iter " << step << "/" 
        << max_newton_iter << ": residual_norm " 
        << residual_norm << " tol: " << newton_tol << std::endl;
    updateSurfaceVertices();
    if (residual_norm < newton_tol || step == max_newton_iter)
    {
        return true;
    }

    T du_norm = 1e10;
    du_norm = lineSearchNewton(residual);
    return false;
}

bool FEM3D::linearSolve(StiffnessMatrix& K, const VectorXT& residual, VectorXT& du)
{
    START_TIMING(LinearSolve)
    Eigen::CholmodSupernodalLLT<StiffnessMatrix, Eigen::Lower> solver;
    // Eigen::CholmodSupernodalLLT<StiffnessMatrix> solver;
    
    T alpha = 1e-6;
    StiffnessMatrix H(K.rows(), K.cols());
    H.setIdentity(); H.diagonal().array() = 1e-10;
    K += H;
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

T FEM3D::computeTotalEnergy()
{
    deformed = undeformed + u;

    T energy = 0.0;
    energy += addElastsicPotential();
    return energy;
}

T FEM3D::computeResidual(VectorXT& residual)
{
    deformed = undeformed + u;
    addElasticForceEntries(residual);
    if (!run_diff_test)
        iterateDirichletDoF([&](int offset, T target)
        {
            residual[offset] = 0;
        });
    return residual.norm();
}

void FEM3D::computeLinearModes(MatrixXT& eigen_vectors, VectorXT& eigen_values)
{
    int nmodes = 10;
    StiffnessMatrix K;
    run_diff_test = true;
    buildSystemMatrix(K);
    run_diff_test = false;
    Spectra::SparseSymShiftSolve<T, Eigen::Lower> op(K);
    
    T shift = -1e-4;
    Spectra::SymEigsShiftSolver<T, 
    Spectra::LARGEST_MAGN, 
    Spectra::SparseSymShiftSolve<T, Eigen::Lower> > 
        eigs(&op, nmodes, 2 * nmodes, shift);

    eigs.init();

    int nconv = eigs.compute();

    if (eigs.info() == Spectra::SUCCESSFUL)
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

void FEM3D::buildSystemMatrix(StiffnessMatrix& K)
{
    deformed = undeformed + u;
    std::vector<Entry> entries;
    
    addElasticHessianEntries(entries);
    
    int n_dof = deformed.rows();
    K.resize(n_dof, n_dof);
    K.setFromTriplets(entries.begin(), entries.end());
    if (!run_diff_test)
        projectDirichletDoFMatrix(K, dirichlet_data);
    
    K.makeCompressed();
}

void FEM3D::projectDirichletDoFMatrix(StiffnessMatrix& A, const std::unordered_map<int, T>& data)
{
    for (auto iter : data)
    {
        A.row(iter.first) *= 0.0;
        A.col(iter.first) *= 0.0;
        A.coeffRef(iter.first, iter.first) = 1.0;
    }
}

T FEM3D::lineSearchNewton(const VectorXT& residual)
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
        std::cout << "\t|du | " << norm << std::endl;
    
    T E0 = computeTotalEnergy();
    // std::cout << "obj: " << E0 << std::endl;
    T alpha = 1.0;
    int cnt = 0;
    VectorXT u_current = u;
    while (true)
    {
        u = u_current + alpha * du;
        T E1 = computeTotalEnergy();
        if (E1 - E0 < 0 || cnt > 10)
        {
            if (cnt > 10)
                std::cout << "cnt > 10" << std::endl;
            break;
        }
        alpha *= 0.5;
        cnt += 1;
    }

    return du.norm();
}

void FEM3D::initializeFromFile(const std::string& filename)
{
    
    MatrixXT V, Vtet; MatrixXi F, Ftet, Ttet; VectorXi tri_flag, tet_flag;
    
    igl::readMSH(filename, Vtet, Ftet, Ttet, tri_flag, tet_flag);

    MatrixXi surface_F;
    igl::boundary_facets(Ttet, surface_F);

    for (int i = 0; i < surface_F.rows(); i++)
    {
        std::swap(surface_F(i, 1), surface_F(i, 2));
    }

    std::unordered_set<int> unique_vtx;
    for (int i = 0; i < surface_F.rows(); i++)
        for (int d = 0; d < 3; d++)
        {
            unique_vtx.insert(surface_F(i, d));
        }
    int n_vtx = unique_vtx.size();
    // std::cout << n_vtx << std::endl;
    
    int cnt = 0;
    for (auto idx : unique_vtx)
    {
        surface_to_tet_node_map[cnt] = idx;
        tet_node_surface_map[idx] = cnt++;

    }
    F.resize(surface_F.rows(), 3);
    for (int i = 0; i < F.rows(); i++)
        for (int d = 0; d < 3; d++)
        {
            F(i, d) = tet_node_surface_map[surface_F(i, d)];
        }
    V.resize(n_vtx, 3);

    cnt = 0;
    for (auto idx : unique_vtx)
        V.row(cnt++) = Vtet.row(idx);
    // igl::writeOBJ("./surface_buggy.obj", V, F);

    TV min_corner = V.colwise().minCoeff();
    TV max_corner = V.colwise().maxCoeff();
    TV center = 0.5 * (min_corner + max_corner);

    T bb_diag = (max_corner - min_corner).norm();
    T scale = 1.0 / bb_diag;
    
    for (int i = 0; i < V.rows(); i++)
    {
        V.row(i) -= center;
    }
    V *= scale;

    for (int i = 0; i < Vtet.rows(); i++)
    {
        Vtet.row(i) -= center;
    }
    Vtet *= scale;

    num_nodes = Vtet.rows();
    undeformed.conservativeResize(num_nodes * 3);
    tbb::parallel_for(0, num_nodes, [&](int i)
    {
        undeformed.segment<3>(i * 3) = Vtet.row(i);
    });
    deformed = undeformed;    
    u = deformed; u.setZero();
    num_ele = Ttet.rows();
    indices.resize(num_ele * 4);
    tbb::parallel_for(0, num_ele, [&](int i)
    {
        indices.segment<4>(i * 4) = Ttet.row(i);
    });

    surface_vertices = V;
    surface_indices = F;
}

void FEM3D::updateSurfaceVertices()
{
    for (const auto& data : surface_to_tet_node_map)
    {
        surface_vertices.row(data.first) = deformed.segment<3>(data.second * 3);
    }
}

void FEM3D::computeBoundingBox(TV& min_corner, TV& max_corner)
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

T FEM3D::computeVolume(const EleNodes& x_undeformed)
{
    TV a = x_undeformed.row(1) - x_undeformed.row(0);
	TV b = x_undeformed.row(2) - x_undeformed.row(0);
	TV c = x_undeformed.row(3) - x_undeformed.row(0);
	T volumeParallelepiped = a.cross(b).dot(c);
	T tetVolume = 1.0 / 6.0 * volumeParallelepiped;
	return tetVolume;
}

T FEM3D::computeInversionFreeStepsize()
{
    Matrix<T, 4, 3> dNdb;
        dNdb << -1.0, -1.0, -1.0, 
            1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, 1.0;
           
    VectorXT step_sizes = VectorXT::Zero(num_ele);

    iterateElementParallel([&](const EleNodes& x_deformed, 
        const EleNodes& x_undeformed, const VtxList& indices, int tet_idx)
    {
        TM dXdb = x_undeformed.transpose() * dNdb;
        TM dxdb = x_deformed.transpose() * dNdb;
        TM A = dxdb * dXdb.inverse();
        T a, b, c, d;
        a = A.determinant();
        b = A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0) + A(0, 0) * A(2, 2) - A(0, 2) * A(2, 0) + A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1);
        c = A.diagonal().sum();
        d = 0.8;

        T t = getSmallestPositiveRealCubicRoot(a, b, c, d);
        if (t < 0 || t > 1) t = 1;
            step_sizes(tet_idx) = t;
    });
    return step_sizes.minCoeff();
}


T FEM3D::addElastsicPotential()
{
    VectorXT energies_neoHookean(num_ele);
    energies_neoHookean.setZero();
    iterateElementSerial([&](const EleNodes& x_deformed, 
        const EleNodes& x_undeformed, const VtxList& indices, int tet_idx)
    {
        T ei = 0.0;
        computeLinearTet3DNeoHookeanEnergy(E, nu, x_deformed, x_undeformed, ei);
        energies_neoHookean[tet_idx] += ei;
        // if (std::isnan(ei))
        // {
        //     std::cout << "NAN " << ei << std::endl;
        // }
    });
    
    return energies_neoHookean.sum();
}

void FEM3D::addElasticForceEntries(VectorXT& residual)
{
    iterateElementSerial([&](const EleNodes& x_deformed, 
        const EleNodes& x_undeformed, const VtxList& indices, int tet_idx)
    {
        // T volume = computeVolume(x_undeformed);
        // if (volume < 1e-8)
        //     return;
        Vector<T, 12> dedx;
        computeLinearTet3DNeoHookeanEnergyGradient(E, nu, x_deformed, x_undeformed, dedx);
        addForceEntry<3>(residual, indices, -dedx);
        // std::cout << dedx.transpose() << std::endl;
    });
}


void FEM3D::addElasticHessianEntries(std::vector<Entry>& entries)
{
    iterateElementSerial([&](const EleNodes& x_deformed, 
        const EleNodes& x_undeformed, const VtxList& indices, int tet_idx)
    {
        T volume = computeVolume(x_undeformed);
        // if (volume < 1e-8)
        //     return;
        Matrix<T, 12, 12> hessian, hessian_ad;
        computeLinearTet3DNeoHookeanEnergyHessian(E, nu, x_deformed, x_undeformed, hessian);
        addHessianEntry<3, 3>(entries, indices, hessian);
    });
}