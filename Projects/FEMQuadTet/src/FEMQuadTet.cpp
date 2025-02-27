#include <igl/readMSH.h>
#include <Eigen/CholmodSupport>
#include "../autodiff/QuadraticTetEnergy.h"
#include "../include/FEMQuadTet.h"

// This function assumes Vtet is an (n x 3) matrix of vertex positions
// and Ttet is an (m x 4) matrix with each row holding indices of the tet's
// vertices. The function produces Vquad (all vertices: original plus midpoints)
// and Tquad (an (m x 10) matrix for the quadratic tets).
void FEMQuadTet::convertToQuadraticTets(const MatrixXT& Vtet,
                            const MatrixXi& Ttet, MatrixXT& Vquad,
                            MatrixXi& Tquad)
{
    // Map a sorted edge (pair of vertex indices) to the new midpoint's index.
    std::map<std::pair<int, int>, int> edgeToMidpoint;
    // Container for the midpoint positions.
    std::vector<Eigen::RowVector3d> midpoints;

    // Number of original vertices (they remain unchanged).
    int numOrigVerts = Vtet.rows();

    // Loop over each tet and process its six edges.
    for (int i = 0; i < Ttet.rows(); ++i)
    {
        // Extract the four vertex indices for this tet.
        int v0 = Ttet(i, 0);
        int v1 = Ttet(i, 1);
        int v2 = Ttet(i, 2);
        int v3 = Ttet(i, 3);

        // List of the six edges; we sort the vertex pair so that each edge is
        // uniquely identified.
        std::pair<int, int> edges[6] = {{std::min(v0, v1), std::max(v0, v1)},
                                        {std::min(v0, v2), std::max(v0, v2)},
                                        {std::min(v0, v3), std::max(v0, v3)},
                                        {std::min(v1, v2), std::max(v1, v2)},
                                        {std::min(v1, v3), std::max(v1, v3)},
                                        {std::min(v2, v3), std::max(v2, v3)}};

        // For each edge, check if the midpoint was already computed.
        for (int j = 0; j < 6; ++j)
        {
            auto edge = edges[j];
            if (edgeToMidpoint.find(edge) == edgeToMidpoint.end())
            {
                // Compute the midpoint position.
                Eigen::RowVector3d midpoint =
                    0.5 * (Vtet.row(edge.first) + Vtet.row(edge.second));
                // The new index is the original vertex count plus the current
                // number of midpoints.
                int newIndex =
                    numOrigVerts + static_cast<int>(midpoints.size());
                edgeToMidpoint[edge] = newIndex;
                midpoints.push_back(midpoint);
            }
        }
    }

    // Create the new vertex matrix by appending midpoints to the original
    // vertices.
    Vquad.resize(numOrigVerts + midpoints.size(), Vtet.cols());
    Vquad.topRows(numOrigVerts) = Vtet;
    for (size_t i = 0; i < midpoints.size(); ++i)
    {
        Vquad.row(numOrigVerts + i) = midpoints[i];
    }

    // Create the new connectivity matrix for quadratic tets.
    // Each row now has 10 nodes: the first 4 are the original vertices,
    // and the next 6 are the midpoints for edges (v0,v1), (v0,v2), (v0,v3),
    // (v1,v2), (v1,v3), (v2,v3) respectively.
    Tquad.resize(Ttet.rows(), 10);
    for (int i = 0; i < Ttet.rows(); ++i)
    {
        int v0 = Ttet(i, 0);
        int v1 = Ttet(i, 1);
        int v2 = Ttet(i, 2);
        int v3 = Ttet(i, 3);

        int m01 = edgeToMidpoint[{std::min(v0, v1), std::max(v0, v1)}];
        int m02 = edgeToMidpoint[{std::min(v0, v2), std::max(v0, v2)}];
        int m03 = edgeToMidpoint[{std::min(v0, v3), std::max(v0, v3)}];
        int m12 = edgeToMidpoint[{std::min(v1, v2), std::max(v1, v2)}];
        int m13 = edgeToMidpoint[{std::min(v1, v3), std::max(v1, v3)}];
        int m23 = edgeToMidpoint[{std::min(v2, v3), std::max(v2, v3)}];

        // Fill the new connectivity row.
        // Ordering: see figure 10.2
        //https://web.archive.org/web/20160303173558/http://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/AFEM.Ch10.d/AFEM.Ch10.pdf
        Tquad.row(i) << v0, v1, v2, v3, m01, m12, m02, m03, m13, m23;
    }
}

void FEMQuadTet::computeBoundingBox(TV& min_corner, TV& max_corner)
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


void FEMQuadTet::initializeFromFile(const std::string& filename) 
{

    MatrixXT V, Vtet; MatrixXi F, Ftet, Ttet; VectorXi tri_flag, tet_flag;
    
    igl::readMSH(filename, Vtet, Ftet, Ttet, tri_flag, tet_flag);

    n_linear_nodes = Vtet.rows();
    MatrixXT Vquad;
    MatrixXi Tquad;
    convertToQuadraticTets(Vtet, Ttet, Vquad, Tquad);

    linear_tet_indices = Ttet;
    // std::cout << "Vquad: " << Vquad.rows() << " " << Vquad.cols() << std::endl;
    // std::cout << "Tquad: " << Tquad.rows() << " " << Tquad.cols() << std::endl;


    num_nodes = Vquad.rows();
    undeformed.conservativeResize(num_nodes * 3);
    tbb::parallel_for(0, num_nodes, [&](int i)
    {
        undeformed.segment<3>(i * 3) = Vquad.row(i);
    });
    deformed = undeformed;    
    u = deformed; u.setZero();
    f = u; 
    for (int i = 0; i < num_nodes; i++)
    {
        f[i*3+1] = -9.8;
    }
    
    num_ele = Tquad.rows();
    indices.resize(num_ele * 10);
    tbb::parallel_for(0, num_ele, [&](int i)
    {
        indices.segment<10>(i * 10) = Tquad.row(i);
    });
    TV min_corner, max_corner;
    computeBoundingBox(min_corner, max_corner);
    std::cout << "BBOX " << min_corner.transpose() << " " << max_corner.transpose() << std::endl;
    
    for (int i = 0; i < num_nodes; i++)
    {
        TV xi = deformed.segment<3>(i * 3);
        if (xi[2] < min_corner[2] + 0.001 * (max_corner[2] - min_corner[2]))
        {
            dirichlet_data[i * 3] = 0.0;
            dirichlet_data[i * 3 + 1] = 0.0;
            dirichlet_data[i * 3 + 2] = 0.0;
        }
    }

    max_newton_iter = 300000;

    E = 1e8;

}

bool FEMQuadTet::advanceOneStep(int step)
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
    
    if (residual_norm < newton_tol || step == max_newton_iter)
    {
        return true;
    }

    T du_norm = 1e10;
    
    du_norm = lineSearchNewton(residual);
    return false;
}

bool FEMQuadTet::linearSolve(StiffnessMatrix& K, const VectorXT& residual, VectorXT& du)
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

T FEMQuadTet::computeTotalEnergy()
{
    deformed = undeformed + u;

    T energy = 0.0;
    energy += addElastsicPotential();
    
    energy -= f.dot(u);
    return energy;
}

T FEMQuadTet::computeResidual(VectorXT& residual)
{
    deformed = undeformed + u;
    addElasticForceEntries(residual);
    
    residual += f;
    if (!run_diff_test)
        iterateDirichletDoF([&](int offset, T target)
        {
            residual[offset] = 0;
        });
    return residual.norm();
}


void FEMQuadTet::buildSystemMatrix(StiffnessMatrix& K)
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

void FEMQuadTet::projectDirichletDoFMatrix(StiffnessMatrix& A, const std::unordered_map<int, T>& data)
{
    for (auto iter : data)
    {
        A.row(iter.first) *= 0.0;
        A.col(iter.first) *= 0.0;
        A.coeffRef(iter.first, iter.first) = 1.0;
    }
}

T FEMQuadTet::lineSearchNewton(const VectorXT& residual)
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
        // std::cout << "E0: " << E0 << " E1 " << E1 << std::endl;
        if (E1 - E0 < 0 || cnt > 10)
        {
            // if (cnt > 10)
            //     std::cout << "cnt > 10" << std::endl;
            break;
        }
        alpha *= 0.5;
        cnt += 1;
    }

    return du.norm();
}


T FEMQuadTet::addElastsicPotential()
{
    VectorXT energies_neoHookean(num_ele);
    energies_neoHookean.setZero();
    iterateElementSerial([&](const EleNodes& x_deformed, 
        const EleNodes& x_undeformed, const VtxList& indices, int tet_idx)
    {
        T ei = 0.0;
        compute3DNeoHookeanQuadraticTetEnergy(E, nu, x_deformed, x_undeformed, ei);
        energies_neoHookean[tet_idx] += ei;
    });
    
    return energies_neoHookean.sum();
}

T FEMQuadTet::elasticPotentialAnalytical()
{
    T lambda = E * nu / (T(1.0) + nu) / (T(1.0) - T(2.0) * nu);
	T mu = E / T(2.0) / (T(1.0) + nu);

    VectorXT energies_neoHookean(num_ele);
    energies_neoHookean.setZero();
    iterateElementSerial([&](const EleNodes& x_deformed, 
        const EleNodes& x_undeformed, const VtxList& indices, int tet_idx)
    {
        T ei = 0.0;
        for (int quadPointIdx = 0; quadPointIdx < 4; quadPointIdx++)
        {
            Matrix<T, 10, 3> dNdb;
            computeTet10BasisJacobian(gauss3DPP2Rule4(quadPointIdx), dNdb);
            TM dXdb = x_undeformed.transpose() * dNdb;
            TM dxdb = x_deformed.transpose() * dNdb;
            TM def_grad = dxdb * dXdb.inverse();
            TM C = def_grad.transpose() * def_grad;

            T J = def_grad.determinant();
            T lnJ = log(J);

            T I1 = C.trace();

            T energy_density = 0.5 * mu * (I1 - 3.0 - 2.0 * lnJ) + lambda * 0.5 * (lnJ*lnJ);

            TM transformation_Jacobian = x_undeformed.transpose()  * dNdb;
		    T det = transformation_Jacobian.determinant();
            T referenceTetVol = 1.0 / 6.0;
            T vol = referenceTetVol * det;
            ei += gauss3DPP2Rule4Weight(quadPointIdx) * energy_density * vol;
        }
        energies_neoHookean[tet_idx] += ei;
    });
    
    return energies_neoHookean.sum();
}

void FEMQuadTet::addElasticForceEntries(VectorXT& residual)
{
    iterateElementSerial([&](const EleNodes& x_deformed, 
        const EleNodes& x_undeformed, const VtxList& indices, int tet_idx)
    {
        
        Vector<T, 30> dedx;
        compute3DNeoHookeanQuadraticTetEnergyGradient(E, nu, x_deformed, x_undeformed, dedx);
        addForceEntry<3>(residual, indices, -dedx);
    });
}


void FEMQuadTet::addElasticHessianEntries(std::vector<Entry>& entries)
{
    iterateElementSerial([&](const EleNodes& x_deformed, 
        const EleNodes& x_undeformed, const VtxList& indices, int tet_idx)
    {
        Matrix<T, 30, 30> hessian;
        compute3DNeoHookeanQuadraticTetEnergyHessian(E, nu, x_deformed, x_undeformed, hessian);
        addHessianEntry<3, 3>(entries, indices, hessian);
    });
}