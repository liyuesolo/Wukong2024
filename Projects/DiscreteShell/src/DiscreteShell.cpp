#include <Eigen/CholmodSupport>
#include <igl/readOBJ.h>
#include <Spectra/SymEigsShiftSolver.h>
#include <Spectra/MatOp/SparseSymShiftSolve.h>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/SparseSymMatProd.h>

#include "../autodiff/CST3DShell.h"
#include "../include/DiscreteShell.h"

bool DiscreteShell::advanceOneStep(int step)
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

bool DiscreteShell::linearSolve(StiffnessMatrix& K, const VectorXT& residual, VectorXT& du)
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

T DiscreteShell::computeTotalEnergy()
{
    deformed = undeformed + u;

    T energy = 0.0;
    addShellEnergy(energy);
    return energy;
}

T DiscreteShell::computeResidual(VectorXT& residual)
{
    deformed = undeformed + u;
    addShellForceEntry(residual);
    if (!run_diff_test)
        iterateDirichletDoF([&](int offset, T target)
        {
            residual[offset] = 0;
        });
    return residual.norm();
}

void DiscreteShell::computeLinearModes(MatrixXT& eigen_vectors, VectorXT& eigen_values)
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

void DiscreteShell::buildSystemMatrix(StiffnessMatrix& K)
{
    deformed = undeformed + u;
    std::vector<Entry> entries;
    addShellHessianEntries(entries);
    
    int n_dof = deformed.rows();
    K.resize(n_dof, n_dof);
    K.setFromTriplets(entries.begin(), entries.end());
    if (!run_diff_test)
        projectDirichletDoFMatrix(K, dirichlet_data);
    
    K.makeCompressed();
}

void DiscreteShell::projectDirichletDoFMatrix(StiffnessMatrix& A, const std::unordered_map<int, T>& data)
{
    for (auto iter : data)
    {
        A.row(iter.first) *= 0.0;
        A.col(iter.first) *= 0.0;
        A.coeffRef(iter.first, iter.first) = 1.0;
    }
}

T DiscreteShell::lineSearchNewton(const VectorXT& residual)
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

void DiscreteShell::initializeFromFile(const std::string& filename)
{
    MatrixXT V; MatrixXi F;
    igl::readOBJ(filename, V, F);
    iglMatrixFatten<T, 3>(V, undeformed);
    iglMatrixFatten<int, 3>(F, faces);
    deformed = undeformed;
    u = VectorXT::Zero(deformed.rows());
    computeRestShape();
    buildHingeStructure();
    updateLameParameters();
}



void DiscreteShell::buildHingeStructure()
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

void DiscreteShell::computeRestShape()
{
    Xinv.resize(nFaces());
    shell_rest_area = VectorXT::Ones(nFaces());
    thickness = VectorXT::Ones(nFaces()) * 0.003;
    
    iterateFaceSerial([&](int i){
        TV E[2]; // edges
        FaceVtx vertices = getFaceVtxUndeformed(i);
        E[0] = vertices.row(1) - vertices.row(0);
        E[1] = vertices.row(2) - vertices.row(0);
        
        TV N = E[0].cross(E[1]);
        T m_A0 = N.norm()*0.5;
        shell_rest_area[i] = m_A0;
        
        // // compute rest
        TV B2D[2];
        B2D[0] = E[0].normalized();
        TV n = B2D[0].cross(E[1]);
        B2D[1] = E[0].cross(n);
        B2D[1] = B2D[1].normalized();
        
        
        Eigen::Matrix2d MatE2D(2, 2);
        MatE2D(0, 0) = E[0].dot(B2D[0]);
        MatE2D(1, 0) = E[0].dot(B2D[1]);
        MatE2D(0, 1) = E[1].dot(B2D[0]);
        MatE2D(1, 1) = E[1].dot(B2D[1]);

        Eigen::Matrix2d Einv = MatE2D.inverse();
        Xinv[i] = Einv;
        // std::cout << Einv << std::endl;
        // std::cout << N.transpose() << std::endl;
        // std::getchar();
    });
}

void DiscreteShell::addShellInplaneEnergy(T& energy)
{
    iterateFaceSerial([&](int face_idx)
    {
        FaceVtx vertices = getFaceVtxDeformed(face_idx);
        FaceVtx undeformed_vertices = getFaceVtxUndeformed(face_idx);

        FaceIdx indices = faces.segment<3>(face_idx * 3);
        

        T m_Einv[2][2];
        m_Einv[0][0] = Xinv[face_idx](0, 0);
        m_Einv[0][1] = Xinv[face_idx](0, 1);
        m_Einv[1][0] = Xinv[face_idx](1, 0);
        m_Einv[1][1] = Xinv[face_idx](1, 1);

        T m_A0 = shell_rest_area[face_idx];
        T m_h = thickness[face_idx];

        std::array<TV, 3> x;
        x[0] = vertices.row(0);
        x[1] = vertices.row(1);
        x[2] = vertices.row(2);
        
        TV x0 = vertices.row(0); TV x1 = vertices.row(1); TV x2 = vertices.row(2);
        TV X0 = undeformed_vertices.row(0); 
        TV X1 = undeformed_vertices.row(1); 
        TV X2 = undeformed_vertices.row(2);

        T k_s = E * m_h / (1.0 - nu * nu);

        energy += compute3DCSTShellEnergy(nu, k_s, x0, x1, x2, X0, X1, X2);

    });
}

void DiscreteShell::addShellBendingEnergy(T& energy)
{
    iterateHingeSerial([&](const HingeIdx& hinge_idx, int hinge_cnt){
        
        HingeVtx deformed_vertices = getHingeVtxDeformed(hinge_idx);
        HingeVtx undeformed_vertices = getHingeVtxUndeformed(hinge_idx);

        std::array<TV, 4> x, X;
        x[0] = deformed_vertices.row(0);
        x[1] = deformed_vertices.row(1);
        x[2] = deformed_vertices.row(2);
        x[3] = deformed_vertices.row(3);

        X[0] = undeformed_vertices.row(0);
        X[1] = undeformed_vertices.row(1);
        X[2] = undeformed_vertices.row(2);
        X[3] = undeformed_vertices.row(3);

        T k_bend = hinge_stiffness[hinge_cnt] * E * std::pow(thickness[0], 3) / (24 * (1.0 - std::pow(nu, 2)));

        TV x0 = deformed_vertices.row(0); TV x1 = deformed_vertices.row(1); TV x2 = deformed_vertices.row(2); TV x3 = deformed_vertices.row(3);
        TV X0 = undeformed_vertices.row(0); TV X1 = undeformed_vertices.row(1); TV X2 = undeformed_vertices.row(2); TV X3 = undeformed_vertices.row(3);

        energy += computeDSBendingEnergy(k_bend, x0, x1, x2, x3, X0, X1, X2, X3);

    });
}

void DiscreteShell::addShellEnergy(T& energy)
{
    T in_plane_energy = 0.0;
    T bending_energy = 0.0;
    addShellInplaneEnergy(in_plane_energy);
    addShellBendingEnergy(bending_energy);    

    energy += in_plane_energy;
    energy += bending_energy;
}

void DiscreteShell::addShellInplaneForceEntry(VectorXT& residual)
{
    iterateFaceSerial([&](int face_idx)
    {
        FaceVtx vertices = getFaceVtxDeformed(face_idx);
        FaceVtx undeformed_vertices = getFaceVtxUndeformed(face_idx);
        FaceIdx indices = faces.segment<3>(face_idx * 3);
        T m_Einv[2][2];
        m_Einv[0][0] = Xinv[face_idx](0, 0);
        m_Einv[0][1] = Xinv[face_idx](0, 1);
        m_Einv[1][0] = Xinv[face_idx](1, 0);
        m_Einv[1][1] = Xinv[face_idx](1, 1);

        T m_A0 = shell_rest_area[face_idx];
        T m_h = thickness[face_idx];

        std::array<TV, 3> x, X;
        x[0] = vertices.row(0);
        x[1] = vertices.row(1);
        x[2] = vertices.row(2);

        X[0] = undeformed_vertices.row(0);
        X[1] = undeformed_vertices.row(1);
        X[2] = undeformed_vertices.row(2);
        

        T k_s = E * m_h / (1.0 - nu * nu);
        TV x0 = vertices.row(0); TV x1 = vertices.row(1); TV x2 = vertices.row(2);
        TV X0 = undeformed_vertices.row(0); TV X1 = undeformed_vertices.row(1); TV X2 = undeformed_vertices.row(2);
        
        Vector<T, 9> dedx;
        compute3DCSTShellEnergyGradient(nu, k_s, x0, x1, x2, X0, X1, X2, dedx);
        
        dedx *= -1.0;
        
        for (int i = 0; i < 3; i++)
        {
            for (int d = 0; d < 3; d++)
            {
                residual[indices[i] * 3 + d] += dedx[i * 3 + d];
            }   
        }
    });
}

void DiscreteShell::addShellBendingForceEntry(VectorXT& residual)
{
    iterateHingeSerial([&](const HingeIdx& hinge_idx, int hinge_cnt){
                
        HingeVtx deformed_vertices = getHingeVtxDeformed(hinge_idx);
        HingeVtx undeformed_vertices = getHingeVtxUndeformed(hinge_idx);

        std::array<TV, 4> x, X;
        x[0] = deformed_vertices.row(0);
        x[1] = deformed_vertices.row(1);
        x[2] = deformed_vertices.row(2);
        x[3] = deformed_vertices.row(3);

        X[0] = undeformed_vertices.row(0);
        X[1] = undeformed_vertices.row(1);
        X[2] = undeformed_vertices.row(2);
        X[3] = undeformed_vertices.row(3);

        T k_bend = hinge_stiffness[hinge_cnt] * E * std::pow(thickness[0], 3) / (24 * (1.0 - std::pow(nu, 2)));

        Vector<T, 12> fx;

        // #include "autodiff/DSBending_fx.mcg"
        
        TV x0 = deformed_vertices.row(0); TV x1 = deformed_vertices.row(1); TV x2 = deformed_vertices.row(2); TV x3 = deformed_vertices.row(3);
        TV X0 = undeformed_vertices.row(0); TV X1 = undeformed_vertices.row(1); TV X2 = undeformed_vertices.row(2); TV X3 = undeformed_vertices.row(3);

        computeDSBendingEnergyGradient(k_bend, x0, x1, x2, x3, X0, X1, X2, X3, fx);
        fx *= -1;

        for (int i = 0; i < 4; i++)
        {
            for (int d = 0; d < 3; d++)
            {
                residual[hinge_idx[i] * 3 + d] += fx[i * 3 + d];
            }
        }
    });
}
void DiscreteShell::addShellForceEntry(VectorXT& residual)
{
    addShellInplaneForceEntry(residual);
    addShellBendingForceEntry(residual);
}

void DiscreteShell::addShellInplaneHessianEntries(std::vector<Entry>& entries)
{
    iterateFaceSerial([&](int face_idx)
    {
        FaceVtx vertices = getFaceVtxDeformed(face_idx);
        FaceVtx undeformed_vertices = getFaceVtxUndeformed(face_idx);

        FaceIdx indices = faces.segment<3>(face_idx * 3);
        T m_Einv[2][2];
        m_Einv[0][0] = Xinv[face_idx](0, 0);
        m_Einv[0][1] = Xinv[face_idx](0, 1);
        m_Einv[1][0] = Xinv[face_idx](1, 0);
        m_Einv[1][1] = Xinv[face_idx](1, 1);

        T m_A0 = shell_rest_area[face_idx];
        T m_h = thickness[face_idx];

        std::array<TV, 3> x, X;
        x[0] = vertices.row(0);
        x[1] = vertices.row(1);
        x[2] = vertices.row(2);

        X[0] = undeformed_vertices.row(0);
        X[1] = undeformed_vertices.row(1);
        X[2] = undeformed_vertices.row(2);

        T k_s = E * m_h / (1.0 - nu * nu);

        TV x0 = vertices.row(0); TV x1 = vertices.row(1); TV x2 = vertices.row(2);
        TV X0 = undeformed_vertices.row(0); TV X1 = undeformed_vertices.row(1); TV X2 = undeformed_vertices.row(2);

        Matrix<T, 9, 9> hess;
        compute3DCSTShellEnergyHessian(nu, k_s, x0, x1, x2, X0, X1, X2, hess);
        // makePD<T, 9>(hess);
        
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                for (int d = 0; d < 3; d++)
                {
                    for (int dd = 0; dd < 3; dd++)
                    {
                        entries.push_back(Entry(indices[i] * 3 + d, 
                            indices[j] * 3 + dd, hess(i * 3 + d, j * 3 + dd)));
                    }   
                } 
            }
        }
    });
}

void DiscreteShell::addShellBendingHessianEntries(std::vector<Entry>& entries)
{
    iterateHingeSerial([&](const HingeIdx& hinge_idx, int hinge_cnt){
        
        HingeVtx deformed_vertices = getHingeVtxDeformed(hinge_idx);
        HingeVtx undeformed_vertices = getHingeVtxUndeformed(hinge_idx);

        std::array<TV, 4> x, X;
        x[0] = deformed_vertices.row(0);
        x[1] = deformed_vertices.row(1);
        x[2] = deformed_vertices.row(2);
        x[3] = deformed_vertices.row(3);

        X[0] = undeformed_vertices.row(0);
        X[1] = undeformed_vertices.row(1);
        X[2] = undeformed_vertices.row(2);
        X[3] = undeformed_vertices.row(3);

        T k_bend = hinge_stiffness[hinge_cnt] * E * std::pow(thickness[0], 3) / (24 * (1.0 - std::pow(nu, 2)));
        
        Matrix<T, 12, 12> hess;
        TV x0 = deformed_vertices.row(0); TV x1 = deformed_vertices.row(1); TV x2 = deformed_vertices.row(2); TV x3 = deformed_vertices.row(3);
        TV X0 = undeformed_vertices.row(0); TV X1 = undeformed_vertices.row(1); TV X2 = undeformed_vertices.row(2); TV X3 = undeformed_vertices.row(3);

        computeDSBendingEnergyHessian(k_bend, x0, x1, x2, x3, X0, X1, X2, X3, hess);
        

        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                for (int d = 0; d < 3; d++)
                {
                    for (int dd = 0; dd < 3; dd++)
                    {
                        entries.push_back(Entry(hinge_idx[i] * 3 + d, 
                            hinge_idx[j] * 3 + dd , hess(i * 3 + d, j * 3 + dd)));
                    }   
                } 
            }
        }
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        

    });
}

void DiscreteShell::addShellHessianEntries(std::vector<Entry>& entries)
{
    addShellInplaneHessianEntries(entries);
    addShellBendingHessianEntries(entries);
}

void DiscreteShell::setHingeStiffness()
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

void DiscreteShell::computeBoundingBox(TV& min_corner, TV& max_corner)
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