#include <Eigen/CholmodSupport>
#include <igl/readOBJ.h>
#include <igl/massmatrix.h>
#include <Spectra/SymEigsShiftSolver.h>
#include <Spectra/MatOp/SparseSymShiftSolve.h>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/SparseSymMatProd.h>

#include "../autodiff/CST3DShell.h"
#include "../include/DiscreteShell.h"

bool DiscreteShell::advanceOneTimeStep()
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

bool DiscreteShell::advanceOneStep(int step)
{
    if (dynamics)
    {
        std::cout << "=================== Time STEP " << step * dt << "s===================" << std::endl;
        bool finished = advanceOneTimeStep();
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

bool DiscreteShell::linearSolve(StiffnessMatrix& K, const VectorXT& residual, VectorXT& du)
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

T DiscreteShell::computeTotalEnergy()
{
    deformed = undeformed + u;

    T energy = 0.0;
    addShellEnergy(energy);
    if (add_gravity)
        addShellGravitionEnergy(energy);
    if (dynamics)
        addInertialEnergy(energy);
    energy -= u.dot(external_force);
    return energy;
}

T DiscreteShell::computeResidual(VectorXT& residual)
{
    deformed = undeformed + u;
    addShellForceEntry(residual);
    if (add_gravity)
        addShellGravitionForceEntry(residual);
    if (dynamics)
        addInertialForceEntry(residual);

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

void DiscreteShell::buildSystemMatrix(StiffnessMatrix& K)
{
    deformed = undeformed + u;
    std::vector<Entry> entries;
    addShellHessianEntries(entries);
    if (add_gravity)
        addShellGravitionHessianEntry(entries);
    if (dynamics)
        addInertialHessianEntries(entries);
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
            // if (cnt > 10)
                // std::cout << "cnt > 10" << " |du| " << norm << " |g| " << residual.norm() << std::endl;
            break;
        }
        alpha *= 0.5;
        cnt += 1;
    }

    return alpha * du.norm();
}

void DiscreteShell::initializeFromFile(const std::string& filename)
{
    MatrixXT V; MatrixXi F;
    igl::readOBJ(filename, V, F);

    TV min_corner = V.colwise().minCoeff();
    TV max_corner = V.colwise().maxCoeff();

    T bb_diag = (max_corner - min_corner).norm();

    V *= 1.0 / bb_diag;

    V *= 0.5;

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
    deformed = undeformed;
    u = VectorXT::Zero(deformed.rows());
    external_force = VectorXT::Zero(deformed.rows());
    external_force[5 * 3 + 2] = -10.0;
    buildHingeStructure();
    dynamics = true;
    add_gravity = true;
    use_consistent_mass_matrix = true;
    // E = 0.0;
    dt = 1e-2;
    simulation_duration = 1000000;
    
    hinge_stiffness.setConstant(10);
    if (dynamics)
    {
        initializeDynamicStates();
    }
    
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

void DiscreteShell::addShellInplaneEnergy(T& energy)
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

void DiscreteShell::addShellBendingEnergy(T& energy)
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

void DiscreteShell::addShellEnergy(T& energy)
{
    T in_plane_energy = 0.0;
    T bending_energy = 0.0;
    addShellInplaneEnergy(in_plane_energy);
    addShellBendingEnergy(bending_energy);    

    energy += in_plane_energy;
    energy += bending_energy;
}

void DiscreteShell::addShellInplaneForceEntries(VectorXT& residual)
{
    iterateFaceSerial([&](int face_idx)
    {
        FaceVtx vertices = getFaceVtxDeformed(face_idx);
        FaceVtx undeformed_vertices = getFaceVtxUndeformed(face_idx);
        FaceIdx indices = faces.segment<3>(face_idx * 3);

        T k_s = E * thickness / (1.0 - nu * nu);
        TV x0 = vertices.row(0); TV x1 = vertices.row(1); TV x2 = vertices.row(2);
        TV X0 = undeformed_vertices.row(0); TV X1 = undeformed_vertices.row(1); TV X2 = undeformed_vertices.row(2);
        
        Vector<T, 9> dedx;
        compute3DCSTShellEnergyGradient(nu, k_s, x0, x1, x2, X0, X1, X2, dedx);
        
        addForceEntry<3>(residual, {indices[0], indices[1], indices[2]}, -dedx);
    });
}

void DiscreteShell::addShellBendingForceEntries(VectorXT& residual)
{
    iterateHingeSerial([&](const HingeIdx& hinge_idx, int hinge_cnt){
                
        HingeVtx deformed_vertices = getHingeVtxDeformed(hinge_idx);
        HingeVtx undeformed_vertices = getHingeVtxUndeformed(hinge_idx);


        T k_bend = hinge_stiffness[hinge_cnt] * E * std::pow(thickness, 3) / (24 * (1.0 - std::pow(nu, 2)));

        Vector<T, 12> dedx;
        
        TV x0 = deformed_vertices.row(0); TV x1 = deformed_vertices.row(1); TV x2 = deformed_vertices.row(2); TV x3 = deformed_vertices.row(3);
        TV X0 = undeformed_vertices.row(0); TV X1 = undeformed_vertices.row(1); TV X2 = undeformed_vertices.row(2); TV X3 = undeformed_vertices.row(3);

        computeDSBendingEnergyGradient(k_bend, x0, x1, x2, x3, X0, X1, X2, X3, dedx);
        addForceEntry<3>(residual, 
            {hinge_idx[0], hinge_idx[1], hinge_idx[2], hinge_idx[3]}, -dedx);
    });
}
void DiscreteShell::addShellForceEntry(VectorXT& residual)
{
    addShellInplaneForceEntries(residual);
    addShellBendingForceEntries(residual);
}

void DiscreteShell::addShellInplaneHessianEntries(std::vector<Entry>& entries)
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

void DiscreteShell::addShellBendingHessianEntries(std::vector<Entry>& entries)
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

void DiscreteShell::addShellGravitionEnergy(T& energy)
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

void DiscreteShell::addShellGravitionForceEntry(VectorXT& residual)
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

        Vector<T, 9> dedx;
        compute3DCSTGravitationalEnergyGradient(density, thickness, gravity, x0, x1, x2, X0, X1, X2, dedx);

        addForceEntry<3>(residual, {indices[0], indices[1], indices[2]}, -dedx);
    });
}

void DiscreteShell::addShellGravitionHessianEntry(std::vector<Entry>& entries)
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

void DiscreteShell::checkTotalGradient(bool perturb)
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

void DiscreteShell::checkTotalGradientScale(bool perturb)
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


void DiscreteShell::checkTotalHessian(bool perturb)
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


void DiscreteShell::checkTotalHessianScale(bool perturb)
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
void DiscreteShell::addInertialEnergy(T& energy)
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
            Vector<T, 9> x_n_plus_1_vec, xn_vec, vn_vec;
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

void DiscreteShell::addInertialForceEntry(VectorXT& residual)
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
            Vector<T, 9> x_n_plus_1_vec, xn_vec, vn_vec;
            for (int i = 0; i < 3; i++)
            {
                x_n_plus_1_vec.segment<3>(i * 3) = vertices.row(i);
                xn_vec.segment<3>(i * 3) = xn.segment<3>(indices[i] * 3);
                vn_vec.segment<3>(i * 3) = vn.segment<3>(indices[i] * 3);
            }
            Vector<T, 9> dedx = mass_mat * (2.0 * x_n_plus_1_vec - 2.0 * (xn_vec + vn_vec * dt)) / (2.0 * std::pow(dt, 2));
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

void DiscreteShell::addInertialHessianEntries(std::vector<Entry>& entries)
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

void DiscreteShell::updateDynamicStates()
{
    vn = (deformed - xn) / dt;
    xn = deformed;
}

void DiscreteShell::initializeDynamicStates()
{
    if (!use_consistent_mass_matrix)
        computeMassMatrix();
    xn = undeformed;
    vn = VectorXT::Zero(undeformed.rows());
}

void DiscreteShell::computeMassMatrix()
{
    MatrixXT V; MatrixXi F;
    vectorToIGLMatrix<T, 3>(undeformed, V);
    vectorToIGLMatrix<int, 3>(faces, F);
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, M);
    mass_diagonal.resize(deformed.rows() / 3);
    mass_diagonal.setZero();
    mass_diagonal = M.diagonal();   
}

void DiscreteShell::computeConsistentMassMatrix(const FaceVtx& vtx_pos, Matrix<T, 9, 9>& mass_mat)
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
