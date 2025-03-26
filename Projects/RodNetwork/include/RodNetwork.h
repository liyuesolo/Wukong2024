#ifndef RODNETWORK_H
#define RODNETWORK_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <tbb/tbb.h>
#include <unordered_map>
#include <utility>

#include "Rod.h"
#include "Util.h"
#include "VecMatDef.h"

class RodNetwork
{
public:
    using VectorXT = Matrix<T, Eigen::Dynamic, 1>;
    using MatrixXT = Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using VectorXi = Vector<int, Eigen::Dynamic>;
    using MatrixXi = Matrix<int, Eigen::Dynamic, Eigen::Dynamic>;
    using VtxList = std::vector<int>;
    using StiffnessMatrix = Eigen::SparseMatrix<T>;
    using Entry = Eigen::Triplet<T>;
    using TV = Vector<T, 3>;
    using TV2 = Vector<T, 2>;
    using TM2 = Matrix<T, 2, 2>;
    using TV3 = Vector<T, 3>;
    using IV = Vector<int, 3>;
    using IV2 = Vector<int, 2>;
    using TM = Matrix<T, 3, 3>;

public:
    T ROD_A = 2.5e-3;
    T ROD_B = 2.5e-3;

    VectorXT deformed_states;
    VectorXT rest_states;
    VectorXT dq;

    std::vector<Rod*> rods;
    std::vector<RodCrossing*> rod_crossings;
    std::unordered_map<int, T> dirichlet_data;

    bool add_stretching = true;
    bool add_bending_and_twisting = true;
    bool add_rigid_joint = true;

    bool verbose = false;
    bool run_diff_test = false;
    int max_newton_iter = 2000;
    int ls_max = 10;
    T newton_tol = 1e-6;

private:
    template <class OP>
    void iterateDirichletDoF(const OP& f)
    {
        for (auto dirichlet : dirichlet_data)
            f(dirichlet.first, dirichlet.second);
    }

    template <int dim_row = 2, int dim_col = 2>
    void addHessianEntry(std::vector<Entry>& triplets,
                         const std::vector<int>& vtx_idx,
                         const MatrixXT& hessian, int shift_row = 0,
                         int shift_col = 0)
    {

        for (int i = 0; i < vtx_idx.size(); i++)
        {
            int dof_i = vtx_idx[i];
            for (int j = 0; j < vtx_idx.size(); j++)
            {
                int dof_j = vtx_idx[j];
                for (int k = 0; k < dim_row; k++)
                    for (int l = 0; l < dim_col; l++)
                        triplets.emplace_back(
                            dof_i * dim_row + k + shift_row,
                            dof_j * dim_col + l + shift_col,
                            hessian(i * dim_row + k, j * dim_col + l));
            }
        }
    }

    template <int dim_row = 2, int dim_col = 2>
    void addJacobianEntry(std::vector<Entry>& triplets,
                          const std::vector<int>& vtx_idx,
                          const std::vector<int>& vtx_idx2,
                          const MatrixXT& jacobian, int shift_row = 0,
                          int shift_col = 0)
    {

        for (int i = 0; i < vtx_idx.size(); i++)
        {
            int dof_i = vtx_idx[i];
            for (int j = 0; j < vtx_idx2.size(); j++)
            {
                int dof_j = vtx_idx2[j];
                for (int k = 0; k < dim_row; k++)
                    for (int l = 0; l < dim_col; l++)
                        triplets.emplace_back(
                            dof_i * dim_row + k + shift_row,
                            dof_j * dim_col + l + shift_col,
                            jacobian(i * dim_row + k, j * dim_col + l));
            }
        }
    }

    template <int dim = 2>
    void addForceEntry(VectorXT& residual, const std::vector<int>& vtx_idx,
                       const VectorXT& gradent, int shift = 0)
    {
        for (int i = 0; i < vtx_idx.size(); i++)
            residual.template segment<dim>(vtx_idx[i] * dim + shift) +=
                gradent.template segment<dim>(i * dim);
    }

public:
    RodNetwork() {}
    ~RodNetwork() {}

    // ================ RodNetwork.cpp ==================
    T computeTotalEnergy();

    T computeResidual(VectorXT& residual);

    void projectDirichletDoFMatrix(StiffnessMatrix& A,
                                   const std::unordered_map<int, T>& data);

    void buildSystemDoFMatrix(StiffnessMatrix& K);

    bool linearSolve(StiffnessMatrix& K, const VectorXT& residual,
                     VectorXT& du);

    T lineSearchNewton(const VectorXT& residual);

    bool advanceOneStep(int step);

    void computeBoundingBox(TV& bottom_left, TV& top_right);

    void resetScene();

    // ================== DerivativeTest.cpp ==================
    void testGradientFD();
    void testGradient2ndOrderTerm();
    void testHessian2ndOrderTerm();

    // ================== Stretching.cpp ==================
    T addStretchingEnergy();
    void addStretchingForce(VectorXT& residual);
    void addStretchingHessian(std::vector<Entry>& entry_K);

    // ================== BendingAndTwisting.cpp ==================
    T addBendingAndTwistingEnergy();
    void addBendingAndTwistingForceEntries(VectorXT& residual);
    void addBendingAndTwistingHessianEntries(std::vector<Entry>& entry_K);

    // ================== Joint.cpp ==================
    T addJointBendingAndTwistingEnergy();
    void addJointBendingAndTwistingForceEntries(VectorXT& residual);
    void addJointBendingAndTwistingHessianEntries(std::vector<Entry>& entry_K);

    // initialize from file
    void initializeFromFile(const std::string& filename);

    // code for initializing rods
    void addAStraightRod(const TV& from, const TV& to, int from_idx, int to_idx,
                         T max_segment_length, int& full_dof_cnt, int& node_cnt,
                         int& rod_cnt);

    void addCrossingPoint(std::vector<TV>& existing_nodes, const TV& point,
                          int& full_dof_cnt, int& node_cnt);
};

#endif
