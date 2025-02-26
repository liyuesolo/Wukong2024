#ifndef FEMQUADTET_H
#define FEMQUADTET_H

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

#include "Timer.h"
#include "Util.h"
#include "VecMatDef.h"

class FEMQuadTet
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

    using EleNodes = Matrix<T, 10, 3>;
    using EleIdx = Vector<int, 10>;

    T E = 1e4;
    T nu = 0.48;

    T lambda, mu;

    int num_nodes = 0;
    int num_ele = 0;

    VectorXi indices;              // tet indices
    VectorXT deformed, undeformed; // tet nodes
    VectorXT u, f;

    int n_linear_nodes = 0;
    MatrixXi linear_tet_indices;
    MatrixXT linear_tet_vertices;

    bool run_diff_test = false;
    int max_newton_iter = 500;
    bool use_Newton = true;
    bool jump_out = true;
    bool verbose = true;
    T newton_tol = 1e-6;

    std::vector<T> residual_norms;
    std::unordered_map<int, T> dirichlet_data;

    template <class OP>
    void iterateDirichletDoF(const OP& f)
    {
        for (auto dirichlet : dirichlet_data)
        {
            f(dirichlet.first, dirichlet.second);
        }
    }

    template <typename OP>
    void iterateElementSerial(const OP& f)
    {
        for (int i = 0; i < int(indices.size() / (10)); i++)
        {
            EleIdx tet_idx = indices.segment<10>(i * (10));
            EleNodes tet_deformed = getEleNodesDeformed(tet_idx);
            EleNodes tet_undeformed = getEleNodesUndeformed(tet_idx);
            f(tet_deformed, tet_undeformed,
              {tet_idx[0], tet_idx[1], tet_idx[2], tet_idx[3], tet_idx[4],
               tet_idx[5], tet_idx[6], tet_idx[7], tet_idx[8], tet_idx[9]},
              i);
        }
    }

    EleNodes getEleNodesDeformed(const EleIdx& nodal_indices)
    {
        if (nodal_indices.size() != 10)
            std::cout << "getEleNodesDeformed() not a tet" << std::endl;
        EleNodes tet_x;
        for (int i = 0; i < 10; i++)
        {
            int idx = nodal_indices[i] * 3;
            if (idx > deformed.rows())
                std::cout << "idx out of bound " << std::endl;
            tet_x.row(i) = deformed.segment<3>(idx);
        }
        return tet_x;
    }

    EleNodes getEleNodesUndeformed(const EleIdx& nodal_indices)
    {
        if (nodal_indices.size() != 10)
            std::cout << "getEleNodesUndeformed() not a tet" << std::endl;
        EleNodes tet_x;
        for (int i = 0; i < 10; i++)
        {
            int idx = nodal_indices[i] * 3;
            if (idx >= undeformed.rows())
                std::cout << "idx out of bound " << std::endl;
            tet_x.row(i) = undeformed.segment<3>(idx);
        }
        return tet_x;
    }

    TV gauss3DPP2Rule4(int idx)
    {
        assert(idx < 4);
        T h1 = (5 + 3 * std::sqrt(5.0)) / 20; // 0.585410....
        T g1 = (5 - std::sqrt(5.0)) / 20; // 0.1381966...
        TV result;
        result << g1, g1, g1;
        if (idx < 3)
        {
            result[idx] = h1;
        }
        return result;
    }

    T gauss3DPP2Rule4Weight(int idx)
    {
        return 0.25;
    }

public:
    bool advanceOneStep(int step);
    bool linearSolve(StiffnessMatrix& K, const VectorXT& residual,
                     VectorXT& du);
    T computeTotalEnergy();
    T computeResidual(VectorXT& residual);
    void buildSystemMatrix(StiffnessMatrix& K);
    T lineSearchNewton(const VectorXT& residual);
    void projectDirichletDoFMatrix(StiffnessMatrix& A,
                                   const std::unordered_map<int, T>& data);

    void initializeFromFile(const std::string& filename);
    void convertToQuadraticTets(const MatrixXT& Vtet, const MatrixXi& Ttet,
                                MatrixXT& Vquad, MatrixXi& Tquad);
    void computeBoundingBox(TV& min_corner, TV& max_corner);
    
    T elasticPotentialAnalytical();
    T addElastsicPotential();
    void addElasticForceEntries(VectorXT& residual);
    void addElasticHessianEntries(std::vector<Entry>& entries);

    FEMQuadTet() {}
    ~FEMQuadTet() {}
};

#endif
