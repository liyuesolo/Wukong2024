#ifndef ROD2D_H
#define ROD2D_H

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

class Rod2D
{
public:
    using VectorXT = Matrix<T, Eigen::Dynamic, 1>;
    using MatrixXT = Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using VectorXi = Vector<int, Eigen::Dynamic>;
    using MatrixXi = Matrix<int, Eigen::Dynamic, Eigen::Dynamic>;
    using VtxList = std::vector<int>;
    using StiffnessMatrix = Eigen::SparseMatrix<T>;
    using Entry = Eigen::Triplet<T>;
    using TV3 = Vector<T, 3>;
    using TV = Vector<T, 2>;
    using TM = Matrix<T, 2, 2>;

    using IV3 = Vector<int, 3>;
    using IV = Vector<int, 2>;
    using TM3 = Matrix<T, 3, 3>;

public:
    bool verbose = false;
    bool dynamics = true;
    VectorXT deformed, undeformed, u;
    VectorXT xn, vn;
    VectorXT external_force;
    bool run_diff_test = false;
    std::unordered_map<int, T> dirichlet_data;
    T newton_tol = 1e-6;
    int max_newton_iter = 100;
    T dt = 1.0 / 30.0;
    int simulation_duration = 1000;
    int ls_max = 10;
    T kb = 1.0;
    T ks = 1000.0;
    // 0 1; 1 2; 2 3; 3 4; 4 5
    VectorXi rod_indices;

    T density = 1.0;

public:
    template <class OP>
    void iterateDirichletDoF(const OP& f)
    {
        for (auto dirichlet : dirichlet_data)
        {
            f(dirichlet.first, dirichlet.second);
        }
    }

    T computeMaterialCurvature2D(const TV& x1, const TV& x2, const TV& x3)
    {
        TV t1 = (x2 - x1).normalized();
        TV t2 = (x3 - x2).normalized();
        return 2.0 * (t1[0] * t2[1] - t1[1] * t2[0]) / (1.0 + t1.dot(t2));
    }

    template <class OP>
    void iterateRodSegments(const OP& f)
    {
        for (int i = 0; i < rod_indices.rows() / 2; i++)
        {
            f(rod_indices[i * 2], rod_indices[i * 2 + 1], i);
        }
    }

    template <class OP>
    void iterateRodBendingSegments(const OP& f)
    {
        for (int i = 0; i < rod_indices.rows() / 2 - 1; i++)
        {
            f(rod_indices[i * 2], rod_indices[i * 2 + 1],
              rod_indices[(i + 1) * 2 + 1], i);
        }
    }

public:
    Rod2D() {}
    ~Rod2D() {}

    void addInertialEnergy(T& energy);
    void addInertialForceEntries(VectorXT& residual);
    void addInertialHessianEntries(std::vector<Entry>& entries);

    void addElasticEnergy(T& energy);
    void addElasticForceEntries(VectorXT& residual);
    void addElasticHessianEntries(std::vector<Entry>& entries);

    void initialize();
    bool advanceOneStep(int step);
    bool advanceOneTimeStep();
    bool linearSolve(StiffnessMatrix& K, const VectorXT& residual,
                     VectorXT& du);
    void buildSystemMatrix(StiffnessMatrix& K);
    T lineSearchNewton(const VectorXT& residual);
    T computeTotalEnergy();
    T computeResidual(VectorXT& residual);
    void projectDirichletDoFMatrix(StiffnessMatrix& A,
                                   const std::unordered_map<int, T>& data);
    void updateDynamicStates();
};

#endif
