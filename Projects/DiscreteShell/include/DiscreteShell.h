#ifndef DISCRETE_SHELL_H
#define DISCRETE_SHELL_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <tbb/tbb.h>
#include <unordered_map>
#include <utility>

#include "Timer.h"
#include "Util.h"
#include "VecMatDef.h"

class DiscreteShell
{
public:
    using VectorXT = Matrix<T, Eigen::Dynamic, 1>;
    using MatrixXT = Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using VectorXi = Vector<int, Eigen::Dynamic>;
    using MatrixXi = Matrix<int, Eigen::Dynamic, Eigen::Dynamic>;
    using VtxList = std::vector<int>;
    using TV = Vector<T, 3>;
    using TV2 = Vector<T, 2>;
    using TM2 = Matrix<T, 2, 2>;
    using TV3 = Vector<T, 3>;
    using IV = Vector<int, 3>;
    using IV2 = Vector<int, 2>;
    using TM = Matrix<T, 3, 3>;
    using StiffnessMatrix = Eigen::SparseMatrix<T>;
    using Entry = Eigen::Triplet<T>;
    using Face = Vector<int, 3>;
    using Triangle = Vector<int, 3>;
    using Edge = Vector<int, 2>;

    using Hinges = Matrix<int, Eigen::Dynamic, 4>;

    using FaceVtx = Matrix<T, 3, 3>;
    using FaceIdx = Vector<int, 3>;
    using HingeIdx = Vector<int, 4>;
    using HingeVtx = Matrix<T, 4, 3>;

public:
    // material parameter setup1
    T density = 1.5e3; // Cotton kg/m^3
    TV gravity = TV(0.0, -9.8, 0.0);
    T E = 1e6;
    T nu = 0.45;

    T lambda, mu;

    Hinges hinges;
    std::vector<TM2> Xinv;
    VectorXT shell_rest_area;
    T thickness = 0.003; // meter
    VectorXi faces;
    std::vector<Triangle> triangles;
    VectorXT rest_area;
    VectorXT undeformed_area;

    VectorXT hinge_stiffness;

    VectorXT deformed, undeformed;
    VectorXT u;
    VectorXT external_force;

    bool run_diff_test = false;
    int max_newton_iter = 500;
    bool use_Newton = true;
    bool add_gravity = false;
    bool jump_out = true;
    bool verbose = false;
    bool use_consistent_mass_matrix = true;
    T newton_tol = 1e-6;

    std::vector<T> residual_norms;
    std::unordered_map<int, T> dirichlet_data;

    bool dynamics = false;
    T dt = 0.01;
    T simulation_duration = 10;
    StiffnessMatrix M;
    VectorXT mass_diagonal;
    VectorXT xn;
    VectorXT vn;
    MatrixXT face_normals;

public:
    template <class OP>
    void iterateDirichletDoF(const OP& f)
    {
        for (auto dirichlet : dirichlet_data)
        {
            f(dirichlet.first, dirichlet.second);
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

    template <int dim = 2>
    void getSubVector(const VectorXT& _vector, const std::vector<int>& vtx_idx,
                      VectorXT& sub_vec, int shift = 0)
    {
        sub_vec.resize(vtx_idx.size() * dim);
        sub_vec.setZero();
        for (int i = 0; i < vtx_idx.size(); i++)
        {
            sub_vec.template segment<dim>(i * dim) =
                _vector.template segment<dim>(vtx_idx[i] * dim + shift);
        }
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

    template <int dim_row = 2, int dim_col = 2>
    void addHessianMatrixEntry(MatrixXT& matrix_global,
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
                    {
                        matrix_global(dof_i * dim_row + k + shift_row,
                                      dof_j * dim_col + l + shift_col) +=
                            hessian(i * dim_row + k, j * dim_col + l);
                    }
            }
        }
    }

    template <int dim_row = 2, int dim_col = 2>
    void addJacobianMatrixEntry(MatrixXT& matrix_global,
                                const std::vector<int>& vtx_idx,
                                const std::vector<int>& vtx_idx2,
                                const MatrixXT& hessian, int shift_row = 0,
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
                    {
                        matrix_global(dof_i * dim_row + k + shift_row,
                                      dof_j * dim_col + l + shift_col) +=
                            hessian(i * dim_row + k, j * dim_col + l);
                    }
            }
        }
    }

    template <int dim0 = 2, int dim1 = 2>
    void addHessianBlock(std::vector<Entry>& triplets,
                         const std::vector<int>& vtx_idx,
                         const MatrixXT& hessian_block, int shift_row = 0,
                         int shift_col = 0)
    {

        int dof_i = vtx_idx[0];
        int dof_j = vtx_idx[1];

        for (int k = 0; k < dim0; k++)
            for (int l = 0; l < dim1; l++)
            {
                triplets.emplace_back(dof_i * dim0 + k + shift_row,
                                      dof_j * dim1 + l + shift_col,
                                      hessian_block(k, l));
            }
    }

    template <int size>
    VectorXT computeHessianBlockEigenValues(const Matrix<T, size, size>& symMtr)
    {
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T, size, size>> eigenSolver(
            symMtr);
        return eigenSolver.eigenvalues();
    }

    std::vector<Entry> entriesFromSparseMatrix(const StiffnessMatrix& A)
    {
        std::vector<Entry> triplets;
        for (int k = 0; k < A.outerSize(); ++k)
            for (StiffnessMatrix::InnerIterator it(A, k); it; ++it)
                triplets.emplace_back(it.row(), it.col(), it.value());
        return triplets;
    }

    // ====================== discrete shell ==========================
    void updateLameParameters()
    {
        lambda = E * nu / (1.0 + nu) / (1.0 - 2.0 * nu);
        mu = E / 2.0 / (1.0 + nu);
    }

    HingeVtx getHingeVtxDeformed(const HingeIdx& hi)
    {
        HingeVtx cellx;
        for (int i = 0; i < 4; i++)
        {
            cellx.row(i) = deformed.segment<3>(hi[i] * 3);
        }
        return cellx;
    }

    HingeVtx getHingeVtxUndeformed(const HingeIdx& hi)
    {
        HingeVtx cellx;
        for (int i = 0; i < 4; i++)
        {
            cellx.row(i) = undeformed.segment<3>(hi[i] * 3);
        }
        return cellx;
    }

    FaceVtx getFaceVtxDeformed(int face)
    {
        FaceVtx cellx;
        FaceIdx nodal_indices = faces.segment<3>(face * 3);
        for (int i = 0; i < 3; i++)
        {
            cellx.row(i) = deformed.segment<3>(nodal_indices[i] * 3);
        }
        return cellx;
    }

    FaceVtx getFaceVtxUndeformed(int face)
    {
        FaceVtx cellx;
        FaceIdx nodal_indices = faces.segment<3>(face * 3);
        for (int i = 0; i < 3; i++)
        {
            cellx.row(i) = undeformed.segment<3>(nodal_indices[i] * 3);
        }
        return cellx;
    }

    template <typename OP>
    void iterateFaceSerial(const OP& f)
    {
        for (int i = 0; i < faces.rows() / 3; i++)
            f(i);
    }

    template <typename OP>
    void iterateNodeSerial(const OP& f)
    {
        for (int i = 0; i < deformed.rows() / 3; i++)
            f(i);
    }

    template <typename OP>
    void iterateFaceParallel(const OP& f)
    {
        tbb::parallel_for(0, int(faces.rows() / 3), [&](int i) { f(i); });
    }

    template <typename OP>
    void iterateTriangleSerial(const OP& f)
    {
        for (int i = 0; i < triangles.size(); i++)
            f(triangles[i], i);
    }

    template <typename OP>
    void iterateHingeSerial(const OP& f)
    {
        for (int i = 0; i < hinges.rows(); i++)
        {
            const Vector<int, 4> nodes = hinges.row(i);
            f(nodes, i);
        }
    }
    void projectDirichletDoFMatrix(StiffnessMatrix& A,
                                   const std::unordered_map<int, T>& data);

    T computeAngle(const TV& x0, const TV& x1, const TV& x2, const TV& x3)
    {
        // we want to measure which angle we would have to rotate a vector from
        // x0 to the edge around the axis (x2 - x1) to get to x3
        //  (this fixes the sign of the angle we are measuring)
        TV n = (x2 - x1).cross(x0 - x1).normalized();
        TV m = (x2 - x1).cross(n).normalized();
        T y = -n.dot(x3 - x1);
        T x = m.dot(x3 - x1);
        return atan2(y, x);
    }

    T computeAngleDiff(const TV& x0, const TV& x1, const TV& x2, const TV& x3,
                       T rest_angle)
    {
        TV n = (x2 - x1).cross(x0 - x1).normalized();
        TV m = (x2 - x1).cross(n).normalized();
        T y = -n.dot(x3 - x1);
        T x = m.dot(x3 - x1);
        // rotate 2d vector s.t. reference angle is at 0....
        T rotAngle = -rest_angle;
        T xr = std::cos(rotAngle) * x - std::sin(rotAngle) * y;
        T yr = std::sin(rotAngle) * x + std::cos(rotAngle) * y;

        return atan2(yr, xr);
    }

    Matrix<T, 2, 3> barycentricJacobian(const TV& a, const TV& b, const TV& c)
    {
        TV v0 = b - a, v1 = c - a;
        // TV v2 = p - a;
        T d00 = v0.dot(v0); // (b - a).dot(b - a) => d d00 / d p = 0
        T d01 = v0.dot(v1); // (b - a).dot(c - a) => d d01 / d p = 0
        T d11 = v1.dot(v1); // (c - a).dot(c - a) => d d11 / dp = 0
        // T d20 = v2.dot(v0); // (p - a).dot(b - a) => d d20 / dp = (b - a)^T
        // T d21 = v2.dot(v1); //  ( p - a).dot(c - a) = > d21 / dp = (c - a)^T
        T denom =
            d00 * d11 - d01 * d01; // => d00 * d11 constant in => drops out, d01
                                   // constant in p => derivative is 0
        // v = (d11 * d20 - d01 * d21) / denom;
        // w = (d00 * d21 - d01 * d20) / denom;
        // u = 1.0f - v - w;
        // TV dvdp = (d11 * dd20 / dp - d01 * d d21 / dp) / denom;
        TV dvdp = (d11 * (b - a) - d01 * (c - a)) / denom;
        TV dwdp = (d00 * (c - a) - d01 * (b - a)) / denom;
        Matrix<T, 2, 3> result;
        result.row(0) = dvdp.transpose();
        result.row(1) = dwdp.transpose();
        return result;
    }

    Matrix<T, 3, 2> compute3DCSTDeformationGradient(const TV& x1Undef,
                                                    const TV& x2Undef,
                                                    const TV& x3Undef,
                                                    const TV& x1, const TV& x2,
                                                    const TV& x3)
    {
        // defGrad = d x / d X
        // X(b) = X * N(b) = X * [ 1 - b1 - b2; b1; b2];
        // b(X) = [X2 - X1, X3 - X1]^-1 [X - X1]
        // x(X) = x * N(b(X)) then take the jacobian of this to get defgrad
        // d x / d X = x * d N / dX = x * dN/db * [X2 - X1, X3 - X1]^-1;
        // however here X means a 2 dimensional vector! so we get a 3x2 defgrad
        // x(X) = x * N(Barycentric(X, X1, X2, X3))
        // defGrad = dx / d X = x * dN/dB * dB /dX
        // that would work except that the defGradient is gonna be 3x3 and in
        // the undef config. something like [1, 0, 0;0,0,0;0,0,1]; and then E =
        // 0.5 * (F^T F - I)  is gonna give a non zero energy in the undef
        // config.
        // instead we choose a 2D coordinate system X* in the undef
        // configuration for which we compute the def grad
        // defGrad = d x / d X*
        // x(X*) = x * N(Barycentric(X(X*)));
        // X(X*) = X1 + t * X*[0] + q *X*[1]

        TV tUndef = (x2Undef - x1Undef).normalized();
        TV e2Undef = (x3Undef - x1Undef);
        TV qUndef = (e2Undef - tUndef * e2Undef.dot(tUndef)).normalized();

        Eigen::Matrix3d x;
        x << x1, x2, x3;

        // N(b) = [1 - b1 - b2, b1, b2]
        Matrix<T, 3, 2> dNdb;
        dNdb << -1.0, -1.0, 1.0, 0.0, 0.0, 1.0;

        Matrix<T, 2, 3> dBdX = barycentricJacobian(x1Undef, x2Undef, x3Undef);
        Matrix<T, 3, 2> dXdXStar;
        dXdXStar << tUndef, qUndef;
        Matrix<T, 3, 2> defGrad =
            x * dNdb * dBdX *
            dXdXStar; // note that this F is not very intuitive it can contain
                      // -1 for undef configuration, but its not a problem as
                      // long as only F^T*F is used
        return defGrad;
    }
    // ================================================================

public:
    void initializeDynamicExampleScene(const std::string& filename);
    void initializeNonManifoldExampleScene(const std::string& filename);

    bool advanceOneStep(int step);
    bool advanceOneTimeStep();
    bool linearSolve(StiffnessMatrix& K, const VectorXT& residual,
                     VectorXT& du);
    T computeTotalEnergy();
    T computeResidual(VectorXT& residual);
    void buildSystemMatrix(StiffnessMatrix& K);
    T lineSearchNewton(const VectorXT& residual);

    void computeLinearModes(MatrixXT& eigen_vectors, VectorXT& eigen_values);
    void initializeFromFile(const std::string& filename);
    void buildHingeStructure();
    int nFaces() { return faces.rows() / 3; }

    void addShellEnergy(T& energy);
    void addShellForceEntry(VectorXT& residual);
    void addShellHessianEntries(std::vector<Entry>& entries);

    // virtual function for different time integration schemes
    //                      different ways of computing the mass matrix
    virtual void addInertialEnergy(T& energy);
    virtual void addInertialForceEntry(VectorXT& residual);
    virtual void addInertialHessianEntries(std::vector<Entry>& entries);
    virtual void updateDynamicStates();
    virtual void initializeDynamicStates();
    virtual void computeMassMatrix();
    virtual void computeConsistentMassMatrix(const FaceVtx& p,
                                             Matrix<T, 9, 9>& mass_mat);

    // stretching and bending energy
    virtual void addShellInplaneEnergy(T& energy);
    virtual void addShellBendingEnergy(T& energy);
    virtual void addShellInplaneForceEntries(VectorXT& residual);
    virtual void addShellBendingForceEntries(VectorXT& residual);
    virtual void addShellInplaneHessianEntries(std::vector<Entry>& entries);
    virtual void addShellBendingHessianEntries(std::vector<Entry>& entries);

    // graviational energy
    void addShellGravitionEnergy(T& energy);
    void addShellGravitionForceEntry(VectorXT& residual);
    void addShellGravitionHessianEntry(std::vector<Entry>& entries);

    void computeBoundingBox(TV& min_corner, TV& max_corner);

    void setHingeStiffness();

    // derivative tests
    void checkTotalGradient(bool perturb = false);
    void checkTotalGradientScale(bool perturb = false);
    void checkTotalHessian(bool perturb = false);
    void checkTotalHessianScale(bool perturb = false);

public:
    DiscreteShell() { updateLameParameters(); }
    ~DiscreteShell() {}
};

#endif