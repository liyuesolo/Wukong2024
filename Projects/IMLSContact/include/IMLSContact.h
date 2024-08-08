#ifndef IMLSContact_H
#define IMLSContact_H

#include <utility>
#include <iostream>
#include <fstream>
#include <Eigen/Geometry>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <tbb/tbb.h>
#include <unordered_map>
#include <unordered_set>
#include <complex>
#include <iomanip>
#include <ipc/ipc.hpp>

#include "VecMatDef.h"
#include "Timer.h"
#include "Util.h"

class IMLSContact
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
    using EleNodes = Matrix<T, 4, 3>;
    using EleIdx = Vector<int, 4>;
    using FaceVtx = Matrix<T, 3, 3>;
    using FaceIdx = Vector<int, 3>;

public:
    // material parameter setup1
    T density = 1.5e3;  //Cotton kg/m^3
    TV gravity = TV(0.0, -9.8, 0.0);    
    T E = 1e6;
    T nu = 0.45;

    T lambda, mu;

    int num_nodes = 0;
    int num_ele = 0;
    
    MatrixXT surface_vertices;
    MatrixXi surface_indices;
    
    VectorXi indices;// tet indices
    VectorXT deformed, undeformed; // tet nodes
    VectorXT u, f;

    bool run_diff_test = false;
    int max_newton_iter = 500;
    bool use_Newton = true;
    bool jump_out = true;
    bool verbose = true;
    T newton_tol = 1e-6;

    std::vector<T> residual_norms;
    std::unordered_map<int, T> dirichlet_data;

    std::unordered_map<int, int> tet_node_surface_map;
    std::unordered_map<int, int> surface_to_tet_node_map;

    bool add_cubic_plane = false;
    T w_plane = 1e3;

    std::vector<std::vector<int>> vtx_tets;
    bool use_VBD = false;
    bool coloring = false;
    VectorXi colors;

    bool use_ipc = false;
    Eigen::MatrixXd ipc_vertices;
    Eigen::MatrixXi ipc_edges, ipc_faces;
    T ipc_barrier_distance = 1e-2;
    T ipc_barrier_weight = 1e3;
    T ipc_min_dis = 1.0;
    T max_barrier_weight = 1e8;
    ipc::CollisionMesh collision_mesh;

    bool dynamics = false;
    T dt = 0.01;
    T simulation_duration = 10;
    StiffnessMatrix M;
    VectorXT mass_diagonal;
    VectorXT xn;
    VectorXT vn;

private:
    template <class OP>
    void iterateDirichletDoF(const OP& f) 
    {
        for (auto dirichlet: dirichlet_data){
            f(dirichlet.first, dirichlet.second);
        } 
    }

    template <typename OP>
    void iterateElementSerial(const OP& f)
    {
        for (int i = 0; i < int(indices.size()/(4)); i++)
        {
            EleIdx tet_idx = indices.segment<4>(i * (4));
            EleNodes tet_deformed = getEleNodesDeformed(tet_idx);
            EleNodes tet_undeformed = getEleNodesUndeformed(tet_idx);
            f(tet_deformed, tet_undeformed, {tet_idx[0], tet_idx[1], tet_idx[2], tet_idx[3]}, i);
        }
    }

    template <typename OP>
    void iterateElementParallel(const OP& f)
    {
        tbb::parallel_for(0, int(indices.size()/(4)), [&](int i)
        {
            EleIdx tet_idx = indices.segment<4>(i * (4));
            EleNodes tet_deformed = getEleNodesDeformed(tet_idx);
            EleNodes tet_undeformed = getEleNodesUndeformed(tet_idx);
            f(tet_deformed, tet_undeformed, {tet_idx[0], tet_idx[1], tet_idx[2], tet_idx[3]}, i);
        });
    }

    EleNodes getEleNodesDeformed(const EleIdx& nodal_indices)
    {
        if (nodal_indices.size() != 4)
            std::cout << "getEleNodesDeformed() not a tet" << std::endl; 
        EleNodes tet_x;
        for (int i = 0; i < 4; i++)
        {
            int idx = nodal_indices[i]*3;
            if (idx > deformed.rows())
                std::cout << "idx out of bound " << std::endl;
            tet_x.row(i) = deformed.segment<3>(idx);
        }
        return tet_x;
    }

    EleNodes getEleNodesUndeformed(const EleIdx& nodal_indices)
    {
        if (nodal_indices.size() != 4)
            std::cout << "getEleNodesUndeformed() not a tet" << std::endl;
        EleNodes tet_x;
        for (int i = 0; i < 4; i++)
        {
            int idx = nodal_indices[i]*3;
            if (idx >= undeformed.rows())
                std::cout << "idx out of bound " << std::endl;
            tet_x.row(i) = undeformed.segment<3>(idx);
        }
        return tet_x;
    }

    EleNodes getEleNodesFromVector(const EleIdx& nodal_indices, const VectorXT& vec)
    {
        EleNodes tet_x;
        for (int i = 0; i < 4; i++)
        {
            int idx = nodal_indices[i]*3;
            if (idx >= vec.rows())
                std::cout << "idx out of bound " << std::endl;
            tet_x.row(i) = vec.segment<3>(idx);
        }
        return tet_x;
    }

    inline T getSmallestPositiveRealQuadRoot(T a, T b, T c, T tol)
    {
        // return negative value if no positive real root is found
        using std::abs;
        using std::sqrt;
        T t;
        if (abs(a) <= tol) {
            if (abs(b) <= tol) // f(x) = c > 0 for all x
                t = -1;
            else
                t = -c / b;
        }
        else {
            T desc = b * b - 4 * a * c;
            if (desc > 0) {
                t = (-b - sqrt(desc)) / (2 * a);
                if (t < 0)
                    t = (-b + sqrt(desc)) / (2 * a);
            }
            else // desv<0 ==> imag
                t = -1;
        }
        return t;
    }

    inline T getSmallestPositiveRealCubicRoot(T a, T b, T c, T d, T tol = 1e-10)
    {
        // return negative value if no positive real root is found
        using std::abs;
        using std::complex;
        using std::pow;
        using std::sqrt;
        T t = -1;
        if (abs(a) <= tol)
            t = getSmallestPositiveRealQuadRoot(b, c, d, tol);
        else {
            complex<T> i(0, 1);
            complex<T> delta0(b * b - 3 * a * c, 0);
            complex<T> delta1(2 * b * b * b - 9 * a * b * c + 27 * a * a * d, 0);
            complex<T> C = pow((delta1 + sqrt(delta1 * delta1 - 4.0 * delta0 * delta0 * delta0)) / 2.0, 1.0 / 3.0);
            if (abs(C) < tol)
                C = pow((delta1 - sqrt(delta1 * delta1 - 4.0 * delta0 * delta0 * delta0)) / 2.0, 1.0 / 3.0);
            complex<T> u2 = (-1.0 + sqrt(3.0) * i) / 2.0;
            complex<T> u3 = (-1.0 - sqrt(3.0) * i) / 2.0;
            complex<T> t1 = (b + C + delta0 / C) / (-3.0 * a);
            complex<T> t2 = (b + u2 * C + delta0 / (u2 * C)) / (-3.0 * a);
            complex<T> t3 = (b + u3 * C + delta0 / (u3 * C)) / (-3.0 * a);
            if ((abs(imag(t1)) < tol) && (real(t1) > 0))
                t = real(t1);
            if ((abs(imag(t2)) < tol) && (real(t2) > 0) && ((real(t2) < t) || (t < 0)))
                t = real(t2);
            if ((abs(imag(t3)) < tol) && (real(t3) > 0) && ((real(t3) < t) || (t < 0)))
                t = real(t3);
        }
        return t;
    }
    
    std::vector<Entry> entriesFromSparseMatrix(const StiffnessMatrix &A) 
    {
        std::vector<Entry> triplets;
        for (int k = 0; k < A.outerSize(); ++k)
            for (StiffnessMatrix::InnerIterator it(A, k); it; ++it)
                triplets.emplace_back(it.row(), it.col(), it.value());
        return triplets;
    }

public:
    bool advanceOneStep(int step);
    bool advanceOneTimeStep();
    bool linearSolve(StiffnessMatrix& K, const VectorXT& residual, VectorXT& du);
    T computeTotalEnergy();
    T computeResidual(VectorXT& residual);
    void buildSystemMatrix(StiffnessMatrix& K);
    T lineSearchNewton(const VectorXT& residual);
    void projectDirichletDoFMatrix(StiffnessMatrix& A, const std::unordered_map<int, T>& data);

    T vertexBlockDescent(const VectorXT& residual);
    void buildVtxTetConnectivity();
    void saveTetsAroundVtx(int vtx_idx);
    T computeGi(int vtx_idx);
    void computeHifi(int vtx_idx, TM& hess, TV& force);
    void graphColoring(const MatrixXT& V, const MatrixXi& TT);

    void computeLinearModes(MatrixXT& eigen_vectors, VectorXT& eigen_values);
    void initializeFromFile(const std::string& filename);
    void initializeSingleKnot();
    void computeBoundingBox(TV& min_corner, TV& max_corner);
    void updateSurfaceVertices();

    // virtual function for different time integration schemes
    //                      different ways of computing the mass matrix
    virtual void addInertialEnergy(T& energy);
    virtual void addInertialForceEntry(VectorXT& residual);
    virtual void addInertialHessianEntries(std::vector<Entry>& entries);
    virtual void updateDynamicStates();
    virtual void initializeDynamicStates();
    virtual void computeMassMatrix();

    T computeInversionFreeStepsize();

    T computeVolume(const EleNodes& x_undeformed);
    void initializeFEMData(const MatrixXT& V, const MatrixXi& F);
    T addElastsicPotential();
    void addElasticForceEntries(VectorXT& residual);
    void addElasticHessianEntries(std::vector<Entry>& entries);

    // Penalty.cpp
    T addCubicPlaneEnergy(T w = 1.0);
    void addCubicPlaneForceEntry(VectorXT& residual, T w = 1.0);
    void addCubicPlaneHessianEntry(std::vector<Entry>& entries, T w = 1.0);

    //DerivativeTest.cpp
    void checkTotalGradient(bool perturb);
    void checkTotalGradientScale(bool perturb);
    void checkTotalHessianScale(bool perturb);

    // IPC.cpp
    void buildIPCRestData();
    T addIPCEnergy();
    void addIPCForceEntries(VectorXT& residual);
    void addIPCHessianEntries(std::vector<Entry>& entries);
    void updateIPCVertices();
    T computeCollisionFreeStepsize(const VectorXT& du);
    void updateBarrierInfo(bool first_step);
public:
    IMLSContact() {}
    ~IMLSContact() {}
};

#endif