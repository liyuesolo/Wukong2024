#ifndef FEM3D_H
#define FEM3D_H

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

#include "VecMatDef.h"
#include "Timer.h"
#include "Util.h"

class FEM3D
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

public:
    T E = 1e3;
    T nu = 0.48;

    T lambda, mu;

    int num_nodes = 0;
    int num_ele = 0;
    
    MatrixXT surface_vertices;
    MatrixXi surface_indices;
    
    VectorXi indices;// tet indices
    VectorXT deformed, undeformed; // tet nodes
    VectorXT u;



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

private:
    template <class OP>
    void iterateDirichletDoF(const OP& f) 
    {
        for (auto dirichlet: dirichlet_data){
            f(dirichlet.first, dirichlet.second);
        } 
    }

    // template<int dim = 2>
    // void addForceEntry(VectorXT& residual, 
    //     const std::vector<int>& vtx_idx, 
    //     const VectorXT& gradent, int shift = 0)
    // {
    //     for (int i = 0; i < vtx_idx.size(); i++)
    //         residual.template segment<dim>(vtx_idx[i] * dim + shift) += gradent.template segment<dim>(i * dim);
    // }

    // template<int dim = 2>
    // void getSubVector(const VectorXT& _vector, 
    //     const std::vector<int>& vtx_idx, 
    //     VectorXT& sub_vec, int shift = 0)
    // {
    //     sub_vec.resize(vtx_idx.size() * dim);
    //     sub_vec.setZero();
    //     for (int i = 0; i < vtx_idx.size(); i++)
    //     {
    //         sub_vec.template segment<dim>(i * dim) = _vector.template segment<dim>(vtx_idx[i] * dim + shift);
    //     }
    // }

    
    // template<int dim_row=2, int dim_col=2>
    // void addHessianEntry(
    //     std::vector<Entry>& triplets,
    //     const std::vector<int>& vtx_idx, 
    //     const MatrixXT& hessian, 
    //     int shift_row = 0, int shift_col=0)
    // {
        
    //     for (int i = 0; i < vtx_idx.size(); i++)
    //     {
    //         int dof_i = vtx_idx[i];
    //         for (int j = 0; j < vtx_idx.size(); j++)
    //         {
    //             int dof_j = vtx_idx[j];
    //             for (int k = 0; k < dim_row; k++)
    //                 for (int l = 0; l < dim_col; l++)
    //                     triplets.emplace_back(
    //                             dof_i * dim_row + k + shift_row, 
    //                             dof_j * dim_col + l + shift_col, 
    //                             hessian(i * dim_row + k, j * dim_col + l)
    //                         );                
    //         }
    //     }
    // }

    // template<int dim_row=2, int dim_col=2>
    // void addJacobianEntry(
    //     std::vector<Entry>& triplets,
    //     const std::vector<int>& vtx_idx,
    //     const std::vector<int>& vtx_idx2, 
    //     const MatrixXT& jacobian, 
    //     int shift_row = 0, int shift_col=0)
    // {
        
    //     for (int i = 0; i < vtx_idx.size(); i++)
    //     {
    //         int dof_i = vtx_idx[i];
    //         for (int j = 0; j < vtx_idx2.size(); j++)
    //         {
    //             int dof_j = vtx_idx2[j];
    //             for (int k = 0; k < dim_row; k++)
    //                 for (int l = 0; l < dim_col; l++)
    //                     triplets.emplace_back(
    //                             dof_i * dim_row + k + shift_row, 
    //                             dof_j * dim_col + l + shift_col, 
    //                             jacobian(i * dim_row + k, j * dim_col + l)
    //                         ); 
    //         }
    //     }
    // }

    // template<int dim_row=2, int dim_col=2>
    // void addHessianMatrixEntry(
    //     MatrixXT& matrix_global,
    //     const std::vector<int>& vtx_idx, 
    //     const MatrixXT& hessian,
    //     int shift_row = 0, int shift_col=0)
    // {
    //     for (int i = 0; i < vtx_idx.size(); i++)
    //     {
    //         int dof_i = vtx_idx[i];
    //         for (int j = 0; j < vtx_idx.size(); j++)
    //         {
    //             int dof_j = vtx_idx[j];
    //             for (int k = 0; k < dim_row; k++)
    //                 for (int l = 0; l < dim_col; l++)
    //                 {
    //                     matrix_global(dof_i * dim_row + k + shift_row, 
    //                     dof_j * dim_col + l + shift_col) 
    //                         += hessian(i * dim_row + k, j * dim_col + l);
    //                 }
    //         }
    //     }
    // }

    // template<int dim_row=2, int dim_col=2>
    // void addJacobianMatrixEntry(
    //     MatrixXT& matrix_global,
    //     const std::vector<int>& vtx_idx, 
    //     const std::vector<int>& vtx_idx2, 
    //     const MatrixXT& hessian,
    //     int shift_row = 0, int shift_col=0)
    // {
    //     for (int i = 0; i < vtx_idx.size(); i++)
    //     {
    //         int dof_i = vtx_idx[i];
    //         for (int j = 0; j < vtx_idx2.size(); j++)
    //         {
    //             int dof_j = vtx_idx2[j];
    //             for (int k = 0; k < dim_row; k++)
    //                 for (int l = 0; l < dim_col; l++)
    //                 {
    //                     matrix_global(dof_i * dim_row + k + shift_row, 
    //                     dof_j * dim_col + l + shift_col) 
    //                         += hessian(i * dim_row + k, j * dim_col + l);
    //                 }
    //         }
    //     }
    // }

    // template<int dim0=2, int dim1=2>
    // void addHessianBlock(
    //     std::vector<Entry>& triplets,
    //     const std::vector<int>& vtx_idx, 
    //     const MatrixXT& hessian_block,
    //     int shift_row = 0, int shift_col=0)
    // {

    //     int dof_i = vtx_idx[0];
    //     int dof_j = vtx_idx[1];
        
    //     for (int k = 0; k < dim0; k++)
    //         for (int l = 0; l < dim1; l++)
    //         {
    //             triplets.emplace_back(dof_i * dim0 + k + shift_row, 
    //                 dof_j * dim1 + l + shift_col, 
    //                 hessian_block(k, l));
    //         }
    // }

    // template<int size>
    // VectorXT computeHessianBlockEigenValues(const Matrix<T, size, size> & symMtr)
    // {
    //     Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T, size, size>> eigenSolver(symMtr);
    //     return eigenSolver.eigenvalues();
    // }

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
    

public:
    bool advanceOneStep(int step);
    bool linearSolve(StiffnessMatrix& K, const VectorXT& residual, VectorXT& du);
    T computeTotalEnergy();
    T computeResidual(VectorXT& residual);
    void buildSystemMatrix(StiffnessMatrix& K);
    T lineSearchNewton(const VectorXT& residual);
    void projectDirichletDoFMatrix(StiffnessMatrix& A, const std::unordered_map<int, T>& data);

    void computeLinearModes(MatrixXT& eigen_vectors, VectorXT& eigen_values);
    void initializeFromFile(const std::string& filename);
    void computeBoundingBox(TV& min_corner, TV& max_corner);
    void updateSurfaceVertices();

    T computeInversionFreeStepsize();

    T computeVolume(const EleNodes& x_undeformed);
    void initializeFEMData(const MatrixXT& V, const MatrixXi& F);
    T addElastsicPotential();
    void addElasticForceEntries(VectorXT& residual);
    void addElasticHessianEntries(std::vector<Entry>& entries);

public:
    FEM3D() {}
    ~FEM3D() {}
};

#endif