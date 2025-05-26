#ifndef IMLSCONTACT_H
#define IMLSCONTACT_H

#include <utility>
#include <iostream>
#include <fstream>
#include <Eigen/Geometry>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <tbb/tbb.h>
#include <unordered_map>
#include <complex>
#include <iomanip>

#include "VecMatDef.h"
#include "Timer.h"
#include "Util.h"

#include <CGAL/Epick_d.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Fuzzy_iso_box.h>
#include <CGAL/Search_traits_d.h>

#include "../include/IMLSdiff.h"
#include "../include/RIMLSdiff.h"
// #include <ceres/gradient_problem_solver.h>
// #include <ceres/first_order_function.h>
// #include <ceres/gradient_problem.h>
#include <algorithm>      // for std::min
#include <boost/math/tools/roots.hpp>


typedef CGAL::Epick_d<CGAL::Dimension_tag<3> > K;
typedef K::Point_d Point_d;
typedef std::tuple<Point_d, int> Point_and_int;
typedef CGAL::Search_traits_d<K>                       Traits_base;
typedef CGAL::Search_traits_adapter<Point_and_int, CGAL::Nth_of_tuple_property_map<0, Point_and_int>, Traits_base>       Traits;
typedef CGAL::Kd_tree<Traits> Tree;
typedef CGAL::Fuzzy_sphere<Traits> Fuzzy_sphere;

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

    using Hinges = Matrix<int, Eigen::Dynamic, 4>;
    
    using FaceVtx = Matrix<T, 3, 3>;
    using FaceIdx = Vector<int, 3>;
    using HingeIdx = Vector<int, 4>;
    using HingeVtx = Matrix<T, 4, 3>;

 
public:
    // material parameter setup1
    T density = 1.5e3;  //Cotton kg/m^3
    TV gravity = TV(0.0, -9.8, 0.0);    
    T E = 1e4;
    T nu = 0.3;

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
    T dt = 0.1;
    T simulation_duration = 10;
    StiffnessMatrix M;
    VectorXT mass_diagonal;
    VectorXT xn;
    VectorXT vn;

    

public:
    template <class OP>
    void iterateDirichletDoF(const OP& f) 
    {
        for (auto dirichlet: dirichlet_data){
            f(dirichlet.first, dirichlet.second);
        } 
    }

    template<int dim = 2>
    void addForceEntry(VectorXT& residual, 
        const std::vector<int>& vtx_idx, 
        const VectorXT& gradent, int shift = 0)
    {
        for (int i = 0; i < vtx_idx.size(); i++)
            residual.template segment<dim>(vtx_idx[i] * dim + shift) += gradent.template segment<dim>(i * dim);
    }

    template<int dim = 2>
    void getSubVector(const VectorXT& _vector, 
        const std::vector<int>& vtx_idx, 
        VectorXT& sub_vec, int shift = 0)
    {
        sub_vec.resize(vtx_idx.size() * dim);
        sub_vec.setZero();
        for (int i = 0; i < vtx_idx.size(); i++)
        {
            sub_vec.template segment<dim>(i * dim) = _vector.template segment<dim>(vtx_idx[i] * dim + shift);
        }
    }

    
    template<int dim_row=2, int dim_col=2>
    void addHessianEntry(
        std::vector<Entry>& triplets,
        const std::vector<int>& vtx_idx, 
        const MatrixXT& hessian, 
        int shift_row = 0, int shift_col=0)
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
                                hessian(i * dim_row + k, j * dim_col + l)
                            );                
            }
        }
    }

    template<int dim_row=2, int dim_col=2>
    void addJacobianEntry(
        std::vector<Entry>& triplets,
        const std::vector<int>& vtx_idx,
        const std::vector<int>& vtx_idx2, 
        const MatrixXT& jacobian, 
        int shift_row = 0, int shift_col=0)
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
                                jacobian(i * dim_row + k, j * dim_col + l)
                            ); 
            }
        }
    }

    template<int dim_row=2, int dim_col=2>
    void addHessianMatrixEntry(
        MatrixXT& matrix_global,
        const std::vector<int>& vtx_idx, 
        const MatrixXT& hessian,
        int shift_row = 0, int shift_col=0)
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
                        dof_j * dim_col + l + shift_col) 
                            += hessian(i * dim_row + k, j * dim_col + l);
                    }
            }
        }
    }

    template<int dim_row=2, int dim_col=2>
    void addJacobianMatrixEntry(
        MatrixXT& matrix_global,
        const std::vector<int>& vtx_idx, 
        const std::vector<int>& vtx_idx2, 
        const MatrixXT& hessian,
        int shift_row = 0, int shift_col=0)
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
                        dof_j * dim_col + l + shift_col) 
                            += hessian(i * dim_row + k, j * dim_col + l);
                    }
            }
        }
    }

    template<int dim0=2, int dim1=2>
    void addHessianBlock(
        std::vector<Entry>& triplets,
        const std::vector<int>& vtx_idx, 
        const MatrixXT& hessian_block,
        int shift_row = 0, int shift_col=0)
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

    template<int size>
    VectorXT computeHessianBlockEigenValues(const Matrix<T, size, size> & symMtr)
    {
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T, size, size>> eigenSolver(symMtr);
        return eigenSolver.eigenvalues();
    }

    std::vector<Entry> entriesFromSparseMatrix(const StiffnessMatrix &A) 
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
            cellx.row(i) = deformed.segment<3>(hi[i]*3);
        }
        return cellx;
    }

    HingeVtx getHingeVtxUndeformed(const HingeIdx& hi)
    {
        HingeVtx cellx;
        for (int i = 0; i < 4; i++)
        {
            cellx.row(i) = undeformed.segment<3>(hi[i]*3);
        }
        return cellx;
    }

    FaceVtx getFaceVtxDeformed(int face)
    {
        FaceVtx cellx;
        FaceIdx nodal_indices = faces.segment<3>(face * 3);
        for (int i = 0; i < 3; i++)
        {
            cellx.row(i) = deformed.segment<3>(nodal_indices[i]*3);
        }
        return cellx;
    }

    FaceVtx getFaceVtxUndeformed(int face)
    {
        FaceVtx cellx;
        FaceIdx nodal_indices = faces.segment<3>(face * 3);
        for (int i = 0; i < 3; i++)
        {
            cellx.row(i) = undeformed.segment<3>(nodal_indices[i]*3);
        }
        return cellx;
    }

    template <typename OP>
    void iterateFaceSerial(const OP& f)
    {
        for (int i = 0; i < faces.rows()/3; i++)
            f(i);
    }

    template <typename OP>
    void iterateNodeSerial(const OP& f)
    {
        for (int i = 0; i < deformed.rows()/3; i++)
            f(i);
    }

    template <typename OP>
    void iterateFaceParallel(const OP& f)
    {
        tbb::parallel_for(0, int(faces.rows()/3), [&](int i)
        {
            f(i);
        });
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
    void projectDirichletDoFMatrix(StiffnessMatrix& A, const std::unordered_map<int, T>& data);
    
    // ================================================================

public:
    
    bool advanceOneStep(int step);
    bool advanceOneTimeStep();
    bool linearSolve(StiffnessMatrix& K, const VectorXT& residual, VectorXT& du);
    T computeTotalEnergy();
    T computeResidual(VectorXT& residual);
    void buildSystemMatrix(StiffnessMatrix& K);
    T lineSearchNewton(const VectorXT& residual);

    void computeLinearModes(MatrixXT& eigen_vectors, VectorXT& eigen_values);
    void initializeFromFile(const std::string& filename);
    void buildHingeStructure();
    int nFaces () { return faces.rows() / 3; }
    
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
    virtual void computeConsistentMassMatrix(const FaceVtx& p, Matrix<T, 9, 9>& mass_mat);
    
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

    // Functions for Robust IMLS Contact
    int dim = 3;
    int num_nodes = 0;
    AScalar sigma_r=0.5, sigma_n= 1;
    std::vector<std::vector<int>> current_indices;
    std::vector<std::pair<int,double>> dist_info;
    Eigen::SparseMatrix<double> geodist_close_matrix;
    std::vector<int> mesh_begin_vertices;
    std::vector<int> is_surface_vertex;
    std::vector<std::vector<int>> nodes_adjacency;
    std::vector<std::vector<int>> master_nodes_adjacency;
    std::vector<std::vector<int>> slave_nodes_adjacency;
    std::vector<int> close_slave_nodes;
    std::vector<int> close_master_nodes;
    Tree accelerator;
    AScalar radius = 0.6;
    AScalar IMLS_pen_scale = 1.;
    bool add_friction = false;
    T max_barrier_weight = 1e8;
    T friction_mu = 0.5;
    T epsv_times_h = 1e-5;
    bool self_contact = false;
    double IMLS_param = 1;
    bool use_ipc = false;
    int num_ipc_vtx = 0;
    T barrier_distance = 1e-5;
    T barrier_weight = 1e7;
    T ipc_min_dis = 1e-6;
    bool BARRIER_ENERGY = true;
    double search_radius = 0.5;

    void BuildAcceleratorTreeSC();
    void IMLS_local_gradient_to_global_gradient_and_multiply(
        VectorXa& local_grad,
        Eigen::VectorXi ids,
        int dim1,
        AScalar scale,
        VectorXa& global_vector,
        std::vector<Eigen::Triplet<double>>& triplets);
    
    void IMLS_local_gradient_to_global_gradient(
        VectorXa& local_grad,
        Eigen::VectorXi ids,
        int dim1,
        Eigen::VectorXd& grad);

    void IMLS_local_hessian_to_global_triplets(
        VectorXa& local_hessian,
        Eigen::VectorXi ids,
        int dim1,
        std::vector<Eigen::Triplet<double>>& triplets);

    void IMLS_local_hessian_matrix_to_global_triplets(
        MatrixXa& local_hessian,
        Eigen::VectorXi ids,
        int dim1,
        std::vector<Eigen::Triplet<double>>& triplets);

    void IMLS_vector_muliplication_to_triplets(
        VectorXa& v1,
        VectorXa& v2,
        AScalar scale,
        std::vector<Eigen::Triplet<double>>& triplets);

    void BuildAcceleratorTreeSCTest();
    void addFastRIMLSSCTestEnergy(T& energy);
    void addFastRIMLSSCEnergy(T& energy);
    void addFastRIMLSSCForceEntries(VectorXT& residual);
    void addFastRIMLSSCHessianEntries(std::vector<Entry>& entries);
    T computeCollisionFreeStepsizeUnsigned(const VectorXT& _u, const VectorXT& du);
    double evaluateUnsignedDistanceSqRIMLS(int i, double t, const VectorXa& du, double ipc_stepsize, double& distance);
    double evaluateUnsignedDistanceSqGradientRIMLS(int i, double t,const VectorXa& du, double ipc_stepsize);

    VectorXa divisionGradient(AScalar& fx, AScalar& gx, VectorXa& dfxdx, VectorXa& dgxdx);
    void divisionHessian(AScalar& fx, AScalar& gx, VectorXa& dfxdx, VectorXa& dgxdx, std::vector<Eigen::Triplet<double>>& d2fxdx2, std::vector<Eigen::Triplet<double>>& d2gxdx2, std::vector<Eigen::Triplet<double>>& out);
    void BuildAcceleratorTree(bool is_master = true);

public:
    IMLSContact() 
    {
        updateLameParameters();
    }
    ~IMLSContact() {}
};

// struct RosenbrockCost : public ceres::FirstOrderFunction {
//     int index;
//     Eigen::VectorXd du;
//     double ipc_stepsize;
//     IMLSContact* solver;

//     virtual ~RosenbrockCost() = default;

//     virtual bool Evaluate(const double* parameters,
//                           double* cost,
//                           double* gradient) const override {
//         const double t = parameters[0];
//         double temp;

//          cost[0] = solver->evaluateUnsignedDistanceSqRIMLS(index, t, du, ipc_stepsize, temp);
//         if (gradient) {
//             gradient[0] = solver->evaluateUnsignedDistanceSqGradientRIMLS(index, t, du, ipc_stepsize);
//         }
//         return true;
//     }

//     virtual int NumParameters() const override { return 1; }
// };

#endif