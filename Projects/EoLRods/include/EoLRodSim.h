#ifndef EOL_ROD_SIM_H
#define EOL_ROD_SIM_H

#include <utility>
#include <iostream>
#include <Eigen/Geometry>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <tbb/tbb.h>
#include "VecMatDef.h"

#include "RestState.h"
#include "Timer.h"
#include "Rod.h"

#define WARP 0
#define WEFT 1


template <int dim>
struct VectorHash
{
    typedef Vector<int, 3> IV;
    size_t operator()(const IV& a) const{
        std::size_t h = 0;
        for (int d = 0; d < dim; ++d) {
            h ^= std::hash<int>{}(a(d)) + 0x9e3779b9 + (h << 6) + (h >> 2); 
        }
        return h;
    }
};

template <int dim>
struct VectorPairHash
{
    typedef Vector<int, 3> IV;
    size_t operator()(const std::pair<IV, IV>& a) const{
        std::size_t h = 0;
        for (int d = 0; d < dim; ++d) {
            h ^= std::hash<int>{}(a.first(d)) + 0x9e3779b9 + (h << 6) + (h >> 2);
            h ^= std::hash<int>{}(a.second(d)) + 0x9e3779b9 + (h << 6) + (h >> 2);
        }
        return h;
    }
};  

class EoLRodSim
{
public:
    using Simulation = EoLRodSim;
    
    using TV = Vector<T, 3>;

    using TV2 = Vector<T, 2>;
    using TV3 = Vector<T, 3>;
    using TVDOF = Vector<T, 3+2>;
    using TVStack = Matrix<T, 3, Eigen::Dynamic>;
    

    using TM = Matrix<T, 3, 3>;
    using TM3 = Matrix<T, 3, 3>;
    using TMDOF = Matrix<T, 3 + 2, 3 + 2>;

    using TV3Stack = Matrix<T, 3, Eigen::Dynamic>;
    using IV3Stack = Matrix<int, 3, Eigen::Dynamic>;
    using IV4Stack = Matrix<int, 4, Eigen::Dynamic>;
    using DOFStack = Matrix<T, 3 + 2, Eigen::Dynamic>;

    using VectorXT = Matrix<T, Eigen::Dynamic, 1>;

    using IV2 = Vector<int, 2>;
    using IV3 = Vector<int, 3>;
    using IV4 = Vector<int, 4>;
    using IV5 = Vector<int, 5>;
    
    using StiffnessMatrix = Eigen::SparseMatrix<T>;
    using Entry = Eigen::Triplet<T>;

    using Offset = Vector<int, 3 + 1>;
    using Range = Vector<T, 2>;
    using Mask = Vector<bool, 3>;
    using Mask2 = Vector<bool, 2>;

    
    int N_PBC_BENDING_ELE = 5;

    int dof = 3 + 2;
    
    DOFStack q, q0;

    VectorXT deformed_states;
    VectorXT rest_states;

    VectorXT dq;

    VectorXT perturb;

    IV3Stack rods;
    IV4Stack connections;
    TV3Stack normal;
    int n_nodes;
    int n_dof;
    int n_rods;
    int n_pb_cons;
    int incremental_steps = 0;
    int max_newton_iter = 500;
    int ls_max = 12;

    IV2 n_rod_uv;
    
    const static int grid_range = 3;
    int final_3;

    T dt = 1;
    T newton_tol = 1e-4;
    T E = 3.5e9; //PLA
    T R = 0.0002;

    T unit = 0.03;
    T visual_R = 0.01;
    T rho = 1;
    T ks = 1.0;  // stretching term
    T kc = 1e2;  //constraint term
    T kn = 1e-3; //
    T kb = 1.0;
    T kb_penalty = 1.0;
    T km = 1e-3; //mass term
    T kx = 1.0; // shearing term
    T k_pbc = 1.0; // perodic BC term
    T k_strain = 1.0;
    T L = 1;
    T ke = 1e-2; // Eulerian DoF penalty
    T kr = 1e3;
    T k_yc = 1.0;

    std::vector<int> dof_offsets;    

    float theta = 0.f;

    T tunnel_u = R * 4.0;
    T tunnel_v = R * 4.0;

    TV gravity = TV::Zero();
    bool verbose = false;

    bool add_rigid_joint = true;
    bool add_stretching = true;
    bool add_bending = true;
    bool add_shearing = true;
    bool add_twisting = true;
    bool add_penalty = false;
    bool add_rotation_penalty = true;
    bool add_regularizor = false;
    bool add_pbc = true;
    bool add_pbc_bending = true;
    bool add_pbc_twisting = true;
    bool add_eularian_reg = true;
    bool disable_sliding = true;
    bool print_force_mag = false;
    bool add_contact_penalty = true;
    bool run_diff_test = false;

    bool new_frame_work = false;

    TVDOF fix_all, fix_eulerian, fix_lagrangian, fix_u, fix_v;
    std::unordered_map<int, std::pair<TVDOF, TVDOF>> dirichlet_data;
    std::unordered_map<int, T> dirichlet_dof;
    
    std::vector<std::vector<Offset>> pbc_bending_pairs;
    std::vector<std::vector<int>> pbc_bending_pairs_rod_id;

    std::vector<std::vector<int>> pbc_bending_bn_pairs;
    std::vector<std::vector<int>> yarns;
    std::unordered_map<int, int> yarn_map;
    
    std::vector<std::pair<int, std::pair<Offset, Offset>>> pbc_pairs;
    std::unordered_map<int, std::pair<std::pair<Offset, Offset>, std::pair<int, int>>> pbc_pairs_reference;

    // pbc_ref[direction] = (node_i, node_j)
    std::vector<std::pair<int, IV2>> pbc_ref;

    std::vector<int> sliding_nodes;
    IV2 slide_over_n_rods = IV2::Zero();

    

    std::vector<std::pair<std::pair<Offset, Offset>, std::pair<TV, T>>> pbc_strain_data;

    std::vector<std::vector<int>> yarn_group;
    std::vector<bool> is_end_nodes;

    std::vector<RestState*> curvature_functions;

    std::vector<Rod*> Rods;
    std::vector<RodCrossing*> rod_crossings;

    StiffnessMatrix W;

    std::vector<std::pair<TV, T>> boundary_spheres;

    std::function<void(EoLRodSim&, int)> incremental_bc;

    // inverse
    std::vector<TV> targets;

    std::vector<T> residual_norms;
public:

    EoLRodSim()
    {
        gravity[1] = -9.8;
        fix_eulerian.setOnes();
        fix_lagrangian.setOnes();
        fix_all.setOnes();
        fix_u.setZero();
        fix_v.setZero();
        fix_v[dof-1] = 1.0;
        fix_u[dof-2] = 1.0;
        fix_lagrangian.template segment<2>(3).setZero();
        fix_eulerian.template segment<3>(0).setZero();

        config();
    }
    ~EoLRodSim() {}

    void config()
    {
        T area = M_PI * R * R;
        ks = E * area;
        kb = E * area * R * R * 0.25;
        kx = E/T(2)/(1.0 + 0.42) * area;

        std::cout << "ks: " << ks << " kb: " << kb << " kx: " << kx << std::endl;
    }
    
    // TODO: use ... operator
    void cout5Nodes(int n0, int n1, int n2, int n3, int n4)
    {
        std::cout << n0 << " " << n1 << " " << n2 << " " << n3 << " " << n4 << std::endl;
    }

    void cout4Nodes(int n0, int n1, int n2, int n3)
    {
        std::cout << n0 << " " << n1 << " " << n2 << " " << n3 << std::endl;
    }

    void cout3Nodes(int n0, int n1, int n2)
    {
        std::cout << n0 << " " << n1 << " " << n2 << std::endl;
    } 

    template <class OP>
    void iteratePBCBendingPairs(const OP& f) {
        for (int i = 0; i < pbc_bending_pairs.size(); i++)
        {
            f(pbc_bending_pairs[i], pbc_bending_pairs_rod_id[i]);
        }
    }

    template <class OP>
    void iterateSlidingNodes(const OP& f) {
        for (int idx : sliding_nodes){
            f(idx);
        } 
    }

    
    template <class OP>
    void iteratePBCBoundaryPairs(const OP& f) {
        for (auto pair : pbc_bending_bn_pairs){
            std::vector<int> node_ids;
            if(std::find(pair.begin(), pair.end(), -1) != pair.end())
                std::cout << "[EoLRodSim.h] -1 in PBC bending pairs" << std::endl;
            for(int i = 0; i < pair.size() - 1; i++)
                node_ids.push_back(pair[i]);
            f(node_ids, pair[pair.size() - 1]);
        } 
    }

    template <class OP>
    void iterateDirichletData(const OP& f) {
        for (auto dirichlet: dirichlet_data){
            f(dirichlet.first, dirichlet.second.first, dirichlet.second.second);
        } 
    }

    template <class OP>
    void iterateDirichletDoF(const OP& f) {
        for (auto dirichlet: dirichlet_dof){
            f(dirichlet.first, dirichlet.second);
        } 
    }

    template <class OP>
    void iteratePBCReferencePairs(const OP& f) {
        for (auto data : pbc_ref){
            f(data.first, data.second(0), data.second(1));
        } 
    }

    template <class OP>
    void iteratePBCStrainData(const OP& f) {
        for (auto data : pbc_strain_data){
            f(data.first.first, data.first.second, data.second.first, data.second.second);
        } 
    }

    template <class OP>
    void iteratePBCPairs(const OP& f) {
        for (auto data : pbc_pairs){
            f(pbc_pairs_reference[data.first].first.first, 
            pbc_pairs_reference[data.first].first.second,
            data.second.first, data.second.second);
        } 
    }

    template <class OP>
    void iteratePBCPairsWithDirection(const OP& f) {
        for (auto data : pbc_pairs){
            f(data.first, data.second.first, data.second.second);
        } 
    }
    
    template <class OP>
    void iterateAllLagrangianDoFs(const OP& f) {
        for (auto& rod : Rods)
        {
            int cnt = 0;
            for (int i = 0; i < rod->indices.size() - 1; i++)
            {
                if (cnt < rod->dof_node_location.size())
                {
                    if (i == rod->dof_node_location[cnt])
                        cnt++;
                    else
                        f(rod->offset_map[rod->indices[i]]);
                }
                if (cnt == rod->dof_node_location.size())
                {
                    f(rod->offset_map[rod->indices[i]]);
                }
            }
        }
        for (auto& crossing : rod_crossings)
        {
            f(Rods[crossing->rods_involved.front()]->offset_map[crossing->node_idx]);
        }
    }

    T computeTotalEnergy(bool verbose = false);

    T computeResidual(Eigen::Ref<VectorXT> residual);
    
    
    void projectDirichletDoFSystemMatrix(StiffnessMatrix& A);

    void projectDirichletDoFMatrix(StiffnessMatrix& A, const std::unordered_map<int, T>& data);
    
    void buildSystemDoFMatrix(StiffnessMatrix& K);

    bool linearSolve(StiffnessMatrix& K, const VectorXT& residual, VectorXT& du);

    T lineSearchNewton(const VectorXT& residual);

    void checkHessianPD(const VectorXT& dq);

    void computeSmallestEigenVector(const StiffnessMatrix& K, Eigen::Ref<VectorXT> eigen_vector);

    bool staticSolve();

    void staticSolveIncremental(int step);

    bool advanceOneStep(int step);

    void setVerbose(bool v) { verbose = v; }
    
    void resetScene();


    void fixCrossing();

    void releaseCrossing()
    {
        for(auto& crossing : rod_crossings)
        {
            int node_idx = crossing->node_idx;
            std::vector<int> rods_involved = crossing->rods_involved;
            for (int rod_idx : rods_involved)
            {
                Offset offset;
                Rods[rod_idx]->getEntry(node_idx, offset);
                dirichlet_dof.erase(Rods[rod_idx]->reduced_map[offset[3]]);
            }
        }
    }

    bool forward(Eigen::Ref<VectorXT> dq);
    void inverse();

    void fixRegion(std::function<bool(const TV&)> inside_region);
    void fixRegionAvoidRod(std::function<bool(const TV&)> inside_region, int rod_idx);
    void fixRegionalDisplacement(std::function<bool(const TV&, TV&, Vector<bool, 3>&)> helper);

    
    void getCrossingPosition(int crossing_idx, TV& x) const
    {
        auto crossing = rod_crossings[crossing_idx];
        auto rod = Rods[crossing->rods_involved.front()];
        rod->x(crossing->node_idx, x);
    }

    void fixCrossingLagrangian(int crossing_idx, TV delta, Mask mask);

    void computeBoundingBox(TV& bottom_left, TV& top_right) const;
    
private:

    
public:
    // Elasticity.cpp
    void setUniaxialStrain(T theta, T s, TV& strain_dir, TV& ortho_dir);
    void setBiaxialStrain(T theta1, T s1, T theta2, T s2, TV& strain_dir, TV& ortho_dir);
    void setBiaxialStrainWeighted(T theta1, T s1, T theta2, T s2, T w);
    // void computeMacroStress(TM& sigma, TV strain_dir);
    void computeDeformationGradientUnitCell();
    void fitDeformationGradientUnitCell();
    

    // Scene.cpp 
    void buildSceneFromUnitPatch(int patch_id);
    void buildPeriodicNetwork(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& C, bool show_rest);
    
    //Visualization.cpp
    void getColorPerYarn(Eigen::MatrixXd& C, int n_rod_per_yarn = 4);
    void getEulerianDisplacement(Eigen::MatrixXd& X, Eigen::MatrixXd& x);
    void showStretching(Eigen::MatrixXd& C);
    void buildMeshFromRodNetwork(Eigen::MatrixXd& V, Eigen::MatrixXi& F, 
        Eigen::Ref<const DOFStack> q_display, Eigen::Ref<const IV3Stack> rods_display,
        Eigen::Ref<const TV3Stack> normal_tile);
    void generateMeshForRendering(Eigen::MatrixXd& V, Eigen::MatrixXi& F, TV shift = TV::Zero(), bool show_rest_shape = true);
    void markSlidingRange(int idx, int dir, int depth, std::vector<bool>& can_slide, int root);
    
    // DerivativeTest.cpp
    
    void derivativeTest();
    void testGradient2ndOrderTerm(Eigen::Ref<VectorXT> dq);
    void testHessian2ndOrderTerm(Eigen::Ref<VectorXT> dq);
    void testGradient(Eigen::Ref<VectorXT> dq);
    void testHessian(Eigen::Ref<VectorXT> dq);

    void checkMaterialPositionDerivatives();


    // ======================== Energy Forces and Hessian Entries ========================
    //                                             so -df/dx
    // Joint.cpp
    void addJointBendingAndTwistingK(std::vector<Entry>& entry_K);
    void addJointBendingAndTwistingForce(Eigen::Ref<VectorXT> residual);
    T addJointBendingAndTwistingEnergy(bool bending = true, bool twisting = true);

    //PBCBendingAndTwisting.cpp
    void buildPBCData(
        std::vector<TV>& data, 
        const std::vector<Offset>& offsets, 
        const std::vector<int>& rod_id,
        std::vector<TV>& dXdu, std::vector<TV>& d2Xdu2,
        bool g = false, bool f = false);

    T add3DPBCBendingAndTwistingEnergy(bool bending = true, bool twisting = true);
    void add3DPBCBendingAndTwistingForce(Eigen::Ref<VectorXT> residual, bool bending = true, bool twisting = true);
    void add3DPBCBendingAndTwistingK(std::vector<Entry>& entry_K, bool bending = true, bool twisting = true);

    // Bending.cpp
    // void addBendingK(std::vector<Eigen::Triplet<T>>& entry_K);  
    // void addBendingForce(Eigen::Ref<VectorXT> residual);
    // T addBendingEnergy();

    // Stretching.cpp
    void addStretchingK(std::vector<Entry>& entry_K);  
    void addStretchingForce(Eigen::Ref<VectorXT> residual);
    T addStretchingEnergy();

    // Shearing.cpp
    // void addShearingK(Eigen::Ref<const DOFStack> q_temp, std::vector<Eigen::Triplet<T>>& entry_K, bool top_right);  
    // void addShearingForce(Eigen::Ref<const DOFStack> q_temp, Eigen::Ref<DOFStack> residual, bool top_right);
    // T addShearingEnergy(Eigen::Ref<const DOFStack> q_temp, bool top_right);
    // void toMapleNodesVector(std::vector<Vector<T, 3>>& x, Eigen::Ref<const DOFStack> q_temp,
    //     std::vector<int>& nodes);

    //PeriodicBC.cpp
    T addPBCEnergy();
    void addPBCForce(Eigen::Ref<VectorXT> residual);
    void addPBCK(std::vector<Eigen::Triplet<T>>& entry_K);  
    
    void buildRotationPenaltyData(
        std::vector<TV2>& data, std::vector<Offset>& offsets,
        std::vector<TV2>& dXdu, std::vector<TV2>& d2Xdu2);


    // EulerianConstraints.cpp    
    void addRegularizingK(std::vector<Eigen::Triplet<T>>& entry_K);
    T addRegularizingEnergy();
    void addRegularizingForce(Eigen::Ref<VectorXT> residual);
    

    // ParallelContact.cpp
    T addParallelContactEnergy();
    void addParallelContactForce(Eigen::Ref<VectorXT> residual);
    void addParallelContactK(std::vector<Entry>& entry_K);
    void parallelContactdfdp(Eigen::Ref<VectorXT> dfdp);

    T add3DBendingAndTwistingEnergy(bool bending = true, bool twisting = true);
    void add3DBendingAndTwistingForce(Eigen::Ref<VectorXT> residual);
    void add3DBendingAndTwistingK(std::vector<Entry>& entry_K);

};

#endif