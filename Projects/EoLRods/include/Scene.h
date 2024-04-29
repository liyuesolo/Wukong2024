#ifndef SCENE_H
#define SCENE_H

#include "EoLRodSim.h"

class EoLRodSim;

class Scene
{
public:
    using TV = Vector<T, 3>;
    using TV2 = Vector<T, 2>;

    using VectorXT = Matrix<T, Eigen::Dynamic, 1>;
    using StiffnessMatrix = Eigen::SparseMatrix<T>;
    using Entry = Eigen::Triplet<T>;

    using Offset = Vector<int, 3 + 1>;
    using Range = Vector<T, 2>;
    using Mask = Vector<bool, 3>;
    using Mask2 = Vector<bool, 2>;

private:
    EoLRodSim& sim;

    VectorXT& deformed_states = sim.deformed_states;
    
public:
    Scene(EoLRodSim& eol_sim) : sim(eol_sim) {}
    ~Scene() {}
    
    // ------------------------------- Scene Setup -------------------------------
    void buildInterlockingSquareScene(int sub_div);
    void buildStraightRodScene(int sub_div);
    void buildGridScene(int sub_div);
    void buildFullScaleSquareScene(int sub_div);

private:

    // ------------------------------- Common Function -------------------------------
    void clearSimData();

    void appendThetaAndJointDoF(std::vector<Entry>& w_entry, 
        int& full_dof_cnt, int& dof_cnt);   
    
    void addStraightYarnCrossNPoints(const TV& from, const TV& to,
        const std::vector<TV>& passing_points, 
        const std::vector<int>& passing_points_id, 
        int sub_div,
        std::vector<TV>& sub_points, std::vector<int>& node_idx,
        std::vector<int>& key_points_location,
        int start, bool pbc = false);

    void markCrossingDoF(
        std::vector<Eigen::Triplet<T>>& w_entry,
        int& dof_cnt);

    void addAStraightRod(const TV& from, const TV& to, 
        const std::vector<TV>& passing_points, 
        const std::vector<int>& passing_points_id, 
        int sub_div,
        int& full_dof_cnt, int& node_cnt, int& rod_cnt);
    
    void addPoint(const TV& point, int& full_dof_cnt, int& node_cnt);

    void addCrossingPoint(std::vector<TV>& existing_nodes, 
        const TV& point, int& full_dof_cnt, int& node_cnt);
};

#endif