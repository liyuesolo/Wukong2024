#ifndef ROD_H
#define ROD_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/Sparse>

#include <iostream>
#include <tbb/tbb.h>

#include "Util.h"
#include "VecMatDef.h"

class Rod;

struct RodCrossing
{
    using TV = Vector<T, 3>;
    using Mask = Vector<bool, 3>;

    int node_idx;
    bool is_fixed;

    std::vector<int> rods_involved;
    std::vector<Vector<T, 2>> undeformed_twist;

    std::unordered_map<int, int> on_rod_idx;
    Vector<T, 3> omega;
    Vector<T, 4> omega_acc;

    Matrix<T, 3, 3> rotation_accumulated;

    int dof_offset;
    int reduced_dof_offset;

    void updateRotation(const Vector<T, 3>& new_omega)
    {
        Matrix<T, 3, 3> R = rotationMatrixFromEulerAngle(
            new_omega[2], new_omega[1], new_omega[0]);
        rotation_accumulated = R * rotation_accumulated;
        omega.setZero();
    }

    RodCrossing(int id, std::vector<int> involved)
        : node_idx(id), rods_involved(involved), is_fixed(false)
    {
        omega_acc.setZero();
        omega_acc[3] = 1.0;
        omega.setZero();
        rotation_accumulated.setIdentity();
    }
};

class Rod
{

public:
    using Entry = Eigen::Triplet<T>;
    using TV = Vector<T, 3>;
    using IV = Vector<int, 3>;
    using VectorXT = Matrix<T, Eigen::Dynamic, 1>;
    using VectorXi = Eigen::VectorXi;
    using Offset = Vector<int, 3>;
    using Mask = Vector<bool, 3>;

    T B[2][2];

    Matrix<T, 2, 2> bending_coeffs;

    // unit m
    T a = 1e-3, b = 1e-3;

    // T E = 3.5e9;
    T E = 1e5;
    T ks;
    T kt;

    int rod_id;
    bool closed;
    int theta_dof_start_offset = 0;
    int theta_reduced_dof_start_offset = 0;
    // full states for all nodes in the system
    VectorXT& full_states;
    VectorXT& rest_states;

    std::vector<TV> reference_frame_us;

    std::vector<TV> rest_tangents;
    std::vector<TV> rest_normals;

    // previous tangent back after for time parallel transport
    //  -> for each Newton iterationw
    std::vector<TV> prev_tangents;

    VectorXT reference_angles;

    VectorXT reference_twist;

    // which node on this rod has been marked as crossing node, ranging from 0
    // to # nodes
    std::vector<int> dof_node_location;

    std::vector<bool> fixed_by_crossing;

    // maps local nodal idx to global dof offset
    std::unordered_map<int, IV> offset_map;

    // reduced map maps global dof to actual dofs
    std::unordered_map<int, int> reduced_map;

    // global indices of nodes on this rod
    std::vector<int> indices;

public:
    // Rod.cpp
    void setupBishopFrame();
    T computeReferenceTwist(const TV& tangent, const TV& prev_tangent,
                            int rod_idx);
    void curvatureBinormal(const TV& t1, const TV& t2, TV& kb);
    void rotateReferenceFrameToLastNewtonStepAndComputeReferenceTwsit();

public:
    // ============================== iterators
    // ===================================

    template <class OP>
    void iterateSegments(const OP& f)
    {
        for (int i = 0; i < indices.size() - 1; i++)
        {
            f(indices[i], indices[i + 1], i);
        }
    }

    template <class OP>
    void iterateSegmentsWithOffset(const OP& f)
    {
        for (int i = 0; i < indices.size() - 1; i++)
        {
            f(indices[i], indices[i + 1], offset_map[indices[i]],
              offset_map[indices[i + 1]], i);
        }
    }

    template <class OP>
    void iterate3NodesWithOffsets(const OP& f)
    {
        int cnt = 0;
        for (int i = 1; i < indices.size() - 1; i++)
        {
            bool is_crossing = false;
            if (cnt < dof_node_location.size())
            {
                if (i ==
                    dof_node_location[cnt]) // || i-1 == dof_node_location[cnt]
                                            // || i + 1 ==
                                            // dof_node_location[cnt])
                {
                    is_crossing = true;
                    cnt++;
                }
            }
            f(indices[i], indices[i + 1], indices[i - 1],
              offset_map[indices[i]], offset_map[indices[i + 1]],
              offset_map[indices[i - 1]], i, is_crossing);
        }
        if (closed)
        {
            bool is_crossing = false;
            if (dof_node_location.size())
                if (dof_node_location[0] ==
                    0) // || dof_node_location[0] == 1 || dof_node_location[0]
                       // == indices.size() - 2)
                    is_crossing = true;
            f(indices.front(), indices[1], indices[indices.size() - 2],
              offset_map[indices.front()], offset_map[indices[1]],
              offset_map[indices[indices.size() - 2]], 0, is_crossing);
        }
    }

    template <class OP>
    void iterate3Nodes(const OP& f)
    {
        int cnt = 0;
        for (int i = 1; i < indices.size() - 1; i++)
            f(indices[i], indices[i + 1], indices[i - 1], i);
    }

    // ============================== Lagrangian Eulerian value helpers
    // =================================== deformed Lagrangian position
    void x(int node_idx, TV& pos)
    {
        pos = full_states.segment<3>(node_idx * 3);
        // pos = TV::Zero();
        // Offset idx = offset_map[node_idx];
        // for (int d = 0; d < 3; d++)
        // {
        //     pos[d] = full_states[idx[d]];
        // }
    }

    // undeformed Lagrangian position
    void X(int node_idx, TV& pos)
    {
        pos = rest_states.segment<3>(node_idx * 3);

        // Offset idx = offset_map[node_idx];
        // for (int d = 0; d < 3; d++)
        // {
        //     pos[d] = rest_states[idx[d]];
        // }
    }

    // ============================== helpers
    // ===================================
    std::pair<int, int> neighbouringCrossingIndex(int current_crossing_location)
    {
        auto iter =
            std::find(dof_node_location.begin(), dof_node_location.end(),
                      current_crossing_location);
        if (iter == dof_node_location.end())
        {
            std::cout << "[Rod.h] invalid crossing node location on rod"
                      << std::endl;
            std::exit(0);
        }
        int location = std::distance(dof_node_location.begin(), iter);
        if (location == 0)
        {
            int left = indices.front();
            int right = -1;
            if (location == dof_node_location.size() - 1)
                right = indices.back();
            else
                right = indices[dof_node_location[location + 1]];
            return std::make_pair(left, right);
        }
        else if (location == dof_node_location.size() - 1)
        {
            int right = indices.back();
            int left = -1;
            if (location == 0)
                left = indices.front();
            else
                left = indices[dof_node_location[location - 1]];
            return std::make_pair(left, right);
        }
        else
        {
            int left = indices[dof_node_location[location - 1]];
            int right = indices[dof_node_location[location + 1]];
            return std::make_pair(left, right);
        }
    }

    int entry(int node_idx) { return 0; }

    int nodeIdx(int node_pos)
    {
        if (node_pos == -1)
        {
            if (closed)
                return indices[indices.size() - 2];
            else
                return -1;
        }

        else if (node_pos == indices.size())
        {
            if (closed)
                return indices.front();
            else
                return -1;
        }

        else
            return indices[node_pos];
    }

    void getEntry(int node_idx, Offset& idx) { idx = offset_map[node_idx]; }

    void getEntryReduced(int node_idx, Offset& idx)
    {
        idx = offset_map[node_idx];
        for (int d = 0; d < 3; d++)
        {
            idx[d] = reduced_map[idx[d]];
        }
    }

    void getEntryByLocation(int node_location, Offset& idx)
    {
        idx = offset_map[indices[node_location]];
    }

    void frontOffset(Offset& idx) { idx = offset_map[indices.front()]; }

    void backOffset(Offset& idx) { idx = offset_map[indices.back()]; }

    void backOffsetReduced(Offset& idx)
    {
        idx = offset_map[indices.back()];
        for (int d = 0; d < 3 + 1; d++)
            idx[d] = reduced_map[idx[d]];
    }

    void frontOffsetReduced(Offset& idx)
    {
        idx = offset_map[indices.front()];
        for (int d = 0; d < 3 + 1; d++)
            idx[d] = reduced_map[idx[d]];
    }

    void frontDoF(Vector<T, 3 + 1>& q)
    {
        q = Vector<T, 3 + 1>::Zero();
        Offset idx = offset_map[indices.front()];
        for (int d = 0; d < 3 + 1; d++)
            q[d] = full_states[idx[d]];
    }

    void backDoF(Vector<T, 3 + 1>& q)
    {
        q = Vector<T, 3 + 1>::Zero();
        Offset idx = offset_map[indices.back()];
        for (int d = 0; d < 3 + 1; d++)
            q[d] = full_states[idx[d]];
    }

    void fixPointLagrangian(int node_idx, TV delta,
                            std::unordered_map<int, T>& dirichlet_data)
    {
        for (int d = 0; d < 3; d++)
            dirichlet_data[reduced_map[offset_map[indices[node_idx]][d]]] =
                delta[d];
    }

    void fixPointLagrangianByID(int node_idx, TV delta, Mask mask,
                                std::unordered_map<int, T>& dirichlet_data)
    {
        for (int d = 0; d < 3; d++)
            if (mask[d])
                dirichlet_data[reduced_map[offset_map[node_idx][d]]] = delta[d];
    }

    void fixEndPointLagrangian(std::unordered_map<int, T>& dirichlet_data)
    {
        for (int d = 0; d < 3; d++)
        {
            dirichlet_data[reduced_map[offset_map[indices.front()][d]]] = 0;
            dirichlet_data[reduced_map[offset_map[indices.back()][d]]] = 0;
        }
    }

    int numSeg() { return indices.size() - 1; }

    void markDoF(std::vector<Entry>& w_entry, int& dof_cnt);

    void initCoeffs()
    {
        T coeff = 0.25 * E * M_PI * a * b;
        B[0][0] = coeff * a * a;
        B[0][1] = 0;
        B[1][0] = 0;
        B[1][1] = coeff * b * b;

        // bending_coeffs << coeff * a * a, 0, 0, coeff * b*b;
        bending_coeffs.setZero();
        bending_coeffs(0, 0) = coeff * a * a;
        bending_coeffs(1, 1) = coeff * b * b;

        ks = E * M_PI * a * b;

        T G = E / T(2) / (1.0 + 0.42);

        // T G = 1e7;

        kt = 0.25 * G * M_PI * a * b * (a * a + b * b);
    }

    bool isFixedNodeForPrinting(int node_idx, int rod_idx);

public:
    Rod(VectorXT& q, VectorXT& q0, int id, T _a, T _b)
        : full_states(q), rest_states(q0), rod_id(id), a(_a), b(_b),
          closed(false)
    {
        initCoeffs();
    }
    Rod(VectorXT& q, VectorXT& q0, int id, bool c, T _a, T _b)
        : full_states(q), rest_states(q0), rod_id(id), closed(c), a(_a), b(_b)
    {
        initCoeffs();
    }
    ~Rod() {}

    // private:
};
#endif