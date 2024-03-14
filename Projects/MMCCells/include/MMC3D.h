#ifndef MMC3D_H
#define MMC3D_H

#include <utility>
#include <iostream>
#include <fstream>
#include <Eigen/Geometry>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <tbb/tbb.h>
#include "VecMatDef.h"

class MMC3D
{
public:
    using VectorXT = Matrix<T, Eigen::Dynamic, 1>;
    using MatrixXT = Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using VectorXi = Vector<int, Eigen::Dynamic>;
    using MatrixXi = Matrix<int, Eigen::Dynamic, Eigen::Dynamic>;
    using StiffnessMatrix = Eigen::SparseMatrix<T>;
    using Entry = Eigen::Triplet<T>;
    using TV = Vector<T, 3>;
    using TM = Matrix<T, 3, 3>; 
    
    VectorXT deformed, undeformed, u;
    int n_components = 1;
    int n_dof_component = 9;
    int p = 6;
    T lambda = 100;
public:

    T levelsetValue(const TV& point)
    {
        T value = 0.0;
        for (int i = 0; i < n_components; i++)
        {
            VectorXT dof_single_component = deformed.segment(i * n_dof_component, n_dof_component);
            T x = dof_single_component[0];
            T y = dof_single_component[1];
            T z = dof_single_component[2];
            T alpha = dof_single_component[3];
            T beta = dof_single_component[4];
            T gamma = dof_single_component[5];
            T L1 = dof_single_component[6];
            T L2 = dof_single_component[7];
            T L3 = dof_single_component[8];

            T ca = std::cos(alpha), sa = std::sin(alpha);
            T cb = std::cos(beta), sb = std::sin(beta);
            T cg = std::cos(gamma), sg = std::sin(gamma);

            TM R;
            R << cb * cg, cb * sg, -sb,
                sa * sb * cg - ca * sg, sa * sb * sg + ca * cg, sa * sb,
                    ca * sb * cg + sa * sg, ca * sb * sg - sa * cg, ca * cb;

            // std::cout << R << std::endl;
            TV prime = R * ( point - TV(x, y, z));

            value += std::exp(lambda * (1.0 - 
                std::pow(std::pow(prime[0]/L1, p) + std::pow(prime[1]/L2, p) + std::pow(prime[2]/L3, p), 1.0/T(p))));
            
        }

        return std::log(value) / lambda;
    }

    void initialize()
    {
        n_components = 2;
        n_dof_component = 9;
        deformed.resize(n_components * n_dof_component);
        deformed.setZero();

        deformed[6] = 1.0;
        deformed[7] = 3.0;
        deformed[8] = 1.0;

        deformed[n_dof_component + 0] = 3.0;

        deformed[n_dof_component + 6] = 1.0;
        deformed[n_dof_component + 7] = 3.0;
        deformed[n_dof_component + 8] = 1.0;
    }
    
public:
    MMC3D() {}
    ~MMC3D() {}
};


#endif