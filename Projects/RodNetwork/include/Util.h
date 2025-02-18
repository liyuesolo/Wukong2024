#ifndef UTIL_H
#define UTIL_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include "VecMatDef.h"
template <class VectorType>
inline void appendToEnd(std::vector<VectorType>& a,
                        const std::vector<VectorType>& b)
{
    a.insert(a.end(), b.begin(), b.end());
}

inline T signedAngle(const Eigen::Vector3d& u, const Eigen::Vector3d& v,
                     const Eigen::Vector3d& n)
{
    Eigen::Vector3d w = u.cross(v);
    T angle = std::atan2(w.norm(), u.dot(v));
    if (n.dot(w) < 0)
        return -angle;
    return angle;
}

inline bool colinear(Eigen::Vector3d a, Eigen::Vector3d b)
{
    if ((a - b).norm() < 1e-2)
        return true;
    if ((a - b).norm() > 1.99)
        return true;
    return false;
}

inline void rotateAxisAngle(Eigen::Vector3d& v, const Eigen::Vector3d& z,
                            const T theta)
{

    if (theta == 0)
        return;

    T c = cos(theta);
    T s = sin(theta);

    v = c * v + s * z.cross(v) + z.dot(v) * (1.0 - c) * z;
}

Eigen::Matrix3d rotationMatrixFromEulerAngle(T angle_z, T angle_y, T angle_x);

Eigen::Vector3d parallelTransport(const Eigen::Vector3d& u,
                                  const Eigen::Vector3d& t0,
                                  const Eigen::Vector3d& t1);
Eigen::Vector3d parallelTransportOrthonormalVector(const Eigen::Vector3d& u,
                                                   const Eigen::Vector3d& t0,
                                                   const Eigen::Vector3d& t1);

bool circleCircleIntersection(const Eigen::Vector3d& x0, T r0,
                              const Eigen::Vector3d& x1, T r1,
                              Eigen::Vector3d& ixn0, Eigen::Vector3d& ixn1);

template <int size>
void projectBlockPD(Eigen::Matrix<T, size, size>& symMtr)
{
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T, size, size>> eigenSolver(
        symMtr);
    if (eigenSolver.eigenvalues()[0] >= 0.0)
    {
        return;
    }
    Eigen::DiagonalMatrix<T, size> D(eigenSolver.eigenvalues());
    int rows = ((size == Eigen::Dynamic) ? symMtr.rows() : size);
    for (int i = 0; i < rows; i++)
    {
        if (D.diagonal()[i] < 0.0)
        {
            D.diagonal()[i] = 0.0;
        }
        else
        {
            break;
        }
    }
    symMtr =
        eigenSolver.eigenvectors() * D * eigenSolver.eigenvectors().transpose();
}
#endif