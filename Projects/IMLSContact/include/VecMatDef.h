#ifndef VECMATDEF_H
#define VECMATDEF_H

#include <Eigen/Core>
#include <vector>
#include <Eigen/Eigen>

typedef double AScalar;
typedef Eigen::Matrix<AScalar, 1, 1> Vector1a;
typedef Eigen::Matrix<AScalar, 2, 1> Vector2a;
typedef Eigen::Matrix<AScalar, 3, 1> Vector3a;
typedef Eigen::Matrix<AScalar, 4, 1> Vector4a;
typedef Eigen::Matrix<AScalar, 6, 1> Vector6a;
typedef Eigen::Matrix<AScalar, 8, 1> Vector8a;
typedef Eigen::Matrix<AScalar, 9, 1> Vector9a;
typedef Eigen::Matrix<AScalar, 10, 1> Vector10a;
typedef Eigen::Matrix<AScalar, 12, 1> Vector12a;
typedef Eigen::Matrix<AScalar, 13, 1> Vector13a;
typedef Eigen::Matrix<AScalar, 14, 1> Vector14a;
typedef Eigen::Matrix<AScalar, 15, 1> Vector15a;
typedef Eigen::Matrix<AScalar, 24, 1> Vector18a;
typedef Eigen::Matrix<AScalar, 24, 1> Vector24a;
typedef Eigen::Matrix<AScalar, Eigen::Dynamic, 1> VectorXa;

typedef Eigen::Matrix<AScalar, 2, 2> Matrix2a;
typedef Eigen::Matrix<AScalar, 3, 3> Matrix3a;
typedef Eigen::Matrix<AScalar, 4, 4> Matrix4a;
typedef Eigen::Matrix<AScalar, 6, 6> Matrix6a;
typedef Eigen::Matrix<AScalar, 8, 8> Matrix8a;
typedef Eigen::Matrix<AScalar, 9, 9> Matrix9a;
typedef Eigen::Matrix<AScalar, 10, 10> Matrix10a;
typedef Eigen::Matrix<AScalar, 12, 12> Matrix12a;
typedef Eigen::Matrix<AScalar, 13, 13> Matrix13a;
typedef Eigen::Matrix<AScalar, 14, 14> Matrix14a;
typedef Eigen::Matrix<AScalar, 15, 15> Matrix15a;
typedef Eigen::Matrix<AScalar, 24, 24> Matrix24a;
//typedef Eigen::Matrix<AScalar, 6, 4> Matrix6x4a;

typedef Eigen::Matrix<AScalar, Eigen::Dynamic, Eigen::Dynamic> MatrixXa;

typedef std::vector<Vector2a, Eigen::aligned_allocator<Vector2a> > Vector2v;
typedef std::vector<Vector3a, Eigen::aligned_allocator<Vector3a> > Vector3v;
typedef std::vector<Vector4a, Eigen::aligned_allocator<Vector4a> > Vector4v;

template <typename T, int dim>
using Vector = Eigen::Matrix<T, dim, 1, 0, dim, 1>;

template <typename T, int n, int m>
using Matrix = Eigen::Matrix<T, n, m, 0, n, m>;


using T = double;


#endif