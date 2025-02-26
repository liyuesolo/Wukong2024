#ifndef UTIL_H
#define UTIL_H

#include <Eigen/Geometry>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <fstream>
#include "VecMatDef.h"

inline bool fileExist (const std::string& name) {
    std::ifstream f(name.c_str());
    return f.good();
}

template <typename _type, int _n_col>
void vectorToIGLMatrix(const Matrix<_type, Eigen::Dynamic, 1>& vec, 
    Matrix<_type, Eigen::Dynamic, Eigen::Dynamic>& mat)
{
    int n_rows = vec.rows() / _n_col;
    mat.resize(n_rows, _n_col);
    for (int i = 0; i < n_rows; i++)
        mat.row(i) = vec.template segment<_n_col>(i * _n_col);
}

template <typename _type, int _n_col>
void iglMatrixFatten(const Matrix<_type, Eigen::Dynamic, Eigen::Dynamic>& mat, 
    Matrix<_type, Eigen::Dynamic, 1>& vec)
{
    int n_rows = mat.rows();
    vec.resize(n_rows * _n_col);
    for (int i = 0; i < n_rows; i++)
        vec.template segment<_n_col>(i * _n_col) = mat.row(i);
}
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