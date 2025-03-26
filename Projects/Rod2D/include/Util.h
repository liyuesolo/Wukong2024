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

template <int dim = 2>
void addForceEntry(Eigen::VectorXd& residual, const std::vector<int>& vtx_idx,
                    const Eigen::VectorXd& gradent, int shift = 0)
{
    for (int i = 0; i < vtx_idx.size(); i++)
        residual.template segment<dim>(vtx_idx[i] * dim + shift) +=
            gradent.template segment<dim>(i * dim);
}

template <int dim = 2>
void getSubVector(const Eigen::VectorXd& _vector, const std::vector<int>& vtx_idx,
                    Eigen::VectorXd& sub_vec, int shift = 0)
{
    sub_vec.resize(vtx_idx.size() * dim);
    sub_vec.setZero();
    for (int i = 0; i < vtx_idx.size(); i++)
    {
        sub_vec.template segment<dim>(i * dim) =
            _vector.template segment<dim>(vtx_idx[i] * dim + shift);
    }
}

template <int dim_row = 2, int dim_col = 2>
void addHessianEntry(std::vector<Eigen::Triplet<T>>& triplets,
                        const std::vector<int>& vtx_idx,
                        const Eigen::MatrixXd& hessian, int shift_row = 0,
                        int shift_col = 0)
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
                        hessian(i * dim_row + k, j * dim_col + l));
        }
    }
}

template <int dim_row = 2, int dim_col = 2>
void addJacobianEntry(std::vector<Eigen::Triplet<T>>& triplets,
                        const std::vector<int>& vtx_idx,
                        const std::vector<int>& vtx_idx2,
                        const Eigen::MatrixXd& jacobian, int shift_row = 0,
                        int shift_col = 0)
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
                        jacobian(i * dim_row + k, j * dim_col + l));
        }
    }
}

template <int dim_row = 2, int dim_col = 2>
void addHessianMatrixEntry(Eigen::MatrixXd& matrix_global,
                            const std::vector<int>& vtx_idx,
                            const Eigen::MatrixXd& hessian, int shift_row = 0,
                            int shift_col = 0)
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
                                    dof_j * dim_col + l + shift_col) +=
                        hessian(i * dim_row + k, j * dim_col + l);
                }
        }
    }
}

template <int dim_row = 2, int dim_col = 2>
void addJacobianMatrixEntry(Eigen::MatrixXd& matrix_global,
                            const std::vector<int>& vtx_idx,
                            const std::vector<int>& vtx_idx2,
                            const Eigen::MatrixXd& hessian, int shift_row = 0,
                            int shift_col = 0)
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
                                    dof_j * dim_col + l + shift_col) +=
                        hessian(i * dim_row + k, j * dim_col + l);
                }
        }
    }
}

template <int dim0 = 2, int dim1 = 2>
void addHessianBlock(std::vector<Eigen::Triplet<T>>& triplets,
                        const std::vector<int>& vtx_idx,
                        const Eigen::MatrixXd& hessian_block, int shift_row = 0,
                        int shift_col = 0)
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

template <int size>
Eigen::VectorXd computeHessianBlockEigenValues(const Matrix<T, size, size>& symMtr)
{
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T, size, size>> eigenSolver(
        symMtr);
    return eigenSolver.eigenvalues();
}



#endif