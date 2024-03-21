#ifndef NEURAL_MATERAL_MODEL_H
#define NEURAL_MATERAL_MODEL_H

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

#include <cppflow/cppflow.h>

class NeuralMaterialModel
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

    cppflow::model& neural_model;

public:
    void queryNetworkDerivatives();

    template<int dim>
    VectorXT valueBatch(const Vector<T, dim>& input, int batch_dim = 1)
    {
        VectorXT value_batch = VectorXT::Zero(batch_dim);
        return value_batch;
    }

    template<int dim>
    VectorXT gradBatch(const Vector<T, dim>& input, int batch_dim = 1)
    {
        VectorXT grad_batch = VectorXT::Zero(batch_dim * dim);
        return grad_batch;
    }

    template<int dim>
    VectorXT hessBatch(const Vector<T, dim>& input, int batch_dim = 1)
    {
        VectorXT hess_batch = VectorXT::Zero(batch_dim * dim * dim);
        return hess_batch;
    }
    
public:
    NeuralMaterialModel(cppflow::model& _model) : neural_model(_model) {}
    ~NeuralMaterialModel() {}
};

#endif