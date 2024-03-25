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

    
    VectorXT valueBatch(const VectorXT& input, int batch_dim)
    {
        int data_dim = input.rows() / batch_dim;
        std::vector<double> nn_input(input.rows());
        for (int i = 0; i < input.rows(); i++)
            nn_input[i] = input[i];
        cppflow::tensor input_tensor(nn_input, {batch_dim, data_dim});
        auto output = neural_model(
            {
                {"serving_default_input_1:0", input_tensor}
            },
            {
                "StatefulPartitionedCall:0", 
                "StatefulPartitionedCall:1",
                "StatefulPartitionedCall:2"
            }
        );
        auto data = output[0].get_data<T>();
        VectorXT value_batch = Eigen::Map<VectorXT>(data.data(), data.size());
        return value_batch;
    }

    
    VectorXT gradBatch(const VectorXT& input, int batch_dim)
    {
        int data_dim = input.rows() / batch_dim;
        std::vector<double> nn_input(input.rows());
        for (int i = 0; i < input.rows(); i++)
            nn_input[i] = input[i];
        cppflow::tensor input_tensor(nn_input, {batch_dim, data_dim});
        auto output = neural_model(
            {
                {"serving_default_input_1:0", input_tensor}
            },
            {
                "StatefulPartitionedCall:0", 
                "StatefulPartitionedCall:1",
                "StatefulPartitionedCall:2"
            }
        );
        auto data = output[1].get_data<T>();
        VectorXT grad_batch = Eigen::Map<VectorXT>(data.data(), data.size());
        return grad_batch;
    }

    
    VectorXT hessBatch(const VectorXT& input, int batch_dim)
    {
        int data_dim = input.rows() / batch_dim;
        std::vector<double> nn_input(input.rows());
        for (int i = 0; i < input.rows(); i++)
            nn_input[i] = input[i];
        cppflow::tensor input_tensor(nn_input, {batch_dim, data_dim});
        auto output = neural_model(
            {
                {"serving_default_input_1:0", input_tensor}
            },
            {
                "StatefulPartitionedCall:0", 
                "StatefulPartitionedCall:1",
                "StatefulPartitionedCall:2"
            }
        );
        auto data = output[2].get_data<T>();
        VectorXT hess_batch = Eigen::Map<VectorXT>(data.data(), data.size());
        
        return hess_batch;
    }

    T value(const VectorXT& input)
    {
        T _value = 0.0;
        std::vector<double> nn_input(input.rows());
        for (int i = 0; i < input.rows(); i++)
            nn_input[i] = input[i];
        cppflow::tensor input_tensor(nn_input, {1, input.rows()});
        auto output = neural_model(
            {
                {"serving_default_input_1:0", input_tensor}
            },
            {
                "StatefulPartitionedCall:0", 
                "StatefulPartitionedCall:1",
                "StatefulPartitionedCall:2"
            }
        );
        _value = output[0].get_data<T>()[0];
        
        return _value;
    }


    Vector<T, 9> grad(const VectorXT& input)
    {
        std::vector<double> nn_input(input.rows());
        for (int i = 0; i < input.rows(); i++)
            nn_input[i] = input[i];
        cppflow::tensor input_tensor(nn_input, {1, input.rows()});
        auto output = neural_model(
            {
                {"serving_default_input_1:0", input_tensor}
            },
            {
                "StatefulPartitionedCall:0", 
                "StatefulPartitionedCall:1",
                "StatefulPartitionedCall:2"
            }
        );
        
        auto data = output[1].get_data<T>();
        VectorXT _grad = Eigen::Map<VectorXT>(data.data(), data.size());
        return _grad;
    }

    
    Matrix<T, 9, 9> hess(const VectorXT& input)
    {
        std::vector<double> nn_input(input.rows());
        for (int i = 0; i < input.rows(); i++)
            nn_input[i] = input[i];
        cppflow::tensor input_tensor(nn_input, {1, input.rows()});
        auto output = neural_model(
            {
                {"serving_default_input_1:0", input_tensor}
            },
            {
                "StatefulPartitionedCall:0", 
                "StatefulPartitionedCall:1",
                "StatefulPartitionedCall:2"
            }
        );
        // this gives a flattened hessian matrix
        auto data = output[2].get_data<T>();
        Matrix<T, 9, 9> _hess = Eigen::Map<Eigen::Matrix<T, 9, 9>,                                                                 
                                        Eigen::Unaligned, 
                                        Eigen::Stride<1, 9>>(data.data());
        
        return _hess;
    }
    
public:
    NeuralMaterialModel(cppflow::model& _model) : neural_model(_model) {}
    ~NeuralMaterialModel() {}
};

#endif