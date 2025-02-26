#ifndef AUTODIFF_H
#define AUTODIFF_H


#include <utility> 
#include <iostream>
#include <fstream>
#include <Eigen/Geometry>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <tbb/tbb.h>
#include <unordered_map>
#include <iomanip>

#include <cppad/cg.hpp>
#include <cppad/cg/support/cppadcg_eigen.hpp>
#include "VecMatDef.h"

using namespace CppAD;
using namespace CppAD::cg;

typedef CG<T> CGD;
typedef AD<CGD> A;

class Autodiff
{
public:
	using VectorXT = Matrix<A, Eigen::Dynamic, 1>;
	using MatrixXT = Matrix<A, Eigen::Dynamic, Eigen::Dynamic>;
	using VectorXi = Vector<int, Eigen::Dynamic>;
	using MatrixXi = Matrix<int, Eigen::Dynamic, Eigen::Dynamic>;
	using VtxList = std::vector<int>;
	using Entry = Eigen::Triplet<T>;
	using TV = Vector<A, 3>;
	using TV2 = Vector<A, 2>;
	using TM2 = Matrix<A, 2, 2>;
	using TV3 = Vector<A, 3>;
	using IV = Vector<int, 3>;
	using IV2 = Vector<int, 2>;
	using TM = Matrix<A, 3, 3>;

public:
	void generateDiscreteShellGravitionalEnergy();

public:
	Autodiff() {} 
	~Autodiff() {} 
};


#endif
