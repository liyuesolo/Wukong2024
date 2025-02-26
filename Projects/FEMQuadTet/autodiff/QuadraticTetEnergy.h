#ifndef QUADRATIC_TET_ENERGY_H
#define QUADRATIC_TET_ENERGY_H

#include "../include/VecMatDef.h"

void computeTet10Basis(const Eigen::Matrix<double,3,1> & xi, Eigen::Matrix<double, 10, 1>& basis);
void computeTet10BasisJacobian(const Eigen::Matrix<double,3,1> & xi, Eigen::Matrix<double, 10, 3>& jacobian);

void compute3DNeoHookeanQuadraticTetEnergy(double E, double nu, const Eigen::Matrix<double,10,3> & x, const Eigen::Matrix<double,10,3> & xUndef, double& energy);
void compute3DNeoHookeanQuadraticTetEnergyGradient(double E, double nu, const Eigen::Matrix<double,10,3> & x, const Eigen::Matrix<double,10,3> & xUndef, Eigen::Matrix<double, 30, 1>& energygradient);
void compute3DNeoHookeanQuadraticTetEnergyHessian(double E, double nu, const Eigen::Matrix<double,10,3> & x, const Eigen::Matrix<double,10,3> & xUndef, Eigen::Matrix<double, 30, 30>& energyhessian);



#endif