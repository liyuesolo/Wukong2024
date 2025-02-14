#ifndef EOL_ROD_STRETCHING_ENERGY_H
#define EOL_ROD_STRETCHING_ENERGY_H

#include "../include/VecMatDef.h"

void computeRodStretchingEnergy(double ks, const Eigen::Matrix<double,3,1> & X1, const Eigen::Matrix<double,3,1> & X2, const Eigen::Matrix<double,3,1> & x1, const Eigen::Matrix<double,3,1> & x2, 
	double& energy);
void computeRodStretchingEnergyGradient(double ks, const Eigen::Matrix<double,3,1> & X1, const Eigen::Matrix<double,3,1> & X2, const Eigen::Matrix<double,3,1> & x1, const Eigen::Matrix<double,3,1> & x2, 
	Eigen::Matrix<double, 6, 1>& energygradient);
void computeRodStretchingEnergyHessian(double ks, const Eigen::Matrix<double,3,1> & X1, const Eigen::Matrix<double,3,1> & X2, const Eigen::Matrix<double,3,1> & x1, const Eigen::Matrix<double,3,1> & x2, 
	Eigen::Matrix<double, 6, 6>& energyhessian);

    
#endif