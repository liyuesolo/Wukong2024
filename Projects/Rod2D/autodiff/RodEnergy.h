#ifndef ROD_ENERGY_H
#define ROD_ENERGY_H

#include "../include/VecMatDef.h"

void computeRod2DStretchEnergy(double stretchStiffness, double undeformedLength, const Eigen::Matrix<double,2,1> & x1, const Eigen::Matrix<double,2,1> & x2, double& energy);
void computeRod2DStretchEnergyGradient(double stretchStiffness, double undeformedLength, const Eigen::Matrix<double,2,1> & x1, const Eigen::Matrix<double,2,1> & x2, Eigen::Matrix<double, 4, 1>& energygradient);
void computeRod2DStretchEnergyHessian(double stretchStiffness, double undeformedLength, const Eigen::Matrix<double,2,1> & x1, const Eigen::Matrix<double,2,1> & x2, Eigen::Matrix<double, 4, 4>& energyhessian);


void computeRod2DBendEnergy(double bendStiffness, double undeformedLength, double undeformedMaterialCurvature, const Eigen::Matrix<double,2,1> & x1, const Eigen::Matrix<double,2,1> & x2, 
	const Eigen::Matrix<double,2,1> & x3, double& energy);
void computeRod2DBendEnergyGradient(double bendStiffness, double undeformedLength, double undeformedMaterialCurvature, const Eigen::Matrix<double,2,1> & x1, const Eigen::Matrix<double,2,1> & x2, 
	const Eigen::Matrix<double,2,1> & x3, Eigen::Matrix<double, 6, 1>& energygradient);
void computeRod2DBendEnergyHessian(double bendStiffness, double undeformedLength, double undeformedMaterialCurvature, const Eigen::Matrix<double,2,1> & x1, const Eigen::Matrix<double,2,1> & x2, 
	const Eigen::Matrix<double,2,1> & x3, Eigen::Matrix<double, 6, 6>& energyhessian);
    
#endif