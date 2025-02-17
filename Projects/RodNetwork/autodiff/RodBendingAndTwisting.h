#ifndef ROD_BENDING_AND_TWISTING_ENERGY_H
#define ROD_BENDING_AND_TWISTING_ENERGY_H

#include "../include/VecMatDef.h"

void computeRodBendingAndTwistEnergy(const Eigen::Matrix<double,2,2> & B, double kt, double undeformedTwist, const Eigen::Matrix<double,3,1> & referenceNormal1, const Eigen::Matrix<double,3,1> & referenceTangent1, 
	const Eigen::Matrix<double,3,1> & referenceNormal2, const Eigen::Matrix<double,3,1> & referenceTangent2, double referenceTwist, const Eigen::Matrix<double,3,1> & x1, const Eigen::Matrix<double,3,1> & x2, 
	const Eigen::Matrix<double,3,1> & x3, const Eigen::Matrix<double,3,1> & X1, const Eigen::Matrix<double,3,1> & X2, const Eigen::Matrix<double,3,1> & X3, double theta1, 
	double theta2, double& energy);
void computeRodBendingAndTwistEnergyGradient(const Eigen::Matrix<double,2,2> & B, double kt, double undeformedTwist, const Eigen::Matrix<double,3,1> & referenceNormal1, const Eigen::Matrix<double,3,1> & referenceTangent1, 
	const Eigen::Matrix<double,3,1> & referenceNormal2, const Eigen::Matrix<double,3,1> & referenceTangent2, double referenceTwist, const Eigen::Matrix<double,3,1> & x1, const Eigen::Matrix<double,3,1> & x2, 
	const Eigen::Matrix<double,3,1> & x3, const Eigen::Matrix<double,3,1> & X1, const Eigen::Matrix<double,3,1> & X2, const Eigen::Matrix<double,3,1> & X3, double theta1, 
	double theta2, Eigen::Matrix<double, 11, 1>& energygradient);
void computeRodBendingAndTwistEnergyHessian(const Eigen::Matrix<double,2,2> & B, double kt, double undeformedTwist, const Eigen::Matrix<double,3,1> & referenceNormal1, const Eigen::Matrix<double,3,1> & referenceTangent1, 
	const Eigen::Matrix<double,3,1> & referenceNormal2, const Eigen::Matrix<double,3,1> & referenceTangent2, double referenceTwist, const Eigen::Matrix<double,3,1> & x1, const Eigen::Matrix<double,3,1> & x2, 
	const Eigen::Matrix<double,3,1> & x3, const Eigen::Matrix<double,3,1> & X1, const Eigen::Matrix<double,3,1> & X2, const Eigen::Matrix<double,3,1> & X3, double theta1, 
	double theta2, Eigen::Matrix<double, 11, 11>& energyhessian);

#endif