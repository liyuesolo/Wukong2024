#ifndef JOINT_BENDING_AND_TWISTING_H
#define JOINT_BENDING_AND_TWISTING_H

#include "../include/VecMatDef.h"

void computeRodEulerAngleBendingAndTwistEnergyRBFirst(const Eigen::Matrix<double,2,2> & bendingStiffness, double twistStiffness, double undeformedTwist, const Eigen::Matrix<double,3,1> & referenceTangent1, const Eigen::Matrix<double,3,1> & referenceNormal1, 
	const Eigen::Matrix<double,3,1> & referenceTangent2, const Eigen::Matrix<double,3,1> & referenceNormal2, double referenceTwist, const Eigen::Matrix<double,3,1> & xRB, const Eigen::Matrix<double,3,1> & x3, 
	const Eigen::Matrix<double,3,1> & XRB, const Eigen::Matrix<double,3,1> & X3, const Eigen::Matrix<double,3,3> & rigidBodyAccumulatingRotation, const Eigen::Matrix<double,3,1> & rigidBodyRotationDelta, double theta2, 
	double& energy);
void computeRodEulerAngleBendingAndTwistEnergyRBFirstGradient(const Eigen::Matrix<double,2,2> & bendingStiffness, double twistStiffness, double undeformedTwist, const Eigen::Matrix<double,3,1> & referenceTangent1, const Eigen::Matrix<double,3,1> & referenceNormal1, 
	const Eigen::Matrix<double,3,1> & referenceTangent2, const Eigen::Matrix<double,3,1> & referenceNormal2, double referenceTwist, const Eigen::Matrix<double,3,1> & xRB, const Eigen::Matrix<double,3,1> & x3, 
	const Eigen::Matrix<double,3,1> & XRB, const Eigen::Matrix<double,3,1> & X3, const Eigen::Matrix<double,3,3> & rigidBodyAccumulatingRotation, const Eigen::Matrix<double,3,1> & rigidBodyRotationDelta, double theta2, 
	Eigen::Matrix<double, 10, 1>& energygradient);
void computeRodEulerAngleBendingAndTwistEnergyRBFirstHessian(const Eigen::Matrix<double,2,2> & bendingStiffness, double twistStiffness, double undeformedTwist, const Eigen::Matrix<double,3,1> & referenceTangent1, const Eigen::Matrix<double,3,1> & referenceNormal1, 
	const Eigen::Matrix<double,3,1> & referenceTangent2, const Eigen::Matrix<double,3,1> & referenceNormal2, double referenceTwist, const Eigen::Matrix<double,3,1> & xRB, const Eigen::Matrix<double,3,1> & x3, 
	const Eigen::Matrix<double,3,1> & XRB, const Eigen::Matrix<double,3,1> & X3, const Eigen::Matrix<double,3,3> & rigidBodyAccumulatingRotation, const Eigen::Matrix<double,3,1> & rigidBodyRotationDelta, double theta2, 
	Eigen::Matrix<double, 10, 10>& energyhessian);

void computeRodEulerAngleBendingAndTwistEnergyRBSecond(const Eigen::Matrix<double,2,2> & bendingStiffness, double twistStiffness, double undeformedTwist, const Eigen::Matrix<double,3,1> & referenceTangent1, const Eigen::Matrix<double,3,1> & referenceNormal1, 
	const Eigen::Matrix<double,3,1> & referenceTangent2, const Eigen::Matrix<double,3,1> & referenceNormal2, double referenceTwist, const Eigen::Matrix<double,3,1> & xRB, const Eigen::Matrix<double,3,1> & x1, 
	const Eigen::Matrix<double,3,1> & XRB, const Eigen::Matrix<double,3,1> & X1, const Eigen::Matrix<double,3,3> & rigidBodyAccumulatingRotation, const Eigen::Matrix<double,3,1> & rigidBodyRotationDelta, double theta1, 
	double& energy);
void computeRodEulerAngleBendingAndTwistEnergyRBSecondGradient(const Eigen::Matrix<double,2,2> & bendingStiffness, double twistStiffness, double undeformedTwist, const Eigen::Matrix<double,3,1> & referenceTangent1, const Eigen::Matrix<double,3,1> & referenceNormal1, 
	const Eigen::Matrix<double,3,1> & referenceTangent2, const Eigen::Matrix<double,3,1> & referenceNormal2, double referenceTwist, const Eigen::Matrix<double,3,1> & xRB, const Eigen::Matrix<double,3,1> & x1, 
	const Eigen::Matrix<double,3,1> & XRB, const Eigen::Matrix<double,3,1> & X1, const Eigen::Matrix<double,3,3> & rigidBodyAccumulatingRotation, const Eigen::Matrix<double,3,1> & rigidBodyRotationDelta, double theta1, 
	Eigen::Matrix<double, 10, 1>& energygradient);
void computeRodEulerAngleBendingAndTwistEnergyRBSecondHessian(const Eigen::Matrix<double,2,2> & bendingStiffness, double twistStiffness, double undeformedTwist, const Eigen::Matrix<double,3,1> & referenceTangent1, const Eigen::Matrix<double,3,1> & referenceNormal1, 
	const Eigen::Matrix<double,3,1> & referenceTangent2, const Eigen::Matrix<double,3,1> & referenceNormal2, double referenceTwist, const Eigen::Matrix<double,3,1> & xRB, const Eigen::Matrix<double,3,1> & x1, 
	const Eigen::Matrix<double,3,1> & XRB, const Eigen::Matrix<double,3,1> & X1, const Eigen::Matrix<double,3,3> & rigidBodyAccumulatingRotation, const Eigen::Matrix<double,3,1> & rigidBodyRotationDelta, double theta1, 
	Eigen::Matrix<double, 10, 10>& energyhessian);


#endif