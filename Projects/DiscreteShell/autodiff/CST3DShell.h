#ifndef CST_3D_SHELL_H
#define CST_3D_SHELL_H

#include "../include/VecMatDef.h"


T compute3DCSTShellEnergy(T poissonsRatio, T stiffness, const Matrix<T,3,1> & x1, const Matrix<T,3,1> & x2, const Matrix<T,3,1> & x3, 
	const Matrix<T,3,1> & x1Undef, const Matrix<T,3,1> & x2Undef, const Matrix<T,3,1> & x3Undef);

void compute3DCSTShellEnergyGradient(T poissonsRatio, T stiffness, const Matrix<T,3,1> & x1, const Matrix<T,3,1> & x2, const Matrix<T,3,1> & x3, 
	const Matrix<T,3,1> & x1Undef, const Matrix<T,3,1> & x2Undef, const Matrix<T,3,1> & x3Undef, Matrix<T, 9, 1>& energygradient);

void compute3DCSTShellEnergyHessian(T poissonsRatio, T stiffness, const Matrix<T,3,1> & x1, const Matrix<T,3,1> & x2, const Matrix<T,3,1> & x3, 
	const Matrix<T,3,1> & x1Undef, const Matrix<T,3,1> & x2Undef, const Matrix<T,3,1> & x3Undef, Matrix<T, 9, 9>& energyhessian);


double computeDSBendingEnergy(double stiffness, const Matrix<double,3,1> & x0, const Matrix<double,3,1> & x1, const Matrix<double,3,1> & x2, const Matrix<double,3,1> & x3, 
	const Matrix<double,3,1> & x0Undef, const Matrix<double,3,1> & x1Undef, const Matrix<double,3,1> & x2Undef, const Matrix<double,3,1> & x3Undef);
void computeDSBendingEnergyGradient(double stiffness, const Matrix<double,3,1> & x0, const Matrix<double,3,1> & x1, const Matrix<double,3,1> & x2, const Matrix<double,3,1> & x3, 
	const Matrix<double,3,1> & x0Undef, const Matrix<double,3,1> & x1Undef, const Matrix<double,3,1> & x2Undef, const Matrix<double,3,1> & x3Undef, Matrix<double, 12, 1>& energygradient);
void computeDSBendingEnergyHessian(double stiffness, const Matrix<double,3,1> & x0, const Matrix<double,3,1> & x1, const Matrix<double,3,1> & x2, const Matrix<double,3,1> & x3, 
	const Matrix<double,3,1> & x0Undef, const Matrix<double,3,1> & x1Undef, const Matrix<double,3,1> & x2Undef, const Matrix<double,3,1> & x3Undef, Matrix<double, 12, 12>& energyhessian);

T compute3DCSTGravitationalEnergy(T undeformedDensity, T thickness, const Matrix<T,3,1> & gravity, const Matrix<T,3,1> & x1, const Matrix<T,3,1> & x2, 
	const Matrix<T,3,1> & x3, const Matrix<T,3,1> & x1Undef, const Matrix<T,3,1> & x2Undef, const Matrix<T,3,1> & x3Undef);
void compute3DCSTGravitationalEnergyGradient(T undeformedDensity, T thickness, const Matrix<T,3,1> & gravity, const Matrix<T,3,1> & x1, const Matrix<T,3,1> & x2, 
	const Matrix<T,3,1> & x3, const Matrix<T,3,1> & x1Undef, const Matrix<T,3,1> & x2Undef, const Matrix<T,3,1> & x3Undef, Matrix<T, 9, 1>& energygradient);
void compute3DCSTGravitationalEnergyHessian(T undeformedDensity, T thickness, const Matrix<T,3,1> & gravity, const Matrix<T,3,1> & x1, const Matrix<T,3,1> & x2, 
	const Matrix<T,3,1> & x3, const Matrix<T,3,1> & x1Undef, const Matrix<T,3,1> & x2Undef, const Matrix<T,3,1> & x3Undef, Matrix<T, 9, 9>& energyhessian);	
#endif