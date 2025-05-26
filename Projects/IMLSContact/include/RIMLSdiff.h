#ifndef RIMLSDIFF_H
#define RIMLSDIFF_H

#include"VecMatDef.h"

 AScalar fx_func_2(const Vector3a& x, const VectorXa& v, const Vector3a& ck, AScalar h, const AScalar f, const Vector3a& grad_f, AScalar sigma_n, AScalar sigma_r);
 VectorXa dfdx_func_2(const Vector3a& x, const VectorXa& v, const Vector3a& ck, AScalar h, const AScalar f, const Vector3a& grad_f, AScalar sigma_n, AScalar sigma_r);
 Eigen::Matrix<AScalar,43,43> d2fdx2_func_2(const Vector3a& x, const VectorXa& v, const Vector3a& ck, AScalar h, const AScalar f, const Vector3a& grad_f, AScalar sigma_n, AScalar sigma_r);

 AScalar gx_func_2(const Vector3a& x, const VectorXa& v, const Vector3a& ck, AScalar h, const AScalar f, const Vector3a& grad_f, AScalar sigma_n, AScalar sigma_r);
 VectorXa dgdx_func_2(const Vector3a& x, const VectorXa& v, const Vector3a& ck, AScalar h, const AScalar f, const Vector3a& grad_f, AScalar sigma_n, AScalar sigma_r);
 Eigen::Matrix<AScalar,43,43> d2gdx2_func_2(const Vector3a& x, const VectorXa& v, const Vector3a& ck, AScalar h, const AScalar f, const Vector3a& grad_f, AScalar sigma_n, AScalar sigma_r);

 AScalar sumGF1_func(const Vector3a& x, const VectorXa& v, const Vector3a& ck, AScalar h);
 VectorXa dsumGF1dx_func(const Vector3a& x, const VectorXa& v, const Vector3a& ck, AScalar h);
 VectorXa d2sumGF1dx2_func(const Vector3a& x, const VectorXa& v, const Vector3a& ck, AScalar h);

 AScalar sumGF2_func(const Vector3a& x, const VectorXa& v, const Vector3a& ck, AScalar h);
 VectorXa dsumGF2dx_func(const Vector3a& x, const VectorXa& v, const Vector3a& ck, AScalar h);
 VectorXa d2sumGF2dx2_func(const Vector3a& x, const VectorXa& v, const Vector3a& ck, AScalar h);

 AScalar sumGF3_func(const Vector3a& x, const VectorXa& v, const Vector3a& ck, AScalar h);
 VectorXa dsumGF3dx_func(const Vector3a& x, const VectorXa& v, const Vector3a& ck, AScalar h);
 VectorXa d2sumGF3dx2_func(const Vector3a& x, const VectorXa& v, const Vector3a& ck, AScalar h);

 AScalar sumGW1_func(const Vector3a& x, const VectorXa& v, const Vector3a& ck, AScalar h);
 VectorXa dsumGW1dx_func(const Vector3a& x, const VectorXa& v, const Vector3a& ck, AScalar h);
 VectorXa d2sumGW1dx2_func(const Vector3a& x, const VectorXa& v, const Vector3a& ck, AScalar h);

 AScalar sumGW2_func(const Vector3a& x, const VectorXa& v, const Vector3a& ck, AScalar h);
 VectorXa dsumGW2dx_func(const Vector3a& x, const VectorXa& v, const Vector3a& ck, AScalar h);
 VectorXa d2sumGW2dx2_func(const Vector3a& x, const VectorXa& v, const Vector3a& ck, AScalar h);

 AScalar sumGW3_func(const Vector3a& x, const VectorXa& v, const Vector3a& ck, AScalar h);
 VectorXa dsumGW3dx_func(const Vector3a& x, const VectorXa& v, const Vector3a& ck, AScalar h);
 VectorXa d2sumGW3dx2_func(const Vector3a& x, const VectorXa& v, const Vector3a& ck, AScalar h);

 AScalar sumN1_func(const Vector3a& x, const VectorXa& v, const Vector3a& ck, AScalar h);
 VectorXa dsumN1dx_func(const Vector3a& x, const VectorXa& v, const Vector3a& ck, AScalar h);
 VectorXa d2sumN1dx2_func(const Vector3a& x, const VectorXa& v, const Vector3a& ck, AScalar h);

 AScalar sumN2_func(const Vector3a& x, const VectorXa& v, const Vector3a& ck, AScalar h);
 VectorXa dsumN2dx_func(const Vector3a& x, const VectorXa& v, const Vector3a& ck, AScalar h);
 VectorXa d2sumN2dx2_func(const Vector3a& x, const VectorXa& v, const Vector3a& ck, AScalar h);

 AScalar sumN3_func(const Vector3a& x, const VectorXa& v, const Vector3a& ck, AScalar h);
 VectorXa dsumN3dx_func(const Vector3a& x, const VectorXa& v, const Vector3a& ck, AScalar h);
 VectorXa d2sumN3dx2_func(const Vector3a& x, const VectorXa& v, const Vector3a& ck, AScalar h);


#endif