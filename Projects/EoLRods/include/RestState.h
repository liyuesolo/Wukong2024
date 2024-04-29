#ifndef CURVATURE_FUNCTION_H
#define CURVATURE_FUNCTION_H

#include <cmath>
#include <iostream>
#include <utility>
#include <vector>

#include "HybridC2Curve.h"
#include "VecMatDef.h"


class RestState
{
public:
    using TV = Vector<T, 3>;
    using TV2 = Vector<T, 2>;

    Vector<T, 3 + 1> starting_point, ending_point;

public:
    RestState(Vector<T, 3 + 1> q0, Vector<T, 3 + 1> q1)
        : starting_point(q0), ending_point(q1) {}
    RestState()
        : starting_point(Vector<T, 3 + 1>::Zero()), ending_point(Vector<T, 3 + 1>::Ones()) {}
    ~RestState() {}

    virtual T value(T u) {return 0;}
    virtual void gradient(T u, T& dedu) { dedu = 0; }
    virtual void hessian(T u, T& de2du2) { de2du2 = 0; }
    
    virtual void getMaterialPos(T u, TV& X, TV& dXdu, TV& d2Xdu2, bool g, bool h) 
    { 
        X = starting_point.template segment<3>(0) + 
            (u - starting_point[3]) / (ending_point[3] - starting_point[3]) * (ending_point.template segment<3>(0) - starting_point.template segment<3>(0));
        dXdu = (ending_point.template segment<3>(0) - starting_point.template segment<3>(0)) / (ending_point[3] - starting_point[3]) ;
        d2Xdu2 = TV::Zero();
    }
};


class LineCurvature : public RestState
{
    using TV = Vector<T, 3>;

public:
    LineCurvature(Vector<T, 3 + 1> q0, 
                Vector<T, 3 + 1> q1) : RestState(q0, q1) {}
    LineCurvature() : RestState(){}
    
    virtual void getMaterialPos(T u, TV& X, TV& dXdu, TV& d2Xdu2, bool g, bool h);
};


class DiscreteHybridCurvature : public RestState
{
public:
    using TV = Vector<T, 3>;
    using TV2 = Vector<T, 2>;
    HybridC2Curve<2>* curve;
    std::vector<T> data_points_discrete_arc_length;

public:
    DiscreteHybridCurvature() : RestState(){}

    DiscreteHybridCurvature(Vector<T, 3 + 1> q0, 
                            Vector<T, 3 + 1> q1) : RestState(q0, q1) {}
    virtual void getMaterialPos(T u, TV& X, TV& dXdu, TV& d2Xdu2, bool g, bool h); 

    void setData(HybridC2Curve<2>* c, const std::vector<T>& v)
    {
        curve = c;
        data_points_discrete_arc_length = v;
    }

};

#endif