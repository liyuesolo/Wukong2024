#ifndef HYBRID_C2_CURVE_H
#define HYBRID_C2_CURVE_H


// reference
// view-source:http://www.cemyuksel.com/research/interpolating_splines/curves.html
// A Class of C2 Interpolating Splines ACM Transactions on Graphics, 39, 5, 2020

#include <utility>
#include <iostream>
#include <Eigen/Geometry>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include "VecMatDef.h"

template<int dim>
struct CurveData
{
    CurveData(Vector<T, dim> c, Vector<T, dim> a1, Vector<T, dim> a2, Vector<T, 2> bound)
    : center(c), axis1(a1), axis2(a2), limits(bound) {}
    Vector<T, dim> center, axis1, axis2;
    Vector<T, 2> limits;
};

template<int dim>
class HybridC2Curve
{
public:
    using TV = Vector<T, dim>;

    int sub_div;

    std::vector<TV> data_points;
    std::vector<TV> points_on_curve;
    std::vector<CurveData<dim>*> curve_data;
 
public:
    HybridC2Curve(int _sub_div) : sub_div(_sub_div) {}
    HybridC2Curve() : sub_div(64) {}
    ~HybridC2Curve() 
    {
        for(auto data : curve_data)
        {
            delete data;
        }
    }

    void getLinearSegments(std::vector<TV>& points)
    {
    
        if (data_points.size() <= 2)
            points = data_points;
        else
        {
            points = points_on_curve;
        }
    }

    void constructIndividualCurves();

    void sampleCurves(std::vector<TV>& points);
    
    void normalizeDataPoints();

    void getPosOnCurve(int curve_idx, T t, TV& pos, bool interpolate);

    void derivativeTestdCdt();

    void getPosOnCurveWithDerivatives(int curve_idx, T t, 
        TV& ci, TV& dci, TV& ddci,
        bool compute_dci, bool compute_ddci,
        bool interpolate);

    // adapted from line 1143 curvePos
    void F(T t, int idx, TV& pos)
    {
        T tt = curve_data[idx]->limits[0] + t * (curve_data[idx]->limits[1] - curve_data[idx]->limits[0]);
        pos = curve_data[idx]->center + 
                curve_data[idx]->axis1 * std::cos(tt) +
                curve_data[idx]->axis2 * std::sin(tt);
    }
    

    void dF(T t, int idx, TV& dpos)
    {
		T theta  = curve_data[idx]->limits[0] + t * (curve_data[idx]->limits[1] - curve_data[idx]->limits[0]);
		T dtheta = (curve_data[idx]->limits[1] - curve_data[idx]->limits[0]);
        
        dpos = (-std::sin(theta) * curve_data[idx]->axis1 + 
                curve_data[idx]->axis2 * std::cos(theta)) * dtheta;
    }

    void ddF(T t, int idx, TV& ddpos)
    {
        T theta  = curve_data[idx]->limits[0] + t * (curve_data[idx]->limits[1] - curve_data[idx]->limits[0]);
		T dtheta = (curve_data[idx]->limits[1] - curve_data[idx]->limits[0]);
        
        ddpos = -dtheta * dtheta * (std::cos(theta) * curve_data[idx]->axis1 + 
                curve_data[idx]->axis2 * std::sin(theta));
    }

private:

    // adapted from source code
    void circularInterpolation(int i, TV& center, TV& axis1, TV& axis2, Vector<T, 3>& limits);

    // adapted from source code
    void ellipticalInterpolation(int i, TV& center, TV& axis1, TV& axis2, Vector<T, 3>& limits);

    // adapted from source code
    void hybridInterpolation(int i, TV& center, TV& axis1, TV& axis2, Vector<T, 3>& limits);

};



#endif