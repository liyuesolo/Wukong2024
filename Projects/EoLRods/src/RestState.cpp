#include "../include/RestState.h"


void DiscreteHybridCurvature::getMaterialPos(T u, TV& X, 
    TV& dXdu, TV& d2Xdu2, bool g, bool h) 
{
    bool debug = false;

    if (u < 0)
    {
        X = this->starting_point.template segment<3>(0);
        dXdu = TV::Zero();
        d2Xdu2 = TV::Zero();
        return;
    }

    if (debug)
        std::cout << "getMaterialPos" << std::endl;
    T dx = 1.0;

    int left_node = 0;
    // find which curve the current points lies on
    for (int i = 0; i < curve->data_points.size(); i++)
    {
        T fraction_data_point = i;

        if( fraction_data_point > u || i == curve->data_points.size() - 1)
        {
            left_node = i - 1;
            break;
        }
    }

    int curve_id = 0;
    
    if(left_node > 1)
        curve_id = left_node - 1;

    T t = u - T(left_node * dx);
    
    bool interpolate = false;
    if ( left_node != 0 && left_node != curve->data_points.size() - 2)
        interpolate = true;
    if (curve->data_points.size() == 3 && u >= 1)
    {
        curve_id += 1;
        interpolate= false;
    }
        
    if (debug)
        std::cout << "t " << t << " u " << u <<  " curve idx " << curve_id << " left node " << left_node << " interpolate " << interpolate << std::endl;
    
    // curve->getPosOnCurve(curve_id, t, X, interpolate);

    TV2 _X, _dXdu, _d2Xdu2;
    curve->getPosOnCurveWithDerivatives(curve_id, t, _X, _dXdu, _d2Xdu2, g, h, interpolate);
    X = TV(_X[0], _X[1], 0);
    dXdu = TV(_dXdu[0], _dXdu[1], 0);
    d2Xdu2 = TV(_d2Xdu2[0], _d2Xdu2[1], 0);
    
    // else if constexpr (dim == 2)
    //     curve->getPosOnCurveWithDerivatives(curve_id, t, X, dXdu, d2Xdu2, g, h, interpolate);

    // T shifted = shift ? t * 0.5 + 0.5 : t * 0.5;
    // X *= 0.03;
    if (debug)
    {
        std::cout << u << " " << X.transpose() << std::endl;
        std::cout << "getMaterialPos done" << std::endl;
        // std::getchar();
    }
}

void LineCurvature::getMaterialPos(T u, TV& X, TV& dXdu, TV& d2Xdu2, bool g, bool h) 
{ 
    T scale = (this->ending_point[3] - this->starting_point[3]);
    X = this->starting_point.template segment<3>(0) + 
        (u - this->starting_point[3]) / scale
         * (this->ending_point.template segment<3>(0) - this->starting_point.template segment<3>(0));
    
    dXdu = (this->ending_point.template segment<3>(0) - this->starting_point.template segment<3>(0)) / scale;
    d2Xdu2 = TV::Zero();
}

