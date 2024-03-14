#include "../include/HybridC2Curve.h"


template<int dim>
void HybridC2Curve<dim>::constructIndividualCurves()
{
    curve_data.clear();
    
    for(int i = 1; i < data_points.size() - 1; i++)
    {
        TV center, axis1, axis2; Vector<T, 3> limits;
        hybridInterpolation(i, center, axis1, axis2, limits);
        curve_data.push_back(new CurveData<dim>(
                                center, axis1, axis2, 
                                Vector<T, 2>(limits[0], limits[2])));
                                
    }
}

template<int dim>
void HybridC2Curve<dim>::getPosOnCurve(int curve_idx, T t, TV& pos, bool interpolate)
{
    if (interpolate)
    {
        T theta = t * M_PI * 0.5;
        TV F0, F1;
        F((theta + M_PI * 0.5)/M_PI, curve_idx, F0);
        F(theta/M_PI, curve_idx+1, F1);
        pos = std::cos(theta) * std::cos(theta) * F0 + std::sin(theta) * std::sin(theta) * F1;
    }
    else
    {
        F(t, curve_idx, pos);
    }
}

template<int dim>
void HybridC2Curve<dim>::normalizeDataPoints()
{
    TV min_point = TV::Ones() * 1e10, max_point = TV::Ones() * -1e10;
    for (auto pt : data_points)
    {
        for (int d = 0; d < dim; d++)
        {
            min_point[d] = std::min(min_point[d], pt[d]);
            max_point[d] = std::min(max_point[d], pt[d]);
        }
    }
    T longest_dis = -1000;
    for (int d = 0; d < dim; d++)
        longest_dis = std::max(longest_dis, max_point[d] - min_point[d]);
    for (TV& pt : data_points)
    {
        pt = (pt - min_point) / longest_dis;
    }
}

template<int dim>
void HybridC2Curve<dim>::sampleCurves(std::vector<TV>& points)
{
    points_on_curve.clear();
    constructIndividualCurves();

    T ti = curve_data.front()->limits[0];
    T ti_plus_1 = 0;
    T dt = (ti_plus_1 - ti) / (sub_div / 2);
    for (int i = 0; i < sub_div / 2; i++)
    {
        T t = (ti + dt * (T)i - curve_data.front()->limits[0])/ (curve_data.front()->limits[1] - curve_data.front()->limits[0]);
        TV F0; F(t, 0, F0);
        points_on_curve.push_back(F0);    
    }
    for (int curve_idx = 0; curve_idx < curve_data.size()-1; curve_idx++)
    {
        T ti = 0;
        T ti_plus_1 = curve_data[curve_idx]->limits[1]; 
        
        T dt = (ti_plus_1 - ti) / (sub_div / 2);
        T dt_next = (0 - curve_data[curve_idx + 1]->limits[0]) / (sub_div / 2); 
        for (int i = 0; i < sub_div / 2; i++)
        {
            T t = (ti + dt * (T)i - curve_data[curve_idx]->limits[0]) / (curve_data[curve_idx]->limits[1] - curve_data[curve_idx]->limits[0]);

            T theta = (dt * T(i)) / (ti_plus_1 - ti) * M_PI * 0.5;
            TV F0, F1;
            F(t, curve_idx, F0);
            T t_next = ( dt_next * T(i)) / (curve_data[curve_idx + 1]->limits[1] - curve_data[curve_idx + 1]->limits[0]);
            F(t_next, curve_idx+1, F1);
            // equation (2)
            TV Ci = std::cos(theta) * std::cos(theta) * F0 + std::sin(theta) * std::sin(theta) * F1;
            points_on_curve.push_back(Ci);
            // std::cout << "curve idx: " << curve_idx << " limit 0: " << curve_data[curve_idx]->limits[0] << " limit 1: " << curve_data[curve_idx]->limits[1] << std::endl;
        }
    }

    ti = 0;
    ti_plus_1 = curve_data.back()->limits[1];
    dt = (ti_plus_1 - ti) / (sub_div / 2);
    for (int i = 0; i < sub_div / 2; i++)
    {
        T t = (ti + dt * (T)i - curve_data.back()->limits[0]) / (curve_data.back()->limits[1] - curve_data.back()->limits[0]);
        TV F1; F(t, curve_data.size() - 1, F1);
        // std::cout << " F1 at last segment " << F1.transpose() << std::endl;
        points_on_curve.push_back(F1);
    }

    points_on_curve.push_back(data_points.back());

    points = points_on_curve;
}

template<int dim>
void HybridC2Curve<dim>::derivativeTestdCdt()
{
    constructIndividualCurves();
    T epsilon = 1e-6;

    //check if the derivative of individual function is correct.
    // passed
    for(int i = 0; i < curve_data.size(); i++)
    {
        for (T t = 0.0; t < 1.0; t += 0.1)
        {
            TV dF0dt, dF1dt, d2F0dt2;
            dF(t, i, dF0dt);
            dF(t + epsilon, i, dF1dt);
            ddF(t, i, d2F0dt2);
            TV F0, F1;
            F(t, i, F0); F(t + epsilon, i, F1);
            for (int d = 0; d < dim; d++)
            {
                std::cout << "dF: " <<  (F1[d] - F0[d]) / epsilon << " " << dF0dt[d] << std::endl;
                std::cout << "ddF " << (dF1dt[d] - dF0dt[d]) / epsilon << " " << d2F0dt2[d] << std::endl;
            }
            std::getchar();
        }
    }

    //check interpolated value

    for(int i = 1; i < curve_data.size() - 1; i++)
    {
        std::cout << "curve " << i << std::endl;
        for (T t = 0.0; t < 1.0; t += 0.1)
        {
            auto Ci = [&](T _t, int i, TV& c)
            {
                T ti = (0 + _t * (curve_data[i]->limits[1] - 0) - curve_data[i]->limits[0]) / (curve_data[i]->limits[1] - curve_data[i]->limits[0]);
                T ti_plus_1 = (_t * (0 - curve_data[i + 1]->limits[0]) ) / (curve_data[i + 1]->limits[1] - curve_data[i + 1]->limits[0]);

                TV F0, F1;
                F(ti, i, F0);
                F(ti_plus_1, i+1, F1);

                T theta = _t * M_PI * 0.5;

                T cos2 = std::cos(theta) * std::cos(theta);
                T sin2 = std::sin(theta) * std::sin(theta);
                c = cos2 * F0 + sin2 * F1;
            };

            auto computedCdt = [&](T _t, int i, TV& dcdt)
            {
                T ti = (0 + _t * (curve_data[i]->limits[1] - 0) - curve_data[i]->limits[0]) / (curve_data[i]->limits[1] - curve_data[i]->limits[0]);
                T ti_plus_1 = (_t * (0 - curve_data[i + 1]->limits[0]) ) / (curve_data[i + 1]->limits[1] - curve_data[i + 1]->limits[0]);
                
                T dtidt = (curve_data[i]->limits[1] - 0) / (curve_data[i]->limits[1] - curve_data[i]->limits[0]);
                T dti_plus_1dt = (0 - curve_data[i + 1]->limits[0]) / (curve_data[i + 1]->limits[1] - curve_data[i + 1]->limits[0]);

                TV F0, F1;
                F(ti, i, F0);
                F(ti_plus_1, i+1, F1);

                T theta = _t * M_PI * 0.5;
                
                TV dF0dt, dF1dt;

                T dtheta_dt = M_PI * 0.5;
                
                T cos2 = std::cos(theta) * std::cos(theta);
                T sin2 = std::sin(theta) * std::sin(theta);


                dF(ti, i, dF0dt); 
                dF(ti_plus_1, i + 1, dF1dt); 
                
                TV dF_dtheta = TV::Zero();
                dF_dtheta += 2.0 * std::cos(theta) * std::sin(theta) * (F1 - F0) * dtheta_dt;
                // std::cout << dF_dtheta.transpose() << std::endl;
                dF_dtheta += cos2 * dF0dt * dtidt ;
                // std::cout << dF_dtheta.transpose() << std::endl;
                dF_dtheta += sin2 * dF1dt * dti_plus_1dt; // write it out for clearity
                // std::cout << dF_dtheta.transpose() << std::endl;
                dcdt = dF_dtheta ;
            };

            auto ddCi = [&](T _t, int i, TV& d2cdt2)
            {
                TV dF0dt, dF1dt;
                TV d2F0dt2, d2F1dt2;
                
                T dtheta_dt = M_PI * 0.5;

                T ti = (0 + _t * (curve_data[i]->limits[1] - 0) - curve_data[i]->limits[0]) / (curve_data[i]->limits[1] - curve_data[i]->limits[0]);
                T ti_plus_1 = (_t * (0 - curve_data[i + 1]->limits[0]) ) / (curve_data[i + 1]->limits[1] - curve_data[i + 1]->limits[0]);
                T dtidt = (curve_data[i]->limits[1] - 0) / (curve_data[i]->limits[1] - curve_data[i]->limits[0]);
                T dti_plus_1dt = (0 - curve_data[i + 1]->limits[0]) / (curve_data[i + 1]->limits[1] - curve_data[i + 1]->limits[0]);

                TV F0, F1;
                F(ti, i, F0);
                F(ti_plus_1, i+1, F1);

                T theta = _t * M_PI * 0.5;

                T cos2 = std::cos(theta) * std::cos(theta);
                T sin2 = std::sin(theta) * std::sin(theta);


                dF(ti, i, dF0dt); 
                dF(ti_plus_1, i + 1, dF1dt); 

                ddF(ti, i, d2F0dt2); 
                ddF(ti_plus_1, i + 1, d2F1dt2); 

                d2cdt2 = TV::Zero();
                d2cdt2 += 2.0 * (cos2 - sin2) * (F1 - F0) * dtheta_dt * dtheta_dt;
                
                d2cdt2 += 4.0 * std::cos(theta) * std::sin(theta) * (dF1dt * dti_plus_1dt - dF0dt * dtidt) * dtheta_dt;

                d2cdt2 += cos2 * d2F0dt2 * dtidt * dtidt + sin2 * d2F1dt2 * dti_plus_1dt * dti_plus_1dt;

            };
            
            TV c0, c1;
            Ci(t, i, c0); Ci(t + epsilon, i, c1);
            
            TV dc0; computedCdt(t, i, dc0);
            // std::cout << "C0 " << c0.transpose() << " C1 " << c1.transpose() << " dc0 " << dc0.transpose() << std::endl;
            TV dc1; computedCdt(t + epsilon, i, dc1);

            TV ddc0; ddCi(t, i, ddc0);
            std::cout << "t = " << t << std::endl;
            for (int d = 0; d < dim; d++)
            {
                std::cout << "dc: " <<  (c1[d] - c0[d]) / epsilon << " " << dc0[d] << std::endl;
                std::cout << "ddc: " <<  (dc1[d] - dc0[d]) / epsilon << " " << ddc0[d] << std::endl;
            }
            std::getchar();
        }
    }

}

// adapted from source code
template<int dim>
void HybridC2Curve<dim>::circularInterpolation(int i, TV& center, TV& axis1, TV& axis2, Vector<T, 3>& limits)
{
    axis1.setZero(); axis2.setZero(); center.setZero(); limits.setZero();

    int j = (i - 1 + data_points.size()) % data_points.size();
    int k = (i + 1) % data_points.size();

    TV vec1 = data_points[i] - data_points[j];
    TV mid1 = data_points[j] + 0.5 * vec1;
    TV vec2 = data_points[k] - data_points[i];
    TV mid2 = data_points[i] + 0.5 * vec2;
    TV dir1, dir2;
    dir1[0] = -vec1[1]; dir1[1] = vec1[0];
    dir2[0] = -vec2[1]; dir2[1] = vec2[0];
    
    T det = dir1[0] * dir2[1] - dir1[1] * dir2[0];

    if (std::abs(det) < 0.001)
    {
        if(vec1[0] * vec1[0] + vec1[1] * vec1[1] >= 0 || data_points.size() <=2 )
        {
            T small_angle = 0.01;
            T s = std::sin(small_angle);
            T l1 = vec1.norm(), l2 = vec2.norm();
            center = data_points[i];
            axis1 = TV::Zero();
            axis2 = vec2 / s;
            limits = Vector<T, 3>(-small_angle * l1 / l2, 0, small_angle);
            return;
        }
        else
            det = 0.001;
    }
    T s = ( dir2[1] * (mid2[0] - mid1[0]) + dir2[0] * (mid1[1] - mid2[1]) ) / det;
    center = mid1 + dir1 * s;
    axis1  = data_points[i] - center;
    axis2[0] = -axis1[1]; axis2[1] = axis1[0];
    T len2   = axis1[0]*axis1[0] + axis1[1]*axis1[1];
    TV toPt2  = data_points[k] - center;
    T limit2 = std::atan2( axis2.dot(toPt2), axis1.dot(toPt2) );
    TV toPt1  = data_points[j] - center;
    T limit1 = std::atan2( axis2.dot(toPt1), axis1.dot(toPt1) );

    if ( limit1 * limit2 > 0 ) 
    {
        if ( std::abs(limit1)<std::abs(limit2) ) limit2 += limit2 > 0 ? -T(2) * M_PI : T(2) * M_PI;
        if ( std::abs(limit1)>std::abs(limit2) ) limit1 += limit1 > 0 ? -T(2) * M_PI : T(2) * M_PI;
    }

    limits = Vector<T, 3>(limit1, 0, limit2);
}

// adapted from source code
template<int dim>
void HybridC2Curve<dim>::ellipticalInterpolation(int i, TV& center, TV& axis1, TV& axis2, Vector<T, 3>& limits)
{
    axis1.setZero(); axis2.setZero(); center.setZero(); limits.setZero();
    int numIter = 16;
    int j = (i - 1 + data_points.size()) % data_points.size();
    int k = (i + 1) % data_points.size();


    TV vec1 = data_points[j] - data_points[i];
    TV vec2 = data_points[k] - data_points[i];

    if ( data_points.size() <= 2 ) 
    {
        T small_angle = 0.01;
        T s = std::sin(small_angle);
        center = data_points[i];
        axis1 = TV::Zero();
        axis2 = vec2 / s;
        limits = Vector<T, 3>(-small_angle, 0, small_angle);
        return;
    }

    T len1 = std::sqrt( vec1[0]*vec1[0] + vec1[1]*vec1[1] );
    T len2 = std::sqrt( vec2[0]*vec2[0] + vec2[1]*vec2[1] );
    T cosa = (vec1[0]*vec2[0] + vec1[1]*vec2[1]) / (len1*len2);
    T maxA = std::acos(cosa);
    T ang  = maxA * 0.5;
    T incA = maxA * 0.25;
    T l1 = len1;
    T l2 = len2;
    if ( len1 < len2 ) { l1=len2; l2=len1; }
    T a, b, c, d;
    for ( int iter=0; iter<numIter; iter++ ) 
    {
        T theta = ang * 0.5;
        a = l1 * std::sin(theta);
        b = l1 * std::cos(theta);
        T beta = maxA - theta;
        c = l2 * std::sin(beta);
        d = l2 * std::cos(beta);
        T v = (1.0-d/b)*(1.0-d/b)+(c*c)/(a*a);	// ellipse equation
        ang += ( v > 1 ) ? incA : -incA;
        incA *= 0.5;
    }

    TV vec, pt2;
    T len;
    if ( len1 < len2 ) 
    {
        vec = vec2;
        len = len2;
        pt2 = data_points[k];
    } 
    else 
    {
        vec = vec1;
        len = len1;
        pt2 = data_points[j];
    }
    TV dir  = vec / len;
    TV perp = TV::Zero();
    perp[0] = -dir[1]; perp[1] = dir[0];
    T cross = vec1[0] * vec2[1] - vec1[1] * vec2[0];
    if ( (len1<len2 && cross>0) || (len1>=len2 && cross<0) ) 
    {
        perp[0] = dir[1]; 
        perp[1] = -dir[0];
    }
        
    T v = b*b/len;
    T h = b*a/len;
    axis1  = -dir * v - perp * h;
    center = data_points[i] - axis1;
    axis2 = pt2 - center;
    T beta   = std::asin(std::min(c/a, T(1)));
    if (len1<len2)
    {
        limits = Vector<T, 3>(-beta, 0, M_PI * 0.5);
    }
    else
    {
        axis2 *= -1;
        limits = Vector<T, 3>(-M_PI * 0.5, 0, beta);
    }
}

// adapted from source code
template<int dim>
void HybridC2Curve<dim>::hybridInterpolation(int i, TV& center, TV& axis1, TV& axis2, Vector<T, 3>& limits)
{
    circularInterpolation(i, center, axis1, axis2, limits);
    T lim0 = limits[0];
    T lim2 = limits[2];
    if (lim2 < lim0)
    {
        T tmp = lim0;
        lim0 = lim2;
        lim2 = tmp;
    }
    if (lim0 < -M_PI * 0.5 || lim2 > M_PI * 0.5)
        ellipticalInterpolation(i, center, axis1, axis2, limits);
}

template<int dim>
void HybridC2Curve<dim>::getPosOnCurveWithDerivatives(int curve_idx, T t, 
        TV& ci, TV& dci, TV& ddci,
        bool compute_dci, bool compute_ddci,
        bool interpolate)
{
    dci = TV::Zero(); ddci = TV::Zero();
    if (interpolate)
    {
        // we keep 0 here to indicate that's pi in the curve
        T ti = (0 + t * (curve_data[curve_idx]->limits[1] - 0) - curve_data[curve_idx]->limits[0]) / (curve_data[curve_idx]->limits[1] - curve_data[curve_idx]->limits[0]);
        T ti_plus_1 = (t * (0 - curve_data[curve_idx + 1]->limits[0]) ) / (curve_data[curve_idx + 1]->limits[1] - curve_data[curve_idx + 1]->limits[0]);
        T dtidt = (curve_data[curve_idx]->limits[1] - 0) / (curve_data[curve_idx]->limits[1] - curve_data[curve_idx]->limits[0]);
        T dti_plus_1dt = (0 - curve_data[curve_idx + 1]->limits[0]) / (curve_data[curve_idx + 1]->limits[1] - curve_data[curve_idx + 1]->limits[0]);
        
        T theta = t * M_PI * 0.5;
        
        TV F0, F1;
        F(ti, curve_idx, F0);
        F(ti_plus_1, curve_idx+1, F1);

        T cos2 = std::cos(theta) * std::cos(theta);
        T sin2 = std::sin(theta) * std::sin(theta);
        ci = cos2 * F0 + sin2 * F1;
        
        TV dF0dt, dF1dt;
        TV d2F0dt2, d2F1dt2;
        T dtheta_dt = M_PI * 0.5;

        if (compute_dci)
        {
            dF(ti, curve_idx, dF0dt); 
            dF(ti_plus_1, curve_idx + 1, dF1dt); 
            
            TV dF_dtheta = TV::Zero();
            dF_dtheta += 2.0 * std::cos(theta) * std::sin(theta) * (F1 - F0)  * dtheta_dt;
            dF_dtheta += cos2 * dF0dt * dtidt ;
            
            dF_dtheta += sin2 * dF1dt * dti_plus_1dt; // write it out for clearity
            
            dci = dF_dtheta;
        }
        if (compute_ddci)
        {
            ddF(ti, curve_idx, d2F0dt2); 
            ddF(ti_plus_1, curve_idx+1, d2F1dt2); 

            TV d2cdt2 = TV::Zero();
            d2cdt2 += 2.0 * (cos2 - sin2) * (F1 - F0)  * dtheta_dt * dtheta_dt;
            d2cdt2 += 4.0 * std::cos(theta) * std::sin(theta) * (dF1dt * dti_plus_1dt - dF0dt * dtidt)  * dtheta_dt;
            d2cdt2 += cos2 * d2F0dt2 * dtidt * dtidt + sin2 * d2F1dt2 * dti_plus_1dt * dti_plus_1dt;

            ddci = d2cdt2;
        }
    }
    else
    {
        T ti;
        T dtidt;
        if (curve_idx == 0)
        {
            ti = (0 + t * (0 - curve_data.front()->limits[0])) / (curve_data.front()->limits[1] - curve_data.front()->limits[0]);
            dtidt= (0 - curve_data.front()->limits[0]) / (curve_data.front()->limits[1] - curve_data.front()->limits[0]);
            F(ti, curve_idx, ci);
            if (compute_dci) dF(ti, curve_idx, dci); 
            if (compute_ddci) ddF(ti, curve_idx, ddci);
        }
        else
        {
            ti = (0 + t * (curve_data.back()->limits[1] - 0) - curve_data.back()->limits[0]) / (curve_data.back()->limits[1] - curve_data.back()->limits[0]);
            dtidt = (curve_data.back()->limits[1] - 0) / (curve_data.back()->limits[1] - curve_data.back()->limits[0]);
            F(ti, curve_data.size() - 1, ci);
            if (compute_dci) dF(ti, curve_data.size() - 1, dci); 
            if (compute_ddci) ddF(ti, curve_data.size() - 1, ddci);
        }
        

        dci *= dtidt;
        ddci *= dtidt * dtidt;
    }
}

template class HybridC2Curve<3>;
template class HybridC2Curve<2>;