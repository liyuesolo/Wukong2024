#ifndef PARAMETERICYARN_H
#define PARAMETERICYARN_H


#include <utility> 
#include <iostream>
#include <fstream>
#include <Eigen/Geometry>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <tbb/tbb.h>
#include <unordered_map>
#include <iomanip>


#include "VecMatDef.h"
#include "Util.h"

class ParametericYarn
{
public:
	using VectorXT = Matrix<T, Eigen::Dynamic, 1>;
	using MatrixXT = Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
	using VectorXi = Vector<int, Eigen::Dynamic>;
	using MatrixXi = Matrix<int, Eigen::Dynamic, Eigen::Dynamic>;
	using VtxList = std::vector<int>;
	using StiffnessMatrix = Eigen::SparseMatrix<T>;
	using Entry = Eigen::Triplet<T>;
	using TV = Vector<T, 3>;
	using TV2 = Vector<T, 2>;
	using TM2 = Matrix<T, 2, 2>;
	using TV3 = Vector<T, 3>;
	using IV = Vector<int, 3>;
	using IV2 = Vector<int, 2>;
	using TM = Matrix<T, 3, 3>;

	struct Vector3
	{
		double x;
		double y;
		double z;
	};

	double sqr( double x )
	{
		return x*x;
	}

	double norm2( struct Vector3 v )
	{
		return sqr(v.x) + sqr(v.y) + sqr(v.z);
	}

	double norm( struct Vector3 v )
	{
		return sqrt( norm2( v ));
	}

	struct Vector3 scale( double c, struct Vector3 v )
	{
		struct Vector3 cv;
		cv.x = c * v.x;
		cv.y = c * v.y;
		cv.z = c * v.z;
		return cv;
	}

	struct Vector3 cross( struct Vector3 u, struct Vector3 v )
	{
		struct Vector3 w;
		w.x = u.y*v.z - u.z*v.y;
		w.y = u.z*v.x - u.x*v.z;
		w.z = u.x*v.y - u.y*v.x;
		return w;
	}

	struct Vector3 yarnCurve( double t, double a, double h, double d )
	{
		struct Vector3 gamma_t;

		gamma_t.x = t + a*sin(2.*t);
		gamma_t.y = h * cos(t);
		gamma_t.z = d * cos(2.*t);

		return gamma_t;
	}

	void frenetFrame( double t, double a, double h, double d,
					struct Vector3* e1, struct Vector3* e2, struct Vector3* e3 )
	{
		double u_t, v_t, x_t, y_t;

		e1->x = 1. + 2.*a*cos(2.*t);
		e1->y = -h*sin(t);
		e1->z = -2.*d*sin(2.*t);

		u_t = norm2(*e1);
		v_t = 2.*h*h*cos(t)*sin(t) + 16.*d*d*cos(2.*t)*sin(2.*t) -
				8.*a*(1. + 2.*a*cos(2.*t))*sin(2.*t);
		x_t = 1./sqrt(u_t);
		y_t = v_t/(2.*pow(u_t,3./2.));

		e2->x = y_t*(-1. - 2.*a*cos(2.*t)) - x_t*4.*a*sin(2.*t);
		e2->y = y_t*h*sin(t) - x_t*h*cos(t);
		e2->z = y_t*2.*d*sin(2.*t) - x_t*4.*d*cos(2.*t);

		*e1 = scale(  x_t, *e1 );
		*e2 = scale( 1./norm(*e2), *e2 );
		*e3 = cross( *e1, *e2 );
	}

	struct Vector3 fiberCurve( double t, double a, double h, double d,
							double r, double omega, double phi )
	{
		struct Vector3 gamma_t;
		struct Vector3 eta_t;
		struct Vector3 e1, e2, e3;
		double theta_t;

		gamma_t = yarnCurve( t, a, h, d );
		frenetFrame( t, a, h, d, &e1, &e2, &e3 );
		theta_t = t*omega - 2.*cos(t) + phi;

		eta_t.x = gamma_t.x + r*( cos(theta_t)*e2.x + sin(theta_t)*e3.x );
		eta_t.y = gamma_t.y + r*( cos(theta_t)*e2.y + sin(theta_t)*e3.y );
		eta_t.z = gamma_t.z + r*( cos(theta_t)*e2.z + sin(theta_t)*e3.z );

		return eta_t;
	}

	VectorXT deformed;

public:
	ParametericYarn() 
	{
		std::vector<TV> points;
		const int nRows = 16;           /* number of rows generated */
		const int nLoops = 12;          /* number of loops in each row */
		const int samplesPerLoop = 64;  /* points sampled per period */
		const int nFibers = 4;          /* number of twisted fibers */
		const double a = 3/2.;          /* loop roundness */
		const double d = 1.;            /* loop depth */
		const double h = 4.;            /* loop height */
		const double rowOffset = h + 1./2.;     /* row spacing */
		const double r = 1./2.;         /* yarn radius */
		const double omega = 5.;        /* fiber twist */
		double t;
		int row, loop, sample;
		int nPoints;
		double dt, y0, t0;
		struct Vector3 gamma_t;
		nPoints = 0;
   		dt = (2.*M_PI)/(double)samplesPerLoop;
		for( row = 0; row < nRows; row++ )
		{
			y0 = rowOffset * (double)row;

			for( loop = 0; loop < nLoops; loop++ )
			{
				t0 = 2.*M_PI * loop;
				for( sample = 0; sample < samplesPerLoop; sample++ )
				{
					t = t0 + dt*(double)sample;

					gamma_t = yarnCurve( t, a, h, d );
					points.emplace_back(gamma_t.x,
						gamma_t.y + y0,
						gamma_t.z);
					
				}
			}
		}
		deformed.resize(points.size() * 3);
		for (int i = 0; i < points.size(); i++)
		{
			deformed.segment<3>(i * 3) = points[i];
		}
		
	} 
	~ParametericYarn() {} 
};


#endif
