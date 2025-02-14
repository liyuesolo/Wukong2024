#include "../include/Util.h"

Eigen::Matrix3d rotationMatrixFromEulerAngle(T angle_z, T angle_y, T angle_x)
{
    Eigen::Matrix3d R, yaw, pitch, roll;
    yaw.setZero(); pitch.setZero(); roll.setZero();
    yaw(0, 0) = cos(angle_z);	yaw(0, 1) = -sin(angle_z);
    yaw(1, 0) = sin(angle_z);	yaw(1, 1) = cos(angle_z);
    yaw(2, 2) = 1.0;
    //y rotation
    pitch(0, 0) = cos(angle_y); pitch(0, 2) = sin(angle_y);
    pitch(1, 1) = 1.0;
    pitch(2, 0) = -sin(angle_y); pitch(2, 2) = cos(angle_y);
    //x rotation
    roll(0, 0) = 1.0;
    roll(1, 1) = cos(angle_x); roll(1, 2) = -sin(angle_x);
    roll(2, 1) = sin(angle_x); roll(2, 2) = cos(angle_x);
    R = yaw * pitch * roll;
    return R;
}

Eigen::Vector3d parallelTransport(const Eigen::Vector3d& u, const Eigen::Vector3d& t0, const Eigen::Vector3d& t1) 
{
  
    Eigen::Vector3d b = t0.cross(t1);


    if(b.norm() < std::numeric_limits<T>::epsilon())
        return u;

    b.normalize();

    Eigen::Vector3d n0 = t0.cross(b).normalized();
    Eigen::Vector3d n1 = t1.cross(b).normalized();

    return u.dot(t0.normalized()) * t1.normalized() + u.dot(n0) * n1 +
            u.dot(b) * b;
}

Eigen::Vector3d parallelTransportOrthonormalVector(const Eigen::Vector3d& u, const Eigen::Vector3d& t0, const Eigen::Vector3d& t1) {
  
    Eigen::Vector3d b = t0.cross(t1);


    if(b.norm() < std::numeric_limits<T>::epsilon())
        return u;

    b.normalize();

    Eigen::Vector3d n0 = t0.cross(b);
    Eigen::Vector3d n1 = t1.cross(b);

    return u.dot(n0) * n1 + u.dot(b) * b;
}

//http://paulbourke.net/geometry/circlesphere/tvoght.c
bool circleCircleIntersection(const Eigen::Vector3d& x0, T r0,
                               const Eigen::Vector3d& x1, T r1,
                               Eigen::Vector3d& ixn0,
                               Eigen::Vector3d& ixn1)
{
    ixn0 = Eigen::Vector3d::Zero();
    ixn1 = Eigen::Vector3d::Zero();
    T a, dx, dy, d, h, rx, ry;
    T x2, y2;

    /* dx and dy are the vertical and horizontal distances between
    * the circle centers.
    */
    dx = x1[0] - x0[0];
    dy = x1[1] - x0[1];

    /* Determine the straight-line distance between the centers. */
    //d = sqrt((dy*dy) + (dx*dx));
    d = std::hypot(dx,dy); // Suggested by Keith Briggs

    /* Check for solvability. */
    if (d > (r0 + r1))
    {
    /* no solution. circles do not intersect. */
    return false;
    }
    if (d < fabs(r0 - r1))
    {
    /* no solution. one circle is contained in the other */
    return false;
    }

    /* 'point 2' is the point where the line through the circle
    * intersection points crosses the line between the circle
    * centers.  
    */

    /* Determine the distance from point 0 to point 2. */
    a = ((r0*r0) - (r1*r1) + (d*d)) / (2.0 * d) ;

    /* Determine the coordinates of point 2. */
    x2 = x0[0] + (dx * a/d);
    y2 = x0[1] + (dy * a/d);

    /* Determine the distance from point 2 to either of the
    * intersection points.
    */
    h = std::sqrt((r0*r0) - (a*a));

    /* Now determine the offsets of the intersection points from
    * point 2.
    */
    rx = -dy * (h/d);
    ry = dx * (h/d);

    /* Determine the absolute intersection points. */
    ixn0[0] = x2 + rx;
    ixn1[0] = x2 - rx;
    ixn0[1] = y2 + ry;
    ixn1[1] = y2 - ry;

    return true;
}