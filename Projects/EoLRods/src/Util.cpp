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