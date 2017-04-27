#include "Line3D.hpp"
#include <iostream>

Line3D::Line3D(const Vector3d& dv, const Vector3d& pt) : m_dvec(dv), m_pt(pt) { }

const bool Line3D::IsPointOnLine(const Vector3d& pt) const {
    double t1 = GetParameterForX(pt.x());
    double t2 = GetParameterForY(pt.y());
    double t3 = GetParameterForZ(pt.z());
    double eps = 0.00001;
    double d1 = std::abs(t1 - t2);
    double d2 = std::abs(t1 - t3);
    double d3 = std::abs(t2 - t3);
    return (d1 < eps && d2 < eps && d3 < eps);
}

double Line3D::GetParameterForX(double x) const {
    return (x - m_pt.x()) / m_dvec.x();
}

double Line3D::GetParameterForY(double y) const {
    return (y - m_pt.y()) / m_dvec.y();
}

double Line3D::GetParameterForZ(double z) const {
    return (z - m_pt.z()) / m_dvec.z();
}
