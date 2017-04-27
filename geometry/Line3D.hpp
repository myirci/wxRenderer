#ifndef LINE3D_HPP
#define LINE3D_HPP

#include "Vector3D.hpp"

class Line3D {
public:
    Line3D(const Vector3d& dv, const Vector3d& pt);
    const bool IsPointOnLine(const Vector3d& pt) const;
    double GetParameterForX(double x) const;
    double GetParameterForY(double y) const;
    double GetParameterForZ(double z) const;
private:
    Vector3d m_dvec;
    Vector3d m_pt;
};

#endif // LINE3D_HPP
