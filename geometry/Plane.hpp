#ifndef PLANE_HPP
#define PLANE_HPP

#include "Vector3D.hpp"

class Plane {
public:
    Plane(const double a, const double b, const double c, const double d);
    Plane(const Vector3d& n, const Vector3d& pt);
    virtual ~Plane();
    bool IsOnPlane(const Vector3d& pt) const;
    friend std::ostream& operator<<(std::ostream& out, const Plane& pl);
protected:
    Vector3d m_normal;
private:
    double* m_coeff;
};

std::ostream& operator<<(std::ostream& out, const Plane& pl);

#endif // PLANE_HPP
