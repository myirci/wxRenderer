#include "Plane.hpp"
#include <iostream>

Plane::Plane(const double a, const double b, const double c, const double d) {
    m_coeff = new double[4];
    m_coeff[0] = a;
    m_coeff[1] = b;
    m_coeff[2] = c;
    m_coeff[3] = d;
    m_normal = Vector3d(a,b,c);
    m_normal.normalize();
}

Plane::Plane(const Vector3d& n, const Vector3d& pt) {
    m_coeff = new double[4];
    m_coeff[0] = n.x();
    m_coeff[1] = n.y();
    m_coeff[2] = n.z();
    m_coeff[3] = -n.dot(pt);
    m_normal = n.normalized();
}

Plane::~Plane() {
    if(m_coeff != nullptr) {
        delete[] m_coeff;
        m_coeff = nullptr;
    }
}

bool Plane::IsOnPlane(const Vector3d& pt) const {
    double val = m_coeff[0]*pt.x() + m_coeff[1]*pt.y() + m_coeff[2]*pt.z() + m_coeff[3];
    return val < 0.0001;
}

std::ostream& operator<<(std::ostream& out, const Plane& pl) {
    out << "normal: " << pl.m_normal << std::endl;
    out << "coeff : " << pl.m_coeff[0] << " "
                      << pl.m_coeff[1] << " "
                      << pl.m_coeff[2] << " "
                      << pl.m_coeff[3] << std::endl;
    return out;
}
