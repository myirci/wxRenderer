#ifndef VECTOR3D_HPP
#define VECTOR3D_HPP

#include <cmath>
#include <ostream>

template <typename T>
class Vector3D {
public:
    Vector3D();
    Vector3D(T x, T y, T z);
    Vector3D(const Vector3D<T>& rhs);
    T x() const;
    T y() const;
    T z() const;
    T& x();
    T& y();
    T& z();
    double norm() const;
    void normalize();
    Vector3D<T> normalized() const;
    double dot(const Vector3D<T>& rhs) const;
    Vector3D<T> cross(const Vector3D<T>& rhs) const;
    Vector3D<T>& operator= (const Vector3D<T>& rhs);
    Vector3D<T>& operator+=(const Vector3D<T>& rhs);
    Vector3D<T>& operator-=(const Vector3D<T>& rhs);
    Vector3D<T>& operator*=(const double s);
    Vector3D<T>& operator/=(const double s);
private:
    T m_x;
    T m_y;
    T m_z;
};

template <typename T>
Vector3D<T>::Vector3D() : m_x(T()), m_y(T()), m_z(T()) { }

template <typename T>
Vector3D<T>::Vector3D(T x, T y, T z) : m_x(x), m_y(y), m_z(z) { }

template <typename T>
Vector3D<T>::Vector3D(const Vector3D<T>& rhs) {
    this->m_x = rhs.x();
    this->m_y = rhs.y();
    this->m_z = rhs.z();
}

template <typename T> T Vector3D<T>::x() const { return m_x; }
template <typename T> T Vector3D<T>::y() const { return m_y; }
template <typename T> T Vector3D<T>::z() const { return m_z; }
template <typename T> T& Vector3D<T>::x() { return m_x; }
template <typename T> T& Vector3D<T>::y() { return m_y; }
template <typename T> T& Vector3D<T>::z() { return m_z; }

template <typename T>
double Vector3D<T>::dot(const Vector3D<T>& rhs) const{
    return m_x*rhs.x() + m_y*rhs.y() + m_z*rhs.z();
}

template <typename T>
double Vector3D<T>::norm() const {
    return std::sqrt(this->dot(*this));
}

template <typename T>
void Vector3D<T>::normalize() {
    *this /= this->norm();
}

template <typename T>
Vector3D<T> Vector3D<T>::normalized() const{
    Vector3D<T> vec(*this);
    vec.normalize();
    return vec;
}

template <typename T>
Vector3D<T> Vector3D<T>::cross(const Vector3D<T>& rhs) const {
    return Vector3D<T>(this->m_y*rhs.z() - this->m_z*rhs.y(),
                       this->m_z*rhs.x() - this->m_x*rhs.z(),
                       this->m_x*rhs.y() - this->m_y*rhs.x());
}

template <typename T>
Vector3D<T>& Vector3D<T>::operator+= (const Vector3D<T>& rhs) {
    this->m_x += rhs.x();
    this->m_y += rhs.y();
    this->m_z += rhs.z();
    return *this;
}

template <typename T>
Vector3D<T>& Vector3D<T>::operator-= (const Vector3D<T>& rhs) {
    this->m_x -= rhs.x();
    this->m_y -= rhs.y();
    this->m_z -= rhs.z();
    return *this;
}

template <typename T>
Vector3D<T>& Vector3D<T>::operator*= (const double s) {
    this->m_x *= s;
    this->m_y *= s;
    this->m_z *= s;
    return *this;
}

template <typename T>
Vector3D<T>& Vector3D<T>::operator/= (const double s) {
    this->m_x /= s;
    this->m_y /= s;
    this->m_z /= s;
    return *this;
}

template <typename T>
Vector3D<T>& Vector3D<T>::operator= (const Vector3D<T>& rhs) {
    this->m_x = rhs.x();
    this->m_y = rhs.y();
    this->m_z = rhs.z();
    return *this;
}

// Non-member functions:
template <typename T>
Vector3D<T> operator+(const Vector3D<T>& v1, const Vector3D<T>& v2) {
    return Vector3D<T>(v1.x() + v2.x(),
                       v1.y() + v2.y(),
                       v1.z() + v2.z());
}

template <typename T>
Vector3D<T> operator-(const Vector3D<T>& v1, const Vector3D<T>& v2) {
    return Vector3D<T>(v1.x() - v2.x(),
                       v1.y() - v2.y(),
                       v1.z() - v2.z());
}

template <typename T>
Vector3D<T> operator*(const Vector3D<T>& v, const double s) {
    return Vector3D<T>(v.x()*s,
                       v.y()*s,
                       v.z()*s);
}

template <typename T>
Vector3D<T> operator*(const double s, const Vector3D<T>& v) {
    return v*s;
}

template <typename T>
std::ostream& operator<<(std::ostream& out, const Vector3D<T>& v) {
    out << v.x() << " " << v.y() << " " << v.z() << std::endl;
    return out;
}

typedef Vector3D<double> Vector3d;
typedef Vector3D<float> Vector3f;

#endif // VECTOR3D_HPP
