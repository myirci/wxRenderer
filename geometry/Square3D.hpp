#ifndef SQUARE3D_HPP
#define SQUARE3D_HPP

#include "Plane.hpp"
#include "Geometry.hpp"

class Square3D : public Plane, public Geometry {
public:
    Square3D(const double a,
             const double b,
             const double c,
             const double d,
             const Vector3d& corner1,
             const Vector3d& corner2);

    Square3D(const Vector3d& n,
             const Vector3d& corner1,
             const Vector3d& corner2);

    void GenerateDataInsert(std::vector<Vector3d>& data,
                            int num_points);

    void GenerateDataFill(std::vector<Vector3d>& data,
                          int num_points);
private:
    Vector3d m_start;
    Vector3d m_vec1;
    Vector3d m_vec2;
    double m_width;
};

#endif // SQUARE3D_HPP
