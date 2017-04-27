#ifndef POINT3D_H
#define POINT3D_H

#include "Geometry.hpp"

class Point3D : public Geometry {
public:
    Point3D(Vector3d p);
    void GenerateDataInsert(std::vector<Vector3d>& data,
                            int num_points);
    void GenerateDataFill(std::vector<Vector3d>& data,
                          int num_points);
private:
    Vector3d pt;
};

#endif // POINT3D_H
