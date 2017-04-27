#include "Point3D.h"

Point3D::Point3D(Vector3d p) : pt(p) { }

void Point3D::GenerateDataFill(std::vector<Vector3d>& data,
                               int num_points) {
    for(int i = 0; i < num_points; ++i) {
        data[i] = pt;
    }
}

void Point3D::GenerateDataInsert(std::vector<Vector3d>& data,
                                 int num_points) {
    for(int i = 0; i < num_points; ++i) {
        data.push_back(pt);
    }
}
