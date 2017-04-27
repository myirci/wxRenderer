#include "LineSegment3D.hpp"
#include <iostream>

LineSegment3D::LineSegment3D(const Vector3d& p0, const Vector3d& p1) : Line3D(p1 - p0, p0),
    pt0(p0), pt1(p1) { }

void LineSegment3D::GenerateDataInsert(std::vector<Vector3d>& data, int num_points) {
    Vector3d v = pt1 - pt0;
    for(int i = 0; i < num_points; ++i) {
        data.push_back(pt0 + ((1.0/(num_points-1))*i)*v);
    }
}

void LineSegment3D::GenerateDataFill(std::vector<Vector3d>& data, int num_points) {
    Vector3d v = pt1 - pt0;
    for(int i = 0; i < num_points; ++i) {
        data[i] = pt0 + ((1.0/(num_points-1))*i)*v;
    }
}

