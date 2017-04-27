#ifndef LINE_SEGMENT_3D_HPP
#define LINE_SEGMENT_3D_HPP

#include "Geometry.hpp"
#include "Line3D.hpp"

class LineSegment3D : public Line3D, public Geometry {
public:
    LineSegment3D(const Vector3d& p0, const Vector3d& p1);
    void GenerateDataInsert(std::vector<Vector3d>& data, int num_points);
    void GenerateDataFill(std::vector<Vector3d>& data, int num_points);
private:
    Vector3d pt0, pt1;
};

#endif // LINE_SEGMENT_3D_HPP


