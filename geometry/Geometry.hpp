#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include <vector>
#include "Vector3D.hpp"

class Geometry {
public:
    Geometry() { }
    virtual ~Geometry() { }
    virtual void GenerateDataInsert(
            std::vector<Vector3d>& data,
            int num_points) = 0;
    virtual void GenerateDataFill(
            std::vector<Vector3d>& data,
            int num_points) = 0;
};

#endif // GEOMETRY_HPP
