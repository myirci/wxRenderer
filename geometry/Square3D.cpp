#include "Square3D.hpp"
#include <iostream>

Square3D::Square3D(const double a,
                         const double b,
                         const double c,
                         const double d,
                         const Vector3d& corner1,
                         const Vector3d& corner2) :
    Plane(a, b, c, d), m_start(corner1) {
    m_vec1 = corner2 - corner1;
    m_width = m_vec1.norm();
    m_vec2 = m_vec1.cross(m_normal);
    m_vec1.normalize();
    m_vec2.normalize();
}

Square3D::Square3D(const Vector3d& n,
                   const Vector3d& corner1,
                   const Vector3d& corner2) :
    Plane(n, corner1), m_start(corner1) {
    m_vec1 = corner2 - corner1;
    m_width = m_vec1.norm();
    m_vec2 = m_vec1.cross(m_normal);
    m_vec1.normalize();
    m_vec2.normalize();
}

void Square3D::GenerateDataInsert(std::vector<Vector3d>& data,
                                  int num_points) {
    unsigned int n = std::floor(std::sqrt(num_points));
    double step = m_width / (n - 1);
    Vector3d point1 = m_start;
    Vector3d point2;
    for(int j = 0; j < n; ++j) {
        for(int i = 0; i < n; ++i) {
            point2 = point1 + i*step*m_vec1;
            data.push_back(point2);
        }
        point1 = m_start + (j+1)*step*m_vec2;
    }

    int remaining = num_points - n*n;
    if(remaining > 0) {
        for(int i = n*n; i < num_points; ++i) {
            data.push_back(data[0]);
        }
    }
}

void Square3D::GenerateDataFill(std::vector<Vector3d>& data,
                                int num_points) {
    unsigned int n = std::floor(std::sqrt(num_points));
    double step = m_width / (n - 1);
    Vector3d point1 = m_start;
    Vector3d point2;
    for(int j = 0; j < n; ++j) {
        for(int i = 0; i < n; ++i) {
            point2 = point1 + i*step*m_vec1;
            data[j*n + i] = point2;
            IsOnPlane(data[j*n + i]);
        }
        point1 = m_start + (j+1)*step*m_vec2;
    }

    int remaining = num_points - n*n;
    if(remaining > 0) {
        for(int i = n*n; i < num_points; ++i) {
            data[i] = data[0];
        }
    }
}
