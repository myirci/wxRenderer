#include "MyRenderer.hpp"
#include "../utility/Utility.hpp"
#include <iostream>
#include <wx/dcclient.h>
#include <Eigen/StdVector>

/* ***********************************************
 * Constructors
 ************************************************/

MyRenderer::MyRenderer(unsigned int w, unsigned int h,
                       unsigned int x0, unsigned int y0,
                       unsigned int vpw, unsigned int vph,
                       float dnval, float dfval) :
    m_screen(w, h, x0, y0, vpw, vph, dnval, dfval),
    m_glprojection(Eigen::Matrix4d::Identity()),
    m_projection(Eigen::Matrix4d::Identity()),
    m_modelview(Eigen::Matrix4d::Identity()) { }

/* ***********************************************
 * Rendering
 ************************************************/

void MyRenderer::render_points(const double* data,
                               const unsigned int num_verts,
                               wxPaintDC& dc) {

    // Combine the vertex data into a matrix
    Eigen::MatrixXd vertices(4, num_verts);
    int k = 0;
    for(int i = 0; i < num_verts; ++i) {
        vertices(0,i) = data[k++];
        vertices(1,i) = data[k++];
        vertices(2,i) = data[k++];
        vertices(3,i) = 1.0;
    }
    // project(vertices, dc);
    project_with_matrices(vertices, dc);
}

void MyRenderer::render_points(const std::vector<Vector3d>& data,
                               wxPaintDC& dc) {

    Eigen::MatrixXd vertices(4, data.size());
    int i = 0;
    for(auto it = data.begin(); it != data.end(); ++it, ++i) {
        vertices(0,i) = it->x();
        vertices(1,i) = it->y();
        vertices(2,i) = it->z();
        vertices(3,i) = 1.0;
    }
    // project(vertices, dc);
    project_with_matrices(vertices, dc);

}

Vector3d MyRenderer::render_point_and_print_info(const Vector3d& pt) const {

    // 1) copy the input point into a 4x1 vector
    Eigen::Vector4d p(pt.x(), pt.y(), pt.z(), 1);
    std::cout << "Point in camera coordinates:\n" << p << std::endl;
    std::cout << "--------------------------------" << std::endl;

    // 2) perform projection
    Eigen::Vector4d projected = m_projection * p;
    std::cout << "Point in projection coordinates:\n" << projected << std::endl;
    std::cout << "--------------------------------" << std::endl;

    // 3) perform clipping coordinates transformation
    Eigen::Vector4d clip = m_glprojection * p;
    std::cout << "Point in clipping coordinates:\n" << clip << std::endl;
    std::cout << "--------------------------------" << std::endl;

    // 4) perform clipping
    if(clip(0) < -clip(3) || clip(0) > clip(3) ||
       clip(1) < -clip(3) || clip(1) > clip(3) ||
       clip(2) < -clip(3) || clip(2) > clip(3)) {
        std::cout << "Point is clipped" << std::endl;
        std::cout << "--------------------------------" << std::endl;
    }

    // 5) perform perspective division
    clip(0) /= clip(3);
    clip(1) /= clip(3);
    clip(2) /= clip(3);
    clip(3) = 1;

    // 6) perform viewport mapping
    Eigen::Vector4d glscreen_coord = m_screen.get_viewport_transform_matrix() * clip;
    std::cout << "Point in OpenGL screen coordinates:\n" << glscreen_coord << std::endl;
    std::cout << "--------------------------------" << std::endl;

    // 7) perform coordinate transformation and depth mapping
    Eigen::Vector4d img_coord = m_screen.get_coordinate_transform_matrix() * glscreen_coord;
    std::cout << "Point in image coordinates:\n" << img_coord << std::endl;
    std::cout << "--------------------------------" << std::endl;

    return Vector3d(img_coord(0), img_coord(1), img_coord(2));
}

Vector3d MyRenderer::render_and_get(const Vector3d& pt) const {
    Eigen::Vector4d p(pt.x(), pt.y(), pt.z(), 1);
    Eigen::Vector4d clip = m_glprojection * p;
    if(clip(0) < -clip(3) || clip(0) > clip(3) ||
       clip(1) < -clip(3) || clip(1) > clip(3) ||
       clip(2) < -clip(3) || clip(2) > clip(3)) {
    }
    clip(0) /= clip(3);
    clip(1) /= clip(3);
    clip(2) /= clip(3);
    clip(3) = 1;
    Eigen::Vector4d glscreen_coord = m_screen.get_viewport_transform_matrix() * clip;
    Eigen::Vector4d img_coord = m_screen.get_coordinate_transform_matrix() * glscreen_coord;
    return Vector3d(img_coord(0),img_coord(1),img_coord(2));
}

void MyRenderer::project(const Eigen::MatrixXd& vertices, wxPaintDC& dc) {
    // Perform projection
    Eigen::MatrixXd projected = m_glprojection * vertices;

    // get the depth buffer
    float* depth_buffer = m_screen.get_depth_buffer();

    int px(0), py(0);
    double x(0.0), y(0.0), z(0.0), w(0.0);
    double half_vpw = static_cast<double>(m_screen.get_viewport_width()) / 2.0;
    double half_vph = static_cast<double>(m_screen.get_viewport_height()) / 2.0;
    float var1 = (m_screen.get_depth_far_val() - m_screen.get_depth_near_val()) / 2.0;
    float var2 = (m_screen.get_depth_far_val() + m_screen.get_depth_near_val()) / 2.0;
    int index = 0;
    int xx = 0;
    int yy = 0;
    // clipping, perspective division and viewport mapping
    for(int i = 0; i < projected.cols(); i++) {
        // these coordinates are called clipping coordinates on which the
        // clipping is applied.
        x = projected.col(i)(0);
        y = projected.col(i)(1);
        z = projected.col(i)(2);
        w = projected.col(i)(3);
        if(x < -w || x > w || y < -w || y > w || z < -w || z > w) { continue; }

        // if clipping test is passed, perform perspective division to find
        // normalized device coorniates
        x /= w;
        y /= w;
        z /= w;

        // after getting the normalized device coordinates, viewport mapping is
        // done to get the pixel coordinates. note that OpenGL screen coordinate
        // frame origin is located at the bottom-left corner of the screen. The
        // calculated px and py are according to the OpenGL coordinate system.
        px = static_cast<int>(half_vpw*(x+1) + m_screen.get_x0());
        py = static_cast<int>(half_vph*(y+1) + m_screen.get_y0());

        // we need to convert from OpenGL screen coordinate system to image coordinate system.
        // In this program, image coordinate system is the wxWidget's windows coordinate system.
        // The origin of this coordinat frame is located at the top left corner of the window.
        // x-coordinate does not change but we need to modify y-coordinate.
        py = m_screen.get_screen_height() - py - 1;

        // perform the depth mapping and fill in the depth buffer
        xx = py - m_screen.get_screen_height() + m_screen.get_viewport_height() + m_screen.get_y0();
        yy = px - m_screen.get_x0();
        index = xx * m_screen.get_viewport_width() + yy;
        depth_buffer[index] = z*var1 + var2;

        // render the point
        if(dc.GetPen().GetWidth() == 1) {
            dc.DrawPoint(px, py);
        }
        else {
            dc.DrawLine(px, py, px, py);
        }
    }
}

void MyRenderer::project_with_matrices(const Eigen::MatrixXd& vertices,
                                       wxPaintDC& dc) {
    // 1) perform glprojection
    Eigen::MatrixXd projected = m_glprojection * vertices;

    // 2) perform clipping and perspective division
    std::vector <Eigen::Vector4d, Eigen::aligned_allocator<Eigen::Vector4d>> ndc;
    for(int i = 0; i < projected.cols(); i++) {
        // these coordinates are called clipping coordinates on which the
        // clipping is applied.
        if(projected.col(i)(0) < -projected.col(i)(3) ||
           projected.col(i)(0) > projected.col(i)(3)  ||
           projected.col(i)(1) < -projected.col(i)(3) ||
           projected.col(i)(1) > projected.col(i)(3)  ||
           projected.col(i)(2) < -projected.col(i)(3) ||
           projected.col(i)(2) > projected.col(i)(3)) { continue; }

        // if clipping test is passed, perform perspective division to find
        // normalized device coorniates
        ndc.push_back(Eigen::Vector4d(projected.col(i)(0) / projected.col(i)(3),
                                      projected.col(i)(1) / projected.col(i)(3),
                                      projected.col(i)(2) / projected.col(i)(3),
                                      1));
    }

    // 3) perform viewport mapping, image coordinate transformation, rendering,
    // and store depth values

    float* depth_buffer = m_screen.get_depth_buffer();
    Eigen::Vector4d img_coord;
    unsigned int index = 0;
    unsigned int x, y;
    int var = m_screen.get_viewport_height() + m_screen.get_y0() - m_screen.get_screen_height();
    for(int i = 0; i < ndc.size(); ++i) {
        img_coord = m_screen.get_coordinate_transform_matrix() *
                    m_screen.get_viewport_transform_matrix() * ndc[i];
        x = static_cast<unsigned int>(round(img_coord(0)));
        y = static_cast<unsigned int>(round(img_coord(1)));
        if(dc.GetPen().GetWidth() == 1) {
            dc.DrawPoint(x, y);
        }
        else {
            dc.DrawLine(x, y, x, y);
        }
        index = (y + var) * m_screen.get_viewport_width() + x - m_screen.get_x0();
        depth_buffer[index] = static_cast<float>(img_coord(2));
    }
}

void MyRenderer::inverse_projection(
        const std::vector<Vector3d>& data) {

    // Step-1: Combine the given vertices into 4 x (data.size) matrix
    Eigen::MatrixXd vertices(4, data.size());
    for(int i = 0; i < data.size(); ++i) {
        vertices(0,i) = data[i].x();
        vertices(1,i) = data[i].y();
        vertices(2,i) = data[i].z();
        vertices(3,i) = 1.0;
    }

    // Step-2:
    // a) Convert from image coordinates to OpenGL screen coordinates.
    // Note that inverse of image coordinate transform matrix is equal to itself.
    // b) Convert from OpenGL screen coordinates to ndc
    // c) Convert from ndc to eye coordinates
    Eigen::MatrixXd eye_coord =
            m_glprojection.inverse() *
            m_screen.get_viewport_transform_matrix().inverse() *
            m_screen.get_coordinate_transform_matrix() * vertices;

    std::cout.precision(16);
    // Step-3:Print the converted cooridnates:
    for(int i = 0; i < eye_coord.cols(); ++i) {
        std::cout << eye_coord.col(i)(0) / eye_coord.col(i)(3) << ", "
                  << eye_coord.col(i)(1) / eye_coord.col(i)(3) << ", "
                  << eye_coord.col(i)(2) / eye_coord.col(i)(3) << std::endl;
    }
}

/* ***********************************************
 * Accessing and modifying the matrices
 ************************************************/

void MyRenderer::set_screen_width_and_height(int w, int h) {
    m_screen.set_screen_size(w, h);
}

void MyRenderer::set_glprojection_matrix(const Eigen::Matrix4d& mat) {
    m_glprojection = mat;
}

void MyRenderer::get_glprojection_matrix(Eigen::Matrix4d& mat) const {
    mat = m_glprojection;
}

void MyRenderer::get_projection_matrix(Eigen::Matrix4d& mat) const {
    mat = m_projection;
}

void MyRenderer::calculate_projection_matrix(double near, double far) {
    m_projection = Eigen::Matrix4d::Zero();
    m_projection(0,0) = near;
    m_projection(1,1) = near;
    m_projection(2,2) = far + near;
    m_projection(2,3) = far*near;
    m_projection(3,2) = -1;
}

void MyRenderer::get_inverse_glprojection_matrix(double fovy, double near,
                                                 double far, Eigen::Matrix4d& mat) {
    mat = Eigen::Matrix4d::Zero();
    mat(0,0) = m_screen.get_viewport_aspect_ratio() * tan(deg2rad(0.5*fovy));
    mat(1,1) = tan(deg2rad(0.5*fovy));
    mat(2,3) = -1;
    mat(3,2) = (near-far)/(2*far*near);
    mat(3,3) = (far+near)/(2*far*near);
}

void MyRenderer::get_inverse_glprojection_matrix(double left, double right,
                                                 double bottom, double top,
                                                 double near, double far,
                                                 Eigen::Matrix4d& mat) {
    mat = Eigen::Matrix4d::Zero();
    mat(0,0) = (right-left) / (2*near);
    mat(0,3) = (right+left) / (2*near);
    mat(1,1) = (top-bottom) / (2*near);
    mat(1,3) = (top+bottom) / (2*near);
    mat(2,3) = -1;
    mat(3,2) = (near-far)/(2*far*near);
    mat(3,3) = (near+far)/(2*far*near);
}

void MyRenderer::get_viewport_size(std::pair<unsigned int, unsigned int>& size) const {
    size.first = m_screen.get_viewport_width();
    size.second = m_screen.get_viewport_height();
}

void MyRenderer::get_viewport_corner(std::pair<unsigned int, unsigned int>& corner) const {
    corner.first = m_screen.get_x0();
    corner.second = m_screen.get_y0();
}

void MyRenderer::get_screen_size(std::pair<unsigned int, unsigned int>& size) const {
    size.first = m_screen.get_screen_width();
    size.second = m_screen.get_screen_height();
}

float* MyRenderer::get_depth_buffer() {
    return m_screen.get_depth_buffer();
}

float MyRenderer::get_depth_near_val() {
    return m_screen.get_depth_near_val();
}

float MyRenderer::get_depth_far_val() {
    return m_screen.get_depth_far_val();
}

/* ***********************************************
 * Printing and testing functions
 ************************************************/

void MyRenderer::print_glprojection_matrix() const {
    std::cout << "glprojection matrix: " << std::endl;
    std::cout << m_glprojection << std::endl;
}

void MyRenderer::print_viewport_transform_matrix() const {
    std::cout << "viewport transform matrix: " << std::endl;
    std::cout << m_screen.get_viewport_transform_matrix() << std::endl;
}

void MyRenderer::print_inverted_viewport_transform_matrix() const {
    std::cout << "inverted viewport transform matrix: " << std::endl;
    if(m_screen.get_viewport_transform_matrix().determinant() != 0) {
        std::cout << m_screen.get_viewport_transform_matrix().inverse() << std::endl;
    }
    else {
        std::cout << "viewport tansform matrix is not invertible" << std::endl;
    }
}

void MyRenderer::print_inverted_glprojection_matrix() const {
    std::cout << "inverted glprojection matrix: " << std::endl;
    if(m_glprojection.determinant() != 0) {
        std::cout << m_glprojection.inverse() << std::endl;
    }
    else {
        std::cout << "glprojection matrix is not invertible" << std::endl;
    }
}

void MyRenderer::print_screen_parameters() const {
    std::cout << "screen parameters: " << std::endl;
    std::cout << m_screen << std::endl;
}

void MyRenderer::print_converted_parameters(double fovy, double near, double far) const {
    double half_fov = deg2rad(0.5*fovy);
    double t = tan(half_fov)*near;
    double b = -t;
    double r = t * m_screen.get_viewport_aspect_ratio();
    double l = -r;
    std::cout << "left: " << l << std::endl;
    std::cout << "right: " << r << std::endl;
    std::cout << "bottom: " << b << std::endl;
    std::cout << "top: " << t << std::endl;
    std::cout << "near: " << near << std::endl;
    std::cout << "far: " << far << std::endl;
}

void MyRenderer::print_converted_parameters(double left, double right,
                                          double bottom, double top,
                                          double near, double far) const {
    if(right != -left || top != -bottom) {
        std::cout << "Viewing volume is not symmetric" << std::endl;
        return;
    }
    double fovy = rad2deg(2*atan(top/near));
    std::cout << "fovy: " << fovy << std::endl;
    std::cout << "aspect: " << m_screen.get_viewport_aspect_ratio() << std::endl;
    std::cout << "near: " << near << std::endl;
    std::cout << "far: " << far << std::endl;
}

void MyRenderer::calculate_near_far_left_right_bottom_top_values_from_glprojection_matrix() const {
    std::cout <<"calculated near, far, bottom, top, left, right values from glprojection matrix: " << std::endl;

    Eigen::Vector4d left = m_glprojection.row(0) + m_glprojection.row(3);
    Eigen::Vector4d right = m_glprojection.row(0) - m_glprojection.row(3);
    Eigen::Vector4d bottom = m_glprojection.row(1) + m_glprojection.row(3);
    Eigen::Vector4d top = m_glprojection.row(1) - m_glprojection.row(3);
    Eigen::Vector4d near = m_glprojection.row(2) + m_glprojection.row(3);
    Eigen::Vector4d far = m_glprojection.row(2) - m_glprojection.row(3);

    std::cout << "near = " << near(3)/near(2) << std::endl;
    std::cout << "far  = " << far(3)/far(2) << std::endl;
    std::cout << "left = " << (left(2)/left(0))*(near(3)/near(2)) << std::endl;
    std::cout << "right = " << (right(2)/right(0))*(near(3)/near(2)) << std::endl;
    std::cout << "bottom = " << (bottom(2)/bottom(1))*(near(3)/near(2)) << std::endl;
    std::cout << "top = " << (top(2)/top(1))*(near(3)/near(2)) << std::endl;
}

/* ***********************************************
 * Implementation of OpenGL API:
 *************************************************/

void MyRenderer::glFrustrum(double left, double right,
                          double bottom, double top,
                          double near, double far) {
    m_glprojection = Eigen::Matrix4d::Zero();
    m_glprojection(0,0) = (2 * near) / (right - left);
    m_glprojection(0,2) = (right + left) / (right - left);
    m_glprojection(1,1) = (2 * near) / (top - bottom);
    m_glprojection(1,2) = (top + bottom) / (top - bottom);
    m_glprojection(2,2) = -(far + near) / (far - near);
    m_glprojection(2,3) = (-2*far*near) / (far - near);
    m_glprojection(3,2) = -1;
    calculate_projection_matrix(far, near);
}

void MyRenderer::gluPerspective_via_glFrustrum(double fovy, double near, double far) {
    double t = tan(deg2rad(0.5*fovy))*near;
    double b = -t;
    double r = t * m_screen.get_viewport_aspect_ratio();
    double l = -r;
    this->glFrustrum(l, r, b, t, near, far);
}

void MyRenderer::gluPerspective(double fovy, double near, double far) {
    m_glprojection = Eigen::Matrix4d::Zero();
    double half_fov = deg2rad(0.5*fovy);
    m_glprojection(0,0) = 1 / (tan(half_fov) * m_screen.get_viewport_aspect_ratio());
    m_glprojection(1,1) = 1 / tan(half_fov);
    m_glprojection(2,2) = -(far + near) / (far - near);
    m_glprojection(2,3) = (-2*far*near) / (far - near);
    m_glprojection(3,2) = -1;
    calculate_projection_matrix(far, near);
}

void MyRenderer::glOrtho(double left, double right,
                         double bottom, double top,
                         double near, double far) {
    m_glprojection = Eigen::Matrix4d::Zero();
    m_glprojection(0,0) = (2) / (right - left);
    m_glprojection(0,3) = -(right + left) / (right - left);
    m_glprojection(1,1) = (2) / (top - bottom);
    m_glprojection(1,3) = -(top + bottom) / (top - bottom);
    m_glprojection(2,2) = (-2) / (far - near);
    m_glprojection(2,3) = -(far + near) / (far - near);
    m_glprojection(3,3) = 1;
    calculate_projection_matrix(far, near);
}

void MyRenderer::glViewport(unsigned int x, unsigned int y,
                            unsigned int width, unsigned int height) {
    m_screen.set_viewport(x, y, width, height);
}

void MyRenderer::glDepthRange(float nearVal, float farVal) {
    m_screen.set_depth_values(nearVal, farVal);
}


