#ifndef MY_RENDERER_HPP
#define MY_RENDERER_HPP

#include "../geometry/Vector3D.hpp"
#include "Screen.hpp"
#include <vector>
#include <Eigen/Dense>
#include <utility>

class wxPaintDC;

class MyRenderer {
public:
    // constructors:
    MyRenderer(unsigned int w, unsigned int h,
               unsigned int x0, unsigned int y0,
               unsigned int vpw, unsigned int vph,
               float dnval = 0.0, float dfval = 1.0);

    // setter and getters
    void set_screen_width_and_height(int w, int h);
    void set_glprojection_matrix(const Eigen::Matrix4d& mat);
    void get_glprojection_matrix(Eigen::Matrix4d& mat) const;
    void get_projection_matrix(Eigen::Matrix4d& mat) const;
    void get_inverse_glprojection_matrix(double fovy,
                                         double near, double far,
                                         Eigen::Matrix4d& mat);
    void get_inverse_glprojection_matrix(double left, double right,
                                         double bottom, double top,
                                         double near, double far,
                                         Eigen::Matrix4d& mat);
    void get_viewport_size(std::pair<unsigned int, unsigned int>& size) const;
    void get_viewport_corner(std::pair<unsigned int, unsigned int>& corner) const;
    void get_screen_size(std::pair<unsigned int, unsigned int>& size) const;
    float* get_depth_buffer();
    float get_depth_near_val();
    float get_depth_far_val();

    // rendering functions:
    Vector3d render_point_and_print_info(const Vector3d& pt) const;
    Vector3d render_and_get(const Vector3d& pt) const;
    void render_points(const double* data, const unsigned int num_verts, wxPaintDC& dc);
    void render_points(const std::vector<Vector3d>& data, wxPaintDC& dc);
    void inverse_projection(const std::vector<Vector3d>& data);

    // OpenGL API implementations
    void glFrustrum(double left, double right,
                    double bottom, double top,
                    double near, double far);
    void gluPerspective(double fovy, double near, double far);
    void gluPerspective_via_glFrustrum(double fovy, double near, double far);
    void glOrtho(double left, double right,
                 double bottom, double top,
                 double near, double far);
    void glViewport(unsigned int x, unsigned int y,
                    unsigned int width, unsigned int height);
    void glDepthRange(float nearVal, float farVal);

    // test and print fucntions
    void print_converted_parameters(double fovy, double near, double far) const;
    void print_converted_parameters(double left, double right,
                                    double bottom, double top,
                                    double near, double far) const;
    void print_glprojection_matrix() const;
    void print_inverted_glprojection_matrix() const;
    void print_viewport_transform_matrix() const;
    void print_inverted_viewport_transform_matrix() const;
    void print_screen_parameters() const;
    void calculate_near_far_left_right_bottom_top_values_from_glprojection_matrix() const;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
private:
    Screen m_screen;
    Eigen::Matrix4d m_glprojection;
    Eigen::Matrix4d m_projection;
    Eigen::Matrix4d m_modelview;
    void project(const Eigen::MatrixXd& vertices, wxPaintDC& dc);
    void project_with_matrices(const Eigen::MatrixXd& vertices, wxPaintDC& dc);
    void calculate_projection_matrix(double near, double far);
};

#endif // MY_RENDERER_HPP
