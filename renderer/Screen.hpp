#ifndef SCREEN_HPP
#define SCREEN_HPP

#include <Eigen/Dense>
#include <ostream>

class Screen {
public:
    Screen(unsigned int w, unsigned int h,
           unsigned int x0, unsigned int y0,
           unsigned int vpw, unsigned int vph,
           float dnval = 0.0, float dfval = 1.0);
    ~Screen();
    unsigned int get_screen_width() const { return m_w; }
    unsigned int get_screen_height() const { return m_h; }
    unsigned int get_viewport_width() const { return m_vpw; }
    unsigned int get_viewport_height() const { return m_vph; }
    unsigned int get_x0() const { return m_x0; }
    unsigned int get_y0() const { return m_y0; }
    float get_depth_near_val() const { return m_dnval; }
    float get_depth_far_val() const { return m_dfval; }
    const Eigen::Matrix4d& get_viewport_transform_matrix() const { return m_vptransform; }
    const Eigen::Matrix4d& get_coordinate_transform_matrix() const { return m_ctransform; }
    float* get_depth_buffer() { return m_dbuffer; }
    double get_viewport_aspect_ratio() const;
    double get_screen_aspect_ratio() const;
    void set_screen_size(unsigned int width, unsigned int heigth);
    void set_viewport(unsigned int x0, unsigned int y0,
                      unsigned int width, unsigned int heigth);
    void set_depth_values(float nval, float fval);
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
private:
    unsigned int m_w, m_h;     // width and height of the screen
    unsigned int m_vpw, m_vph; // viewport width and height
    unsigned int m_x0, m_y0;   // vivewport rectangle lower-left coordinates
    float m_dnval;             // near clipping plane depth val
    float m_dfval;             // far clipping plane depth val
    Eigen::Matrix4d m_vptransform;  // viewport transformation matrix
    Eigen::Matrix4d m_ctransform;   // coordinate transform matrix (OpenGL to WxWidgets)
    float* m_dbuffer;
    unsigned int m_bsize;
    void update_viewport_transformation_matrix();
    void update_coordinate_transformation_matrix();
    void update_depth_map_coordinate_transformation_matrix();
    void reset_depth_buffer();
    friend std::ostream& operator << (std::ostream& out, const Screen& scr);
};

std::ostream& operator << (std::ostream& out, const Screen& scr);

#endif // SCREEN_HPP
