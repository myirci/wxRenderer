#include "Screen.hpp"
#include <iostream>

Screen::Screen(unsigned int w, unsigned int h,
               unsigned int x0, unsigned int y0,
               unsigned int vpw, unsigned int vph,
               float dnval, float dfval) :
    m_w(w),
    m_h(h),
    m_x0(x0),
    m_y0(y0),
    m_vpw(vpw),
    m_vph(vph),
    m_dnval(dnval),
    m_dfval(dfval),
    m_bsize(0),
    m_dbuffer(nullptr) {

    update_viewport_transformation_matrix();
    update_coordinate_transformation_matrix();
    reset_depth_buffer();
}

Screen::~Screen() {
    if(m_dbuffer != nullptr) {
        delete[] m_dbuffer;
        m_dbuffer = nullptr;
    }
}

void Screen::reset_depth_buffer() {
    if(m_bsize == m_vpw * m_vph) { return; }

    m_bsize = m_vpw * m_vph;
    if(m_dbuffer != nullptr) {
        delete[] m_dbuffer;
        m_dbuffer = nullptr;
    }

    m_dbuffer = new float[m_bsize];
    for(int i = 0; i < m_bsize; ++i) {
        m_dbuffer[i] = m_dfval;
    }
}

double Screen::get_viewport_aspect_ratio() const {
    return static_cast<double>(m_vpw) / static_cast<double>(m_vph);
}

double Screen::get_screen_aspect_ratio() const {
    return static_cast<double>(m_w) / static_cast<double>(m_h);
}

void Screen::set_screen_size(unsigned int width,
                             unsigned int height) {
    double rw = static_cast<double>(m_vpw) / static_cast<double> (m_w);
    double rh = static_cast<double>(m_vph) / static_cast<double> (m_h);
    double rx0 = static_cast<double>(width) / static_cast<double> (m_w);
    double ry0 = static_cast<double>(height) / static_cast<double> (m_h);
    m_w = width;
    m_h = height;
    m_vpw = static_cast<unsigned int>(round(m_w * rw));
    m_vph = static_cast<unsigned int>(round(m_h * rh));
    m_x0 = static_cast<unsigned int>(round(m_x0 * rx0));
    m_y0 = static_cast<unsigned int>(round(m_y0 * ry0));
    update_coordinate_transformation_matrix();
    update_viewport_transformation_matrix();
    reset_depth_buffer();
}

void Screen::set_viewport(unsigned int x0,
                          unsigned int y0,
                          unsigned int width,
                          unsigned int heigth) {
    m_vpw = width;
    m_vph = heigth;
    m_x0 = x0;
    m_y0 = y0;
    update_viewport_transformation_matrix();
    reset_depth_buffer();
}

void Screen::set_depth_values(float nval, float fval) {
    m_dnval = nval;
    m_dfval = fval;
}

void Screen::update_viewport_transformation_matrix() {
    m_vptransform = Eigen::Matrix4d::Identity();
    m_vptransform(0,0) = m_vpw/2.0;
    m_vptransform(0,3) = m_x0 + m_vpw/2.0;
    m_vptransform(1,1) = m_vph/2.0;
    m_vptransform(1,3) = m_y0 + m_vph/2.0;
    m_vptransform(2,2) = (m_dfval - m_dnval)/2.0;
    m_vptransform(2,3) = (m_dfval + m_dnval)/2.0;
}

void Screen::update_coordinate_transformation_matrix() {
    m_ctransform = Eigen::Matrix4d::Identity();
    m_ctransform(1,1) = -1;
    m_ctransform(1,3) = m_h - 1;
}

std::ostream& operator << (std::ostream& out, const Screen& scr) {
    out << "window width: " << scr.m_w << std::endl;
    out << "window height: " << scr.m_h << std::endl;
    out << "viewport width: " << scr.m_vpw << std::endl;
    out << "viewport height: " << scr.m_vph << std::endl;
    out << "x0: " << scr.m_x0 << std::endl;
    out << "y0: " << scr.m_y0 << std::endl;
    out << "depth far val: " << scr.m_dfval << std::endl;
    out << "depth near val: " << scr.m_dnval << std::endl;
    return out;
}
