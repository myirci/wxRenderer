#include "MyCanvas.hpp"
#include "MyFrame.hpp"
#include "../renderer/MyRenderer.hpp"
#include "../geometry/LineSegment3D.hpp"
#include "../geometry/Point3D.h"
#include "../geometry/Square3D.hpp"
#include "../data/teapot.hpp"
#include "../utility/Utility.hpp"

BEGIN_EVENT_TABLE(MyCanvas, wxScrolledWindow)
EVT_PAINT(MyCanvas::OnPaint)
EVT_LEFT_DOWN(MyCanvas::OnMouseLeftClick)
EVT_RIGHT_DOWN(MyCanvas::OnMouseRightClick)
EVT_MOTION(MyCanvas::OnMouseMove)
EVT_KEY_DOWN(MyCanvas::OnKeyDown)
EVT_SIZE(MyCanvas::OnResize)
END_EVENT_TABLE()

MyCanvas::MyCanvas(MyFrame* parent) :
    wxScrolledWindow(parent,
                     wxID_ANY,
                     wxDefaultPosition,
                     wxDefaultSize,
                     wxHSCROLL | wxVSCROLL | wxFULL_REPAINT_ON_RESIZE),
    m_parent(parent),
    m_rtype(render_type::geometry),
    m_mindepth(0.0), m_maxdepth(0.0) {
    SetBackgroundColour(wxColour(*wxWHITE));
    m_myPens[0] = new wxPen(*wxBLUE, 1, wxPENSTYLE_DOT_DASH);
    m_myPens[1] = new wxPen(*wxBLACK, 3, wxPENSTYLE_SOLID);
    m_myPens[2] = new wxPen(*wxBLACK, 1, wxPENSTYLE_SOLID);

    m_near = 1;
    m_far = 100;
    m_fovy = 40;
    InitializeRenderer();
}

void MyCanvas::InitializeRenderer() {
    wxSize size = m_parent->GetClientSize();
    m_renderer = new MyRenderer(size.x, size.y, 25, 25, 800, 600);
    // m_renderer->glOrtho(-100,100,-100,100,1,100);
    m_renderer->gluPerspective(m_fovy, m_near, m_far);
    // m_renderer->glFrustrum(-1.0, 1.0, -1.0, 1.0, 2.0, 50.0);
    // m_renderer->print_converted_parameters(40.0, 1, 100);
}

MyCanvas::~MyCanvas() {
    if(m_renderer) {
        delete m_renderer;
        m_renderer = nullptr;
    }
    if(!m_geometry.empty()) {
        for(auto it = m_geometry.begin(); it != m_geometry.end(); ++it) {
            delete *it;
            *it = nullptr;
        }
        m_geometry.clear();
    }
}

void MyCanvas::InverseProjection(const std::vector<Vector3d>& data) {
    m_renderer->inverse_projection(data);
}

void MyCanvas::DisplacementTestConstantDepth() {
    // 1: project the first point
    Vector3d pt1(1.5, -3.4, -10);
    Vector3d pt1_prj = m_renderer->render_point_and_print_info(pt1);
    std::cout << "********************************" << std::endl;

    // 2: project the second point. oberve that the two projected points
    //    have the same depth.
    Vector3d pt2(1.3, -2.7, -10);
    Vector3d pt2_prj = m_renderer->render_point_and_print_info(pt2);
    std::cout << "********************************" << std::endl;

    // 3: print the displacement vectors
    Vector3d d1 = pt2_prj - pt1_prj;
    std::cout << "Direct: displacement vector in image coordinates:\n"
              <<  d1 << std::endl;
    Vector3d d2 = pt2 - pt1;
    std::cout << "Direct: displacement vector in camera coordinates:\n"
              << d2 << std::endl;

    std::pair<unsigned int, unsigned int> size;
    m_renderer->get_viewport_size(size);

    // 4: calculated displacement in the projection coordinates
    double k1 = (2 * m_near * tan(deg2rad(0.5 * m_fovy))) / size.second;
    double disp_prj_x = k1 * d1.x();
    double disp_prj_y = k1 * (-d1.y());

    // 5: depth in the projection coordinates:
    double zp = -(m_far + m_near) - (m_far * m_near) / pt1.z();

    // 6: calculated displacement in camera coordinates:
    double k2 = m_far / (m_far + m_near + zp);
    double disp_eye_x = k2 * disp_prj_x;
    double disp_eye_y = k2 * disp_prj_y;

    std::cout << "Calculated displacement in camera coordinates:\n"
              << disp_eye_x << ", " << disp_eye_y << ", " << "0 "<< std::endl;

    std::cout << "Calcuated eye coordinate: \n"
              << pt1.x() + disp_eye_x << " "
              << pt1.y() + disp_eye_y << " "
              << pt1.z() << std::endl;
}

void MyCanvas::DisplacementTestCoplanarPoints() {
    // 1: project a point which is on the plane (1.6, -0.5, 1, 33.7)
    double l1 = 1.6;
    double l2 = -0.5;
    double l4 = 33.7;

    Vector3d pt1(-4.7665, -5.45518, -28.8012);
    Vector3d pt1_prj = m_renderer->render_point_and_print_info(pt1);

    std::cout << "********************************" << std::endl;

    // 2: project a second point which lies on the same plane
    Vector3d pt2(-0.688672, -11.6069, -38.4016);
    Vector3d pt2_prj = m_renderer->render_point_and_print_info(pt2);
    std::cout << "********************************" << std::endl;

    // "3: print the displacement vectors
    Vector3d d1 = pt2_prj - pt1_prj;
    std::cout << "Displacement vector in image coordinates:\n"
              <<  d1 << std::endl;
    Vector3d d2 = pt2 - pt1;
    std::cout << "Displacement vector in camera coordinates:\n"
              << d2 << std::endl;

    std::pair<unsigned int, unsigned int> size;
    m_renderer->get_viewport_size(size);

    // 4: displacement in the projection coordinates
    double k1 = (2 * tan(deg2rad(0.5*m_fovy))) / size.second;
    double disp_prj_x = k1 * d1.x();
    double disp_prj_y = -k1 * d1.y();

    // 5: x, y coordinates and depth of the first point in the projection coordinates:
    double xp = -pt1.x() / pt1.z();
    double yp = -pt1.y() / pt1.z();
    double zp = -(m_far + m_near) - (m_far * m_near) / pt1.z();

    double xpp = xp + disp_prj_x;
    double ypp = yp + disp_prj_y;
    double zpp = zp - m_far*(l1/l4)*disp_prj_x - m_far*(l2/l4)*disp_prj_y;

    Vector3d calculated_eye_coord;
    calculated_eye_coord.x() = func(m_far, m_near, xpp, zpp);
    calculated_eye_coord.y() = func(m_far, m_near, ypp, zpp);
    calculated_eye_coord.z() = func(m_far, m_near, -1, zpp);

    std::cout << "Calcuated eye coordinate: \n"
              << calculated_eye_coord.x() << " "
              << calculated_eye_coord.y() << " "
              << calculated_eye_coord.z() << std::endl;
}

double MyCanvas::func(const double f, const double n, const double x, const double y) const {
    return (f*x) / (y+f+n);
}

void MyCanvas::DisplacementTestPointsOnA3DLine() {
    // Project a point which is on a given 3D Line
    Vector3d pt1(-1.512820512820513, -1.615384615384615, -10.74615384615385);
    Vector3d pt1_prj = m_renderer->render_point_and_print_info(pt1);

    std::cout << "********************************" << std::endl;

    // Project a second point which lies on the same 3D line
    Vector3d pt2(6.153846153846155, 2.384615384615385, -22.14615384615385);
    Vector3d pt2_prj = m_renderer->render_point_and_print_info(pt2);
    std::cout << "********************************" << std::endl;

    // direction vector for the 3D line
    Vector3d dirvec = pt2 - pt1;
    dirvec.normalize();

    // print the displacement vectors
    Vector3d d1 = pt2_prj - pt1_prj;
    std::cout << "Displacement vector in image coordinates:\n"
              <<  d1 << std::endl;
    Vector3d d2 = pt2 - pt1;
    std::cout << "Displacement vector in camera coordinates:\n"
              << d2 << std::endl;

    std::pair<unsigned int, unsigned int> size;
    m_renderer->get_viewport_size(size);

    // displacement in the projection coordinates
    double k1 = (2 * tan(deg2rad(0.5*40))) / size.second;
    double disp_prj_x = k1 * d1.x();
    double disp_prj_y = -k1 * d1.y();

    // x, y coordinates and depth of the first point in the projection coordinates:
    double xp = -pt1.x() / pt1.z();
    double yp = -pt1.y() / pt1.z();
    double zp = -101-100/pt1.z();

    // line equation for displacement on projection space
    double a = pt1.z()*dirvec.y() - pt1.y()*dirvec.z();
    double b = pt1.x()*dirvec.z() - pt1.z()*dirvec.x();
    double c = pt1.x()*1*dirvec.y() + pt1.x()*dirvec.z()*yp
            - pt1.y()*1*dirvec.x() - pt1.y()*dirvec.z()*xp
            - pt1.z()*dirvec.x()*yp + pt1.z()*dirvec.y()*xp;

    std::cout << "Line equation in the projection space" << std::endl;

}

void MyCanvas::PrintTestInfo() {
    std::cout.precision(20);
    m_renderer->print_glprojection_matrix();
    m_renderer->print_inverted_glprojection_matrix();
    m_renderer->print_viewport_transform_matrix();
    m_renderer->print_inverted_viewport_transform_matrix();
    m_renderer->print_screen_parameters();
    m_renderer->calculate_near_far_left_right_bottom_top_values_from_glprojection_matrix();
}

void MyCanvas::InsertLineSegment(const Vector3d& pt1,
                                 const Vector3d& pt2) {
    m_geometry.push_back(new LineSegment3D(pt1, pt2));
    Refresh();
}

void MyCanvas::InsertSquare(const Vector3d& normal,
                            const Vector3d& corner1,
                            const Vector3d& corner2) {
    m_geometry.push_back(new Square3D(normal, corner1, corner2));
    Refresh();
}

void MyCanvas::InsertPoints(const std::vector<Vector3d>& data) {
    for(int i = 0; i < data.size(); ++i) {
        m_geometry.push_back(new Point3D(data[i]));
    }
    Refresh();
}

void MyCanvas::ToggleRender(render_type rtype) {
    m_rtype = rtype;
    if(m_rtype == render_type::depth_image) {
        UpdateMinMaxDepthValues();
    }
    Refresh();
}

void MyCanvas::OnPaint(wxPaintEvent& event) {
    wxPaintDC dc(this);
    DrawMidpointLines(dc);
    if(m_rtype == render_type::geometry) {
        dc.SetPen(*m_myPens[1]);
        RenderGeometry(dc);
    }
    else if(m_rtype == render_type::depth_image) {
        RenderDepthImage(dc);
    }
    else if(m_rtype == render_type::teapot) {
        dc.SetPen(*m_myPens[2]);
        RenderTeapot(dc);
    }

    // RenderTest(dc);
}

void MyCanvas::DrawMidpointLines(wxPaintDC& dc) {
    dc.SetPen(*m_myPens[0]);
    std::pair<unsigned int, unsigned int> vpsize, vpcorner;
    m_renderer->get_viewport_size(vpsize);
    m_renderer->get_viewport_corner(vpcorner);
    ToImageCoordinates(vpcorner);
    dc.DrawLine(vpcorner.first, vpcorner.second,
                vpcorner.first + vpsize.first, vpcorner.second);
    dc.DrawLine(vpcorner.first, vpcorner.second - vpsize.second,
                vpcorner.first + vpsize.first, vpcorner.second - vpsize.second);
    dc.DrawLine(vpcorner.first, vpcorner.second,
                vpcorner.first, vpcorner.second - vpsize.second);
    dc.DrawLine(vpcorner.first + vpsize.first, vpcorner.second,
                vpcorner.first + vpsize.first, vpcorner.second - vpsize.second);
    dc.DrawLine(vpcorner.first, vpcorner.second - vpsize.second/2,
                vpcorner.first + vpsize.first, vpcorner.second - vpsize.second/2);
    dc.DrawLine(vpcorner.first + vpsize.first/2, vpcorner.second,
                vpcorner.first + vpsize.first/2, vpcorner.second - vpsize.second);
}

void MyCanvas::RenderGeometry(wxPaintDC& dc) {
    std::vector<Vector3d> data(40);
    std::vector<Vector3d> ptdata(1);
    for(auto it = m_geometry.begin(); it != m_geometry.end(); ++it) {
        if(dynamic_cast<Point3D*>(*it)) {
            (*it)->GenerateDataFill(ptdata, 1);
            m_renderer->render_points(ptdata, dc);
        }
        else {
            (*it)->GenerateDataFill(data, 40);
            m_renderer->render_points(data, dc);
        }
    }
}

void MyCanvas::RenderDepthImage(wxPaintDC& dc) {

    // set the pen
    wxPen pen(*wxBLACK, 5, wxPENSTYLE_SOLID);
    dc.SetPen(pen);

    // get the size info
    std::pair<unsigned int, unsigned int> vp_size, scr_size, vp_corner;
    m_renderer->get_viewport_size(vp_size);
    m_renderer->get_screen_size(scr_size);
    m_renderer->get_viewport_corner(vp_corner); // in OpenGL screen coordinates

    // get the depth buffer
    float* dbuffer = m_renderer->get_depth_buffer();

    float diff = m_maxdepth - m_mindepth;
    unsigned int limit = vp_size.first * vp_size.second;
    unsigned int val = 0;
    unsigned int x, y;

    for(int i = 0; i < limit; ++i) {
        // normalize the depth value
        val = static_cast<unsigned int>((255*(dbuffer[i] - m_mindepth))/diff);
        if(val != 255) {
            pen.SetColour(val, val, val);
            dc.SetPen(pen);
            x = (i % vp_size.first) + vp_corner.first;
            y = (i / vp_size.first) + scr_size.second - vp_size.second - vp_corner.second;
            dc.DrawLine(x, y, x, y);
        }
    }
}

void MyCanvas::RenderTeapot(wxPaintDC& dc) {
    m_renderer->render_points(teapot, num_verts_teapot, dc);
}

void MyCanvas::UpdateMinMaxDepthValues() {
    float* dbuffer = m_renderer->get_depth_buffer();
    std::pair<unsigned int, unsigned int> size;
    m_renderer->get_viewport_size(size);
    unsigned int limit = size.first * size.second;
    m_mindepth = m_renderer->get_depth_far_val();
    m_maxdepth = m_renderer->get_depth_near_val();
    for(int i = 0; i < limit; ++i) {
        if(dbuffer[i] > m_maxdepth) {
            m_maxdepth = dbuffer[i];
        }
        if(dbuffer[i] < m_mindepth) {
            m_mindepth = dbuffer[i];
        }
    }
}

void MyCanvas::RenderTest(wxPaintDC& dc) {
    dc.SetPen(*m_myPens[0]);
    Vector3d p0(-20, -5, -1.1);
    Vector3d p1(30, 7, -35.3);
    Vector3d ip0 = m_renderer->render_and_get(p0);
    Vector3d ip1 = m_renderer->render_and_get(p1);
    dc.DrawLine(ip0.x(), ip0.y(), ip1.x(), ip1.y());
}

void MyCanvas::OnResize(wxSizeEvent& event) {
    // Note that event.GetSize() == m_parent->GetClientSize()
    m_renderer->set_screen_width_and_height(event.GetSize().GetWidth(),
                                            event.GetSize().GetHeight());
}

void MyCanvas::OnMouseMove(wxMouseEvent& event) {
    wxString str;
    str << "WxWidgets Scr-Coord: (" << event.GetX() << ", " << event.GetY() << ")";
    m_parent->SetStatusText(str, 0);

    std::pair<unsigned int, unsigned int> mouse_pos(event.GetX(), event.GetY());
    ToImageCoordinates(mouse_pos);
    str = wxEmptyString;
    str << "OpenGL Scr-Coord: (" << mouse_pos.first << ", " << mouse_pos.second << ")";
    m_parent->SetStatusText(str, 1);
    Refresh();
}

void MyCanvas::OnMouseLeftClick(wxMouseEvent& event) { }

void MyCanvas::OnMouseRightClick(wxMouseEvent& event) { }

void MyCanvas::OnKeyDown(wxKeyEvent &event) {
    if(event.GetKeyCode() == WXK_ESCAPE) { }
}

void MyCanvas::ToImageCoordinates(std::pair<unsigned int, unsigned int>& pt) {
    wxSize size = m_parent->GetClientSize();
    pt.second = size.GetHeight() - pt.second - 1;
}

