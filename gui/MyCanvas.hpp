#ifndef MYCANVAS_HPP
#define MYCANVAS_HPP

#include "../geometry/Geometry.hpp"
#include <wx/scrolwin.h>
#include <vector>

class MyFrame;
class wxPaintDC;
class MyRenderer;

enum class render_type : unsigned char {
    geometry,
    depth_image,
    teapot
};

class MyCanvas : public wxScrolledWindow {
public:
    MyCanvas(MyFrame* parent);
    ~MyCanvas();
    void InsertPoints(const std::vector<Vector3d>& data);
    void InsertLineSegment(const Vector3d& pt1,
                           const Vector3d& pt2);
    void InsertSquare(const Vector3d& normal,
                      const Vector3d& corner1,
                      const Vector3d& corner2);
    void ToggleRender(render_type rtype);
    void PrintTestInfo();
    void InverseProjection(const std::vector<Vector3d>& data);
    // Tests for regularization
    void DisplacementTestConstantDepth();
    void DisplacementTestCoplanarPoints();
    void DisplacementTestPointsOnA3DLine();
private:
    MyFrame* m_parent;
    MyRenderer* m_renderer;
    std::vector<Geometry*> m_geometry;
    render_type m_rtype;
    wxPen* m_myPens[3];
    float m_mindepth;
    float m_maxdepth;
    double m_near;
    double m_far;
    double m_fovy;
    void OnPaint(wxPaintEvent& event);
    void OnMouseLeftClick(wxMouseEvent& event);
    void OnMouseRightClick(wxMouseEvent& event);
    void OnMouseMove(wxMouseEvent& event);
    void OnKeyDown(wxKeyEvent& event);
    void OnResize(wxSizeEvent& event);
    void RenderGeometry(wxPaintDC& dc);
    void RenderDepthImage(wxPaintDC& dc);
    void RenderTeapot(wxPaintDC& dc);
    void RenderTest(wxPaintDC& dc);
    void UpdateMinMaxDepthValues();
    inline void InitializeRenderer();
    inline void DrawMidpointLines(wxPaintDC& dc);
    inline double func(const double f, const double n, const double x, const double y) const;
    void ToImageCoordinates(std::pair<unsigned int, unsigned int>& pt);
    DECLARE_EVENT_TABLE()
};

#endif // MYCANVAS_HPP
