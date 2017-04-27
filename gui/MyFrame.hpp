#include <wx/wx.h>

#define wxID_GEOMETRY_RENDER_LINE                   wxID_HIGHEST + 1
#define wxID_GEOMETRY_RENDER_POINT                  wxID_HIGHEST + 2
#define wxID_GEOMETRY_RENDER_SQUARE                 wxID_HIGHEST + 3
#define wxID_TEST_PRINT_MATRICES                    wxID_HIGHEST + 4
#define wxID_RENDER_DEPTH_IMAGE                     wxID_HIGHEST + 5
#define wxID_RENDER_GEOMETRY                        wxID_HIGHEST + 6
#define wxID_RENDER_TEAPOT                          wxID_HIGHEST + 7
#define wxID_PROJECTION_INVERSE_PROJECTION          wxID_HIGHEST + 8
#define wxID_TEST_DISPLACEMENT_CONSTANT_DEPTH       wxID_HIGHEST + 9
#define wxID_TEST_DISPLACEMENT_COPLANAR_POINTS      wxID_HIGHEST + 10
#define wxID_TEST_DISPLACEMENT_POINTS_ON_A_LINE     wxID_HIGHEST + 11

class MyCanvas;

class MyFrame : public wxFrame {
public:
    MyFrame(const wxString& title);
private:
    MyCanvas* m_canvas;
    void OnQuit(wxCommandEvent& event);
    void OnGeometryInsertLine(wxCommandEvent& event);
    void OnGeometryInsertPoint(wxCommandEvent& event);
    void OnGeometryInsertSquare(wxCommandEvent& event);
    void OnPrintTestInfo(wxCommandEvent& event);
    void OnDisplacementTestConstantDepth(wxCommandEvent& event);
    void OnDisplacementTestCoplanarPoints(wxCommandEvent& event);
    void OnDisplacementTestPointsOnA3DLine(wxCommandEvent& event);
    void OnToggleRender(wxCommandEvent& event);
    void OnInverseProjection(wxCommandEvent& event);
    DECLARE_EVENT_TABLE()
};		
