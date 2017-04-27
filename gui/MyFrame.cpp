#include "MyFrame.hpp"
#include "MyCanvas.hpp"

BEGIN_EVENT_TABLE(MyFrame, wxFrame)
    EVT_MENU(wxID_EXIT, MyFrame::OnQuit)
    EVT_MENU(wxID_TEST_PRINT_MATRICES, MyFrame::OnPrintTestInfo)
    EVT_MENU(wxID_TEST_DISPLACEMENT_CONSTANT_DEPTH, MyFrame::OnDisplacementTestConstantDepth)
    EVT_MENU(wxID_TEST_DISPLACEMENT_COPLANAR_POINTS, MyFrame::OnDisplacementTestCoplanarPoints)
    EVT_MENU(wxID_TEST_DISPLACEMENT_POINTS_ON_A_LINE, MyFrame::OnDisplacementTestPointsOnA3DLine)
    EVT_MENU(wxID_GEOMETRY_RENDER_LINE, MyFrame::OnGeometryInsertLine)
    EVT_MENU(wxID_GEOMETRY_RENDER_POINT, MyFrame::OnGeometryInsertPoint)
    EVT_MENU(wxID_GEOMETRY_RENDER_SQUARE, MyFrame::OnGeometryInsertSquare)
    EVT_MENU(wxID_PROJECTION_INVERSE_PROJECTION, MyFrame::OnInverseProjection)
    EVT_MENU(wxID_RENDER_GEOMETRY, MyFrame::OnToggleRender)
    EVT_MENU(wxID_RENDER_DEPTH_IMAGE, MyFrame::OnToggleRender)
    EVT_MENU(wxID_RENDER_TEAPOT, MyFrame::OnToggleRender)
END_EVENT_TABLE()

MyFrame::MyFrame(const wxString& title) :
    wxFrame(NULL, wxID_ANY, title, wxDefaultPosition, wxSize(850, 650)) {
    m_canvas = new MyCanvas(this);

    wxMenuBar* menubar = new wxMenuBar;
    wxMenu* file = new wxMenu;
    file->Append(wxID_EXIT, wxT("&Quit"));
    menubar->Append(file, wxT("&File"));

    wxMenu* render = new wxMenu;
    render->AppendRadioItem(wxID_RENDER_GEOMETRY, wxT("Geometry"));
    render->AppendRadioItem(wxID_RENDER_DEPTH_IMAGE, wxT("Depth Image"));
    render->AppendRadioItem(wxID_RENDER_TEAPOT, wxT("Teapot"));
    menubar->Append(render, "&Render");

    wxMenu* geometry = new wxMenu;
    geometry->Append(wxID_GEOMETRY_RENDER_POINT, wxT("Insert 3D Point"));
    geometry->Append(wxID_GEOMETRY_RENDER_LINE, wxT("Insert 3D Line Segment"));
    geometry->Append(wxID_GEOMETRY_RENDER_SQUARE, wxT("Insert 3D Square"));
    menubar->Append(geometry, wxT("&Geometry"));

    wxMenu* projection = new wxMenu;
    projection->Append(wxID_PROJECTION_INVERSE_PROJECTION, wxT("Inverse Project"));
    menubar->Append(projection, wxT("&Projection"));

    wxMenu* test = new wxMenu;
    test->Append(wxID_TEST_PRINT_MATRICES, wxT("Print Test Info"));
    wxMenu* disp = new wxMenu;
    disp->Append(wxID_TEST_DISPLACEMENT_CONSTANT_DEPTH, wxT("Displacement Test - Constant Depth"));
    disp->Append(wxID_TEST_DISPLACEMENT_COPLANAR_POINTS, wxT("Displacement Test - Coplanar Points"));
    disp->Append(wxID_TEST_DISPLACEMENT_POINTS_ON_A_LINE, wxT("Displacement Test - Points on a 3D Line"));
    test->AppendSubMenu(disp, wxT("Displacement"));
    menubar->Append(test, wxT("&Test"));
    SetMenuBar(menubar);

    CreateStatusBar(2);
    int widths[2] = {-1, -1};
    SetStatusWidths(2, widths);

    SetClientSize(850, 650);
	this->Centre();
}

void MyFrame::OnQuit(wxCommandEvent& event) {
    Close();
}

#include "../geometry/Plane.hpp"

void MyFrame::OnGeometryInsertPoint(wxCommandEvent& event) {  
    std::vector<Vector3d> data;
    data.push_back(Vector3d(-1.512820512820513, -1.615384615384615, -10.74615384615385));
    data.push_back(Vector3d(6.153846153846155, 2.384615384615385, -22.14615384615385));
    m_canvas->InsertPoints(data);
}

void MyFrame::OnGeometryInsertLine(wxCommandEvent& event) {
    m_canvas->InsertLineSegment(Vector3d(-20, -5, -1.1),
                                Vector3d(30, 7, -35.3));
}

void MyFrame::OnGeometryInsertSquare(wxCommandEvent& event) {
     m_canvas->InsertSquare(Vector3d(1.6,-0.5,1),
                            Vector3d(-2,1,-30),
                            Vector3d(3,-3,-40));
}

void MyFrame::OnPrintTestInfo(wxCommandEvent& event) {
    m_canvas->PrintTestInfo();
}

void MyFrame::OnDisplacementTestConstantDepth(wxCommandEvent& event) {
    m_canvas->DisplacementTestConstantDepth();
}

void MyFrame::OnDisplacementTestCoplanarPoints(wxCommandEvent& event) {
    m_canvas->DisplacementTestCoplanarPoints();
}

void MyFrame::OnDisplacementTestPointsOnA3DLine(wxCommandEvent& event) {
    m_canvas->DisplacementTestPointsOnA3DLine();
}

void MyFrame::OnToggleRender(wxCommandEvent& event) {
    if(event.GetId() == wxID_RENDER_GEOMETRY) {
        m_canvas->ToggleRender(render_type::geometry);
    }
    else if(event.GetId() == wxID_RENDER_DEPTH_IMAGE) {
        m_canvas->ToggleRender(render_type::depth_image);
    } else if(event.GetId() == wxID_RENDER_TEAPOT) {
        m_canvas->ToggleRender(render_type::teapot);
    }
}

void MyFrame::OnInverseProjection(wxCommandEvent& event) {
    std::vector<Vector3d> data;
    data.push_back(Vector3d(523, 579, 0.909091));
    m_canvas->InverseProjection(data);
}
