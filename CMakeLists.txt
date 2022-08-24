CMAKE_MINIMUM_REQUIRED(VERSION 3.20)

PROJECT(wxRenderer)

SET(CMAKE_BUILD_TYPE Debug)

FIND_PACKAGE(wxWidgets COMPONENTS core base REQUIRED)
INCLUDE(${wxWidgets_USE_FILE})

# INCLUDE_DIRECTORIES(/usr/include/eigen3)
FIND_PACKAGE(Eigen3 3.3 REQUIRED NO_MODULE)

SET(SOURCE gui/MyApp.cpp
           gui/MyFrame.cpp
           gui/MyCanvas.cpp
           renderer/MyRenderer.cpp
           renderer/Screen.cpp
           utility/Utility.cpp
           geometry/Geometry.hpp
           geometry/LineSegment3D.cpp
           geometry/Point3D.cpp
           geometry/Plane.cpp
           geometry/Line3D.cpp
           geometry/Square3D.cpp
           geometry/Vector3D.hpp
           data/teapot.hpp)
             
ADD_EXECUTABLE(${PROJECT_NAME} WIN32 ${SOURCE})
TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${wxWidgets_LIBRARIES} Eigen3::Eigen)
