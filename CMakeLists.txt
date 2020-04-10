CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

PROJECT(wxRenderer)

SET(CMAKE_BUILD_TYPE Debug)

ADD_DEFINITIONS(-std=c++17)



FIND_PACKAGE(wxWidgets REQUIRED core)
IF(wxWidgets_FOUND)
	INCLUDE(${wxWidgets_USE_FILE})
ENDIF(wxWidgets_FOUND)

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
             
ADD_EXECUTABLE(${PROJECT_NAME} ${SOURCE})
TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${wxWidgets_LIBRARIES} Eigen3::Eigen)
