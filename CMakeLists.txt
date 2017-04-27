cmake_minimum_required(VERSION 2.8)

project(wxRenderer)

set(CMAKE_BUILD_TYPE Debug)

add_definitions(-std=c++11)

include_directories(/usr/include/eigen3)

find_package(wxWidgets REQUIRED core)
if(wxWidgets_FOUND)
	include(${wxWidgets_USE_FILE})
endif(wxWidgets_FOUND)

set(SOURCE gui/MyApp.cpp
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
             
add_executable(MyRenderer ${SOURCE})
target_link_libraries(MyRenderer ${wxWidgets_LIBRARIES})
