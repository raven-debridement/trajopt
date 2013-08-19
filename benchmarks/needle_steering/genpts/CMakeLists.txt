cmake_minimum_required(VERSION 2.8)
 
PROJECT(genpts)
 
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
 
add_executable(genpts MACOSX_BUNDLE genpts)
 
target_link_libraries(genpts ${VTK_LIBRARIES})
