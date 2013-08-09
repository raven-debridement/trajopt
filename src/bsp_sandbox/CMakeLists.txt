find_package(Qt4 REQUIRED)

include(${QT_USE_FILE})
add_definitions(${QT_DEFINITIONS})
add_definitions(-DDATA_DIR="${CMAKE_SOURCE_DIR}/data")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")

include_directories(
  SYSTEM
  ${BULLET_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
  ${Eigen_INCLUDE_DIRS}
  ${OpenRAVE_INCLUDE_DIRS}
  ${JSON_INCLUDE_DIRS}
)

#SET(toy_qt_headers toy.hpp)
#QT4_WRAP_CPP(toy_qt_headers_moc ${toy_qt_headers})
#add_executable(toy toy.cpp geometry_2d.cpp ${toy_qt_headers_moc})
#target_link_libraries(toy bsp bsp_qt trajopt ${QT_LIBRARIES})

#SET(arm_qt_headers arm.hpp)
#QT4_WRAP_CPP(arm_qt_headers_moc ${arm_qt_headers})
#add_executable(arm arm.cpp geometry_2d.cpp ${arm_qt_headers_moc})
#target_link_libraries(arm bsp bsp_qt trajopt ${QT_LIBRARIES} ${Boost_LIBRARIES})

#add_executable(geometry_2d_test geometry_2d_test.cpp geometry_2d.cpp)
#target_link_libraries(geometry_2d_test bsp trajopt)

#SET(car_qt_headers car.hpp)
#QT4_WRAP_CPP(car_qt_headers_moc ${car_qt_headers})
#add_executable(car car.cpp geometry_2d.cpp ${car_qt_headers_moc})
#target_link_libraries(car bsp bsp_qt trajopt ${QT_LIBRARIES})

#SET(car_static_qt_headers car_static.hpp)
#QT4_WRAP_CPP(car_static_qt_headers_moc ${car_static_qt_headers})
#add_executable(car_static car_static.cpp geometry_2d.cpp ${car_static_qt_headers_moc})
#target_link_libraries(car_static bsp bsp_qt trajopt ${QT_LIBRARIES})

#SET(needle_qt_headers needle.hpp)
#QT4_WRAP_CPP(needle_qt_headers_moc ${needle_qt_headers})
#add_executable(needle needle.cpp geometry_2d.cpp ${needle_qt_headers_moc})
#target_link_libraries(needle bsp bsp_qt trajopt ${QT_LIBRARIES})

add_executable(barrett_robot barrett_robot.cpp)
target_link_libraries(barrett_robot bsp trajopt)

#add_executable(four_links_robot four_links_robot.cpp)
#target_link_libraries(four_links_robot bsp trajopt)

#add_executable(point_robot point_robot.cpp)
#target_link_libraries(point_robot bsp trajopt)

add_executable(raven raven.cpp)
target_link_libraries(raven bsp trajopt)

#add_library(ravenbsp raven.cpp)
#target_link_libraries(ravenbsp bsp trajopt)

#include_directories(${PYTHON_NUMPY_INCLUDE_DIR})
#boost_python_module(cravenbsppy ravenpy.cpp)
#target_link_libraries(cravenbsppy ravenbsp bsp trajopt ${Boost_SYSTEM_LIBRARY} ${OpenRAVE_BOTH_LIBRARIES})
