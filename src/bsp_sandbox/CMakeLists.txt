find_package(Qt4 REQUIRED)

include(${QT_USE_FILE})
ADD_DEFINITIONS(${QT_DEFINITIONS})

include_directories(
  SYSTEM
  ${BULLET_INCLUDE_DIRS}  
  ${Boost_INCLUDE_DIRS}
  ${Eigen_INCLUDE_DIRS}
  ${OpenRAVE_INCLUDE_DIRS}
  ${JSON_INCLUDE_DIRS}
)

SET(toy_qt_headers toy.hpp)
QT4_WRAP_CPP(toy_qt_headers_moc ${toy_qt_headers})
add_executable(toy toy.cpp ${toy_qt_headers_moc})
target_link_libraries(toy bsp trajopt ${QT_LIBRARIES})

SET(arm_qt_headers arm.hpp)
QT4_WRAP_CPP(arm_qt_headers_moc ${arm_qt_headers})
add_executable(arm arm.cpp geometry_2d.cpp ${arm_qt_headers_moc})
target_link_libraries(arm bsp trajopt ${QT_LIBRARIES} ${Boost_LIBRARIES})

add_executable(geometry_2d_test geometry_2d_test.cpp geometry_2d.cpp)
target_link_libraries(geometry_2d_test bsp trajopt)

SET(car_qt_headers car.hpp)
QT4_WRAP_CPP(car_qt_headers_moc ${car_qt_headers})
add_executable(car car.cpp geometry_2d.cpp ${car_qt_headers_moc})
target_link_libraries(car bsp trajopt ${QT_LIBRARIES})
