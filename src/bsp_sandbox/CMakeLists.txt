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

#SET(car_qt_headers car.hpp)
#QT4_WRAP_CPP(car_qt_headers_moc ${car_qt_headers})
#add_executable(car car.cpp ${car_qt_headers_moc})
#target_link_libraries(car bsp trajopt ${QT_LIBRARIES})
