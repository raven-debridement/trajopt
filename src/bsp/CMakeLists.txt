set(RBDL_INCLUDE_DIRS "${CMAKE_SOURCE_DIR}/ext/rbdl/src" "${CMAKE_BINARY_DIR}/ext/rbdl/src")

find_package(Qt4 REQUIRED)
find_package(Qwt REQUIRED)

include(${QT_USE_FILE})
ADD_DEFINITIONS(${QT_DEFINITIONS})
include_directories(
  SYSTEM
  ${BULLET_INCLUDE_DIRS}  
  ${Boost_INCLUDE_DIRS}
  ${Eigen_INCLUDE_DIRS}
  ${OpenRAVE_INCLUDE_DIRS}
  ${JSON_INCLUDE_DIRS}
  ${RBDL_INCLUDE_DIRS}
  ${QWT_INCLUDE_DIR}
)

SET(toy_qt_HEADERS toy_optimizer_task.hpp toy_plotter.hpp)
QT4_WRAP_CPP(toy_qt_HEADERS_MOC ${toy_qt_HEADERS})
add_executable(toy toy_main.cpp toy.cpp costs.cpp optimizers.cpp utils.cpp toy_optimizer_task.cpp toy_plotter.cpp ${toy_qt_HEADERS_MOC})
target_link_libraries(toy trajopt ${QT_LIBRARIES} ${QWT_LIBRARY})
