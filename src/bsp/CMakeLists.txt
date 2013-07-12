find_package(Qt4 REQUIRED)

include(${QT_USE_FILE})
ADD_DEFINITIONS(${QT_DEFINITIONS})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
include_directories(
  SYSTEM
  ${BULLET_INCLUDE_DIRS}  
  ${Boost_INCLUDE_DIRS}
  ${Eigen_INCLUDE_DIRS}
  ${OpenRAVE_INCLUDE_DIRS}
  ${JSON_INCLUDE_DIRS}
)

set(BSP_SOURCE_FILES
  bsp_optimizer_task.cpp
  bsp_qt_plotter.cpp
  costs.cpp
  optimizers.cpp
  utils.cpp
)

SET(QT_HEADERS bsp_optimizer_task.hpp bsp_qt_plotter.hpp)
QT4_WRAP_CPP(QT_HEADERS_MOC ${QT_HEADERS})

add_library(bsp ${BSP_SOURCE_FILES} ${QT_HEADERS_MOC})
target_link_libraries(bsp trajopt ${QT_LIBRARIES})
