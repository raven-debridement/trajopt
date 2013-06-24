set(RBDL_INCLUDE_DIRS "${CMAKE_SOURCE_DIR}/ext/rbdl/src" "${CMAKE_BINARY_DIR}/ext/rbdl/src")

include_directories(
  SYSTEM
  ${BULLET_INCLUDE_DIRS}  
  ${Boost_INCLUDE_DIRS}
  ${Eigen_INCLUDE_DIRS}
  ${OpenRAVE_INCLUDE_DIRS}
  ${JSON_INCLUDE_DIRS}
  ${RBDL_INCLUDE_DIRS}
)

add_executable(toy toy_main.cpp toy.cpp costs.cpp optimizers.cpp utils.cpp gnuplot_i.cc)
target_link_libraries(toy trajopt)
