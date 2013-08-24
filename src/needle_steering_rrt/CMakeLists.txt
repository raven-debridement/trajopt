add_definitions(-DDATA_DIR="${CMAKE_SOURCE_DIR}/data")

include_directories(
  SYSTEM
  ${BULLET_INCLUDE_DIRS}  
  ${Boost_INCLUDE_DIRS}
  ${Eigen_INCLUDE_DIRS}
  ${OpenRAVE_INCLUDE_DIRS}
  ${JSON_INCLUDE_DIRS}
  #${RBDL_INCLUDE_DIRS}
)

add_executable(needle_steering_rrt main.cpp planner.cpp interface.cpp planner2D.cpp planner3D.cpp ../trajopt/bullet_unity.cpp)
#target_link_libraries(needle_steering_rrt trajopt)
