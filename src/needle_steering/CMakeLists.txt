add_definitions(-DDATA_DIR="${CMAKE_SOURCE_DIR}/data")

#set(RBDL_INCLUDE_DIRS "${CMAKE_SOURCE_DIR}/ext/rbdl/src" "${CMAKE_BINARY_DIR}/ext/rbdl/src")

include_directories(
  SYSTEM
  ${BULLET_INCLUDE_DIRS}  
  ${Boost_INCLUDE_DIRS}
  ${Eigen_INCLUDE_DIRS}
  ${OpenRAVE_INCLUDE_DIRS}
  ${JSON_INCLUDE_DIRS}
  #${RBDL_INCLUDE_DIRS}
)

add_executable(needle_steering needle_steering.cpp constraints.cpp costs.cpp local_configuration.cpp needle_problem_helper.cpp traj_plotter.cpp utils.cpp ../trajopt/bullet_unity.cpp needle_problem_instance.cpp needle_problem_planner.cpp)
target_link_libraries(needle_steering trajopt)
