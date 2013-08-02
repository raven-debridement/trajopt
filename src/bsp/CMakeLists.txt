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
  costs.cpp
  optimizers.cpp
  utils.cpp
  collision/belief_collision.cpp
  collision/belief_collision_checker.cpp
  collision/belief_continuous_collision_collector.cpp
  collision/belief_discrete_collision_collector.cpp
  collision/sigma_hull_shape.cpp
  collision/sigma_hull_cast_shape.cpp
  collision/utils.cpp
  ../trajopt/bullet_unity.cpp
)

add_library(bsp ${BSP_SOURCE_FILES})
target_link_libraries(bsp trajopt)
