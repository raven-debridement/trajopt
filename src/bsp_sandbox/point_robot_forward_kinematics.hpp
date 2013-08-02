#pragma once

#include "point_robot.hpp"

namespace PointRobotBSP {
  Vector3d forward_kinematics(const StateT& state) {
    Vector3d eetrans;
    eetrans[0] = state[0];
    eetrans[1] = state[1];
    eetrans[2] = 0;
    return eetrans;
  }
}
