#pragma once

#include "bsp/common.hpp"

namespace BSPCollision {
  struct BeliefCollision : public Collision {
    struct SigmaHullInfo {
      vector<float> alpha;
      vector<int> instance_ind;
    } mi[2];//, mi1;
    BeliefCollision(const KinBody::Link* linkA, const KinBody::Link* linkB, const OR::Vector& ptA, const OR::Vector& ptB, const OR::Vector& normalB2A, double distance, float weight=1, float time=0);
  };
}
