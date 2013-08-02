#pragma once

#include "bsp/common.hpp"

namespace BSPCollision {

  class BeliefCollisionEvaluator : public CollisionEvaluator {
  public:
    virtual void GetCollisionsCached(const DblVec& x, vector<BeliefCollision>& collisions) = 0;
    virtual void CalcCollisions(const DblVec& x, vector<BeliefCollision>& collisions) = 0;
    virtual void CustomPlot(const DblVec& x, std::vector<OR::GraphHandlePtr>& handles) = 0;
  };
}
