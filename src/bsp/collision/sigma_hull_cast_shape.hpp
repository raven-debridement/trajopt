#pragma once

#include "bsp/common.hpp"

namespace BSPCollision {
  class SigmaHullCastShape : public SigmaHullShape {
  public:
    //btConvexShape* m_shape;
    vector<btTransform> m_t0i; // T_0_i = T_w_0_0^-1 * T_w_0_i
    vector<btTransform> m_t1i; // T_1_i = T_w_0_0^-1 * T_w_1_i
    SigmaHullCastShape(btConvexShape* shape, const vector<btTransform>& t0i, const vector<btTransform>& t1i);
    virtual const char* getName() const;
  };
}
