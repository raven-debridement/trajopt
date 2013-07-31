#include "bsp/collision/sigma_hull_cast_shape.hpp"

namespace BSPCollision {

  SigmaHullCastShape::SigmaHullCastShape(btConvexShape* shape, const vector<btTransform>& t0i, const vector<btTransform>& t1i) : SigmaHullShape(shape, concat(t0i, t1i)), m_t0i(t0i), m_t1i(t1i) {}

  const char* SigmaHullCastShape::getName() const {return "SigmaHullCastShape";}

}
