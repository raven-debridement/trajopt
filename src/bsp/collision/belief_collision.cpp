#include "bsp/collision/belief_collision.hpp"

namespace BSPCollision {
  BeliefCollision::BeliefCollision(const KinBody::Link* linkA, const KinBody::Link* linkB, const OR::Vector& ptA, const OR::Vector& ptB, const OR::Vector& normalB2A, double distance, float weight, float time) :
    Collision(linkA, linkB, ptA, ptB, normalB2A, distance, weight, time) {}
}
