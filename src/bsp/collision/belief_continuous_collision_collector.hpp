#pragma once

#include "bsp/common.hpp"
#include "bsp/collision/belief_collision.hpp"

namespace BSPCollision {

 struct BeliefContinuousCollisionCollector : public btCollisionWorld::ContactResultCallback {
    std::vector<BeliefCollision>& m_collisions;
    const CollisionObjectWrapper* m_cow;
    BulletCollisionChecker* m_cc;
    BeliefContinuousCollisionCollector(vector<BeliefCollision>& collisions, CollisionObjectWrapper* cow, BulletCollisionChecker* cc);
    bool needsCollision(btBroadphaseProxy* proxy0) const;
    virtual btScalar addSingleResult(btManifoldPoint& cp,
        const btCollisionObjectWrapper* colObj0Wrap,int partId0,int index0,
        const btCollisionObjectWrapper* colObj1Wrap,int partId1,int index1);
  };
}
