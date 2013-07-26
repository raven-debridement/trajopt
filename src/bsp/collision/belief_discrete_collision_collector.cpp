#include "bsp/collision/utils.hpp"
#include "bsp/collision/belief_discrete_collision_collector.hpp"
#include "bsp/collision/sigma_hull_shape.hpp"

namespace BSPCollision {

  inline const KinBody::Link* getLink(const btCollisionObject* o) {
    return static_cast<const CollisionObjectWrapper*>(o)->m_link;
  }

  BeliefDiscreteCollisionCollector::BeliefDiscreteCollisionCollector(vector<BeliefCollision>& collisions, CollisionObjectWrapper* cow, BulletCollisionChecker* cc) :
      m_collisions(collisions), m_cow(cow), m_cc(cc) {}

  bool BeliefDiscreteCollisionCollector::needsCollision(btBroadphaseProxy* proxy0) const {
    return (proxy0->m_collisionFilterGroup & m_collisionFilterMask)
        && (m_collisionFilterGroup & proxy0->m_collisionFilterMask)
        && m_cc->CanCollide(m_cow, static_cast<CollisionObjectWrapper*>(proxy0->m_clientObject));
  }
  
  btScalar BeliefDiscreteCollisionCollector::addSingleResult(btManifoldPoint& cp,
      const btCollisionObjectWrapper* colObj0Wrap,int partId0,int index0,
      const btCollisionObjectWrapper* colObj1Wrap,int partId1,int index1) {
    if (cp.m_distance1 > m_cc->GetContactDistance()) return 0;
    const KinBody::Link* linkA = getLink(colObj0Wrap->getCollisionObject());
    const KinBody::Link* linkB = getLink(colObj1Wrap->getCollisionObject());
    m_collisions.push_back(BeliefCollision(linkA, linkB, toOR(cp.m_positionWorldOnA), toOR(cp.m_positionWorldOnB),
        toOR(cp.m_normalWorldOnB), cp.m_distance1));
    LOG_DEBUG("collide %s-%s", linkA->GetName().c_str(), linkB->GetName().c_str());
    bool shapeIsFirst =  (colObj0Wrap->getCollisionObject() == m_cow);
    btVector3 normalWorldFromShape = -(shapeIsFirst ? 1 : -1) * cp.m_normalWorldOnB;
    const SigmaHullShape* shape = dynamic_cast<const SigmaHullShape*>(m_cow->getCollisionShape());
    assert(!!shape);

    vector<float> sup;
    vector<btVector3> ptWorld;
    compute_points_and_supports(shape->m_shape, shape->m_ti, normalWorldFromShape, m_cow, &sup, &ptWorld);

    vector<float> sups;
    vector<btVector3> max_ptWorlds;
    vector<int> instance_inds;
    compute_max_support_points(sup, ptWorld, &sups, &max_ptWorlds, &instance_inds);

    assert(max_ptWorlds.size()>0);
    assert(max_ptWorlds.size()<4);

    const btVector3& ptOnShape = shapeIsFirst ? cp.m_positionWorldOnA : cp.m_positionWorldOnB;
    BeliefCollision& collision = m_collisions.back();
    computeSupportingWeights(max_ptWorlds, ptOnShape, collision.mi[0].alpha);

    collision.mi[0].instance_ind = instance_inds;
    return 1;
  }
}
