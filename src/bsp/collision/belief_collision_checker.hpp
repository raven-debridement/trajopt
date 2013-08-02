#pragma once

#include "bsp/common.hpp"
#include "bsp/collision/belief_collision.hpp"

namespace BSPCollision {

  class BeliefCollisionChecker;

  typedef boost::shared_ptr<BeliefCollisionChecker> BeliefCollisionCheckerPtr;

  class BeliefCollisionChecker : public BulletCollisionChecker {
  public:

	  BeliefCollisionChecker(OR::EnvironmentBaseConstPtr env);

    static BeliefCollisionCheckerPtr GetOrCreate(OR::EnvironmentBase& env);

    void RenderCollisionShape(btCollisionShape* shape, const btTransform& tf,
        OpenRAVE::EnvironmentBase& env, vector<OpenRAVE::GraphHandlePtr>& handles, OR::RaveVector<float> color = OR::RaveVector<float>(1,1,1,.1)) const;

    void PlotSigmaHull(btCollisionShape* shape, const vector<btTransform>& tfi,
		    CollisionObjectWrapper* cow, vector<OpenRAVE::GraphHandlePtr>& handles, OR::RaveVector<float> color);

    void compute_multi_tf(Configuration& rad, const vector<KinBody::LinkPtr>& links, const vector<DblVec>& multi_joints, vector<vector<btTransform> >* output_multi_tf);

    void PlotSigmaHull(Configuration& rad, const vector<KinBody::LinkPtr>& links,
        const vector<DblVec>& multi_joints, vector<OpenRAVE::GraphHandlePtr>& handles);

    void PlotSigmaHullCast(btCollisionShape* shape, const vector<btTransform>& tf0i, const vector<btTransform>& tf1i,
		    CollisionObjectWrapper* cow, vector<OpenRAVE::GraphHandlePtr>& handles, OR::RaveVector<float> color);

    void PlotSigmaHullCast(Configuration& rad, const vector<KinBody::LinkPtr>& links,
        const vector<DblVec>& multi_joints0, const vector<DblVec>& multi_joints1, vector<OpenRAVE::GraphHandlePtr>& handles);
    
    virtual void SigmaHullVsAll(Configuration& rad, const vector<KinBody::LinkPtr>& links,
        const vector<DblVec>& multi_joints, vector<BeliefCollision>& collisions);

    void CheckShapeSigmaHull(btCollisionShape* shape, const vector<btTransform>& tfi,
			  CollisionObjectWrapper* cow, btCollisionWorld* world, vector<BeliefCollision>& collisions);

    virtual void SigmaHullCastVsAll(Configuration& rad, const vector<KinBody::LinkPtr>& links,
        const vector<DblVec>& multi_joints0, const vector<DblVec>& multi_joints1, vector<BeliefCollision>& collisions);

    void CheckShapeSigmaHullCast(btCollisionShape* shape, const vector<btTransform>& tf0i, const vector<btTransform>& tf1i,
			  CollisionObjectWrapper* cow, btCollisionWorld* world, vector<BeliefCollision>& collisions);

  };
}
