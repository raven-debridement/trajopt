#include "bsp/collision/utils.hpp"
#include "bsp/collision/belief_collision_checker.hpp"
#include "bsp/collision/sigma_hull_shape.hpp"
#include "bsp/collision/sigma_hull_cast_shape.hpp"
#include "bsp/collision/belief_discrete_collision_collector.hpp"
#include "bsp/collision/belief_continuous_collision_collector.hpp"

namespace BSPCollision {

  BeliefCollisionChecker::BeliefCollisionChecker(OR::EnvironmentBaseConstPtr env) : BulletCollisionChecker(env) { }

  BeliefCollisionCheckerPtr BeliefCollisionChecker::GetOrCreate(OR::EnvironmentBase& env) {
    UserDataPtr ud = GetUserData(env, "trajopt_cc");
    if (!ud) {
      LOG_INFO("creating bullet belief collision checker for environment");
      BeliefCollisionCheckerPtr checker(new BeliefCollisionChecker(env.shared_from_this()));
      ud = checker;
      SetUserData(env, "trajopt_cc", ud);
    }
    else {
      LOG_DEBUG("already have a belief collision checker for this environment");
    }
    return boost::dynamic_pointer_cast<BeliefCollisionChecker>(ud);

  }

  void BeliefCollisionChecker::RenderCollisionShape(btCollisionShape* shape, const btTransform& tf,
      OpenRAVE::EnvironmentBase& env, vector<OpenRAVE::GraphHandlePtr>& handles, OR::RaveVector<float> color) const {
    switch (shape->getShapeType()) {
    case COMPOUND_SHAPE_PROXYTYPE:
    case CONVEX_HULL_SHAPE_PROXYTYPE:
      break;
    default:
      if (shape->getShapeType() <= CUSTOM_CONVEX_SHAPE_TYPE) {
        btConvexShape* convex = dynamic_cast<btConvexShape*>(shape);
        btShapeHull* hull = new btShapeHull(convex);
        hull->buildHull(convex->getMargin());
        int num_triangles = hull->numTriangles();
        const unsigned int* indices = hull->getIndexPointer();
        const btVector3* vertices = hull->getVertexPointer();
        btVector3 tf_vertices[hull->numVertices()];
        for (int i=0; i<hull->numVertices(); i++) tf_vertices[i] = tf * vertices[i];

        handles.push_back(OSGViewer::GetOrCreate(env.shared_from_this())->drawtrimesh((float*)tf_vertices, 16, (int*) indices, num_triangles, color));
        return;
      }
    }
    BulletCollisionChecker::RenderCollisionShape(shape, tf, env, handles);
  }

  void BeliefCollisionChecker::PlotSigmaHull(btCollisionShape* shape, const vector<btTransform>& tfi,
      CollisionObjectWrapper* cow, vector<OpenRAVE::GraphHandlePtr>& handles, OR::RaveVector<float> color) {
    if (btConvexShape* convex = dynamic_cast<btConvexShape*>(shape)) {
      vector<btTransform> t0i(tfi.size());
      // transform all the points with respect to the first transform
      for (int i=0; i<tfi.size(); i++) t0i[i] = tfi[0].inverseTimes(tfi[i]);
      SigmaHullShape* shape = new SigmaHullShape(convex, t0i);
      RenderCollisionShape(shape, tfi[0], *boost::const_pointer_cast<OpenRAVE::EnvironmentBase>(m_env), handles, color);
      // ADDED
      SetTransparency(handles.back(), 0.2);
      // END ADDED
      delete shape;
    } else if (btCompoundShape* compound = dynamic_cast<btCompoundShape*>(shape)) {
      for (int child_ind = 0; child_ind < compound->getNumChildShapes(); ++child_ind) {
        vector<btTransform> tfi_child(tfi.size());
        for (int i=0; i<tfi.size(); i++) tfi_child[i] = tfi[i]*compound->getChildTransform(child_ind);
        PlotSigmaHull(compound->getChildShape(child_ind), tfi_child, cow, handles, color);
      }
    } else {
      throw std::runtime_error("I can only plot convex shapes and compound shapes made of convex shapes");
    }
  }

  void BeliefCollisionChecker::compute_multi_tf(Configuration& rad, const vector<KinBody::LinkPtr>& links, const vector<DblVec>& multi_joints, vector<vector<btTransform> >* output_multi_tf) {
    assert (output_multi_tf != nullptr);
    output_multi_tf->clear();
    int nlinks = links.size();
    for (int i = 0; i < nlinks; ++i) {
      output_multi_tf->push_back(vector<btTransform>(multi_joints.size()));
    }
    for (int i_multi=0; i_multi<multi_joints.size(); i_multi++) {
      rad.SetDOFValues(multi_joints[i_multi]);
      for (int i_link=0; i_link < nlinks; ++i_link) {
        (*output_multi_tf)[i_link][i_multi] = toBt(links[i_link]->GetTransform());
      }
    }
  }

  void BeliefCollisionChecker::PlotSigmaHull(Configuration& rad, const vector<KinBody::LinkPtr>& links,
      const vector<DblVec>& multi_joints, vector<OpenRAVE::GraphHandlePtr>& handles) {
    rad.Save();
    int nlinks = links.size();
    vector<vector<btTransform> > multi_tf;
    compute_multi_tf(rad, links, multi_joints, &multi_tf);

    for (int i_link=0; i_link < nlinks; ++i_link) {
      assert(m_link2cow[links[i_link].get()] != NULL);
      CollisionObjectWrapper* cow = m_link2cow[links[i_link].get()];
      vector<btTransform>& tfi = multi_tf[i_link];
      float color_param = ((float)i_link)/((float)(nlinks-1));
      PlotSigmaHull(cow->getCollisionShape(), tfi, cow, handles, OR::RaveVector<float>(color_param, 1.0-color_param, 0));
    }
  }

  void BeliefCollisionChecker::PlotSigmaHullCast(btCollisionShape* shape, const vector<btTransform>& tf0i, const vector<btTransform>& tf1i,
      CollisionObjectWrapper* cow, vector<OpenRAVE::GraphHandlePtr>& handles, OR::RaveVector<float> color) {
    if (btConvexShape* convex = dynamic_cast<btConvexShape*>(shape)) {
      vector<btTransform> t0i(tf0i.size());
      vector<btTransform> t1i(tf1i.size());
      // transform all the points with respect to the first transform
      for (int i=0; i<tf0i.size(); i++) t0i[i] = tf0i[0].inverseTimes(tf0i[i]);
      for (int i=0; i<tf1i.size(); i++) t1i[i] = tf0i[0].inverseTimes(tf1i[i]);
      SigmaHullCastShape* shape = new SigmaHullCastShape(convex, t0i, t1i);
      RenderCollisionShape(shape, tf0i[0], *boost::const_pointer_cast<OpenRAVE::EnvironmentBase>(m_env), handles, color);
      SetTransparency(handles.back(), 0.2);
      delete shape;
    } else if (btCompoundShape* compound = dynamic_cast<btCompoundShape*>(shape)) {
      for (int child_ind = 0; child_ind < compound->getNumChildShapes(); ++child_ind) {
        vector<btTransform> tf0i_child(tf0i.size());
        vector<btTransform> tf1i_child(tf1i.size());
        for (int i=0; i<tf0i.size(); i++) tf0i_child[i] = tf0i[i]*compound->getChildTransform(child_ind);
        for (int i=0; i<tf1i.size(); i++) tf1i_child[i] = tf1i[i]*compound->getChildTransform(child_ind);
        PlotSigmaHullCast(compound->getChildShape(child_ind), tf0i_child, tf1i_child, cow, handles, color);
      }
    } else {
      throw std::runtime_error("I can only plot convex shapes and compound shapes made of convex shapes");
    }
  }

  void BeliefCollisionChecker::PlotSigmaHullCast(Configuration& rad, const vector<KinBody::LinkPtr>& links,
      const vector<DblVec>& multi_joints0, const vector<DblVec>& multi_joints1, vector<OpenRAVE::GraphHandlePtr>& handles) {
    rad.Save();
    int nlinks = links.size();
    vector<vector<btTransform> > multi_tf0, multi_tf1;
    compute_multi_tf(rad, links, multi_joints0, &multi_tf0);
    compute_multi_tf(rad, links, multi_joints1, &multi_tf1);
    
    for (int i_link=0; i_link < nlinks; ++i_link) {
      assert(m_link2cow[links[i_link].get()] != NULL);
      CollisionObjectWrapper* cow = m_link2cow[links[i_link].get()];
      vector<btTransform>& tf0i = multi_tf0[i_link];
      vector<btTransform>& tf1i = multi_tf1[i_link];
      float color_param = ((float)i_link)/((float)(nlinks-1));
      PlotSigmaHullCast(cow->getCollisionShape(), tf0i, tf1i, cow, handles, OR::RaveVector<float>(color_param, 1.0-color_param, 0));
    }
  }

  
  void BeliefCollisionChecker::SigmaHullVsAll(Configuration& rad, const vector<KinBody::LinkPtr>& links,
      const vector<DblVec>& multi_joints, vector<BeliefCollision>& collisions) {
    Configuration::SaverPtr saver = rad.Save();
    int nlinks = links.size();
    vector<vector<btTransform> > multi_tf;
    compute_multi_tf(rad, links, multi_joints, &multi_tf);

    rad.SetDOFValues(multi_joints[0]); // is this necessary? => yes it is
    UpdateBulletFromRave();
    m_world->updateAabbs();

    for (int i_link=0; i_link < nlinks; ++i_link) {
      assert(m_link2cow[links[i_link].get()] != NULL);
      CollisionObjectWrapper* cow = m_link2cow[links[i_link].get()];
      CheckShapeSigmaHull(cow->getCollisionShape(), multi_tf[i_link], cow, m_world, collisions);
    }

    LOG_DEBUG("SigmaHullVsAll checked %i links and found %i collisions\n", (int)links.size(), (int)collisions.size());
  }

  void BeliefCollisionChecker::CheckShapeSigmaHull(btCollisionShape* shape, const vector<btTransform>& tfi,
      CollisionObjectWrapper* cow, btCollisionWorld* world, vector<BeliefCollision>& collisions) {
    if (btConvexShape* convex = dynamic_cast<btConvexShape*>(shape)) {
      vector<btTransform> t0i(tfi.size());
      // transform all the points with respect to the first transform
      for (int i=0; i<tfi.size(); i++) t0i[i] = tfi[0].inverseTimes(tfi[i]);
      SigmaHullShape* shape = new SigmaHullShape(convex, t0i);
      CollisionObjectWrapper* obj = new CollisionObjectWrapper(cow->m_link);
      obj->setCollisionShape(shape);
      obj->setWorldTransform(tfi[0]);
      obj->m_index = cow->m_index;
      BeliefDiscreteCollisionCollector cc(collisions, obj, this);
      cc.m_collisionFilterMask = KinBodyFilter;
      cc.m_collisionFilterGroup = RobotFilter;
      world->contactTest(obj, cc);

      delete obj;
      delete shape;
    }
    else if (btCompoundShape* compound = dynamic_cast<btCompoundShape*>(shape)) {
      for (int child_ind = 0; child_ind < compound->getNumChildShapes(); ++child_ind) {
        vector<btTransform> tfi_child(tfi.size());
        for (int i=0; i<tfi.size(); i++) tfi_child[i] = tfi[i]*compound->getChildTransform(child_ind);
        CheckShapeSigmaHull(compound->getChildShape(child_ind), tfi_child, cow, world, collisions);
      }
    }
    else {
      throw std::runtime_error("I can only collision check convex shapes and compound shapes made of convex shapes");
    }
  }

  void BeliefCollisionChecker::SigmaHullCastVsAll(Configuration& rad, const vector<KinBody::LinkPtr>& links,
      const vector<DblVec>& multi_joints0, const vector<DblVec>& multi_joints1, vector<BeliefCollision>& collisions) {
    Configuration::SaverPtr saver = rad.Save();
    int nlinks = links.size();
    vector<vector<btTransform> > multi_tf0, multi_tf1;
    compute_multi_tf(rad, links, multi_joints0, &multi_tf0);
    compute_multi_tf(rad, links, multi_joints1, &multi_tf1);

    rad.SetDOFValues(multi_joints0[0]); // is this necessary?
    UpdateBulletFromRave();
    m_world->updateAabbs();

    for (int i_link=0; i_link < nlinks; ++i_link) {
      assert(m_link2cow[links[i_link].get()] != NULL);
      CollisionObjectWrapper* cow = m_link2cow[links[i_link].get()];
      CheckShapeSigmaHullCast(cow->getCollisionShape(), multi_tf0[i_link], multi_tf1[i_link], cow, m_world, collisions);
    }

    LOG_DEBUG("SigmaHullCastVsAll checked %i links and found %i collisions\n", (int)links.size(), (int)collisions.size());
  }

  void BeliefCollisionChecker::CheckShapeSigmaHullCast(btCollisionShape* shape, const vector<btTransform>& tf0i, const vector<btTransform>& tf1i,
      CollisionObjectWrapper* cow, btCollisionWorld* world, vector<BeliefCollision>& collisions) {
    if (btConvexShape* convex = dynamic_cast<btConvexShape*>(shape)) {
      vector<btTransform> t0i(tf0i.size());
      vector<btTransform> t1i(tf1i.size());
      // transform all the points with respect to the first transform
      for (int i=0; i<tf0i.size(); i++) t0i[i] = tf0i[0].inverseTimes(tf0i[i]);
      for (int i=0; i<tf1i.size(); i++) t1i[i] = tf0i[0].inverseTimes(tf1i[i]);
      SigmaHullCastShape* shape = new SigmaHullCastShape(convex, t0i, t1i);
      CollisionObjectWrapper* obj = new CollisionObjectWrapper(cow->m_link);
      obj->setCollisionShape(shape);
      obj->setWorldTransform(tf0i[0]);
      obj->m_index = cow->m_index;
      BeliefContinuousCollisionCollector cc(collisions, obj, this);
      cc.m_collisionFilterMask = KinBodyFilter;
      cc.m_collisionFilterGroup = RobotFilter;
      world->contactTest(obj, cc);

      delete obj;
      delete shape;
    }
    else if (btCompoundShape* compound = dynamic_cast<btCompoundShape*>(shape)) {
      for (int child_ind = 0; child_ind < compound->getNumChildShapes(); ++child_ind) {
        vector<btTransform> tf0i_child(tf0i.size());
        vector<btTransform> tf1i_child(tf1i.size());
        for (int i=0; i<tf0i.size(); i++) tf0i_child[i] = tf0i[i]*compound->getChildTransform(child_ind);
        for (int i=0; i<tf1i.size(); i++) tf1i_child[i] = tf1i[i]*compound->getChildTransform(child_ind);
        CheckShapeSigmaHullCast(compound->getChildShape(child_ind), tf0i_child, tf1i_child, cow, world, collisions);
      }
    }
    else {
      throw std::runtime_error("I can only continuous collision check convex shapes and compound shapes made of convex shapes");
    }
  }
}
