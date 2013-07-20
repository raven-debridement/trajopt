#pragma once

#include "common.hpp"
#include "trajopt/bullet_collision_checker.hpp"
#include "sco/expr_vec_ops.hpp"

namespace {

  btVector3 toBt(const OR::Vector& v){
    return btVector3(v[0], v[1], v[2]);
  }

  OR::Vector toOR(const btVector3& v) {
    return OR::Vector(v.x(), v.y(), v.z());
  }

  btQuaternion toBtQuat(const OR::Vector& q) {
    return btQuaternion(q[1], q[2], q[3], q[0]);
  }

  btTransform toBt(const OR::Transform& t){
    return btTransform(toBtQuat(t.rot), toBt(t.trans));
  }

  Matrix4d toMatrix(const btTransform& b) {
    btMatrix3x3 rot = b.getBasis();
    btVector3 origin = b.getOrigin();
    Matrix4d ret;
    ret << rot[0][0], rot[0][1], rot[0][2], origin[0],
           rot[1][0], rot[1][1], rot[1][2], origin[1],
           rot[2][0], rot[2][1], rot[2][2], origin[2],
                   0,         0,         0,         1;
    return ret;
  }

  Vector3d toVector(const btVector3& v) {
    return Vector3d(v.x(), v.y(), v.z());
  }
}

namespace BSPCollision {

  inline const KinBody::Link* getLink(const btCollisionObject* o) {
    return static_cast<const CollisionObjectWrapper*>(o)->m_link;
  }

  class MultiCastHullShape : public btConvexShape {
  public:
    btConvexShape* m_shape;
    vector<btTransform> m_t0i; // T_0_i = T_w_0^-1 * T_w_i
    MultiCastHullShape(btConvexShape* shape, const vector<btTransform>& t0i) : m_shape(shape), m_t0i(t0i) {
      m_shapeType = CUSTOM_CONVEX_SHAPE_TYPE;
    }
    btVector3   localGetSupportingVertex(const btVector3& vec)const {
      vector<btVector3> svi (m_t0i.size());
      double max_vec_dot_sv = -INFINITY;
      int max_ind = -1;
      m_shape->localGetSupportingVertex(vec);
      for (int i=0; i<m_t0i.size(); i++) {
        svi[i] = m_t0i[i]*m_shape->localGetSupportingVertex(vec*m_t0i[i].getBasis());
        double vec_dot_sv = vec.dot(svi[i]);
        if (vec_dot_sv > max_vec_dot_sv) {
          max_vec_dot_sv = vec_dot_sv;
          max_ind = i;
        }
      }
      assert(max_vec_dot_sv != -INFINITY);
      assert(max_ind != -1);
      return svi[max_ind];
    }
    //notice that the vectors should be unit length
    void    batchedUnitVectorGetSupportingVertexWithoutMargin(const btVector3* vectors,btVector3* supportVerticesOut,int numVectors) const {
      throw std::runtime_error("not implemented");
    }
    ///getAabb's default implementation is brute force, expected derived classes to implement a fast dedicated version
    void getAabb(const btTransform& t_w0,btVector3& aabbMin,btVector3& aabbMax) const {
      m_shape->getAabb(t_w0, aabbMin, aabbMax);
      btVector3 min_i, max_i;
      for (int i=0; i<m_t0i.size(); i++) {
        m_shape->getAabb(t_w0*m_t0i[i], min_i, max_i );
        aabbMin.setMin(min_i);
        aabbMax.setMax(max_i);
      }
    }
    virtual void getAabbSlow(const btTransform& t,btVector3& aabbMin,btVector3& aabbMax) const {
      throw std::runtime_error("shouldn't happen");
    }
    virtual void    setLocalScaling(const btVector3& scaling) {}
    virtual const btVector3& getLocalScaling() const {
      static btVector3 out(1,1,1);
      return out;
    }
    virtual void    setMargin(btScalar margin) {}
    virtual btScalar    getMargin() const {return 0;}
    virtual int     getNumPreferredPenetrationDirections() const {return 0;}
    virtual void    getPreferredPenetrationDirection(int index, btVector3& penetrationVector) const {throw std::runtime_error("not implemented");}
    virtual void calculateLocalInertia(btScalar, btVector3&) const {throw std::runtime_error("not implemented");}
    virtual const char* getName() const {return "MultiCastHullShape";}
    virtual btVector3 localGetSupportingVertexWithoutMargin(const btVector3& v) const {return localGetSupportingVertex(v);}
  };

  struct BeliefCollision : public Collision {
    struct MultiCastInfo {
      vector<float> alpha;
      vector<int> instance_ind;
    } mi;
    BeliefCollision(const KinBody::Link* linkA, const KinBody::Link* linkB, const OR::Vector& ptA, const OR::Vector& ptB, const OR::Vector& normalB2A, double distance, float weight=1, float time=0) :
      Collision(linkA, linkB, ptA, ptB, normalB2A, distance, weight, time) {}

  };

  struct MultiCastCollisionCollector : public btCollisionWorld::ContactResultCallback {//CollisionCollector {
    std::vector<BeliefCollision>& m_collisions;
    const CollisionObjectWrapper* m_cow;
    BulletCollisionChecker* m_cc;
    MultiCastCollisionCollector(vector<BeliefCollision>& collisions, CollisionObjectWrapper* cow, BulletCollisionChecker* cc) :
        m_collisions(collisions), m_cow(cow), m_cc(cc) {}
    bool needsCollision(btBroadphaseProxy* proxy0) const {
      return (proxy0->m_collisionFilterGroup & m_collisionFilterMask)
          && (m_collisionFilterGroup & proxy0->m_collisionFilterMask)
          && m_cc->CanCollide(m_cow, static_cast<CollisionObjectWrapper*>(proxy0->m_clientObject));
    }
    btVector3 barycentricCoordinates(const btVector3& a, const btVector3& b, const btVector3& c, const btVector3& p) {
      btVector3 n = (b-a).cross(c-a);
      btVector3 na = (c-b).cross(p-b);
      btVector3 nb = (a-c).cross(p-c);
      btVector3 nc = (b-a).cross(p-a);
      float n_length2_inv = 1.0/n.length2();
      return btVector3(n.dot(na)*n_length2_inv, n.dot(nb)*n_length2_inv, n.dot(nc)*n_length2_inv);
    }
    void computeSupportingWeights(const vector<btVector3>& v, const btVector3& p, vector<float>& alpha) {
      alpha.resize(v.size());
      switch ( v.size() )
      {
      case 1:
      {
        alpha[0] = 1;
        break;
      }
      case 2:
      {
        float l0p = (p-v[0]).length();
        float l01 = (v[1]-v[0]).length();
        alpha[0] = l01 > 0  ?  fmin(l0p/l01, 1) : .5;
        alpha[1] = 1 - alpha[0];
        break;
      }
      case 3:
      {
        btVector3 bary = barycentricCoordinates(v[0], v[1], v[2], p);
        alpha[0] = bary[0];
        alpha[1] = bary[1];
        alpha[2] = bary[2];
        break;
      }
      default:
        throw runtime_error("Unsupported case for computeSupportingWeights");
      }
    }
    virtual btScalar addSingleResult(btManifoldPoint& cp,
        const btCollisionObjectWrapper* colObj0Wrap,int partId0,int index0,
        const btCollisionObjectWrapper* colObj1Wrap,int partId1,int index1) {
		  if (cp.m_distance1 > m_cc->GetContactDistance()) return 0;
      const KinBody::Link* linkA = getLink(colObj0Wrap->getCollisionObject());
		  const KinBody::Link* linkB = getLink(colObj1Wrap->getCollisionObject());
      m_collisions.push_back(BeliefCollision(linkA, linkB, toOR(cp.m_positionWorldOnA), toOR(cp.m_positionWorldOnB),
				  toOR(cp.m_normalWorldOnB), cp.m_distance1));
		  LOG_DEBUG("collide %s-%s", linkA->GetName().c_str(), linkB->GetName().c_str());
      bool castShapeIsFirst =  (colObj0Wrap->getCollisionObject() == m_cow);
      btVector3 normalWorldFromCast = -(castShapeIsFirst ? 1 : -1) * cp.m_normalWorldOnB;
      const MultiCastHullShape* shape = dynamic_cast<const MultiCastHullShape*>(m_cow->getCollisionShape());
      assert(!!shape);
      vector<float> sup(shape->m_t0i.size());
      vector<btVector3> ptWorld(shape->m_t0i.size());
      for (int i=0; i<sup.size(); i++) {
        btTransform tfWorld = m_cow->getWorldTransform() * shape->m_t0i[i];
        btVector3 normalLocal = normalWorldFromCast * tfWorld.getBasis();
        ptWorld[i] = tfWorld * shape->localGetSupportingVertex(normalLocal);
        sup[i] = normalWorldFromCast.dot(ptWorld[i]);
      }

      const float SUPPORT_FUNC_TOLERANCE = 1e-5;
      const float COLINEARITY_TOLERANCE = 1e-5;
      float max_sup = *max_element(sup.begin(), sup.end());
      vector<float> sups;
      vector<btVector3> max_ptWorlds;
      vector<int> instance_inds;
      for (int i=0; i<sup.size(); i++) {
        if (max_sup-sup[i] < SUPPORT_FUNC_TOLERANCE) {
          int j;
          for (j=0; j<max_ptWorlds.size(); j++)
            if ((max_ptWorlds[j] - ptWorld[i]).length2() < COLINEARITY_TOLERANCE) break;
          if (j==max_ptWorlds.size()) { // if this ptWorld[i] is not already in the max_ptWorlds
            sups.push_back(sup[i]);
            max_ptWorlds.push_back(ptWorld[i]);
            instance_inds.push_back(i);
          }
        }
      }

      assert(max_ptWorlds.size()>0);
      assert(max_ptWorlds.size()<4);

      const btVector3& ptOnCast = castShapeIsFirst ? cp.m_positionWorldOnA : cp.m_positionWorldOnB;
      BeliefCollision& collision = m_collisions.back();
      computeSupportingWeights(max_ptWorlds, ptOnCast, collision.mi.alpha);


      collision.mi.instance_ind = instance_inds;
      return 1;
    }
  };

  class BeliefCollisionChecker : public BulletCollisionChecker {
  public:
    typedef boost::shared_ptr<BeliefCollisionChecker> BeliefCollisionCheckerPtr;

	  BeliefCollisionChecker(OR::EnvironmentBaseConstPtr env) : BulletCollisionChecker(env) { }

    static BeliefCollisionCheckerPtr GetOrCreate(OR::EnvironmentBase& env) {
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

    void RenderCollisionShape(btCollisionShape* shape, const btTransform& tf,
        OpenRAVE::EnvironmentBase& env, vector<OpenRAVE::GraphHandlePtr>& handles, OR::RaveVector<float> color = OR::RaveVector<float>(1,1,1,.1)) const {
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

    void PlotCastHull(btCollisionShape* shape, const vector<btTransform>& tfi,
		    CollisionObjectWrapper* cow, vector<OpenRAVE::GraphHandlePtr>& handles, OR::RaveVector<float> color) {
	    if (btConvexShape* convex = dynamic_cast<btConvexShape*>(shape)) {
        vector<btTransform> t0i(tfi.size());
        // transform all the points with respect to the first transform
        for (int i=0; i<tfi.size(); i++) t0i[i] = tfi[0].inverseTimes(tfi[i]);
        MultiCastHullShape* shape = new MultiCastHullShape(convex, t0i);
        RenderCollisionShape(shape, tfi[0], *boost::const_pointer_cast<OpenRAVE::EnvironmentBase>(m_env), handles, color);
        SetTransparency(handles.back(), 0.2);
        delete shape;
      } else if (btCompoundShape* compound = dynamic_cast<btCompoundShape*>(shape)) {
        for (int child_ind = 0; child_ind < compound->getNumChildShapes(); ++child_ind) {
          vector<btTransform> tfi_child(tfi.size());
          for (int i=0; i<tfi.size(); i++) tfi_child[i] = tfi[i]*compound->getChildTransform(child_ind);
          PlotCastHull(compound->getChildShape(child_ind), tfi_child, cow, handles, color);
        }
      } else {
		    throw std::runtime_error("I can only plot convex shapes and compound shapes made of convex shapes");
      }
    }
    void PlotCastHull(Configuration& rad, const vector<KinBody::LinkPtr>& links,
        const vector<DblVec>& multi_joints, vector<OpenRAVE::GraphHandlePtr>& handles) {
      rad.Save();
      int nlinks = links.size();
      vector<vector<btTransform> > multi_tf(nlinks, vector<btTransform>(multi_joints.size())); // multi_tf[i_link][i_multi]
      for (int i_multi=0; i_multi<multi_joints.size(); i_multi++) {
        rad.SetDOFValues(multi_joints[i_multi]);
        for (int i_link=0; i_link < nlinks; ++i_link) {
          multi_tf[i_link][i_multi] = toBt(links[i_link]->GetTransform());
        }
      }

      for (int i_link=0; i_link < nlinks; ++i_link) {
        assert(m_link2cow[links[i_link].get()] != NULL);
        CollisionObjectWrapper* cow = m_link2cow[links[i_link].get()];
        vector<btTransform>& tfi = multi_tf[i_link];
        float color_param = ((float)i_link)/((float)(nlinks-1));
        PlotCastHull(cow->getCollisionShape(), multi_tf[i_link], cow, handles, OR::RaveVector<float>(color_param, 1.0-color_param, 0));
      }
    }
    virtual void MultiCastVsAll(Configuration& rad, const vector<KinBody::LinkPtr>& links,
        const vector<DblVec>& multi_joints, vector<BeliefCollision>& collisions) {
      rad.Save();
      int nlinks = links.size();
      vector<vector<btTransform> > multi_tf(nlinks, vector<btTransform>(multi_joints.size())); // multi_tf[i_link][i_multi]
      for (int i_multi=0; i_multi<multi_joints.size(); i_multi++) {
        rad.SetDOFValues(multi_joints[i_multi]);
        for (int i_link=0; i_link < nlinks; ++i_link) {
          multi_tf[i_link][i_multi] = toBt(links[i_link]->GetTransform());
        }
      }
      rad.SetDOFValues(multi_joints[0]); // is this necessary?
      UpdateBulletFromRave();
      m_world->updateAabbs();

      for (int i_link=0; i_link < nlinks; ++i_link) {
        assert(m_link2cow[links[i_link].get()] != NULL);
        CollisionObjectWrapper* cow = m_link2cow[links[i_link].get()];
        CheckShapeMultiCast(cow->getCollisionShape(), multi_tf[i_link], cow, m_world, collisions);
      }

      LOG_DEBUG("MultiCastVsAll checked %i links and found %i collisions\n", (int)links.size(), (int)collisions.size());
    }

    void CheckShapeMultiCast(btCollisionShape* shape, const vector<btTransform>& tfi,
			  CollisionObjectWrapper* cow, btCollisionWorld* world, vector<BeliefCollision>& collisions) {
      if (btConvexShape* convex = dynamic_cast<btConvexShape*>(shape)) {
        vector<btTransform> t0i(tfi.size());
        // transform all the points with respect to the first transform
        for (int i=0; i<tfi.size(); i++) t0i[i] = tfi[0].inverseTimes(tfi[i]);
        MultiCastHullShape* shape = new MultiCastHullShape(convex, t0i);
        CollisionObjectWrapper* obj = new CollisionObjectWrapper(cow->m_link);
        obj->setCollisionShape(shape);
        obj->setWorldTransform(tfi[0]);
        obj->m_index = cow->m_index;
        MultiCastCollisionCollector cc(collisions, obj, this);
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
          CheckShapeMultiCast(compound->getChildShape(child_ind), tfi_child, cow, world, collisions);
        }
      }
      else {
        throw std::runtime_error("I can only continuous collision check convex shapes and compound shapes made of convex shapes");
      }
    }

  };

  typedef boost::shared_ptr<BeliefCollisionChecker> BeliefCollisionCheckerPtr;

  template< class BeliefFuncT >
  class BeliefSingleTimestepCollisionEvaluator : public CollisionEvaluator {
  public:
    typedef typename BeliefFuncT::Ptr BeliefFuncPtr;
    BeliefSingleTimestepCollisionEvaluator(ConfigurationPtr rad, const VarVector& vars, BeliefFuncPtr belief_func, OR::KinBody::LinkPtr endeffector) :
        m_env(rad->GetEnv()),
        m_cc(BeliefCollisionChecker::GetOrCreate(*m_env)),
        m_rad(rad),
        m_vars(vars),
        m_link2ind(),
        m_links(),
        m_filterMask(-1),
        belief_func(belief_func),
        endeffector(endeffector) {
      vector<KinBody::LinkPtr> links;
      vector<int> inds;
      rad->GetAffectedLinks(m_links, true, inds);
      for (int i=0; i < m_links.size(); ++i) {
        m_link2ind[m_links[i].get()] = inds[i];
      }
    }

    virtual void CalcCollisions(const DblVec& x, vector<Collision>& collisions) {
      throw std::runtime_error("not implemented");
    }

    virtual VarVector GetVars() {
      return m_vars;
    }

    void CalcCollisions(const DblVec& x, vector<BeliefCollision>& collisions) {
      typename BeliefFuncT::BeliefT belief = getVec(x, this->m_vars);
      vector<DblVec> dofvals;
      belief_func->sigma_points(belief, &dofvals);
      m_cc->MultiCastVsAll(*m_rad, m_links, dofvals, collisions);
    }

    virtual void GetCollisionsCached(const DblVec& x, vector<BeliefCollision>& collisions) {
      double key = vecSum(x);
      vector<BeliefCollision>* it = m_cache.get(key);
      if (it != NULL) {
        RAVELOG_DEBUG("using cached collision check\n");
        collisions = *it;
      }
      else {
        RAVELOG_DEBUG("not using cached collision check\n");
        CalcCollisions(x, collisions);
        m_cache.put(key, collisions);
      }
    }

    virtual void CollisionsToDistanceExpressions(const vector<BeliefCollision>& collisions, Configuration& rad,
        const Link2Int& link2ind, const VarVector& theta_vars, const DblVec& theta_vals, vector<AffExpr>& exprs, bool isTimestep1) {
      exprs.clear();
      exprs.reserve(collisions.size());
      typename BeliefFuncT::StateT state;
      typename BeliefFuncT::BeliefT belief = toVectorXd(theta_vals);
      belief_func->extract_state(belief, &state);
      rad.SetDOFValues(toDblVec(state));
      BOOST_FOREACH(const BeliefCollision& col, collisions) {
        Link2Int::const_iterator itA = link2ind.find(col.linkA);
        Link2Int::const_iterator itB = link2ind.find(col.linkB);
        AffExpr dist;
        bool linkAFound = itA != link2ind.end();
        bool linkBFound = itB != link2ind.end();
        for (int i=0; i<col.mi.alpha.size(); i++) {
          typename BeliefFuncT::SigmaPointsGradT grad;
          belief_func->sigma_points_grad(belief, col.mi.instance_ind[i], &grad);

          AffExpr dist_a(col.distance);
          if (linkAFound) {
	          MatrixXd pos_jac = rad.PositionJacobian(endeffector->GetIndex(), col.ptA);
            VectorXd dist_grad = toVector3d(col.normalB2A).transpose()*pos_jac*grad;
            exprInc(dist_a, varDot(dist_grad, theta_vars));
            exprInc(dist_a, -dist_grad.dot(toVectorXd(theta_vals)));
          }
          if (linkBFound) {
	          MatrixXd pos_jac = rad.PositionJacobian(endeffector->GetIndex(), col.ptB);
            VectorXd dist_grad = -toVector3d(col.normalB2A).transpose()*pos_jac*grad;
            exprInc(dist_a, varDot(dist_grad, theta_vars));
            exprInc(dist_a, -dist_grad.dot(toVectorXd(theta_vals)));
          }
          if (linkAFound || linkBFound) {
            exprScale(dist_a, col.mi.alpha[i]);
            exprInc(dist, dist_a);
          }
        }
        if (dist.constant!=0 || dist.coeffs.size()!=0 || dist.vars.size()!=0) {
          exprs.push_back(dist);
        }
      }
      RAVELOG_DEBUG("%i distance expressions\n", exprs.size());
    }

    virtual void CalcDistExpressions(const DblVec& x, vector<AffExpr>& exprs) {
      vector<BeliefCollision> collisions;
      GetCollisionsCached(x, collisions);
      DblVec theta = getDblVec(x, m_vars);
      CollisionsToDistanceExpressions(collisions, *m_rad, m_link2ind, m_vars, theta, exprs, false);
    }
    virtual void CalcDists(const DblVec& x, DblVec& dists) {
      vector<BeliefCollision> collisions;
      GetCollisionsCached(x, collisions);
      vector<Collision> raw_collisions;
      for (int i = 0; i < collisions.size(); ++i) {
        raw_collisions.push_back(collisions[i]);
      }
      CollisionsToDistances(raw_collisions, m_link2ind, dists);
    }
    virtual void CustomPlot(const DblVec& x, std::vector<OR::GraphHandlePtr>& handles) {
      typename BeliefFuncT::BeliefT belief = getVec(x, this->m_vars);
      vector<DblVec> dofvals;
      belief_func->sigma_points(belief, &dofvals);
      m_cc->PlotCastHull(*m_rad, m_links, dofvals, handles);
    }
  protected:
    OR::EnvironmentBasePtr m_env;
    BeliefCollisionCheckerPtr m_cc;
    ConfigurationPtr m_rad;
    VarVector m_vars;
    Link2Int m_link2ind;
    vector<OR::KinBody::LinkPtr> m_links;
    short m_filterMask;
    BeliefFuncPtr belief_func;
    OR::KinBody::LinkPtr endeffector;
    Cache<double, vector<BeliefCollision>, 3> m_cache;
  };
}
