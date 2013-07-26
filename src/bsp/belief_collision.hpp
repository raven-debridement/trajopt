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
    int vsize = v.size();
    // if the three vertices are on the same line, only use the first two
    if (vsize == 3 && ((v[1] - v[0]).cross(v[2] - v[0])).length2() < BSP::eps) {
      vsize = 2;
    }

    
      
    switch ( vsize )
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
      //if (isnan(bary[0]) || isnan(bary[1]) || isnan(bary[2])) {
      //  cout << "calculating bary centric coordinates" << endl;
      //  cout << "three points: " << endl;
      //  cout << toVector(v[0]).transpose() << endl;
      //  cout << toVector(v[1]).transpose() << endl;
      //  cout << toVector(v[2]).transpose() << endl;
      //  cout << "for: " << endl;
      //  cout << toVector(p).transpose() << endl;
      //  throw runtime_error("nan!");
      //}
        
      alpha[0] = bary[0];
      alpha[1] = bary[1];
      alpha[2] = bary[2];
      break;
    }
    default: {
      for (auto& pt : v) {
        cout << toVector(pt).transpose() << endl;
      }
      cout << "point: " << toVector(p).transpose() << endl;
      throw runtime_error("Unsupported case for computeSupportingWeights: v.size() = " + std::to_string(v.size()));
    }
    }
  }

  class SigmaHullShape : public btConvexShape {
  public:
    btConvexShape* m_shape;
    vector<btTransform> m_ti; // T_i = T_w_0^-1 * T_w_i
    SigmaHullShape(btConvexShape* shape, const vector<btTransform>& ti) : m_shape(shape), m_ti(ti) {
      m_shapeType = CUSTOM_CONVEX_SHAPE_TYPE;
    }
    virtual btVector3   localGetSupportingVertex(const btVector3& vec)const {
      vector<btVector3> svi (m_ti.size());
      double max_vec_dot_sv = -INFINITY;
      int max_ind = -1;
      m_shape->localGetSupportingVertex(vec);
      for (int i=0; i<m_ti.size(); i++) {
        svi[i] = m_ti[i]*m_shape->localGetSupportingVertex(vec*m_ti[i].getBasis());
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
    virtual void getAabb(const btTransform& t_w0,btVector3& aabbMin,btVector3& aabbMax) const {
      m_shape->getAabb(t_w0, aabbMin, aabbMax);
      btVector3 min_i, max_i;
      for (int i=0; i<m_ti.size(); i++) {
        m_shape->getAabb(t_w0*m_ti[i], min_i, max_i );
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
    virtual const char* getName() const {return "SigmaHullShape";}
    virtual btVector3 localGetSupportingVertexWithoutMargin(const btVector3& v) const {return localGetSupportingVertex(v);}
  };

  class SigmaHullCastShape : public SigmaHullShape {
  public:
    //btConvexShape* m_shape;
    vector<btTransform> m_t0i; // T_0_i = T_w_0_0^-1 * T_w_0_i
    vector<btTransform> m_t1i; // T_1_i = T_w_0_0^-1 * T_w_1_i
    SigmaHullCastShape(btConvexShape* shape, const vector<btTransform>& t0i, const vector<btTransform>& t1i) : SigmaHullShape(shape, concat(t0i, t1i)), m_t0i(t0i), m_t1i(t1i) {}
    virtual const char* getName() const {return "SigmaHullCastShape";}

    void get_support_timestep0(const btVector3& vec, float* output_sup, btVector3* output_pt) {
      assert (output_sup != nullptr);
      assert (output_pt != nullptr);
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
      *output_sup = max_vec_dot_sv;
      *output_pt = svi[max_ind];

    }

    void get_support_timestep1(const btVector3& vec, float* output_sup, btVector3* output_pt) {
      assert (output_sup != nullptr);
      assert (output_pt != nullptr);

      vector<btVector3> svi (m_t1i.size());
      double max_vec_dot_sv = -INFINITY;
      int max_ind = -1;
      m_shape->localGetSupportingVertex(vec);
      for (int i=0; i<m_t1i.size(); i++) {
        svi[i] = m_t1i[i]*m_shape->localGetSupportingVertex(vec*m_t1i[i].getBasis());
        double vec_dot_sv = vec.dot(svi[i]);
        if (vec_dot_sv > max_vec_dot_sv) {
          max_vec_dot_sv = vec_dot_sv;
          max_ind = i;
        }
      }
      assert(max_vec_dot_sv != -INFINITY);
      assert(max_ind != -1);
      *output_sup = max_vec_dot_sv;
      *output_pt = svi[max_ind];

    }
  };

  struct BeliefCollision : public Collision {
    struct SigmaHullInfo {
      vector<float> alpha;
      vector<int> instance_ind;
    } mi[2];//, mi1;
    BeliefCollision(const KinBody::Link* linkA, const KinBody::Link* linkB, const OR::Vector& ptA, const OR::Vector& ptB, const OR::Vector& normalB2A, double distance, float weight=1, float time=0) :
      Collision(linkA, linkB, ptA, ptB, normalB2A, distance, weight, time) {}
  };

  void compute_points_and_supports(const btConvexShape* shape, const vector<btTransform>& ts, const btVector3& normalWorldFromShape, const CollisionObjectWrapper* cow, vector<float>* output_sup, vector<btVector3>* output_ptWorld) {
    assert (output_sup != nullptr);
    assert (output_ptWorld != nullptr);

    output_sup->resize(ts.size());
    output_ptWorld->resize(ts.size());

    //cout << "computing points and supports" << endl;
    for (int i = 0; i < ts.size(); ++i) {
      btTransform tfWorld = cow->getWorldTransform() * ts[i];
      //cout << "world transform: " << endl << toMatrix(tfWorld) << endl << endl;
      btVector3 normalLocal = normalWorldFromShape * tfWorld.getBasis();
      //cout << "normal local: " << toVector(normalLocal).transpose() << endl;
      (*output_ptWorld)[i] = tfWorld * shape->localGetSupportingVertex(normalLocal);
      (*output_sup)[i] = normalWorldFromShape.dot((*output_ptWorld)[i]);
    }
  }

  void compute_max_support_points(const vector<float>& sup, const vector<btVector3>& ptWorld, vector<float>* output_sups, vector<btVector3>* output_max_ptWorlds, vector<int>* output_instance_inds) {
    assert (output_sups != nullptr);
    assert (output_max_ptWorlds != nullptr);
    assert (output_instance_inds != nullptr);

    output_sups->clear();
    output_max_ptWorlds->clear();
    output_instance_inds->clear();

    const float SUPPORT_FUNC_TOLERANCE = 1e-5;
    const float COLINEARITY_TOLERANCE = 1e-5;
    float max_sup = *max_element(sup.begin(), sup.end());
    
    for (int i = 0; i < sup.size(); ++i) {
      if (max_sup-sup[i] < SUPPORT_FUNC_TOLERANCE) {
        int j;
        for (j = 0; j < output_max_ptWorlds->size(); ++j) {
          if (((*output_max_ptWorlds)[j] - ptWorld[i]).length2() < COLINEARITY_TOLERANCE) break;
        }
        if (j == output_max_ptWorlds->size()) { // if this ptWorld[i] is not already in the max_ptWorlds
          output_sups->push_back(sup[i]);
          output_max_ptWorlds->push_back(ptWorld[i]);
          output_instance_inds->push_back(i);
        }
      }
    }
  }

  struct BeliefDiscreteCollisionCollector : public btCollisionWorld::ContactResultCallback {
    std::vector<BeliefCollision>& m_collisions;
    const CollisionObjectWrapper* m_cow;
    BulletCollisionChecker* m_cc;
    BeliefDiscreteCollisionCollector(vector<BeliefCollision>& collisions, CollisionObjectWrapper* cow, BulletCollisionChecker* cc) :
        m_collisions(collisions), m_cow(cow), m_cc(cc) {}
    bool needsCollision(btBroadphaseProxy* proxy0) const {
      return (proxy0->m_collisionFilterGroup & m_collisionFilterMask)
          && (m_collisionFilterGroup & proxy0->m_collisionFilterMask)
          && m_cc->CanCollide(m_cow, static_cast<CollisionObjectWrapper*>(proxy0->m_clientObject));
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
      bool shapeIsFirst =  (colObj0Wrap->getCollisionObject() == m_cow);
      btVector3 normalWorldFromShape = -(shapeIsFirst ? 1 : -1) * cp.m_normalWorldOnB;
      const SigmaHullShape* shape = dynamic_cast<const SigmaHullShape*>(m_cow->getCollisionShape());
      assert(!!shape);

      vector<float> sup;
      vector<btVector3> ptWorld;
      compute_points_and_supports(shape->m_shape, shape->m_ti, normalWorldFromShape, m_cow, &sup, &ptWorld);
      //compute_points_and_supports(shape, shape->m_ti, normalWorldFromShape, m_cow, &sup, &ptWorld);

      vector<float> sups;
      vector<btVector3> max_ptWorlds;
      vector<int> instance_inds;
      compute_max_support_points(sup, ptWorld, &sups, &max_ptWorlds, &instance_inds);

      assert(max_ptWorlds.size()>0);
      assert(max_ptWorlds.size()<4);

      const btVector3& ptOnShape = shapeIsFirst ? cp.m_positionWorldOnA : cp.m_positionWorldOnB;
      BeliefCollision& collision = m_collisions.back();
      computeSupportingWeights(max_ptWorlds, ptOnShape, collision.mi[0].alpha);

      ///cout << "alpha size: " << collision.mi[0].alpha.size() << endl;
      ///for (auto& i : collision.mi[0].alpha) {
      ///  cout << "alpha: " << i << endl;
      ///}

      collision.mi[0].instance_ind = instance_inds;
      return 1;
    }
  };

  struct BeliefContinuousCollisionCollector : public btCollisionWorld::ContactResultCallback {
    std::vector<BeliefCollision>& m_collisions;
    const CollisionObjectWrapper* m_cow;
    BulletCollisionChecker* m_cc;
    BeliefContinuousCollisionCollector(vector<BeliefCollision>& collisions, CollisionObjectWrapper* cow, BulletCollisionChecker* cc) :
        m_collisions(collisions), m_cow(cow), m_cc(cc) {}
    bool needsCollision(btBroadphaseProxy* proxy0) const {
      return (proxy0->m_collisionFilterGroup & m_collisionFilterMask)
          && (m_collisionFilterGroup & proxy0->m_collisionFilterMask)
          && m_cc->CanCollide(m_cow, static_cast<CollisionObjectWrapper*>(proxy0->m_clientObject));
    }

    virtual btScalar addSingleResult(btManifoldPoint& cp,
        const btCollisionObjectWrapper* colObj0Wrap,int partId0,int index0,
        const btCollisionObjectWrapper* colObj1Wrap,int partId1,int index1) {
		  if (cp.m_distance1 > m_cc->GetContactDistance()) return 0;
      const KinBody::Link* linkA = getLink(colObj0Wrap->getCollisionObject());
		  const KinBody::Link* linkB = getLink(colObj1Wrap->getCollisionObject());

      m_collisions.push_back(BeliefCollision(linkA, linkB, toOR(cp.m_positionWorldOnA), toOR(cp.m_positionWorldOnB),
				  toOR(cp.m_normalWorldOnB), cp.m_distance1));
		  LOG_INFO("collide %s-%s", linkA->GetName().c_str(), linkB->GetName().c_str());
      //cout << "m position world on A: " << toVector(cp.m_positionWorldOnA).transpose() << endl;
      //cout << "m position world on B: " << toVector(cp.m_positionWorldOnB).transpose() << endl;
      bool castShapeIsFirst =  (colObj0Wrap->getCollisionObject() == m_cow);
      //if ((colObj0Wrap->getCollisionObject() == m_cow)) {
      //  cout << "cast shape is first" << endl;
      //} else if ((colObj1Wrap->getCollisionObject() == m_cow)) {
      //  cout << "cast shape is second" << endl;
      //} else {
      //  cout << "ERROR! CAST SHAPE IS NEITHER!!" << endl;
      //}

      btVector3 normalWorldFromCast = -(castShapeIsFirst ? 1 : -1) * cp.m_normalWorldOnB;
      //cout << "normal world from cast: " << toVector(normalWorldFromCast).transpose() << endl;
      const SigmaHullCastShape* shape = dynamic_cast<const SigmaHullCastShape*>(m_cow->getCollisionShape());
      assert(!!shape);

      BeliefCollision& collision = m_collisions.back();

      if (castShapeIsFirst) {
        swap(collision.ptA, collision.ptB);
        swap(collision.linkA, collision.linkB);
        collision.normalB2A *= -1;
      }

      //cout << "collision ptA: " << toVector(toBt(collision.ptA)).transpose() << endl;
      //cout << "collision ptB: " << toVector(toBt(collision.ptB)).transpose() << endl;

      vector<float> sup0, sup1;
      vector<btVector3> ptWorld0, ptWorld1;

      //cout << "shape transforms 0: " << endl;
      //for (auto& t : shape->m_t0i) {
      //  cout << toMatrix(t) << endl << endl;
      //}
      //cout << "shape transforms 1: " << endl;
      //for (auto& t : shape->m_t1i) {
      //  cout << toMatrix(t) << endl << endl;
      //}
      compute_points_and_supports(shape->m_shape, shape->m_t0i, normalWorldFromCast, m_cow, &sup0, &ptWorld0);
      compute_points_and_supports(shape->m_shape, shape->m_t1i, normalWorldFromCast, m_cow, &sup1, &ptWorld1);
      SigmaHullShape* shape0 = new SigmaHullShape(shape->m_shape, shape->m_t0i);
      SigmaHullShape* shape1 = new SigmaHullShape(shape->m_shape, shape->m_t1i);
      //vector<btTransform> trans_m_t1i;
      //for (int i = 0; i < shape->m_t1i.size(); ++i) {
      //  trans_m_t1i.push_back(shape->m_t1i[0].inverseTimes(shape->m_t1i[i]));
      //}
      //compute_points_and_supports(shape0, shape->m_t0i, normalWorldFromCast, m_cow, &sup0, &ptWorld0);
      //compute_points_and_supports(shape1, trans_m_t1i, normalWorldFromCast, m_cow, &sup1, &ptWorld1);

      btVector3 ptOnShape0 = shape0->localGetSupportingVertex(normalWorldFromCast);
      btVector3 ptOnShape1 = shape1->localGetSupportingVertex(normalWorldFromCast);
      //compute_points_and_supports(shape1, shape->m_t1i, normalWorldFromCast, m_cow, &sup1, &ptWorld1);

      delete shape0;
      delete shape1;

      //cout << "pt worlds 0: " << endl;
      //for (auto& pt : ptWorld0) {
      //  cout << toVector(pt).transpose() << endl;
      //}
      //cout << "pt worlds 1: " << endl;
      //for (auto& pt : ptWorld1) {
      //  cout << toVector(pt).transpose() << endl;
      //}

      vector<float> sups0, sups1;
      vector<btVector3> max_ptWorlds0, max_ptWorlds1;
      vector<int> instance_inds0, instance_inds1;
      compute_max_support_points(sup0, ptWorld0, &sups0, &max_ptWorlds0, &instance_inds0);
      compute_max_support_points(sup1, ptWorld1, &sups1, &max_ptWorlds1, &instance_inds1);

      assert(max_ptWorlds0.size()>0);
      assert(max_ptWorlds0.size()<4);
      assert(max_ptWorlds1.size()>0);
      assert(max_ptWorlds1.size()<4);

      const btVector3& ptOnCast = castShapeIsFirst ? cp.m_positionWorldOnA : cp.m_positionWorldOnB;
      
      computeSupportingWeights(max_ptWorlds0, ptOnShape0, collision.mi[0].alpha);
      computeSupportingWeights(max_ptWorlds1, ptOnShape1, collision.mi[1].alpha);
      //cout << "alpha size 0: " << collision.mi[0].alpha.size() << endl;
      //for (auto& i : collision.mi[0].alpha) {
      //  cout << "alpha 0: " << i << endl;
      //}
      //cout << "alpha size 1: " << collision.mi[1].alpha.size() << endl;
      //for (auto& i : collision.mi[1].alpha) {
      //  cout << "alpha 1: " << i << endl;
      //}
      //computeSupportingWeights(max_ptWorlds0, ptOnCast, collision.mi[0].alpha);
      //computeSupportingWeights(max_ptWorlds1, ptOnCast, collision.mi[1].alpha);

      collision.mi[0].instance_ind = instance_inds0;
      collision.mi[1].instance_ind = instance_inds1;

      const float SUPPORT_FUNC_TOLERANCE = .01;
      
      // TODO: this section is _definitely_ problematic. think hard about the math

      if (sups0[0] - sups1[0]> SUPPORT_FUNC_TOLERANCE) {
        collision.time = 0;
        collision.cctype = CCType_Time0;
      }
      else if (sups1[0] - sups0[0]> SUPPORT_FUNC_TOLERANCE) {
        collision.time = 1;
        collision.cctype = CCType_Time1;
      }
      else {
        float l0c = (ptOnCast - max_ptWorlds0[0]).length(), 
              l1c = (ptOnCast - max_ptWorlds1[0]).length();
        //cout << "max point worlds 0: " << endl;
        //for (auto& pt : max_ptWorlds0) {
        //  cout << toVector(pt).transpose() << endl;
        //}
        //cout << "max point worlds 1: " << endl;
        //for (auto& pt : max_ptWorlds1) {
        //  cout << toVector(pt).transpose() << endl;
        //}
        collision.ptB = toOR(max_ptWorlds0[0]);
        collision.ptB1 = toOR(max_ptWorlds1[0]);
        collision.cctype = CCType_Between;
        const float LENGTH_TOLERANCE = .001;
        if ( l0c + l1c < LENGTH_TOLERANCE) {
          collision.time = .5;
        }
        else {
          collision.time = l0c/(l0c + l1c); 
        }
      }
      //cout << "collision ptA: " << toVector(toBt(collision.ptA)).transpose() << endl;
      //cout << "collision ptB: " << toVector(toBt(collision.ptB)).transpose() << endl;
      //cout << "collision ptB1: " << toVector(toBt(collision.ptB1)).transpose() << endl;
      //cout << "collision time: " << collision.time << endl;
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

    void PlotSigmaHull(btCollisionShape* shape, const vector<btTransform>& tfi,
		    CollisionObjectWrapper* cow, vector<OpenRAVE::GraphHandlePtr>& handles, OR::RaveVector<float> color) {
	    if (btConvexShape* convex = dynamic_cast<btConvexShape*>(shape)) {
        vector<btTransform> t0i(tfi.size());
        // transform all the points with respect to the first transform
        for (int i=0; i<tfi.size(); i++) t0i[i] = tfi[0].inverseTimes(tfi[i]);
        SigmaHullShape* shape = new SigmaHullShape(convex, t0i);
        RenderCollisionShape(shape, tfi[0], *boost::const_pointer_cast<OpenRAVE::EnvironmentBase>(m_env), handles, color);
        SetTransparency(handles.back(), 0.2);
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

    void compute_multi_tf(Configuration& rad, const vector<KinBody::LinkPtr>& links, const vector<DblVec>& multi_joints, vector<vector<btTransform> >* output_multi_tf) {
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

    void PlotSigmaHull(Configuration& rad, const vector<KinBody::LinkPtr>& links,
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

    void PlotSigmaHullCast(btCollisionShape* shape, const vector<btTransform>& tf0i, const vector<btTransform>& tf1i,
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

    void PlotSigmaHullCast(Configuration& rad, const vector<KinBody::LinkPtr>& links,
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

    
    virtual void SigmaHullVsAll(Configuration& rad, const vector<KinBody::LinkPtr>& links,
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

    void CheckShapeSigmaHull(btCollisionShape* shape, const vector<btTransform>& tfi,
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

    virtual void SigmaHullCastVsAll(Configuration& rad, const vector<KinBody::LinkPtr>& links,
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

    void CheckShapeSigmaHullCast(btCollisionShape* shape, const vector<btTransform>& tf0i, const vector<btTransform>& tf1i,
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

  };

  typedef boost::shared_ptr<BeliefCollisionChecker> BeliefCollisionCheckerPtr;

  class BeliefCollisionEvaluator : public CollisionEvaluator {
  public:
    virtual void GetCollisionsCached(const DblVec& x, vector<BeliefCollision>& collisions) = 0;
    virtual void CalcCollisions(const DblVec& x, vector<BeliefCollision>& collisions) = 0;
    virtual void CustomPlot(const DblVec& x, std::vector<OR::GraphHandlePtr>& handles) = 0;
  };

  template< class BeliefFuncT >
  class BeliefDiscreteCollisionEvaluator : public BeliefCollisionEvaluator {
  public:
    typedef typename BeliefFuncT::Ptr BeliefFuncPtr;
    BeliefDiscreteCollisionEvaluator(ConfigurationPtr rad, const VarVector& vars, BeliefFuncPtr belief_func, OR::KinBody::LinkPtr endeffector) :
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
      m_cc->SigmaHullVsAll(*m_rad, m_links, dofvals, collisions);
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
        for (int i=0; i<col.mi[isTimestep1].alpha.size(); i++) {
          typename BeliefFuncT::SigmaPointsGradT grad;
          belief_func->sigma_points_grad(belief, col.mi[isTimestep1].instance_ind[i], &grad);

          AffExpr dist_a(col.distance);
          if (linkAFound) {
	          //MatrixXd pos_jac = rad.PositionJacobian(endeffector->GetIndex(), col.ptA);
	          MatrixXd pos_jac = rad.PositionJacobian(itA->second, col.ptA);
            VectorXd dist_grad = toVector3d(col.normalB2A).transpose()*pos_jac*grad;
            //cout << "dist grad A: " << dist_grad.transpose() << endl;
            exprInc(dist_a, varDot(dist_grad, theta_vars));
            exprInc(dist_a, -dist_grad.dot(toVectorXd(theta_vals)));
          }
          if (linkBFound) {
	          MatrixXd pos_jac = rad.PositionJacobian(itB->second, (isTimestep1 && (col.cctype == CCType_Between)) ? col.ptB1 : col.ptB);
            VectorXd dist_grad = -toVector3d(col.normalB2A).transpose()*pos_jac*grad;
            //cout << "dist grad B: " << dist_grad.transpose() << endl;
            exprInc(dist_a, varDot(dist_grad, theta_vars));
            exprInc(dist_a, -dist_grad.dot(toVectorXd(theta_vals)));
          }
          if (linkAFound || linkBFound) {
            exprScale(dist_a, col.mi[isTimestep1].alpha[i]);
            //cout << "alpha: " << col.mi[isTimestep1].alpha[i] << endl;
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
      this->CollisionsToDistances(raw_collisions, m_link2ind, dists);
    }
    virtual void CustomPlot(const DblVec& x, std::vector<OR::GraphHandlePtr>& handles) {
      typename BeliefFuncT::BeliefT belief = getVec(x, this->m_vars);
      vector<DblVec> dofvals;
      belief_func->sigma_points(belief, &dofvals);
      m_cc->PlotSigmaHull(*m_rad, m_links, dofvals, handles);
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

  template< class BeliefFuncT >
  class BeliefContinuousCollisionEvaluator : public BeliefDiscreteCollisionEvaluator<BeliefFuncT> {
  public:
    typedef typename BeliefFuncT::Ptr BeliefFuncPtr;
    BeliefContinuousCollisionEvaluator(ConfigurationPtr rad, const VarVector& vars0, const VarVector& vars1, BeliefFuncPtr belief_func, OR::KinBody::LinkPtr endeffector) :
        BeliefDiscreteCollisionEvaluator<BeliefFuncT>(rad, concat(vars0, vars1), belief_func, endeffector),
        m_vars0(vars0),
        m_vars1(vars1) {}

    void CalcCollisions(const DblVec& x, vector<BeliefCollision>& collisions) {
      typename BeliefFuncT::BeliefT belief0 = getVec(x, this->m_vars0);
      typename BeliefFuncT::BeliefT belief1 = getVec(x, this->m_vars1);
      vector<DblVec> dofvals0, dofvals1;
      this->belief_func->sigma_points(belief0, &dofvals0);
      this->belief_func->sigma_points(belief1, &dofvals1);
      this->m_cc->SigmaHullCastVsAll(*this->m_rad, this->m_links, dofvals0, dofvals1, collisions);
    }

    virtual void CollisionsToDistanceExpressions(const vector<BeliefCollision>& collisions, Configuration& rad,
        const Link2Int& link2ind, const VarVector& theta_vars0, const VarVector& theta_vars1, const DblVec& theta_vals0, const DblVec& theta_vals1, vector<AffExpr>& exprs, bool isTimestep1) {
      vector<AffExpr> exprs0, exprs1;
      BeliefDiscreteCollisionEvaluator<BeliefFuncT>::CollisionsToDistanceExpressions(collisions, rad, link2ind, theta_vars0, theta_vals0, exprs0, false);
      BeliefDiscreteCollisionEvaluator<BeliefFuncT>::CollisionsToDistanceExpressions(collisions, rad, link2ind, theta_vars1, theta_vals1, exprs1, true);

      exprs.resize(exprs0.size());
      for (int i=0; i < exprs0.size(); ++i) {
        exprScale(exprs0[i], (1-collisions[i].time));
        exprScale(exprs1[i], collisions[i].time);
        exprs[i] = AffExpr(0);
        exprInc(exprs[i], exprs0[i]);
        exprInc(exprs[i], exprs1[i]);
        cleanupAff(exprs[i]);
      }
    }

    virtual void CalcDistExpressions(const DblVec& x, vector<AffExpr>& exprs) {
      vector<BeliefCollision> collisions;
      this->GetCollisionsCached(x, collisions);
      DblVec theta0 = getDblVec(x, m_vars0);
      DblVec theta1 = getDblVec(x, m_vars1);
      CollisionsToDistanceExpressions(collisions, *this->m_rad, this->m_link2ind, m_vars0, m_vars1, theta0, theta1, exprs, false);
    }

    virtual void CustomPlot(const DblVec& x, std::vector<OR::GraphHandlePtr>& handles) {
      typename BeliefFuncT::BeliefT belief0 = getVec(x, this->m_vars0);
      typename BeliefFuncT::BeliefT belief1 = getVec(x, this->m_vars1);
      vector<DblVec> dofvals0, dofvals1;
      this->belief_func->sigma_points(belief0, &dofvals0);
      this->belief_func->sigma_points(belief1, &dofvals1);
      this->m_cc->PlotSigmaHullCast(*this->m_rad, this->m_links, dofvals0, dofvals1, handles);
    }
  protected:
    VarVector m_vars0, m_vars1;
  };

}
