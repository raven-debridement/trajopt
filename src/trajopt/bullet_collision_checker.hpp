#pragma once
#include "trajopt/collision_checker.hpp"
#include <btBulletCollisionCommon.h>
#include <BulletCollision/CollisionShapes/btShapeHull.h>
#include <BulletCollision/CollisionDispatch/btConvexConvexAlgorithm.h>
#include <openrave-core.h>
#include "utils/eigen_conversions.hpp"
#include <boost/foreach.hpp>
#include <vector>
#include <iostream>
#include <LinearMath/btConvexHull.h>
#include <utils/stl_to_string.hpp>
#include "utils/logging.hpp"
#include "openrave_userdata_utils.hpp"

using namespace util;
using namespace std;
using namespace trajopt;
using namespace OpenRAVE;

namespace trajopt {

class TRAJOPT_API CollisionObjectWrapper : public btCollisionObject {
public:
  CollisionObjectWrapper(KinBody::Link* link) : m_link(link), m_index(-1) {}
  vector<boost::shared_ptr<void> > m_data;
  KinBody::Link* m_link;
  int m_index; // index into collision matrix
  template<class T>
  void manage(T* t) { // manage memory of this object
    m_data.push_back(boost::shared_ptr<T>(t));
  }
  template<class T>
  void manage(boost::shared_ptr<T> t) {
    m_data.push_back(t);
  }
};
typedef CollisionObjectWrapper COW;
typedef boost::shared_ptr<CollisionObjectWrapper> COWPtr;

class TRAJOPT_API BulletCollisionChecker : public CollisionChecker {
public:
  BulletCollisionChecker(OR::EnvironmentBaseConstPtr env);
  ~BulletCollisionChecker();

  ///////// public interface /////////
  virtual void SetContactDistance(float distance);
  virtual double GetContactDistance() {return m_contactDistance;}
  virtual void PlotCollisionGeometry(vector<OpenRAVE::GraphHandlePtr>& handles);
  virtual void ExcludeCollisionPair(const KinBody::Link& link0, const KinBody::Link& link1) {
    m_excludedPairs.insert(LinkPair(&link0, &link1));
    COW *cow0 = GetCow(&link0), *cow1 = GetCow(&link1);
    if (cow0 && cow1) m_allowedCollisionMatrix(cow0->m_index, cow1->m_index) = 0;
  }
  virtual void IncludeCollisionPair(const KinBody::Link& link0, const KinBody::Link& link1) {
    m_excludedPairs.erase(LinkPair(&link0, &link1));
    COW *cow0 = GetCow(&link0), *cow1 = GetCow(&link1);
    if (cow0 && cow1) m_allowedCollisionMatrix(cow0->m_index, cow1->m_index) = 1;
  }
  // collision checking
  virtual void AllVsAll(vector<Collision>& collisions);
  virtual void LinksVsAll(const vector<KinBody::LinkPtr>& links, vector<Collision>& collisions, short filterMask);
  virtual void LinkVsAll(const KinBody::Link& link, vector<Collision>& collisions, short filterMask);
  virtual void ContinuousCheckTrajectory(const TrajArray& traj, Configuration& rad, vector<Collision>&);
  virtual void CastVsAll(Configuration& rad, const vector<KinBody::LinkPtr>& links, const DblVec& startjoints, const DblVec& endjoints, vector<Collision>& collisions);
  ////
  ///////

  CollisionObjectWrapper* GetCow(const KinBody::Link* link) {
    Link2Cow::iterator it = m_link2cow.find(link);
    return (it == m_link2cow.end()) ? 0 : it->second;
  }
  void SetCow(const KinBody::Link* link, COW* cow) {m_link2cow[link] = cow;}
  void LinkVsAll_NoUpdate(const KinBody::Link& link, vector<Collision>& collisions, short filterMask);
  void UpdateBulletFromRave();
  void AddKinBody(const OR::KinBodyPtr& body);
  void RemoveKinBody(const OR::KinBodyPtr& body);
  void AddAndRemoveBodies(const vector<OR::KinBodyPtr>& curVec, const vector<OR::KinBodyPtr>& prevVec, vector<KinBodyPtr>& addedBodies);
  bool CanCollide(const CollisionObjectWrapper* cow0, const CollisionObjectWrapper* cow1) {
    return m_allowedCollisionMatrix(cow0->m_index, cow1->m_index);
  }
  void SetLinkIndices();
  void UpdateAllowedCollisionMatrix();
  void CheckShapeCast(btCollisionShape* shape, const btTransform& tf0, const btTransform& tf1,
      CollisionObjectWrapper* cow, btCollisionWorld* world, vector<Collision>& collisions);

protected:
  btCollisionWorld* m_world;
  btBroadphaseInterface* m_broadphase;
  btCollisionDispatcher* m_dispatcher;
  btCollisionConfiguration* m_coll_config;
  typedef map<const OR::KinBody::Link*, CollisionObjectWrapper*> Link2Cow;
  Link2Cow m_link2cow;
  double m_contactDistance;
  vector<KinBodyPtr> m_prevbodies;
  typedef std::pair<const KinBody::Link*, const KinBody::Link*> LinkPair;
  set< LinkPair > m_excludedPairs;
  Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> m_allowedCollisionMatrix;
  btCollisionShape* createShapePrimitive(OR::KinBody::Link::GeometryPtr geom, bool useTrimesh, CollisionObjectWrapper* cow) const;
  COWPtr CollisionObjectFromLink(OR::KinBody::LinkPtr link, bool useTrimesh) const;
  virtual void RenderCollisionShape(btCollisionShape* shape, const btTransform& tf,
    OpenRAVE::EnvironmentBase& env, vector<OpenRAVE::GraphHandlePtr>& handles) const;
  void ContinuousCheckShape(btCollisionShape* shape, const vector<btTransform>& transforms,
    KinBody::Link* link, btCollisionWorld* world, vector<Collision>& collisions) const;

};

struct TRAJOPT_API CollisionCollector : public btCollisionWorld::ContactResultCallback {
  std::vector<Collision>& m_collisions;
  const CollisionObjectWrapper* m_cow;
  BulletCollisionChecker* m_cc;

  CollisionCollector(vector<Collision>& collisions, CollisionObjectWrapper* cow, BulletCollisionChecker* cc);
  virtual btScalar addSingleResult(btManifoldPoint& cp,
      const btCollisionObjectWrapper* colObj0Wrap,int partId0,int index0,
      const btCollisionObjectWrapper* colObj1Wrap,int partId1,int index1);
  bool needsCollision(btBroadphaseProxy* proxy0) const;
};

struct KinBodyCollisionData;
typedef boost::shared_ptr<KinBodyCollisionData> CDPtr;
struct KinBodyCollisionData : public OpenRAVE::UserData {
  OpenRAVE::KinBodyWeakPtr body;
  std::vector<KinBody::Link*> links;
  std::vector<COWPtr> cows;
  KinBodyCollisionData(OR::KinBodyPtr _body) : body(_body) {}
};

template <typename T>
void SetDifferences(const vector<T>& A, const vector<T>& B, vector<T>& AMinusB, vector<T>& BMinusA) {
  set<T> Aset, Bset;
  AMinusB.clear();
  BMinusA.clear();
  BOOST_FOREACH(const T& a, A) {
    Aset.insert(a);
  }
  BOOST_FOREACH(const T& b, B) {
    Bset.insert(b);
  }
  BOOST_FOREACH(const T& a, A) {
    if (Bset.count(a) == 0) AMinusB.push_back(a);
  }
  BOOST_FOREACH(const T& b, B) {
    if (Aset.count(b) == 0) BMinusA.push_back(b);
  }
}

struct CastHullShape : public btConvexShape {
public:
  btConvexShape* m_shape;
  btTransform m_t01, m_t10; // T_0_1 = T_w_0^-1 * T_w_1
  CastHullShape(btConvexShape* shape, const btTransform& t01) : m_shape(shape), m_t01(t01) {
    m_shapeType = CUSTOM_CONVEX_SHAPE_TYPE;



  }
  btVector3   localGetSupportingVertex(const btVector3& vec)const {
    btVector3 sv0 = m_shape->localGetSupportingVertex(vec);
    btVector3 sv1 = m_t01*m_shape->localGetSupportingVertex(vec*m_t01.getBasis());
    return (vec.dot(sv0) > vec.dot(sv1)) ? sv0 : sv1;
  }
#if 0
  void project(const btTransform& trans, const btVector3& dir, btScalar& min, btScalar& max) const {
    m_children[0]->project(trans, dir, min, max);
    for (int i=1; i < m_children.size(); ++i) {
      btScalar newmin, newmax;
      m_children[i]->project(trans, dir, newmin, newmax);
      btSetMin(min, newmin);
      btSetMax(max, newmax);
    }
  }
#endif

  //notice that the vectors should be unit length
  void    batchedUnitVectorGetSupportingVertexWithoutMargin(const btVector3* vectors,btVector3* supportVerticesOut,int numVectors) const {
    throw std::runtime_error("not implemented");
  }

  ///getAabb's default implementation is brute force, expected derived classes to implement a fast dedicated version
  void getAabb(const btTransform& t_w0,btVector3& aabbMin,btVector3& aabbMax) const {
    m_shape->getAabb(t_w0, aabbMin, aabbMax);
    btVector3 min1, max1;
    m_shape->getAabb(t_w0*m_t01, min1, max1 );
    aabbMin.setMin(min1);
    aabbMax.setMax(max1);
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
  virtual const char* getName() const {return "CastHull";}
  virtual btVector3 localGetSupportingVertexWithoutMargin(const btVector3& v) const {return localGetSupportingVertex(v);}

  void calculateContactTime(Collision& col) {    
    // float support0 = localGetSupportingVertex(col.)
  }

};


struct CastCollisionCollector : public CollisionCollector {
  CastCollisionCollector(vector<Collision>& collisions, CollisionObjectWrapper* cow, BulletCollisionChecker* cc) :
    CollisionCollector(collisions, cow, cc) {}  
  virtual btScalar addSingleResult(btManifoldPoint& cp,
      const btCollisionObjectWrapper* colObj0Wrap,int partId0,int index0,
      const btCollisionObjectWrapper* colObj1Wrap,int partId1,int index1);
  virtual void GetAverageSupport(const btConvexShape* shape, const btVector3& localNormal, float& outsupport, btVector3& outpt) const;
};
}
