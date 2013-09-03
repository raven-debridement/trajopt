#include "trajopt/collision_checker.hpp"
#include <btBulletCollisionCommon.h>
#include <BulletCollision/CollisionShapes/btShapeHull.h>
#include <BulletCollision/CollisionDispatch/btConvexConvexAlgorithm.h>

#include "BulletCollision/CollisionShapes/btConvexHullShape.h"
#include "BulletCollision/NarrowPhaseCollision/btGjkEpaPenetrationDepthSolver.h"
#include "BulletCollision/NarrowPhaseCollision/btGjkPairDetector.h"
#include "BulletCollision/NarrowPhaseCollision/btPointCollector.h"
#include "BulletCollision/NarrowPhaseCollision/btVoronoiSimplexSolver.h"
#include "BulletCollision/NarrowPhaseCollision/btConvexPenetrationDepthSolver.h"

#include <openrave-core.h>
#include "utils/eigen_conversions.hpp"
#include <boost/foreach.hpp>
#include <vector>
#include <iostream>
#include <LinearMath/btConvexHull.h>
#include <utils/stl_to_string.hpp>
#include "utils/logging.hpp"
#include "openrave_userdata_utils.hpp"
#include "bullet_collision_checker.hpp"

using namespace util;
using namespace std;
using namespace trajopt;
using namespace OpenRAVE;

namespace {

#define METERS 
// there's some scale-dependent parameters. By convention I'll put METERS to mark it
const float MARGIN = 0;

#if 1
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
ostream &operator<<(ostream &stream, const btVector3& v) {
  stream << v.x() << " " << v.y() << " " << v.z();
  return stream;
}
ostream &operator<<(ostream &stream, const btQuaternion& v) {
  stream << v.w() << " " << v.x() << " " << v.y() << " " << v.z();
  return stream;
}
ostream &operator<<(ostream &stream, const btTransform& v) {
  stream << v.getOrigin() << " " << v.getRotation();
  return stream;
}
#pragma GCC diagnostic pop
#endif



inline const KinBody::Link* getLink(const btCollisionObject* o) {
  return static_cast<const CollisionObjectWrapper*>(o)->m_link;
}


// only used for AllVsAll
void nearCallback(btBroadphasePair& collisionPair,
    btCollisionDispatcher& dispatcher, const btDispatcherInfo& dispatchInfo) {
  BulletCollisionChecker* cc = static_cast<BulletCollisionChecker*>(dispatcher.m_userData);
  if ( cc->CanCollide(static_cast<CollisionObjectWrapper*>(collisionPair.m_pProxy0->m_clientObject),
                      static_cast<CollisionObjectWrapper*>(collisionPair.m_pProxy1->m_clientObject)))
    dispatcher.defaultNearCallback(collisionPair, dispatcher, dispatchInfo);
}

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

bool isIdentity(const OpenRAVE::Transform& T) {
  float e = 1e-6;
  return
      fabs(T.trans.x) < e &&
      fabs(T.trans.y) < e &&
      fabs(T.trans.z) < e &&
      fabs(T.rot[0]-1) < e &&
      fabs(T.rot[1]) < e &&
      fabs(T.rot[2]) < e &&
      fabs(T.rot[3]) < e;
}

}

namespace trajopt {




btCollisionShape* BulletCollisionChecker::createShapePrimitive(OR::KinBody::Link::GeometryPtr geom, bool useTrimesh, CollisionObjectWrapper* cow) const {

  btCollisionShape* subshape=0;

#if OPENRAVE_VERSION_MINOR <= 8
    #define GT_Box KinBody::Link::GEOMPROPERTIES::GeomBox 
    #define GT_Sphere KinBody::Link::GEOMPROPERTIES::GeomSphere 
    #define GT_Cylinder KinBody::Link::GEOMPROPERTIES::GeomCylinder 
    #define GT_TriMesh KinBody::Link::GEOMPROPERTIES::GeomTrimesh 
    #define TriMesh KinBody::Link::TRIMESH
#endif

  switch (geom->GetType()) {
  case OpenRAVE::GT_Box:
    subshape = new btBoxShape(toBt(geom->GetBoxExtents()));
    break;
  case OpenRAVE::GT_Sphere:
    subshape = new btSphereShape(geom->GetSphereRadius());
    break;
  case OpenRAVE::GT_Cylinder:
    // cylinder axis aligned to Y
  {
    float r = geom->GetCylinderRadius(), h = geom->GetCylinderHeight() / 2;
    subshape = new btCylinderShapeZ(btVector3(r, r, h));
    break;
  }
  case OpenRAVE::GT_TriMesh: {
    const OpenRAVE::TriMesh &mesh = geom->GetCollisionMesh();
    assert(mesh.indices.size() >= 3);
    boost::shared_ptr<btTriangleMesh> ptrimesh(new btTriangleMesh());

    for (size_t i = 0; i < mesh.indices.size(); i += 3) {
      ptrimesh->addTriangle(toBt(mesh.vertices[mesh.indices[i]]), toBt(mesh.vertices[mesh.indices[i + 1]]),
              toBt(mesh.vertices[mesh.indices[i + 2]]));
    }

    if (useTrimesh) {
      subshape = new btBvhTriangleMeshShape(ptrimesh.get(), true);
      cow->manage(ptrimesh);
    }
    else { // CONVEX HULL
      btConvexTriangleMeshShape convexTrimesh(ptrimesh.get());
      convexTrimesh.setMargin(MARGIN); // margin: hull padding
      //Create a hull shape to approximate Trimesh

      bool useShapeHull;

      btShapeHull shapeHull(&convexTrimesh);
      if (mesh.vertices.size() >= 50) {
        bool success = shapeHull.buildHull(-666); // note: margin argument not used
        if (!success) LOG_WARN("shapehull convex hull failed! falling back to original vertices");
        useShapeHull = success;
      }
      else {
        useShapeHull = false;
      }

      btConvexHullShape *convexShape = new btConvexHullShape();
      subshape = convexShape;
      if (useShapeHull) {
        for (int i = 0; i < shapeHull.numVertices(); ++i)
          convexShape->addPoint(shapeHull.getVertexPointer()[i]);
        break;
      }
      else {
        for (int i = 0; i < mesh.vertices.size(); ++i)
          convexShape->addPoint(toBt(mesh.vertices[i]));
        break;
      }
      
    }     
  }
  default:
    assert(0 && "unrecognized collision shape type");
    break;
  }
  return subshape;
}


COWPtr BulletCollisionChecker::CollisionObjectFromLink(OR::KinBody::LinkPtr link, bool useTrimesh) const {
  LOG_DEBUG("creating bt collision object from from %s",link->GetName().c_str());

  const std::vector<boost::shared_ptr<OpenRAVE::KinBody::Link::Geometry> > & geometries=link->GetGeometries();

  if (geometries.empty()) return COWPtr();

  COWPtr cow(new CollisionObjectWrapper(link.get()));

  if ((link->GetGeometries().size() == 1) && isIdentity(link->GetGeometry(0)->GetTransform())) {
    LOG_DEBUG("using identity for %s", link->GetName().c_str());
    btCollisionShape* shape = createShapePrimitive(link->GetGeometry(0), useTrimesh, cow.get());
    shape->setMargin(MARGIN);
    cow->manage(shape);
    cow->setCollisionShape(shape);

  }
  else {
    LOG_DEBUG("NOT using identity for %s", link->GetName().c_str());
    btCompoundShape* compound = new btCompoundShape(/*dynamicAABBtree=*/false);
    cow->manage(compound);
    compound->setMargin(MARGIN); //margin: compound. seems to have no effect when positive but has an effect when negative
    cow->setCollisionShape(compound);

    BOOST_FOREACH(const boost::shared_ptr<OpenRAVE::KinBody::Link::Geometry>& geom, geometries) {

      btCollisionShape* subshape = createShapePrimitive(geom, useTrimesh, cow.get());
      if (subshape != NULL) {
        cow->manage(subshape);
        subshape->setMargin(MARGIN);
        btTransform geomTrans = toBt(geom->GetTransform());
        compound->addChildShape(geomTrans, subshape);
      }
    }

  }

  cow->setWorldTransform(toBt(link->GetTransform()));

  return cow;
}



void BulletCollisionChecker::RenderCollisionShape(btCollisionShape* shape, const btTransform& tf,
    OpenRAVE::EnvironmentBase& env, vector<OpenRAVE::GraphHandlePtr>& handles) const {

  typedef map<btCollisionShape*, HullResult > Shape2Inds;
  Shape2Inds gHullCache;

  switch (shape->getShapeType()) {
  case COMPOUND_SHAPE_PROXYTYPE: {
    btCompoundShape* compound = static_cast<btCompoundShape*>(shape);
    for (int i = 0; i < compound->getNumChildShapes(); ++i) {
      RenderCollisionShape(compound->getChildShape(i),
          tf * compound->getChildTransform(i), env, handles);
    }
    break;
  }
  case CONVEX_HULL_SHAPE_PROXYTYPE: {
    btConvexHullShape* convex = static_cast<btConvexHullShape*>(shape);

    Shape2Inds::iterator it = gHullCache.find(convex);

    btAlignedObjectArray<unsigned int> inds;
    HullResult hr;
    if ( it != gHullCache.end() )
      hr = it->second;
    else {

      HullDesc hd;
      hd.mFlags = QF_TRIANGLES;
      hd.mVcount = convex->getNumPoints();
      hd.mVertices = convex->getPoints();
      hd.mVertexStride = sizeof(btVector3);
      HullLibrary hl;

      if (hl.CreateConvexHull(hd, hr) == QE_FAIL) {
        LOG_ERROR("convex hull computation failed on shape with %i vertices", convex->getNumPoints());
        hr.mNumFaces = 0;
      }
      else {
      }
      gHullCache[convex] = hr;
    }

    if (hr.mNumFaces > 0) {
      vector<btVector3> tverts(hr.mNumOutputVertices);
      for (int i=0; i < tverts.size(); ++i) tverts[i] = tf * hr.m_OutputVertices[i];


      handles.push_back(env.drawtrimesh((float*)&tverts[0], 16,
          (int*) &hr.m_Indices[0], hr.mNumFaces, OR::RaveVector<float>(1,1,1,.1)));
    }
    break;


  }

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

      handles.push_back(env.drawtrimesh((float*)tf_vertices, 16, (int*) indices, num_triangles, OR::RaveVector<float>(1,1,1,1)));
    } else {
      LOG_INFO("not rendering shape of type %i", shape->getShapeType());
    }
    break;
  }
}

BulletCollisionChecker::BulletCollisionChecker(OR::EnvironmentBaseConstPtr env) :
  CollisionChecker(env) {
  m_coll_config = new btDefaultCollisionConfiguration();
  m_dispatcher = new btCollisionDispatcher(m_coll_config);
  m_broadphase = new btDbvtBroadphase();
  m_world = new btCollisionWorld(m_dispatcher, m_broadphase, m_coll_config);
  m_dispatcher->registerCollisionCreateFunc(BOX_SHAPE_PROXYTYPE,BOX_SHAPE_PROXYTYPE,
      m_coll_config->getCollisionAlgorithmCreateFunc(CONVEX_SHAPE_PROXYTYPE, CONVEX_SHAPE_PROXYTYPE));
  m_dispatcher->setNearCallback(&nearCallback);
  m_dispatcher->m_userData = this;
  SetContactDistance(.05);
  UpdateBulletFromRave();
}

BulletCollisionChecker::~BulletCollisionChecker() {
  delete m_world;
  delete m_broadphase;
  delete m_dispatcher;
  delete m_coll_config;
}


void BulletCollisionChecker::SetContactDistance(float dist) {
  LOG_DEBUG("setting contact distance to %.2f", dist);
  m_contactDistance = dist;
  SHAPE_EXPANSION = btVector3(1,1,1)*dist;
  gContactBreakingThreshold = 2.001*dist; // wtf. when I set it to 2.0 there are no contacts with distance > 0
  btCollisionObjectArray& objs = m_world->getCollisionObjectArray();
  for (int i=0; i < objs.size(); ++i) {
    objs[i]->setContactProcessingThreshold(dist);
  }
  btCollisionDispatcher* dispatcher = static_cast<btCollisionDispatcher*>(m_world->getDispatcher());
  dispatcher->setDispatcherFlags(dispatcher->getDispatcherFlags() & ~btCollisionDispatcher::CD_USE_RELATIVE_CONTACT_BREAKING_THRESHOLD);
}


void BulletCollisionChecker::AllVsAll(vector<Collision>& collisions) {
  LOG_WARN("WARNING: AllVsAll seems to be broken! (since a8f8da01)");
  UpdateBulletFromRave();
  LOG_DEBUG("AllVsAll");
  m_world->performDiscreteCollisionDetection();
  int numManifolds = m_dispatcher->getNumManifolds();
  LOG_DEBUG("number of manifolds: %i", numManifolds);
  for (int i = 0; i < numManifolds; ++i) {
    btPersistentManifold* contactManifold = m_dispatcher->getManifoldByIndexInternal(i);
    int numContacts = contactManifold->getNumContacts();
    LOG_DEBUG("number of contacts in manifold %i: %i", i, numContacts);
    const CollisionObjectWrapper* objA = static_cast<const CollisionObjectWrapper*>(contactManifold->getBody0());
    const CollisionObjectWrapper* objB = static_cast<const CollisionObjectWrapper*>(contactManifold->getBody1());
    for (int j = 0; j < numContacts; ++j) {
      btManifoldPoint& pt = contactManifold->getContactPoint(j);
//      stringstream ss; ss << pt.m_localPointA << " | " << pt.m_localPointB;
//      LOG_DEBUG("local pts: %s\n",ss.str().c_str());
      // adjustContactPoint(pt, objA, objB);

      const KinBody::Link* bodyA = objA->m_link;
      const KinBody::Link* bodyB = objB->m_link;

      if (CanCollide(objA, objB)) {
        collisions.push_back(Collision(bodyA, bodyB, toOR(pt.getPositionWorldOnA()), toOR(pt.getPositionWorldOnB()),
            toOR(pt.m_normalWorldOnB), pt.m_distance1, 1./numContacts));
      }
      else {
        LOG_DEBUG("ignoring collision between %s and %s", bodyA->GetName().c_str(), bodyB->GetName().c_str());
        assert(0 && "this shouldn't happen because we're filtering at narrowphase");
      }
      LOG_DEBUG("%s - %s collided", bodyA->GetName().c_str(), bodyB->GetName().c_str());
    }
    // caching helps performance, but for optimization the cost should not be history-dependent
    contactManifold->clearManifold();
  }
}

void BulletCollisionChecker::LinksVsAll(const vector<KinBody::LinkPtr>& links, vector<Collision>& collisions, short filterMask) {
//  AllVsAll(collisions);
//  return;

  UpdateBulletFromRave();
  m_world->updateAabbs();
  
  for (int i=0; i < links.size(); ++i) {
    LinkVsAll_NoUpdate(*links[i], collisions, filterMask);
  }
}


void BulletCollisionChecker::LinkVsAll(const KinBody::Link& link, vector<Collision>& collisions, short filterMask) {
  UpdateBulletFromRave();
  LinkVsAll_NoUpdate(link, collisions, filterMask);
}

void BulletCollisionChecker::LinkVsAll_NoUpdate(const KinBody::Link& link, vector<Collision>& collisions, short filterMask) {
  if (link.GetGeometries().empty()) return;
  CollisionObjectWrapper* cow = GetCow(&link);
  CollisionCollector cc(collisions, cow, this);
  cc.m_collisionFilterMask = filterMask;
  m_world->contactTest(cow, cc);
}

bool BulletCollisionChecker::SegmentVsAll(const btVector3& pt0, const btVector3& pt1) {
  btCollisionWorld::ClosestRayResultCallback RayCallback(pt0, pt1);
  m_world->rayTest(pt0, pt1, RayCallback);
  return RayCallback.hasHit();
}



void BulletCollisionChecker::AddKinBody(const OR::KinBodyPtr& body) {
  CDPtr cd(new KinBodyCollisionData(body));

  int filterGroup = body->IsRobot() ? RobotFilter : KinBodyFilter;
  const vector<OR::KinBody::LinkPtr> links = body->GetLinks();

  trajopt::SetUserData(*body, "bt", cd);
  
  bool useTrimesh = trajopt::GetUserData(*body, "bt_use_trimesh");
  BOOST_FOREACH(const OR::KinBody::LinkPtr& link, links) {
    if (link->GetGeometries().size() > 0) {
      COWPtr new_cow = CollisionObjectFromLink(link, useTrimesh); 
      if (new_cow) {
        SetCow(link.get(), new_cow.get());
        m_world->addCollisionObject(new_cow.get(), filterGroup);
        new_cow->setContactProcessingThreshold(m_contactDistance);
        LOG_DEBUG("added collision object for  link %s", link->GetName().c_str());
        cd->links.push_back(link.get());
        cd->cows.push_back(new_cow);
      }
      else {
        LOG_WARN("ignoring link %s", link->GetName().c_str());
      }
    }
  }

}
void BulletCollisionChecker::RemoveKinBody(const OR::KinBodyPtr& body) {
  LOG_DEBUG("removing %s", body->GetName().c_str());
  BOOST_FOREACH(const OR::KinBody::LinkPtr& link, body->GetLinks()) {
    CollisionObjectWrapper* cow = GetCow(link.get());
    if (cow) {
      m_world->removeCollisionObject(cow);
      m_link2cow.erase(link.get());      
    }
    
  }
  trajopt::RemoveUserData(*body, "bt");
}


void BulletCollisionChecker::AddCastHullShape(Configuration& rad0, Configuration& rad1, const vector<KinBody::LinkPtr>& links, const DblVec& startjoints, const DblVec endjoints) {
  // Almost copied from CastVsAll
  Configuration::SaverPtr saver = rad0.Save();
  rad0.SetDOFValues(startjoints);
  int nlinks = links.size();
  vector<btTransform> tbefore(nlinks), tafter(nlinks);
  for (int i=0; i < nlinks; ++i) {
    tbefore[i] = toBt(links[i]->GetTransform());
  }
  rad1.SetDOFValues(endjoints);
  for (int i=0; i < nlinks; ++i) {
    tafter[i] = toBt(links[i]->GetTransform());
  }
  rad0.SetDOFValues(startjoints);
  bool useTrimesh = trajopt::GetUserData(*links[0]->GetParent(), "bt_use_trimesh");
  CDPtr cd = boost::static_pointer_cast<KinBodyCollisionData>(trajopt::GetUserData(*links[0]->GetParent(), "bt"));
  for (int i = 0; i < cd->cows.size(); ++i) {
    m_managed_cows.push_back(cd->cows[i]);//CowPtr(GetCow(links[i].get())));
  }
  for (int i = 0; i < nlinks; ++i) {
    if (links[i]->GetGeometries().size() > 0) {
      COWPtr cow = CollisionObjectFromLink(links[i], useTrimesh); 
      cow->m_dynamic = false;
      //cow->manage(boost::shared_ptr<btCollisionShape>(cow->getCollisionShape()));
      //assert(m_link2cow[links[i].get()] != NULL);
      //CollisionObjectWrapper* cow = m_link2cow[links[i].get()];
      m_managed_cows.push_back(cow);
      AddCastHullShape(cow->getCollisionShape(), tbefore[i], tafter[i], cow.get(), m_world);
    }
  }
}

void BulletCollisionChecker::AddCastHullShape(btCollisionShape* shape, const btTransform& tf0, const btTransform& tf1,
    CollisionObjectWrapper* cow, btCollisionWorld* world) {
  if (btConvexShape* convex = dynamic_cast<btConvexShape*>(shape)) {
    //boost::shared_ptr<btConvexShape> convex_ptr(convex);
    boost::shared_ptr<CastHullShape> shape(new CastHullShape(convex, tf0.inverseTimes(tf1)));
    COWPtr obj(new CollisionObjectWrapper(cow->m_link));
    obj->setCollisionShape(shape.get());
    obj->setWorldTransform(tf0);
    obj->m_index = cow->m_index;
    obj->manage(shape);
    //obj->manage(convex_ptr);
    world->addCollisionObject(obj.get(), KinBodyFilter);
    obj->setContactProcessingThreshold(m_contactDistance);
    m_managed_cows.push_back(obj);
  } else if (btCompoundShape* compound = dynamic_cast<btCompoundShape*>(shape)) {
    for (int i = 0; i < compound->getNumChildShapes(); ++i) {
      AddCastHullShape(compound->getChildShape(i), tf0*compound->getChildTransform(i), tf1*compound->getChildTransform(i), cow, world);
    }
  } else {
    throw std::runtime_error("I can only add cast hull of convex shapes and compound shapes made of convex shapes");
  }
}

void BulletCollisionChecker::AddAndRemoveBodies(const vector<KinBodyPtr>& curVec, const vector<KinBodyPtr>& prevVec, vector<KinBodyPtr>& toAdd) {
  vector<KinBodyPtr> toRemove;
  SetDifferences(curVec, prevVec, toAdd, toRemove);
  BOOST_FOREACH(const KinBodyPtr& body, toAdd) {
    //assert(!trajopt::GetUserData(*body, "bt"));
    AddKinBody(body);
  }
  BOOST_FOREACH(const KinBodyPtr& body, toRemove) {
    RemoveKinBody(body);
  }
  SetLinkIndices();
}

void BulletCollisionChecker::SetLinkIndices() {
  btCollisionObjectArray& objs = m_world->getCollisionObjectArray();
  for (int i=0; i < objs.size(); ++i) {
    CollisionObjectWrapper* cow = static_cast<CollisionObjectWrapper*>(objs[i]);
    cow->m_index = i;
  }
  m_allowedCollisionMatrix.resize(objs.size(), objs.size());
  m_allowedCollisionMatrix.setOnes();
}

void BulletCollisionChecker::UpdateAllowedCollisionMatrix() {
  BOOST_FOREACH(const LinkPair& pair, m_excludedPairs) {
    const KinBody::Link* linkA = pair.first;
    const KinBody::Link* linkB = pair.second;
    const CollisionObjectWrapper* cowA = GetCow(linkA);
    const CollisionObjectWrapper* cowB = GetCow(linkB);
    if (cowA != NULL && cowB != NULL) {
      m_allowedCollisionMatrix(cowA->m_index, cowB->m_index) = 0;
      m_allowedCollisionMatrix(cowB->m_index, cowA->m_index) = 0;
    }
  }
}

void BulletCollisionChecker::UpdateBulletFromRave() {
  vector<OR::KinBodyPtr> bodies, addedBodies;
  m_env->GetBodies(bodies);

  if (bodies.size() != m_prevbodies.size() || !std::equal(bodies.begin(), bodies.end(), m_prevbodies.begin())) {
    LOG_DEBUG("need to add and remove stuff");
    AddAndRemoveBodies(bodies, m_prevbodies, addedBodies);
    m_prevbodies=bodies;
    float contactDistanceOld = GetContactDistance();
    SetContactDistance(.1 METERS);
    BOOST_FOREACH(const KinBodyPtr& body, addedBodies) {
      IgnoreZeroStateSelfCollisions(body);
    }
    SetContactDistance(contactDistanceOld);
    UpdateAllowedCollisionMatrix();
  }
  else {
    LOG_DEBUG("don't need to add or remove stuff");
  }

  btCollisionObjectArray& objs = m_world->getCollisionObjectArray();
  LOG_DEBUG("%i objects in bullet world", objs.size());
  for (int i=0; i < objs.size(); ++i) {
    CollisionObjectWrapper* cow = static_cast<CollisionObjectWrapper*>(objs[i]);
    if (cow->isDynamic()) {
      cow->setWorldTransform(toBt(cow->m_link->GetTransform()));
    }
  }

}


void BulletCollisionChecker::PlotCollisionGeometry(vector<OpenRAVE::GraphHandlePtr>& handles) {
  UpdateBulletFromRave();
  btCollisionObjectArray& objs = m_world->getCollisionObjectArray();
  LOG_DEBUG("%i objects in bullet world", objs.size());
  for (int i=0; i < objs.size(); ++i) {
    RenderCollisionShape(objs[i]->getCollisionShape(), objs[i]->getWorldTransform(), *boost::const_pointer_cast<OpenRAVE::EnvironmentBase>(m_env), handles);
  }
}

CollisionCollector::CollisionCollector(vector<Collision>& collisions, CollisionObjectWrapper* cow, BulletCollisionChecker* cc) :
  m_collisions(collisions), m_cow(cow), m_cc(cc) {}

btScalar CollisionCollector::addSingleResult(btManifoldPoint& cp,
    const btCollisionObjectWrapper* colObj0Wrap,int partId0,int index0,
    const btCollisionObjectWrapper* colObj1Wrap,int partId1,int index1) {
  if (cp.m_distance1 > m_cc->GetContactDistance()) return 0;
  const KinBody::Link* linkA = getLink(colObj0Wrap->getCollisionObject());
  const KinBody::Link* linkB = getLink(colObj1Wrap->getCollisionObject());
  m_collisions.push_back(Collision(linkA, linkB, toOR(cp.m_positionWorldOnA), toOR(cp.m_positionWorldOnB),
      toOR(cp.m_normalWorldOnB), cp.m_distance1));
  LOG_DEBUG("CollisionCollector: adding collision %s-%s (%.4f)", linkA->GetName().c_str(), linkB->GetName().c_str(), cp.m_distance1);
  return 1;
}
bool CollisionCollector::needsCollision(btBroadphaseProxy* proxy0) const {
  return (proxy0->m_collisionFilterGroup & m_collisionFilterMask)
      && (m_collisionFilterGroup & proxy0->m_collisionFilterMask)
      && m_cc->CanCollide(m_cow, static_cast<CollisionObjectWrapper*>(proxy0->m_clientObject));
}

}




////////// Continuous collisions ////////////////////////

namespace {

vector<btTransform> rightMultiplyAll(const vector<btTransform>& xs, const btTransform& y) {
  vector<btTransform> out(xs.size());
  for (int i=0; i < xs.size(); ++i) out[i] = xs[i]*y;
  return out;
}


}

namespace trajopt {

void BulletCollisionChecker::ContinuousCheckShape(btCollisionShape* shape, const vector<btTransform>& transforms,
    KinBody::Link* link, btCollisionWorld* world, vector<Collision>& collisions) const {
  if (btConvexShape* convex = dynamic_cast<btConvexShape*>(shape)) {
    for (int i=0; i < transforms.size()-1; ++i) {
      btCollisionWorld::ClosestConvexResultCallback ccc(btVector3(NAN, NAN, NAN), btVector3(NAN, NAN, NAN));
      ccc.m_collisionFilterMask = KinBodyFilter;
      world->convexSweepTest(convex, transforms[i], transforms[i+1], ccc, 0);
      if (ccc.hasHit()) {
        collisions.push_back(Collision(link, getLink(ccc.m_hitCollisionObject),
            toOR(ccc.m_hitPointWorld), toOR(ccc.m_hitPointWorld), toOR(ccc.m_hitNormalWorld), 0, 1, i+ccc.m_closestHitFraction));
      }
    }
  }
  else if (btCompoundShape* compound = dynamic_cast<btCompoundShape*>(shape)) {
    for (int i = 0; i < compound->getNumChildShapes(); ++i) {
      ContinuousCheckShape(compound->getChildShape(i), rightMultiplyAll(transforms, compound->getChildTransform(i)),  link, world, collisions);
    }
  }
  else {
    throw std::runtime_error("I can only continuous collision check convex shapes and compound shapes made of convex shapes");
  }

}


void BulletCollisionChecker::ContinuousCheckTrajectory(const TrajArray& traj, Configuration& rad, vector<Collision>& collisions) {
  UpdateBulletFromRave();
  m_world->updateAabbs();

  // first calculate transforms of all the relevant links at each step
  vector<KinBody::LinkPtr> links;
  vector<int> link_inds;
  rad.GetAffectedLinks(links, true, link_inds);


  // don't need to remove them anymore because now I only check collisions
  // against KinBodyFilter stuff
  // remove them, because we can't check moving stuff against each other
  vector<CollisionObjectWrapper*> cows;
  BOOST_FOREACH(KinBody::LinkPtr& link, links) {
    CollisionObjectWrapper* cow = GetCow(link.get());
    assert(cow != NULL);
    cows.push_back(cow);
#if 0
    m_world->removeCollisionObject(cow);
#endif
  }


  typedef vector<btTransform> TransformVec;
  vector<TransformVec> link2transforms(links.size(), TransformVec(traj.rows()));
  Configuration::SaverPtr save = rad.Save();

  for (int iStep=0; iStep < traj.rows(); ++iStep) {
    rad.SetDOFValues(toDblVec(traj.row(iStep)));
    for (int iLink = 0; iLink < links.size(); ++iLink) {
      link2transforms[iLink][iStep] = toBt(links[iLink]->GetTransform());
    }
  }

  for (int iLink = 0; iLink < links.size(); ++iLink) {
    ContinuousCheckShape(cows[iLink]->getCollisionShape(), link2transforms[iLink], links[iLink].get(), m_world, collisions);
  }

#if 0
  // add them back
  BOOST_FOREACH(CollisionObjectWrapper* cow, cows) {
    m_world->addCollisionObject(cow);
  }
#endif
}

#if 0
class CompoundHullShape : public btConvexShape {
  std::vector<btConvexHullShape*> m_children;
  btVector3   localGetSupportingVertex(const btVector3& vec)const {
    btVector3 sv = m_children[0]->localGetSupportingVertex(vec);
    float support = sv.dot(vec);
    for (int i=1; i < m_children.size(); ++i) {
      btVector3 newsv = m_children[i]->localGetSupportingVertex(vec);
      float newsupport = vec.dot(newsv);
      if (newsupport > support) {
        support = newsupport;
        sv = newsv;
      }
    }
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
  void getAabb(const btTransform& t,btVector3& aabbMin,btVector3& aabbMax) const {
    m_children[0]->getAabb(t, aabbMin, aabbMax);
    for (int i=1; i < m_children.size(); ++i) {
      btVector3 newmin, newmax;
      m_children[i]->getAabb(t, newmin, newmax);
      aabbMin.setMin(newmin);
      aabbMax.setMax(newmax);
    }
  }

  virtual void getAabbSlow(const btTransform& t,btVector3& aabbMin,btVector3& aabbMax) const {
    throw std::runtime_error("shouldn't happen");
  }

  virtual void    setLocalScaling(const btVector3& scaling) {}
  virtual const btVector3& getLocalScaling() const {return btVector3(1,1,1);}

  virtual void    setMargin(btScalar margin) {}
  virtual btScalar    getMargin() const {return 0;}

  virtual int     getNumPreferredPenetrationDirections() const {return 0;}
  virtual void    getPreferredPenetrationDirection(int index, btVector3& penetrationVector) const=0;


};
#endif

void GetAverageSupport(const btConvexShape* shape, const btVector3& localNormal, float& outsupport, btVector3& outpt) {
  btVector3 ptSum(0,0,0);
  float ptCount = 0;
  float maxSupport=-1000;
  const float EPSILON = 1e-3;
  const btPolyhedralConvexShape* pshape = dynamic_cast<const btPolyhedralConvexShape*>(shape);
  if (pshape) {
    int nPts = pshape->getNumVertices();

    for (int i=0; i < nPts; ++i) {
      btVector3 pt;
      pshape->getVertex(i, pt);
      float sup  = pt.dot(localNormal);
      if (sup > maxSupport + EPSILON) {
        ptCount=1;
        ptSum = pt;
        maxSupport = sup;
      }
      else if (sup < maxSupport - EPSILON) {
      }
      else {
        ptCount += 1;
        ptSum += pt;
      }
    }
    outsupport = maxSupport;
    outpt = ptSum / ptCount;
  }
  else  {
    outpt = shape->localGetSupportingVertexWithoutMargin(localNormal);
    outsupport = localNormal.dot(outpt);
  }
}

btScalar CastCollisionCollector::addSingleResult(btManifoldPoint& cp,
    const btCollisionObjectWrapper* colObj0Wrap,int partId0,int index0,
    const btCollisionObjectWrapper* colObj1Wrap,int partId1,int index1) {      
      float retval = CollisionCollector::addSingleResult(cp, colObj0Wrap,partId0,index0, colObj1Wrap,partId1,index1); // call base class func
      if (retval == 1) { // if contact was added
        bool castShapeIsFirst =  (colObj0Wrap->getCollisionObject() == m_cow);
        btVector3 normalWorldFromCast = -(castShapeIsFirst ? 1 : -1) * cp.m_normalWorldOnB;
        const CastHullShape* shape = dynamic_cast<const CastHullShape*>((castShapeIsFirst ? colObj0Wrap : colObj1Wrap)->getCollisionObject()->getCollisionShape());
        assert(!!shape);
        btTransform tfWorld0 = m_cow->getWorldTransform();
        btTransform tfWorld1 = m_cow->getWorldTransform() * shape->m_t01;
        btVector3 normalLocal0 = normalWorldFromCast * tfWorld0.getBasis();
        btVector3 normalLocal1 = normalWorldFromCast * tfWorld1.getBasis();

        Collision& col = m_collisions.back();
        const float SUPPORT_FUNC_TOLERANCE = .01 METERS;

        if (castShapeIsFirst) {
          swap(col.ptA, col.ptB);
          swap(col.linkA, col.linkB);
          col.normalB2A *= -1;
        }

#if 0
        btVector3 ptWorld0 = tfWorld0*shape->m_shape->localGetSupportingVertex(normalLocal0);
        btVector3 ptWorld1 = tfWorld1*shape->m_shape->localGetSupportingVertex(normalLocal1);
#else
        btVector3 ptLocal0;
        float localsup0;
        GetAverageSupport(shape->m_shape, normalLocal0, localsup0, ptLocal0);
        btVector3 ptWorld0 = tfWorld0 * ptLocal0;
        btVector3 ptLocal1;
        float localsup1;
        GetAverageSupport(shape->m_shape, normalLocal1, localsup1, ptLocal1);
        btVector3 ptWorld1 = tfWorld1 * ptLocal1;



#endif
        float sup0 = normalWorldFromCast.dot(ptWorld0);
        float sup1 = normalWorldFromCast.dot(ptWorld1);



        // TODO: this section is potentially problematic. think hard about the math
        if (sup0 - sup1 > SUPPORT_FUNC_TOLERANCE) {
          col.time = 0;
          col.cctype = CCType_Time0;
        }
        else if (sup1 - sup0 > SUPPORT_FUNC_TOLERANCE) {
          col.time = 1;
          col.cctype = CCType_Time1;
        }
        else {
          const btVector3& ptOnCast = castShapeIsFirst ? cp.m_positionWorldOnA : cp.m_positionWorldOnB;
          float l0c = (ptOnCast - ptWorld0).length(), 
                l1c = (ptOnCast - ptWorld1).length();

          col.ptB = toOR(ptWorld0);
          col.ptB1 = toOR(ptWorld1);
          col.cctype = CCType_Between;

          const float LENGTH_TOLERANCE = .001 METERS;

          if ( l0c + l1c < LENGTH_TOLERANCE) {

            col.time = .5;
          }
          else {
            col.time = l0c/(l0c + l1c); 
          }

        }
          
      }
      return retval;          
}


void BulletCollisionChecker::CheckShapeCast(btCollisionShape* shape, const btTransform& tf0, const btTransform& tf1,
    CollisionObjectWrapper* cow, btCollisionWorld* world, vector<Collision>& collisions) {
  if (btConvexShape* convex = dynamic_cast<btConvexShape*>(shape)) {
    CastHullShape* shape = new CastHullShape(convex, tf0.inverseTimes(tf1));
    CollisionObjectWrapper* obj = new CollisionObjectWrapper(cow->m_link);
    obj->setCollisionShape(shape);
    obj->setWorldTransform(tf0);
    obj->m_index = cow->m_index;
    CastCollisionCollector cc(collisions, obj, this);
    cc.m_collisionFilterMask = KinBodyFilter;
    // cc.m_collisionFilterGroup = cow->m_collisionFilterGroup;
    world->contactTest(obj, cc);
    delete obj;
    delete shape;
  }
  else if (btCompoundShape* compound = dynamic_cast<btCompoundShape*>(shape)) {
    for (int i = 0; i < compound->getNumChildShapes(); ++i) {
      CheckShapeCast(compound->getChildShape(i), tf0*compound->getChildTransform(i), tf1*compound->getChildTransform(i), cow, world, collisions);
    }
  }
  else {
    throw std::runtime_error("I can only continuous collision check convex shapes and compound shapes made of convex shapes");
  }

}

void BulletCollisionChecker::CastVsAll(Configuration& rad0, Configuration& rad1, const vector<KinBody::LinkPtr>& links,
    const DblVec& startjoints, const DblVec& endjoints, vector<Collision>& collisions) {
  Configuration::SaverPtr saver = rad0.Save();
  rad0.SetDOFValues(startjoints);
  int nlinks = links.size();
  vector<btTransform> tbefore(nlinks), tafter(nlinks);
  for (int i=0; i < nlinks; ++i) {
    tbefore[i] = toBt(links[i]->GetTransform());
  }
  rad1.SetDOFValues(endjoints);
  for (int i=0; i < nlinks; ++i) {
    tafter[i] = toBt(links[i]->GetTransform());
  }
  rad0.SetDOFValues(startjoints);
  UpdateBulletFromRave();
  m_world->updateAabbs();

  for (int i=0; i < nlinks; ++i) {
    assert(m_link2cow[links[i].get()] != NULL);
    CollisionObjectWrapper* cow = m_link2cow[links[i].get()];
    CheckShapeCast(cow->getCollisionShape(), tbefore[i], tafter[i], cow, m_world, collisions);
  }
  LOG_DEBUG("CastVsAll checked %li links and found %li collisions", links.size(), collisions.size());
}


void BulletCollisionChecker::CastVsCastGJKDistance(CastHullShape* shape0, const btTransform& tf0, CastHullShape* shape1, const btTransform& tf1, vector<Collision>& collisions) {
  btGjkEpaPenetrationDepthSolver epa;
  btVoronoiSimplexSolver sGjkSimplexSolver;
  btGjkPairDetector convexConvex(shape0, shape1,&sGjkSimplexSolver,&epa);

  btPointCollector gjkOutput;
  btGjkPairDetector::ClosestPointInput input;
  input.m_transformA = tf0;
  input.m_transformB = tf1;

  convexConvex.getClosestPoints(input, gjkOutput, 0);

  if (gjkOutput.m_hasResult && GetContactDistance() >= gjkOutput.m_distance) {
    Collision collision(NULL, NULL, toOR(gjkOutput.m_pointInWorld), toOR(gjkOutput.m_pointInWorld + gjkOutput.m_normalOnBInWorld*gjkOutput.m_distance),
        toOR(gjkOutput.m_normalOnBInWorld), gjkOutput.m_distance);

    btVector3 ptOn0 = gjkOutput.m_pointInWorld;
    btVector3 ptOn1 = gjkOutput.m_pointInWorld + gjkOutput.m_normalOnBInWorld*gjkOutput.m_distance;
    btVector3 normalWorldFrom0 = -gjkOutput.m_normalOnBInWorld;
    btVector3 normalWorldFrom1 = gjkOutput.m_normalOnBInWorld;
    btTransform tfWorld00 = tf0;
    btTransform tfWorld01 = tf0 * shape0->m_t01;
    btTransform tfWorld10 = tf1;
    btTransform tfWorld11 = tf1 * shape1->m_t01;
    btVector3 normalLocal00 = normalWorldFrom0 * tfWorld00.getBasis();
    btVector3 normalLocal01 = normalWorldFrom0 * tfWorld01.getBasis();
    btVector3 normalLocal10 = normalWorldFrom1 * tfWorld10.getBasis();
    btVector3 normalLocal11 = normalWorldFrom1 * tfWorld11.getBasis();

    const float SUPPORT_FUNC_TOLERANCE = .01 METERS;

    btVector3 ptLocal00;
    float localsup00;
    GetAverageSupport(shape0->m_shape, normalLocal00, localsup00, ptLocal00);
    btVector3 ptWorld00 = tfWorld00 * ptLocal00;
    btVector3 ptLocal01;
    float localsup01;
    GetAverageSupport(shape0->m_shape, normalLocal01, localsup01, ptLocal01);
    btVector3 ptWorld01 = tfWorld01 * ptLocal01;

    btVector3 ptLocal10;
    float localsup10;
    GetAverageSupport(shape1->m_shape, normalLocal10, localsup10, ptLocal10);
    btVector3 ptWorld10 = tfWorld10 * ptLocal10;
    btVector3 ptLocal11;
    float localsup11;
    GetAverageSupport(shape1->m_shape, normalLocal11, localsup11, ptLocal11);
    btVector3 ptWorld11 = tfWorld11 * ptLocal11;

    float sup00 = normalWorldFrom0.dot(ptWorld00);
    float sup01 = normalWorldFrom0.dot(ptWorld01);
    float sup10 = normalWorldFrom1.dot(ptWorld10);
    float sup11 = normalWorldFrom1.dot(ptWorld11);

    // TODO: this section is potentially problematic. think hard about the math
    if (sup00 - sup01 > SUPPORT_FUNC_TOLERANCE) {
      collision.timeA = 0;
      collision.cctypeA = CCType_Time0;
    }
    else if (sup01 - sup00 > SUPPORT_FUNC_TOLERANCE) {
      collision.timeA = 1;
      collision.cctypeA = CCType_Time1;
    }
    else {
      float l0c = (ptOn0 - ptWorld00).length(), 
            l1c = (ptOn0 - ptWorld01).length();
      collision.pt00 = toOR(ptWorld00);
      collision.pt01 = toOR(ptWorld01);
      collision.cctypeA = CCType_Between;
      const float LENGTH_TOLERANCE = .001 METERS;
      if ( l0c + l1c < LENGTH_TOLERANCE) {
        collision.timeA = .5;
      }
      else {
        collision.timeA = l0c/(l0c + l1c); 
      }
    }

    if (sup10 - sup11 > SUPPORT_FUNC_TOLERANCE) {
      collision.timeB = 0;
      collision.cctypeB = CCType_Time0;
    }
    else if (sup11 - sup10 > SUPPORT_FUNC_TOLERANCE) {
      collision.timeB = 1;
      collision.cctypeA = CCType_Time1;
    }
    else {
      float l0c = (ptOn1 - ptWorld10).length(), 
            l1c = (ptOn1 - ptWorld11).length();
      collision.pt10 = toOR(ptWorld10);
      collision.pt11 = toOR(ptWorld11);
      collision.cctypeB = CCType_Between;
      const float LENGTH_TOLERANCE = .001 METERS;
      if ( l0c + l1c < LENGTH_TOLERANCE) {
        collision.timeB = .5;
      }
      else {
        collision.timeB = l0c/(l0c + l1c); 
      }
    }

    collisions.push_back(collision);
  }
}

void createCastHullShape(btCollisionShape* shape, const btTransform& tf0,
    const btTransform& tf1, CollisionObjectWrapper* cow, vector<CollisionObjectWrapper*>& objs) {
  if (btConvexShape* convex = dynamic_cast<btConvexShape*>(shape)) {
    btTransform t01 = tf0.inverseTimes(tf1);
    CastHullShape* shape = new CastHullShape(convex, t01);
    CollisionObjectWrapper* obj = new CollisionObjectWrapper(cow->m_link);
    obj->setCollisionShape(shape);
    obj->setWorldTransform(tf0);
    obj->m_index = cow->m_index;
    objs.push_back(obj);
  }
  else if (btCompoundShape* compound = dynamic_cast<btCompoundShape*>(shape)) {
    for (int child_ind = 0; child_ind < compound->getNumChildShapes(); ++child_ind) {
      btTransform tf0_child = tf0*compound->getChildTransform(child_ind);
      btTransform tf1_child = tf1*compound->getChildTransform(child_ind);
      createCastHullShape(compound->getChildShape(child_ind), tf0_child, tf1_child, cow, objs);
    }
  }
  else {
    throw std::runtime_error("I can only create MultiCastHullShape made of convex shapes and compound shapes made of convex shapes");
  }
}

void BulletCollisionChecker::CastVsCast(Configuration& rad00, Configuration& rad01, Configuration& rad10, Configuration& rad11, const vector<KinBody::LinkPtr>& links0, const vector<KinBody::LinkPtr>& links1, const DblVec& startjoints0, const DblVec& endjoints0, const DblVec& startjoints1, const DblVec& endjoints1, vector<Collision>& collisions) {

  Configuration::SaverPtr saver = rad00.Save();
  rad00.SetDOFValues(startjoints0);
  int nlinks0 = links0.size();
  vector<btTransform> tbefore0(nlinks0), tafter0(nlinks0);
  for (int i=0; i < nlinks0; ++i) {
    tbefore0[i] = toBt(links0[i]->GetTransform());
  }
  rad01.SetDOFValues(endjoints0);
  for (int i=0; i < nlinks0; ++i) {
    tafter0[i] = toBt(links0[i]->GetTransform());
  }

  Configuration::SaverPtr saver2 = rad10.Save();
  rad10.SetDOFValues(startjoints1);
  int nlinks1 = links1.size();
  vector<btTransform> tbefore1(nlinks1), tafter1(nlinks1);
  for (int i=0; i < nlinks1; ++i) {
    tbefore1[i] = toBt(links1[i]->GetTransform());
  }
  rad11.SetDOFValues(endjoints1);
  for (int i=0; i < nlinks1; ++i) {
    tafter1[i] = toBt(links1[i]->GetTransform());
  }

  rad00.SetDOFValues(startjoints0);
  rad10.SetDOFValues(startjoints1);
  UpdateBulletFromRave();
  m_world->updateAabbs();

  for (int i = 0; i < nlinks0; ++i) {
    for (int j = 0; j < nlinks1; ++j) {
      CheckShapeCastVsCast(
          links0[i], tbefore0[i], tafter0[i],
          links1[j], tbefore1[j], tafter1[j],
          collisions);
    }
  }

  LOG_DEBUG("CastVsCast checked %i links vs %i links and found %i collisions", (int)links0.size(), (int)links1.size(), (int)collisions.size());

}

void BulletCollisionChecker::CheckShapeCastVsCast(
    KinBody::LinkPtr link0, const btTransform& tf00, const btTransform& tf01,
    KinBody::LinkPtr link1, const btTransform& tf10, const btTransform& tf11,
    vector<Collision>& collisions) {


  vector<CollisionObjectWrapper*> objs0, objs1;

  assert(m_link2cow[link0.get()] != NULL);
  CollisionObjectWrapper* cow0 = m_link2cow[link0.get()];
  createCastHullShape(cow0->getCollisionShape(), tf00, tf01, cow0, objs0);

  assert(m_link2cow[link1.get()] != NULL);
  CollisionObjectWrapper* cow1 = m_link2cow[link1.get()];
  createCastHullShape(cow1->getCollisionShape(), tf10, tf11, cow1, objs1);

  for (int i=0; i<objs0.size(); i++) {
    for (int j=0; j<objs1.size(); j++) {
      CastHullShape* shape0 = dynamic_cast<CastHullShape*>(objs0[i]->getCollisionShape());
      CastHullShape* shape1 = dynamic_cast<CastHullShape*>(objs1[j]->getCollisionShape());
      assert(!!shape0);
      assert(!!shape1);
      CastVsCastGJKDistance(shape0, objs0[i]->getWorldTransform(), shape1, objs1[j]->getWorldTransform(), collisions);
      for (int i=0; i<collisions.size(); i++) {
        collisions[i].linkA = link0.get();
        collisions[i].linkB = link1.get();
      }
    }
  }
  for (int i=0; i<objs0.size(); i++) {
    delete objs0[i]->getCollisionShape();
    delete objs0[i];
  }
  for (int i=0; i<objs0.size(); i++) {
    delete objs1[i]->getCollisionShape();
    delete objs1[i];
  }
}

}


namespace trajopt {



CollisionCheckerPtr CreateCollisionChecker(OR::EnvironmentBaseConstPtr env) {
  CollisionCheckerPtr checker(new BulletCollisionChecker(env));
  return checker;
}
}
