#include "bsp/collision/sigma_hull_shape.hpp"

namespace BSPCollision {

  SigmaHullShape::SigmaHullShape(btConvexShape* shape, const vector<btTransform>& ti) : m_shape(shape), m_ti(ti) {
    m_shapeType = CUSTOM_CONVEX_SHAPE_TYPE;
  }

  btVector3 SigmaHullShape::localGetSupportingVertex(const btVector3& vec)const {
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
  void SigmaHullShape::batchedUnitVectorGetSupportingVertexWithoutMargin(const btVector3* vectors,btVector3* supportVerticesOut,int numVectors) const {
    throw std::runtime_error("not implemented");
  }

  ///getAabb's default implementation is brute force, expected derived classes to implement a fast dedicated version
  void SigmaHullShape::getAabb(const btTransform& t_w0,btVector3& aabbMin,btVector3& aabbMax) const {
    m_shape->getAabb(t_w0, aabbMin, aabbMax);
    btVector3 min_i, max_i;
    for (int i=0; i<m_ti.size(); i++) {
      m_shape->getAabb(t_w0*m_ti[i], min_i, max_i );
      aabbMin.setMin(min_i);
      aabbMax.setMax(max_i);
    }
  }

  void SigmaHullShape::getAabbSlow(const btTransform& t,btVector3& aabbMin,btVector3& aabbMax) const {
    throw std::runtime_error("shouldn't happen");
  }

  void SigmaHullShape::setLocalScaling(const btVector3& scaling) {}

  const btVector3& SigmaHullShape::getLocalScaling() const {
    static btVector3 out(1,1,1);
    return out;
  }

  void SigmaHullShape::setMargin(btScalar margin) {}

  btScalar SigmaHullShape::getMargin() const {return 0;}

  int SigmaHullShape::getNumPreferredPenetrationDirections() const {return 0;}

  void SigmaHullShape::getPreferredPenetrationDirection(int index, btVector3& penetrationVector) const {throw std::runtime_error("not implemented");}

  void SigmaHullShape::calculateLocalInertia(btScalar, btVector3&) const {throw std::runtime_error("not implemented");}

  const char* SigmaHullShape::getName() const {return "SigmaHullShape";}

  btVector3 SigmaHullShape::localGetSupportingVertexWithoutMargin(const btVector3& v) const {return localGetSupportingVertex(v);}

}
