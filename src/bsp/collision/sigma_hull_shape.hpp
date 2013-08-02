#pragma once

#include "bsp/common.hpp"

namespace BSPCollision {
  class SigmaHullShape : public btConvexShape {
  public:
    btConvexShape* m_shape;
    vector<btTransform> m_ti; // T_i = T_w_0^-1 * T_w_i
    SigmaHullShape(btConvexShape* shape, const vector<btTransform>& ti);
    virtual btVector3   localGetSupportingVertex(const btVector3& vec)const;
    void    batchedUnitVectorGetSupportingVertexWithoutMargin(const btVector3* vectors,btVector3* supportVerticesOut,int numVectors) const;
    virtual void getAabb(const btTransform& t_w0,btVector3& aabbMin,btVector3& aabbMax) const;
    virtual void getAabbSlow(const btTransform& t,btVector3& aabbMin,btVector3& aabbMax) const;
    virtual void    setLocalScaling(const btVector3& scaling);
    virtual const btVector3& getLocalScaling() const;
    virtual void    setMargin(btScalar margin);
    virtual btScalar    getMargin() const;
    virtual int     getNumPreferredPenetrationDirections() const;
    virtual void    getPreferredPenetrationDirection(int index, btVector3& penetrationVector) const;
    virtual void calculateLocalInertia(btScalar, btVector3&) const;
    virtual const char* getName() const;
    virtual btVector3 localGetSupportingVertexWithoutMargin(const btVector3& v) const;
  };
}
