#pragma once

#include "bsp/common.hpp"

namespace BSPCollision {

  btVector3 toBt(const OR::Vector& v);

  OR::Vector toOR(const btVector3& v);

  btQuaternion toBtQuat(const OR::Vector& q);

  btTransform toBt(const OR::Transform& t);

  Matrix4d toMatrix(const btTransform& b);

  Vector3d toVector(const btVector3& v);

  btVector3 barycentricCoordinates(const btVector3& a, const btVector3& b, const btVector3& c, const btVector3& p);

  void computeSupportingWeights(const vector<btVector3>& v, const btVector3& p, vector<float>& alpha);

  void compute_points_and_supports(const btConvexShape* shape, const vector<btTransform>& ts, const btVector3& normalWorldFromShape, const CollisionObjectWrapper* cow, vector<float>* output_sup, vector<btVector3>* output_ptWorld);

  void compute_max_support_points(const vector<float>& sup, const vector<btVector3>& ptWorld, vector<float>* output_sups, vector<btVector3>* output_max_ptWorlds, vector<int>* output_instance_inds);
}
