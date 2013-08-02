#include "bsp/collision/utils.hpp"

namespace BSPCollision {

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
}
