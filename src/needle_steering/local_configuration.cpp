#include "needle_steering.hpp"

namespace Needle {

  LocalConfiguration::LocalConfiguration(KinBodyPtr body, const Matrix4d& pose) :
    body(body), pose(pose) {}
      
  LocalConfiguration::LocalConfiguration(KinBodyPtr body) :
    body(body) {}

  void LocalConfiguration::SetDOFValues(const DblVec& dofs) {
    assert(dofs.size() == 6);
    Vector6d x;
    for (int i = 0; i < dofs.size(); ++i) {
      x[i] = dofs[i];
    }
    OpenRAVE::Transform T = matrixToTransform(pose * expUp(x));
    body->SetTransform(T);
  }

  void LocalConfiguration::GetDOFLimits(DblVec& lower, DblVec& upper) const {
    lower = DblVec(6, -INFINITY);
    upper = DblVec(6, INFINITY);
  }

  DblVec LocalConfiguration::GetDOFValues() {
    DblVec out(6);
    OpenRAVE::Transform T = body->GetTransform();
    out[0] = T.trans.x;
    out[1] = T.trans.y;
    out[2] = T.trans.z;
    OpenRAVE::Vector rot = OpenRAVE::geometry::axisAngleFromQuat(T.rot);
    out[3] = rot.x;
    out[4] = rot.y;
    out[5] = rot.z;
    return out;
  }

  int LocalConfiguration::GetDOF() const {
    return 6;
  }

  OpenRAVE::EnvironmentBasePtr LocalConfiguration::GetEnv() {
    return body->GetEnv();
  }

  DblMatrix LocalConfiguration::PositionJacobian(int link_ind, const OpenRAVE::Vector& pt) const {
    MatrixXd out(3, 6);
    out.leftCols(3) = Matrix3d::Identity();
    assert(link_ind == 0);
    KinBody::LinkPtr link = body->GetLinks()[link_ind];
    OpenRAVE::Vector dr = pt - link->GetTransform().trans;
    double matdata[9] = { 0, dr[2], -dr[1], -dr[2], 0, dr[0], dr[1], -dr[0], 0 };
    out.rightCols(3) = Eigen::Map<MatrixXd>(matdata, 3, 3);
    return out;
  }

  DblMatrix LocalConfiguration::RotationJacobian(int link_ind) const {
    PRINT_AND_THROW("not implemented");
  }

  bool LocalConfiguration::DoesAffect(const KinBody::Link& link) {
    const vector<KinBody::LinkPtr>& links = body->GetLinks();
    for (int i=0; i < links.size(); ++i) {
      if (links[i].get() == &link) return true;
    }
    return false;
  }

  std::vector<KinBody::LinkPtr> LocalConfiguration::GetAffectedLinks() {
    return body->GetLinks();
  }

  void LocalConfiguration::GetAffectedLinks(std::vector<KinBody::LinkPtr>& links,
      bool only_with_geom, vector<int>& link_inds) {
    links = GetAffectedLinks();
    link_inds.resize(links.size());
    for (int i = 0; i < links.size(); ++i)
      link_inds.push_back(links[i]->GetIndex());
  }

  DblVec LocalConfiguration::RandomDOFValues() {
    return toDblVec(Vector6d::Random());
  }

  vector<OpenRAVE::KinBodyPtr> LocalConfiguration::GetBodies() {
    return singleton(body);
  }

}
