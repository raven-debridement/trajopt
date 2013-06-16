#include "needle_steering.hpp"
#include "sco/expr_ops.hpp"
#include "sco/modeling_utils.hpp"
#include "sco/modeling.hpp"
#include "osgviewer/osgviewer.hpp"
#include "trajopt/collision_checker.hpp"
#include "trajopt/collision_terms.hpp"
#include "trajopt/common.hpp"
#include "trajopt/plot_callback.hpp"
#include "trajopt/problem_description.hpp"
#include "trajopt/rave_utils.hpp"
#include "trajopt/trajectory_costs.hpp"
#include "utils/clock.hpp"
#include "utils/config.hpp"
#include "utils/eigen_conversions.hpp"
#include "utils/stl_to_string.hpp"
#include <boost/assign.hpp>
#include <boost/foreach.hpp>
#include <ctime>
#include <openrave-core.h>
#include <openrave/openrave.h>

using namespace trajopt;
using namespace std;
using namespace OpenRAVE;
using namespace util;
using namespace boost::assign;
using namespace Eigen;


namespace Needle {

  template<typename T, size_t N>
  T* end(T (&ra)[N]) {
    return ra + N;
  }

  inline double bound_inf(double result, double bound) {
    return min(max(result, -bound), bound);
  }

  void AddVarArrays(OptProb& prob, int rows, const vector<int>& cols, const vector<double>& lbs, const vector<double>& ubs, const vector<string>& name_prefix, const vector<VarArray*>& newvars) {
    int n_arr = name_prefix.size();
    assert(n_arr == newvars.size());

    vector<MatrixXi> index(n_arr);
    for (int i=0; i < n_arr; ++i) {
      newvars[i]->resize(rows, cols[i]);
      index[i].resize(rows, cols[i]);
    }

    vector<string> names;
    vector<double> all_lbs;
    vector<double> all_ubs;
    int var_idx = prob.getNumVars();
    for (int i=0; i < rows; ++i) {
      for (int k=0; k < n_arr; ++k) {
        for (int j=0; j < cols[k]; ++j) {
          index[k](i,j) = var_idx;
          names.push_back( (boost::format("%s_%i_%i")%name_prefix[k]%i%j).str() );
          all_lbs.push_back(lbs[k]);
          all_ubs.push_back(ubs[k]);
          ++var_idx;
        }
      }
    }
    prob.createVariables(names, all_lbs, all_ubs); // note that w,r, are both unbounded

    const vector<Var>& vars = prob.getVars();
    for (int k=0; k < n_arr; ++k) {
      for (int i=0; i < rows; ++i) {
        for (int j=0; j < cols[k]; ++j) {
          (*newvars[k])(i,j) = vars[index[k](i,j)];
        }
      }
    }
  }

  void AddVarArrays(OptProb& prob, int rows, const vector<int>& cols, const vector<string>& name_prefix, const vector<VarArray*>& newvars) {
    vector<double> lbs(newvars.size(), -INFINITY);
    vector<double> ubs(newvars.size(), INFINITY);
    AddVarArrays(prob, rows, cols, lbs, ubs, name_prefix, newvars);
  }

  void AddVarArray(OptProb& prob, int rows, int cols, double lb, double ub, const string& name_prefix, VarArray& newvars) {
    vector<VarArray*> arrs(1, &newvars);
    vector<string> prefixes(1, name_prefix);
    vector<int> colss(1, cols);
    vector<double> lbs(1, lb);
    vector<double> ubs(1, ub);
    AddVarArrays(prob, rows, colss, lbs, ubs, prefixes, arrs);
  }

  void AddVarArray(OptProb& prob, int rows, int cols, const string& name_prefix, VarArray& newvars) {
    AddVarArray(prob, rows, cols, -INFINITY, INFINITY, name_prefix, newvars);
  }

  Matrix3d rotMat(const Vector3d& x) {
    Matrix3d out;
    out << 0, -x(2), x(1),
           x(2), 0, -x(0),
           -x(1), x(0), 0;
    return out;
  }

  Vector3d rotVec(const Matrix3d& X) {
    Vector3d out;
    out << X(2, 1), X(0, 2), X(1, 0);
    return out;
  }

  Matrix3d expA(const Vector3d& w) {
    double theta = w.norm();
    if (fabs(theta) < 1e-10) {
      return Matrix3d::Identity();
    } else {
      Matrix3d w_hat = rotMat(w);
      return Matrix3d::Identity() + w_hat / (theta*theta) * (1 - cos(theta)) + w_hat*w_hat / (theta*theta*theta) * (theta - sin(theta));
    }
  }

  Matrix3d logInvA(const Vector3d& w) {
    double theta = w.norm();
    Matrix3d w_hat = rotMat(w);
    if (fabs(theta) < 1e-8) {
      return Matrix3d::Identity();
    }
    return Matrix3d::Identity() - 0.5*w_hat + (2*sin(theta) - theta*(1 + cos(theta))) / (2 * theta*theta * sin(theta)) * w_hat*w_hat;
  }

  Matrix3d expRot(const Vector3d& x) {
    double rr = x.squaredNorm();
    if (fabs(rr) < 1e-10) {
      return Matrix3d::Identity();
    } else {
      double r = sqrt(rr);
      return rotMat(x * (sin(r) / r)) + Matrix3d::Identity() * cos(r) + (x*x.transpose()) * ((1 - cos(r)) / rr);
    }
  }

  Vector3d logRot(const Matrix3d& X) {
    // Using the old implementation since it seems more robust in practice
    Vector3d x;
    x << X(2, 1) - X(1, 2),
         X(0, 2) - X(2, 0),
         X(1, 0) - X(0, 1);
    double r = x.norm();
    double t = X(0, 0) + X(1, 1) + X(2, 2) - 1;

    if (fabs(r) < 1e-8) {
      return Vector3d::Zero();
    } else {
      return x * (atan2(r, t) / r);
    }
  }

  Matrix4d expUp(const VectorXd& x) {
    assert(x.size() == 6);
    Matrix4d X = Matrix4d::Identity();
    X.block<3, 3>(0, 0) = expRot(x.tail<3>());
    X.block<3, 1>(0, 3) = expA(x.tail<3>()) * x.head<3>();
    X(3, 3) = 1;
    return X;
  }

  VectorXd logDown(const Matrix4d& X) {
    VectorXd x(6);
    x.tail<3>() = logRot(X.block<3, 3>(0, 0));
    x.head<3>() = (expA(x.tail<3>())).inverse() * X.block<3, 1>(0, 3);
    return x;
  }

  OpenRAVE::Transform matrixToTransform(const Matrix4d& X) {
    OpenRAVE::TransformMatrix M;
    M.trans.x = X(0, 3);
    M.trans.y = X(1, 3);
    M.trans.z = X(2, 3);
    for (int row = 0; row < 3; ++row) {
      for (int col = 0; col < 3; ++col) {
        M.m[row*4+col] = X(row, col);
      }
    }
    return OpenRAVE::Transform(M);
  }

  OpenRAVE::Transform vecToTransform(const VectorXd& x) {
    OpenRAVE::Transform T;
    OpenRAVE::Vector trans(x[0], x[1], x[2]);
    OpenRAVE::Vector rot(x[3], x[4], x[5]);
    T.trans = trans;
    T.rot = OpenRAVE::geometry::quatFromAxisAngle(rot);
    return T;
  }

  SpeedCost::SpeedCost(const Var& var, double coeff) : Cost("Speed"), var_(var), coeff_(coeff) {
    exprInc(expr_, exprMult(var, coeff));
  }
  double SpeedCost::value(const vector<double>& xvec) {
    double speed = getVec(xvec, singleton<Var>(var_))[0];
    return speed * coeff_;
  }
  ConvexObjectivePtr SpeedCost::convex(const vector<double>& xvec, Model* model) {
    ConvexObjectivePtr out(new ConvexObjective(model));
    out->addAffExpr(expr_);
    return out;
  }

  RotationCost::RotationCost(const VarVector& vars, double coeff) : Cost("Rotation"), vars_(vars), coeff_(coeff) {
    for (int i = 0; i < vars.size(); ++i) {
      exprInc(expr_, exprMult(exprSquare(vars[i]), coeff));
    }
  }
  double RotationCost::value(const vector<double>& xvec) {
    VectorXd vals = getVec(xvec, vars_);
    return vals.array().square().sum() * coeff_;
  }
  ConvexObjectivePtr RotationCost::convex(const vector<double>& xvec, Model* model) {
    ConvexObjectivePtr out(new ConvexObjective(model));
    out->addQuadExpr(expr_);
    return out;
  }

  LocalConfiguration::LocalConfiguration(KinBodyPtr body, const Matrix4d& pose) :
    body_(body), pose_(pose) {}
      
  LocalConfiguration::LocalConfiguration(KinBodyPtr body) :
    body_(body) {}

  void LocalConfiguration::SetDOFValues(const DblVec& dofs) {
    VectorXd x(dofs.size());
    for (int i = 0; i < dofs.size(); ++i) {
      x[i] = dofs[i];
    }
    OpenRAVE::Transform T = matrixToTransform(pose_ * expUp(x));
    body_->SetTransform(T);
  }

  void LocalConfiguration::GetDOFLimits(DblVec& lower, DblVec& upper) const {
    lower = DblVec(6, -INFINITY);
    upper = DblVec(6, INFINITY);
  }

  DblVec LocalConfiguration::GetDOFValues() {
    DblVec out(6);
    OpenRAVE::Transform T = body_->GetTransform();
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
    return body_->GetEnv();
  }

  DblMatrix LocalConfiguration::PositionJacobian(int link_ind, const OpenRAVE::Vector& pt) const {
    MatrixXd out(3, 6);
    out.leftCols(3) = Matrix3d::Identity();
    assert(link_ind == 0);
    KinBody::LinkPtr link = body_->GetLinks()[link_ind];
    OpenRAVE::Vector dr = pt - link->GetTransform().trans;
    double matdata[9] = { 0, dr[2], -dr[1], -dr[2], 0, dr[0], dr[1], -dr[0], 0 };
    out.rightCols(3) = Eigen::Map<MatrixXd>(matdata, 3, 3);
    return out;
  }

  DblMatrix LocalConfiguration::RotationJacobian(int link_ind) const {
    PRINT_AND_THROW("not implemented");
  }

  bool LocalConfiguration::DoesAffect(const KinBody::Link& link) {
    const vector<KinBody::LinkPtr>& links = body_->GetLinks();
    for (int i=0; i < links.size(); ++i) {
      if (links[i].get() == &link) return true;
    }
    return false;
  }

  std::vector<KinBody::LinkPtr> LocalConfiguration::GetAffectedLinks() {
    return body_->GetLinks();
  }

  void LocalConfiguration::GetAffectedLinks(std::vector<KinBody::LinkPtr>& links,
      bool only_with_geom, vector<int>& link_inds) {
    links = GetAffectedLinks();
    link_inds.resize(links.size());
    for (int i = 0; i < links.size(); ++i)
      link_inds.push_back(links[i]->GetIndex());
  }

  DblVec LocalConfiguration::RandomDOFValues() {
    return toDblVec(VectorXd::Random(6));
  }

  vector<OpenRAVE::KinBodyPtr> LocalConfiguration::GetBodies() {
    return singleton(body_);
  }

  PositionError::PositionError(LocalConfigurationPtr cfg, const VectorXd& target_pos) : cfg(cfg), target_pos(target_pos), body(cfg->GetBodies()[0]) {}

  VectorXd PositionError::operator()(const VectorXd& a) const {
    Matrix4d X = cfg->pose_ * expUp(a);
    VectorXd x(6);
    x.head<3>() = X.block<3, 1>(0, 3) - target_pos.head<3>();
    x.tail<3>() = Vector3d::Zero();
    return logDown((cfg->pose_ * expUp(a)).inverse() * expUp(target_pos));
  }

  ControlError::ControlError(LocalConfigurationPtr cfg0, LocalConfigurationPtr cfg1, NeedleProblemHelperPtr helper) : cfg0(cfg0), cfg1(cfg1), body(cfg0->GetBodies()[0]), helper(helper) {}

  VectorXd ControlError::operator()(const VectorXd& a) const {
    Matrix4d pose1 = cfg0->pose_ * expUp(a.topRows(6));
    Matrix4d pose2 = cfg1->pose_ * expUp(a.middleRows(6,6));
    double phi = a(12), Delta = a(13);
    double curvature;
    double radius;
    switch (helper->curvature_formulation) {
      case NeedleProblemHelper::UseCurvature:
        switch (helper->curvature_constraint) {
          case NeedleProblemHelper::ConstantRadius:
            curvature = 1.0 / helper->r_min;
            break;
          case NeedleProblemHelper::BoundedRadius:
            curvature = a(14);
            break;
          SWITCH_DEFAULT;
        }
        break;
      case NeedleProblemHelper::UseRadius:
        switch (helper->curvature_constraint) {
          case NeedleProblemHelper::ConstantRadius:
            radius = helper->r_min;
            break;
          case NeedleProblemHelper::BoundedRadius:
            radius = a(14);
            break;
          SWITCH_DEFAULT;
        }
        break;
      SWITCH_DEFAULT;
    }
    switch (helper->formulation) {
      case NeedleProblemHelper::Form1:
      case NeedleProblemHelper::Form3: {
        switch (helper->curvature_constraint) {
          case NeedleProblemHelper::UseCurvature:
            return logDown(helper->TransformPose(pose1, phi, Delta, curvature).inverse() * pose2);
          case NeedleProblemHelper::UseRadius:
            return logDown(helper->TransformPose(pose1, phi, Delta, radius).inverse() * pose2);
          SWITCH_DEFAULT;
        }
      }
      case NeedleProblemHelper::Form2: {
        VectorXd trans1 = pose1.block<3, 1>(0, 3);
        VectorXd trans2 = pose2.block<3, 1>(0, 3);
        Matrix3d R1 = pose1.block<3, 3>(0, 0);
        Vector3d xyz = R1.transpose() * (trans2 - trans1);
        double x = xyz[0], y = xyz[1], z = xyz[2];
        Vector3d err;// err <<
        //  phi - atan2(x, -y),
        //  bound_inf(radius - (x*x + y*y + z*z) / (2 * sqrt(x*x + y*y)), 100), // sketchy lol
        //  Delta - radius * atan2(z, radius - sqrt(x*x + y*y));
        return err;
      }
      SWITCH_DEFAULT;
    }
  }

  int ControlError::outputSize() const {
    switch (helper->formulation) {
      case NeedleProblemHelper::Form1:
        return 6;
      case NeedleProblemHelper::Form2:
        return 3;
      case NeedleProblemHelper::Form3:
        return 6;
      SWITCH_DEFAULT;
    }
  }

  void NeedleProblemHelper::ConfigureProblem(OptProb& prob) {
    CreateVariables(prob);
    InitLocalConfigurations(robot, prob);
    InitTrajectory(prob);
    prob.addCost(CostPtr(new RotationCost(phivars.col(0), coeff_rotation)));
    prob.addCost(CostPtr(new SpeedCost(Deltavar, coeff_speed)));
    AddStartConstraint(prob);
    AddGoalConstraint(prob);
    AddControlConstraint(prob);
    AddCollisionConstraint(prob);
  }

  void NeedleProblemHelper::InitOptimizeVariables(BasicTrustRegionSQP& opt) {
    DblVec initVec;
    // Initialize twistvars
    for (int i = 0; i <= T; ++i) {
      for (int j = 0; j < n_dof; ++j) {
        initVec.push_back(0.);
      }
    }
    // Initialize phivars
    for (int i = 0; i < T; ++i) {
      initVec.push_back(0.);
    }
    // Initialize Delta
    initVec.push_back(Delta_lb);
    // Initialize time frame radii
    if (curvature_constraint == BoundedRadius) {
      for (int i = 0; i < T; ++i) {
        switch (curvature_formulation) {
          case UseCurvature:
            initVec.push_back(1.0 / r_min);
            break;
          case UseRadius:
            initVec.push_back(r_min);
            break;
          SWITCH_DEFAULT;
        }
      }
    }
    opt.initialize(initVec);
  }

  Matrix4d NeedleProblemHelper::TransformPose(const Matrix4d& pose, double phi, double Delta, double curvature_or_radius) const {
    double theta;
    switch (curvature_formulation) {
      case NeedleProblemHelper::UseCurvature:
        theta = Delta * curvature_or_radius;
        break;
      case NeedleProblemHelper::UseRadius:
        theta = Delta / curvature_or_radius;
        break;
      SWITCH_DEFAULT;
    }
    switch (formulation) {
      case NeedleProblemHelper::Form1:
      case NeedleProblemHelper::Form2: {
        VectorXd w(6); w << 0, 0, 0, 0, 0, phi;
        VectorXd v(6); v << 0, 0, Delta, theta, 0, 0;
        return pose * expUp(w) * expUp(v);
      }
      case NeedleProblemHelper::Form3: {
        VectorXd w(6); w << 0, 0, Delta, theta, 0, phi;
        return pose * expUp(w);
      }
      SWITCH_DEFAULT;
    }
  }

  double NeedleProblemHelper::GetPhi(const DblVec& x, int i) const {
    return x[phivars.row(i)[0].var_rep->index];
  }

  double NeedleProblemHelper::GetDelta(const DblVec& x, int i) const {
    return x[Deltavar.var_rep->index];
  }
  
  double NeedleProblemHelper::GetCurvature(const DblVec& x, int i) const {
    assert (curvature_formulation == UseCurvature);
    switch (curvature_constraint) {
      case ConstantRadius:
        return 1.0 / r_min;
      case BoundedRadius:
        return x[curvaturevars.row(i)[0].var_rep->index];
      SWITCH_DEFAULT;
    }
  }

  double NeedleProblemHelper::GetRadius(const DblVec& x, int i) const {
    assert (curvature_formulation == UseRadius);
    switch (curvature_constraint) {
      case ConstantRadius:
        return r_min;
      case BoundedRadius:
        return x[radiusvars.row(i)[0].var_rep->index];
      SWITCH_DEFAULT;
    }
  }

  #ifdef NEEDLE_TEST
  void NeedleProblemHelper::checkAlignment(DblVec& x) {
    double diff = 0.;
    for (int i = 0; i < T; ++i) {
      double phi = GetPhi(x, i);
      double Delta = GetDelta(x, i);
      switch (curvature_formulation) {
        case UseCurvature:
          double curvature = GetCurvature(x, i);
          diff += (local_configs[i+1]->pose_ - TransformPose(local_configs[i]->pose_, phi, Delta, curvature)).norm();
          break;
        case UseRadius:
          double radius = GetRadius(x, i);
          diff += (local_configs[i+1]->pose_ - TransformPose(local_configs[i]->pose_, phi, Delta, radius)).norm();
          break;
        SWITCH_DEFAULT;
      }
    }
  }
  #endif

  void NeedleProblemHelper::OptimizerCallback(OptProb*, DblVec& x) {
    
    switch (method) {
      case Colocation: {
        MatrixXd twistvals = getTraj(x, twistvars);
        for (int i = 0; i < local_configs.size(); ++i) {
          local_configs[i]->pose_ = local_configs[i]->pose_ * expUp(twistvals.row(i));
        }
        #ifdef NEEDLE_TEST
        checkAlignment(x);
        #endif
        setVec(x, twistvars.m_data, DblVec(twistvars.size(), 0));
        break;
      }
      case Shooting: {
        // execute the control input to set local configuration poses
        local_configs[0]->pose_ = expUp(start);
        for (int i = 0; i < T; ++i) {
          double phi = GetPhi(x, i);
          double Delta = GetDelta(x, i);
          switch (curvature_formulation) {
            case UseCurvature: {
              double curvature = GetCurvature(x, i);
              local_configs[i+1]->pose_ = TransformPose(local_configs[i]->pose_, phi, Delta, curvature);
              break;
            }
            case UseRadius: {
              double radius = GetRadius(x, i);
              local_configs[i+1]->pose_ = TransformPose(local_configs[i]->pose_, phi, Delta, radius);
              break;
            }
            SWITCH_DEFAULT;
          }
        }
        setVec(x, twistvars.m_data, DblVec(twistvars.size(), 0));
        break;
      }
      SWITCH_DEFAULT;
    }
  }

  void NeedleProblemHelper::ConfigureOptimizer(BasicTrustRegionSQP& opt) {
    InitOptimizeVariables(opt);
    opt.addCallback(boost::bind(&Needle::NeedleProblemHelper::OptimizerCallback, this, _1, _2));
  }

  void NeedleProblemHelper::CreateVariables(OptProb& prob) {
    // Time frame varies from 0 to T instead of from 0 to T-1
    AddVarArray(prob, T+1, n_dof, "twist", twistvars);
    AddVarArray(prob, T, 1, -PI, PI, "phi", phivars);
    Delta_lb = (goal.topRows(3) - start.topRows(3)).norm() / T / r_min;
    Deltavar = prob.createVariables(singleton<string>("Delta"), singleton<double>(Delta_lb),singleton<double>(INFINITY))[0];
    // Only the twist variables are incremental (i.e. their trust regions should be around zero)
    prob.setIncremental(twistvars.flatten());
    if (curvature_constraint == BoundedRadius) {
      //AddVarArray(prob, T, 1, r_min, INFINITY, "radius", radiusvars);
      switch (curvature_formulation) {
        case UseCurvature:
          AddVarArray(prob, T, 1, 0.01, 1. / r_min, "curvature", curvaturevars);
          break;
        case UseRadius:
          AddVarArray(prob, T, 1, r_min, 100.0, "radius", radiusvars);
          break;
        SWITCH_DEFAULT;
      }
    }
  }

  void NeedleProblemHelper::InitLocalConfigurations(const KinBodyPtr robot, OptProb& prob) {
    for (int i = 0; i <= T; ++i) {
      local_configs.push_back(LocalConfigurationPtr(new LocalConfiguration(robot)));
    }
  }

  void NeedleProblemHelper::InitTrajectory(OptProb& prob) {
    MatrixXd initTraj(T+1, n_dof);
    for (int idof = 0; idof < n_dof; ++idof) {
      initTraj.col(idof) = VectorXd::LinSpaced(T+1, start[idof], goal[idof]);
    }
    for (int i = 0; i <= T; ++i) {
      local_configs[i]->pose_ = expUp(initTraj.row(i));
    }
  }

  void NeedleProblemHelper::AddStartConstraint(OptProb& prob) {
    VarVector vars = twistvars.row(0);
    VectorOfVectorPtr f(new Needle::PositionError(local_configs[0], start));
    VectorXd coeffs = VectorXd::Ones(6);
    prob.addConstraint(ConstraintPtr(new ConstraintFromFunc(f, vars, coeffs, EQ, "entry")));
  }

  void NeedleProblemHelper::AddGoalConstraint(OptProb& prob) {
    VarVector vars = twistvars.row(T);
    VectorOfVectorPtr f(new Needle::PositionError(local_configs[T], goal));
    VectorXd coeffs(n_dof); coeffs << 1., 1., 1., 0., 0., 0.;
    prob.addConstraint(ConstraintPtr(new ConstraintFromFunc(f, vars, coeffs, EQ, "goal")));
  }

  void NeedleProblemHelper::AddControlConstraint(OptProb& prob) {
    for (int i = 0; i < T; ++i) {
      VarVector vars = concat(concat(twistvars.row(i), twistvars.row(i+1)), phivars.row(i));
      vars.push_back(Deltavar);
      if (curvature_constraint == BoundedRadius) {
        switch (curvature_formulation) {
          case UseCurvature:
            vars = concat(vars, curvaturevars.row(i));
            break;
          case UseRadius:
            vars = concat(vars, radiusvars.row(i));
            break;
          SWITCH_DEFAULT;
        }
      }
      VectorOfVectorPtr f(new Needle::ControlError(local_configs[i], local_configs[i+1], NeedleProblemHelperPtr(this)));
      VectorXd coeffs = VectorXd::Ones(boost::static_pointer_cast<Needle::ControlError>(f)->outputSize());
      prob.addConstraint(ConstraintPtr(new ConstraintFromFunc(f, vars, coeffs, EQ, (boost::format("control%i")%i).str())));
    }
  }

  void NeedleProblemHelper::AddCollisionConstraint(OptProb& prob) {
    //Str2Dbl tag2dist_pen(collision_dist_pen), tag2coeff(collision_coeff);
    //for (int i = 0; i < ignored_kinbody_names.size(); ++i) {
    //  tag2coeff.insert( std::pair<string, double>(ignored_kinbody_names[i], 0.0) );
    //}
    //for (int i = 0; i <= T; ++i) {
    //  prob.addConstraint(ConstraintPtr(new CollisionTaggedConstraint(tag2dist_pen, tag2coeff, local_configs[i], twistvars.row(i))));
    //}
  }

  TrajPlotter::TrajPlotter(const vector<LocalConfigurationPtr>& local_configs, const VarArray& vars) : local_configs(local_configs), vars(vars) {
    viewer = OSGViewer::GetOrCreate(local_configs[0]->GetEnv());
  }

  void TrajPlotter::OptimizerCallback(OptProb*, DblVec& x, NeedleProblemHelperPtr helper) {
    vector<GraphHandlePtr> handles;
    vector<KinBodyPtr> bodies = local_configs[0]->GetBodies();
    MatrixXd vals = getTraj(x, vars);
    for (int i=0; i < vals.rows(); ++i) {
      local_configs[i]->SetDOFValues(toDblVec(vals.row(i)));
      BOOST_FOREACH(const KinBodyPtr& body, bodies) {
        handles.push_back(viewer->PlotKinBody(body));
        SetTransparency(handles.back(), .35);
      }
    }
    handles.push_back(viewer->PlotKinBody(GetGoalKinBody(helper))); 
    SetTransparency(handles.back(), 1);
    viewer->Idle();
  }

  KinBodyPtr TrajPlotter::GetGoalKinBody(NeedleProblemHelperPtr helper) {
    OpenRAVE::Transform T;
    T.trans.x = helper->goal[0];
    T.trans.y = helper->goal[1];
    T.trans.z = helper->goal[2];
    OpenRAVE::Vector rot(helper->goal[3], helper->goal[4], helper->goal[5]);
    T.rot = OpenRAVE::geometry::quatFromAxisAngle(rot);
    helper->robot->SetTransform(T);
    return helper->robot;
  }

  void TrajPlotter::PlotBothTrajectories(OptProbPtr prob, const BasicTrustRegionSQP& opt, NeedleProblemHelperPtr helper) {
    DblVec x = prob->getModel()->getVarValues(prob->getModel()->getVars());
    vector<GraphHandlePtr> handles;
    KinBodyPtr robot = helper->robot;
    // plot real trajectory
    Matrix4d current_pose = expUp(helper->start);
    robot->SetTransform(Needle::matrixToTransform(current_pose));
    handles.push_back(viewer->PlotKinBody(robot));
    SetTransparency(handles.back(), .5);
    for (int i = 0; i < helper->T; ++i) {
      double phi = helper->GetPhi(x, i);
      double Delta = helper->GetDelta(x, i);
      switch (helper->curvature_formulation) {
        case NeedleProblemHelper::UseCurvature: {
          double curvature = helper->GetCurvature(x, i);
          current_pose = helper->TransformPose(current_pose, phi, Delta, curvature);
          break;
        }
        case NeedleProblemHelper::UseRadius: {
          double radius = helper->GetRadius(x, i);
          current_pose = helper->TransformPose(current_pose, phi, Delta, radius);
          break;
        }
        SWITCH_DEFAULT;
      }
      robot->SetTransform(Needle::matrixToTransform(current_pose));
      handles.push_back(viewer->PlotKinBody(robot));
      SetTransparency(handles.back(), .5);
    }
    // plot ideal trajectory
    MatrixXd vals = getTraj(x, vars);
    for (int i=0; i < vals.rows(); ++i) {
      local_configs[i]->SetDOFValues(toDblVec(vals.row(i)));
      handles.push_back(viewer->PlotKinBody(robot));
      SetTransparency(handles.back(), .35);
    }
    viewer->Idle();
  }
}

void printVector(VectorXd x) {
  for (int i = 0; i < x.size(); ++i) {
    cout << x(i);
    if (i < x.size() - 1) cout << ", ";
  }
}

int main(int argc, char** argv)
{

  bool plotting=false, verbose=false;
  bool plot_final_result=false;
  double env_transparency = 0.5;

  int T = 10;
  double r_min = 3.5;//2;
  int n_dof = 6;

  int formulation = Needle::NeedleProblemHelper::Form1;
  int curvature_constraint = Needle::NeedleProblemHelper::ConstantRadius;
  int method = Needle::NeedleProblemHelper::Shooting;
  int curvature_formulation = Needle::NeedleProblemHelper::UseRadius;

  double improve_ratio_threshold = 0.1;//0.25;
  double trust_shrink_ratio = 0.9;//0.7;
  double trust_expand_ratio = 1.3;//1.2;
  
  double start_vec_array[] = {-12.82092, 6.80976, 0.06844, 0, 0, 0};
  double goal_vec_array[] = {-3.21932, 6.87362, -1.21877, 0, 0, 0};

  vector<double> start_vec(start_vec_array, start_vec_array + n_dof);
  vector<double> goal_vec(goal_vec_array, goal_vec_array + n_dof);

  {
    Config config;
    config.add(new Parameter<bool>("plotting", &plotting, "plotting"));
    config.add(new Parameter<bool>("plot_final_result", &plot_final_result, "plot_final_result"));
    config.add(new Parameter<bool>("verbose", &verbose, "verbose"));
    config.add(new Parameter<double>("env_transparency", &env_transparency, "env_transparency"));
    config.add(new Parameter<int>("T", &T, "T"));
    config.add(new Parameter<int>("formulation", &formulation, "formulation"));
    config.add(new Parameter<int>("curvature_constraint", &curvature_constraint, "curvature_constraint"));
    config.add(new Parameter<int>("method", &method, "method"));
    config.add(new Parameter<int>("curvature_formulation", &curvature_formulation, "curvature_formulation"));
    config.add(new Parameter<double>("r_min", &r_min, "r_min"));
    config.add(new Parameter<double>("improve_ratio_threshold", &improve_ratio_threshold, "improve_ratio_threshold"));
    config.add(new Parameter<double>("trust_shrink_ratio", &trust_shrink_ratio, "trust_shrink_ratio"));
    config.add(new Parameter<double>("trust_expand_ratio", &trust_expand_ratio, "trust_expand_ratio"));
    config.add(new Parameter< vector<double> >("s", &start_vec, "s"));
    config.add(new Parameter< vector<double> >("g", &goal_vec, "g"));
    CommandParser parser(config);
    parser.read(argc, argv);
  }

  double coeff_rotation = 1.;
  double coeff_speed = 1.;

  RaveInitialize(false, verbose ? Level_Debug : Level_Info);
  EnvironmentBasePtr env = RaveCreateEnvironment();
  env->StopSimulation();
  OSGViewerPtr viewer = OSGViewer::GetOrCreate(env);
  assert(viewer);

  env->Load(string(DATA_DIR) + "/prostate.env.xml");//needleprob.env.xml");
  viewer->SetAllTransparency(env_transparency);
  RobotBasePtr robot = GetRobot(*env);

  VectorXd start(n_dof); for (int i = 0; i < n_dof; ++i) start[i] = start_vec[i];
  VectorXd goal(n_dof); for (int i = 0; i < n_dof; ++i) goal[i] = goal_vec[i];

  const char *ignored_kinbody_c_strs[] = { "KinBodyProstate", "KinBodyDermis", "KinBodyEpidermis", "KinBodyHypodermis" };
  vector<string> ignored_kinbody_names(ignored_kinbody_c_strs, Needle::end(ignored_kinbody_c_strs));

  OptProbPtr prob(new OptProb());

  Needle::NeedleProblemHelperPtr helper(new Needle::NeedleProblemHelper());
  helper->start = start;
  helper->goal = goal;
  helper->coeff_rotation = coeff_rotation;
  helper->coeff_speed = coeff_speed;
  helper->T = T;
  helper->r_min = r_min;
  helper->n_dof = n_dof;
  helper->ignored_kinbody_names = ignored_kinbody_names;
  helper->collision_dist_pen = 0.025;
  helper->collision_coeff = 20;
  helper->formulation = formulation;
  helper->curvature_constraint = curvature_constraint;
  helper->method = method;
  helper->curvature_formulation = curvature_formulation;
  helper->robot = robot;
  helper->ConfigureProblem(*prob);

  BasicTrustRegionSQP opt(prob);
  opt.max_iter_ = 500;    
  opt.improve_ratio_threshold_ = improve_ratio_threshold;
  opt.trust_shrink_ratio_ = trust_shrink_ratio;
  opt.trust_expand_ratio_ = trust_expand_ratio;

  helper->ConfigureOptimizer(opt);

  boost::shared_ptr<Needle::TrajPlotter> plotter;
  plotter.reset(new Needle::TrajPlotter(helper->local_configs, helper->twistvars));
  if (plotting) {
    opt.addCallback(boost::bind(&Needle::TrajPlotter::OptimizerCallback, boost::ref(plotter), _1, _2, helper));
  }

  opt.optimize();
  
  if (plot_final_result) plotter->PlotBothTrajectories(prob, opt, helper);

  RaveDestroy();


}
