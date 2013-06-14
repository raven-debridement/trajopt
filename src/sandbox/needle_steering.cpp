#include "o3.hpp"
#include "sco/expr_ops.hpp"
#include "sco/modeling_utils.hpp"
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
  if (theta < 1e-8) {
    return Matrix3d::Identity();
  }
  return Matrix3d::Identity() - 0.5*w_hat + (2*sin(theta) - theta*(1 + cos(theta))) / (2 * theta*theta * sin(theta)) * w_hat*w_hat;
}

Matrix3d expRot(const Vector3d& x) {
  double theta = x.norm();
  if (fabs(theta) < 1e-8) {
    return Matrix3d::Identity();
  } else {
    Vector3d w = x / theta; 
    Matrix3d w_hat = rotMat(w);
    return Matrix3d::Identity() + w_hat * sin(theta) + w_hat*w_hat * (1 - cos(theta));
  }
}

Vector3d logRot(const Matrix3d& X) {
  //double theta = acos(0.5 * (X.trace() - 1));
  //if (fabs(sin(theta)) < 1e-8) {
  //  return Vector3d::Zero();
  //}
  //Matrix3d w_hat = (X - X.transpose()) / (2 * sin(theta));
  //return theta * rotVec(w_hat);

  // Using the old implementation since it seems more robust in practice
  Vector3d x;
  x << X(2, 1) - X(1, 2),
       X(0, 2) - X(2, 0),
       X(1, 0) - X(0, 1);
  double r = x.norm();
  double t = X(0, 0) + X(1, 1) + X(2, 2) - 1;

  if (r == 0) {
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
  x.head<3>() = logInvA(x.tail<3>()) * X.block<3, 1>(0, 3);
  return x;
}

void printVector(VectorXd x) {
  for (int i = 0; i < x.size(); ++i) {
    cout << x(i);
    if (i < x.size() - 1) cout << ", ";
  }
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

struct NeedleError : public VectorOfVector {
  ConfigurationPtr cfg0, cfg1;
  double radius;
  KinBodyPtr body;
  NeedleError(ConfigurationPtr cfg0, ConfigurationPtr cfg1, double radius) : cfg0(cfg0), cfg1(cfg1), radius(radius), body(cfg0->GetBodies()[0]) {}
  VectorXd operator()(const VectorXd& a) const {
    cfg0->SetDOFValues(toDblVec(a.topRows(6)));
    OR::Transform Tw0 = body->GetTransform();
    cfg1->SetDOFValues(toDblVec(a.middleRows(6,6)));
    OR::Transform Tw1 = body->GetTransform();
    double theta = a(12);
    OR::Transform Ttarg0(OR::geometry::quatFromAxisAngle(OR::Vector(0,0,1), theta),
        OR::Vector(radius*sin(theta), radius*(1-cos(theta)),0));
        
    OR::Transform Ttargw = Tw0 * Ttarg0;
    OR::Vector position_errA = Ttargw.trans - Tw1.trans;
    OR::Vector ori_err = Ttargw * OR::Vector(1,0,0) - Tw1 * OR::Vector(1,0,0);
    return concat(toVector3d(position_errA), toVector3d(ori_err));

  }
};




struct TrajPlotter {
  vector<IncrementalRBPtr> rbs;
  VarArray vars;
  OSGViewerPtr viewer;
  TrajPlotter(const vector<IncrementalRBPtr>& rbs, const VarArray& vars);
  void OptimizerCallback(OptProb*, DblVec& x);
  void PlotBothTrajectories(OptProbPtr prob, const BasicTrustRegionSQP& opt, const VectorXd& start);
};

TrajPlotter::TrajPlotter(const vector<IncrementalRBPtr>& rbs, const VarArray& vars) : rbs(rbs), vars(vars) {
  viewer = OSGViewer::GetOrCreate(rbs[0]->GetEnv());
}

void TrajPlotter::OptimizerCallback(OptProb*, DblVec& x) {
  vector<GraphHandlePtr> handles;
  vector<KinBodyPtr> bodies = rbs[0]->GetBodies();
  MatrixXd traj = getTraj(x,vars);
  for (int i=0; i < traj.rows(); ++i) {
    rbs[i]->SetDOFValues(toDblVec(traj.row(i)));
    BOOST_FOREACH(const KinBodyPtr& body, bodies) {
      handles.push_back(viewer->PlotKinBody(body));
      SetTransparency(handles.back(), .35);
    }
  }
  viewer->Idle();
}  

void TrajPlotter::PlotBothTrajectories(OptProbPtr prob, const BasicTrustRegionSQP& opt, const VectorXd& start) {
  DblVec x = prob->getModel()->getVarValues(prob->getModel()->getVars());
  MatrixXd traj = getTraj(x,vars);
  vector<GraphHandlePtr> handles;
  KinBodyPtr robot = rbs[0]->GetBodies()[0];
  // plot real trajectory
  vector<Matrix4d> poses;
  vector<VectorXd> twists;
  for (int i = 0; i < rbs.size(); ++i) {
    rbs[i]->SetDOFValues(toDblVec(traj.row(i)));
    OR::TransformMatrix M(rbs[i]->m_body->GetTransform());
    Matrix4d pose;
    pose << M.m[0], M.m[1], M.m[2], M.trans[0],
            M.m[4], M.m[5], M.m[6], M.trans[1],
            M.m[8], M.m[9], M.m[10], M.trans[2],
            0,      0,      0,                1;
    poses.push_back(pose);
  }
  for (int i = 0; i < poses.size() - 1; ++i) {
    twists.push_back(logDown(poses[i].inverse() * poses[i+1]));
  }
  Matrix4d current_pose = expUp(start);
  robot->SetTransform(matrixToTransform(current_pose));
  handles.push_back(viewer->PlotKinBody(robot));
  SetTransparency(handles.back(), .5);
  for (int i = 0; i < twists.size(); ++i) {
    twists[i][1] = twists[i][2] = twists[i][4] = 0;
    current_pose = current_pose * expUp(twists[i]);
    robot->SetTransform(matrixToTransform(current_pose));
    handles.push_back(viewer->PlotKinBody(robot));
    SetTransparency(handles.back(), .5);
  }
  // plot ideal trajectory
  MatrixXd vals = getTraj(x, vars);
  for (int i=0; i < vals.rows(); ++i) {
    rbs[i]->SetDOFValues(toDblVec(vals.row(i)));
    handles.push_back(viewer->PlotKinBody(robot));
    SetTransparency(handles.back(), .35);
  }
  viewer->Idle();
}




int main(int argc, char** argv)
{
  bool plotting=false, verbose=false;
  double env_transparency = 0.5;
  int n_steps = 25;
  double turning_radius = 2.85; // turning radius for needle
  int n_dof = 6;
  double improve_ratio_threshold = 0.1;
  double trust_shrink_ratio = 0.9;
  double trust_expand_ratio = 1.3;
  
  double start_vec_array[] = {0, 0, 0, 0, 0, 0};
  double goal_vec_array[] = {6, 0, 6, 0, 0, 0};

  vector<double> start_vec(start_vec_array, start_vec_array + n_dof);
  vector<double> goal_vec(goal_vec_array, goal_vec_array + n_dof);

  {
    Config config;
    config.add(new Parameter<bool>("plotting", &plotting, "plotting"));
    config.add(new Parameter<bool>("verbose", &verbose, "verbose"));
    config.add(new Parameter<double>("env_transparency", &env_transparency, "env_transparency"));
    config.add(new Parameter<int>("n_steps", &n_steps, "n_steps"));
    config.add(new Parameter<double>("turning_radius", &turning_radius, "turning_radius"));
    config.add(new Parameter<double>("improve_ratio_threshold", &improve_ratio_threshold, "improve_ratio_threshold"));
    config.add(new Parameter<double>("trust_shrink_ratio", &trust_shrink_ratio, "trust_shrink_ratio"));
    config.add(new Parameter<double>("trust_expand_ratio", &trust_expand_ratio, "trust_expand_ratio"));
    config.add(new Parameter< vector<double> >("s", &start_vec, "s"));
    config.add(new Parameter< vector<double> >("g", &goal_vec, "g"));
    CommandParser parser(config);
    parser.read(argc, argv);
  }

  RaveInitialize(false, verbose ? Level_Debug : Level_Info);
  EnvironmentBasePtr env = RaveCreateEnvironment();
  env->StopSimulation();
  OSGViewerPtr viewer = OSGViewer::GetOrCreate(env);
  assert(viewer);

  env->Load(string(DATA_DIR) + "/prostate.env.xml");//needleprob.env.xml");
  viewer->SetAllTransparency(env_transparency);
  RobotBasePtr robot = GetRobot(*env);
  RobotAndDOFPtr rad(new RobotAndDOF(robot, vector<int>(), 11, OR::Vector(0,0,1)));


  OptProbPtr prob(new OptProb());
  VarArray trajvars;
  AddVarArray(*prob, n_steps, n_dof, "j", trajvars);


  O3Helper helper(robot, trajvars.block(0,3,n_steps,3));
  helper.ConfigureProblem(*prob);
  
  VectorXd start(n_dof); for (int i = 0; i < n_dof; ++i) start[i] = start_vec[i];
  VectorXd goal(n_dof); for (int i = 0; i < n_dof; ++i) goal[i] = goal_vec[i];

  VectorXd vel_coeffs = VectorXd::Ones(3);
  prob->addCost(CostPtr(new JointVelCost(trajvars.block(0,0,n_steps, 3), vel_coeffs)));

  helper.AddAngVelCosts(*prob, vel_coeffs[0]);

  double dtheta_lb = (goal.topRows(3) - start.topRows(3)).norm() / (n_steps-1)/turning_radius;
  Var dthetavar = prob->createVariables(singleton<string>("speed"), singleton<double>(dtheta_lb),singleton<double>(INFINITY))[0];

  Str2Dbl tag2dist_pen(0.025), tag2coeff(20);

  tag2coeff.insert( std::pair<string, double>("KinBodyProstate", 0.0) );
  tag2coeff.insert( std::pair<string, double>("KinBodyDermis", 0.0) );
  tag2coeff.insert( std::pair<string, double>("KinBodyEpidermis", 0.0) );
  tag2coeff.insert( std::pair<string, double>("KinBodyHypodermis", 0.0) );

  for (int i=0; i < n_steps-1; ++i) {
    VarVector vars0 = trajvars.row(i), vars1 = trajvars.row(i+1);
    VectorOfVectorPtr f(new NeedleError(helper.m_rbs[i], helper.m_rbs[i+1], turning_radius));
    VectorXd coeffs = VectorXd::Ones(6);
    VarVector vars = concat(vars0, vars1);  
    vars.push_back(dthetavar);
    prob->addConstraint(ConstraintPtr(new ConstraintFromFunc(f, vars, coeffs, EQ, (boost::format("needle%i")%i).str())));
    prob->addConstraint(ConstraintPtr(new CollisionTaggedConstraint(tag2dist_pen, tag2coeff, helper.m_rbs[i], vars0)));//, vars1)));
  }


  for (int j=0; j < n_dof; ++j) {
    prob->addLinearConstraint(exprSub(AffExpr(trajvars(0,j)), start[j]), EQ);
  }
  for (int j=0; j < 3; ++j) { // NO orientation constraint
    prob->addLinearConstraint(exprSub(AffExpr(trajvars(n_steps-1,j)), goal[j]), EQ);
  }
  

  BasicTrustRegionSQP opt(prob);
  helper.ConfigureOptimizer(opt);
  opt.max_iter_ = 500;    
  opt.improve_ratio_threshold_ = improve_ratio_threshold;
  opt.trust_shrink_ratio_ = trust_shrink_ratio;
  opt.trust_expand_ratio_ = trust_expand_ratio;

  boost::shared_ptr<TrajPlotter> plotter;
  plotter.reset(new TrajPlotter(helper.m_rbs, trajvars));
  if (plotting) {
    opt.addCallback(boost::bind(&TrajPlotter::OptimizerCallback, boost::ref(plotter), _1, _2));
  }

  MatrixXd initTraj(n_steps, n_dof);  
  for (int idof = 0; idof < n_dof; ++idof) {
    initTraj.col(idof) = VectorXd::LinSpaced(n_steps, start[idof], goal[idof]);
  }
  DblVec initVec = trajToDblVec(initTraj);
  initVec.push_back(dtheta_lb);
  opt.initialize(initVec);
  
  opt.optimize();

  vector<Matrix4d> poses;

  for (int i = 0; i < helper.m_rbs.size(); ++i) {
    helper.m_rbs[i]->SetDOFValues(prob->getModel()->getVarValues(trajvars.row(i)));
    OR::TransformMatrix M(helper.m_rbs[i]->m_body->GetTransform());
    Matrix4d pose;
    pose << M.m[0], M.m[1], M.m[2], M.trans[0],
            M.m[4], M.m[5], M.m[6], M.trans[1],
            M.m[8], M.m[9], M.m[10], M.trans[2],
            0,      0,      0,                1;

    //OR::Vector trans = T.trans;
    //OR::Vector rot = OpenRAVE::geometry::axisAngleFromQuat(T.rot);
    //VectorXd x(6); x << trans[0], trans[1], trans[2], rot[0], rot[1], rot[2];
    poses.push_back(pose);//expUp(x));
  }

  for (int i = 0; i < poses.size() - 1; ++i) {
    cout << "twist at time " << i << ": ";
    printVector(logDown(poses[i].inverse() * poses[i+1]));
    cout << endl;
  }

  plotter->PlotBothTrajectories(prob, opt, start);
  
  RaveDestroy();

}
