#pragma once
#include "trajopt/common.hpp"
#include "osgviewer/osgviewer.hpp"
#include "robot_and_dof.hpp"
#include <time.h>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/enable_shared_from_this.hpp>
#include <Eigen/Dense>
#include <Eigen/Geometry>

namespace trajopt {

class BeliefRobotAndDOF : public RobotAndDOF, public boost::enable_shared_from_this<BeliefRobotAndDOF> {
private:
  boost::variate_generator<boost::mt19937, boost::normal_distribution<> > generator;
public:
  BeliefRobotAndDOF(OpenRAVE::RobotBasePtr _robot, const IntVec& _joint_inds, int _affinedofs=0, const OR::Vector _rotationaxis=OR::Vector());

  Eigen::MatrixXd GetDynNoise();
  inline int GetQSize() { return GetDynNoise().cols(); }

  Eigen::MatrixXd GetObsNoise();
  inline int GetRSize() { return GetObsNoise().cols(); }
  inline int GetObsSize() { return GetObsNoise().rows(); }

  Eigen::VectorXd Observe(Eigen::VectorXd dofs, Eigen::VectorXd r);
  Eigen::VectorXd Dynamics(Eigen::VectorXd dofs, Eigen::VectorXd u, Eigen::VectorXd q);
  Eigen::VectorXd BeliefDynamics(Eigen::VectorXd theta, Eigen::VectorXd u);
  Eigen::VectorXd beta(Eigen::VectorXd theta0, Eigen::VectorXd u0);
  Eigen::VectorXd VectorXdRand(int size, double sigma=1);

  void composeBelief(Eigen::VectorXd x, Eigen::MatrixXd V, VectorXd& theta);
  void decomposeBelief(Eigen::VectorXd theta, VectorXd& x, Eigen::MatrixXd& V);
  void ekfUpdate(Eigen::VectorXd z0, Eigen::VectorXd u0, Eigen::VectorXd theta0, VectorXd& theta);
  void ekfUpdate(Eigen::VectorXd z0, Eigen::VectorXd u0, Eigen::VectorXd xest0, Eigen::MatrixXd Vest0, VectorXd& xest, Eigen::MatrixXd& Vest);

  Eigen::VectorXd sqrt_sigma;
};
typedef boost::shared_ptr<BeliefRobotAndDOF> BeliefRobotAndDOFPtr;

class BeliefDynamicsConstraint: public Constraint {
public:
	BeliefDynamicsConstraint(const VarVector& theta0_vars,	const VarVector& theta1_vars, const VarVector& u_vars, BeliefRobotAndDOFPtr brad);
  vector<double> value(const vector<double>& x);
  ConvexConstraintsPtr convex(const vector<double>& x, Model* model);
  ConstraintType type() {return type_;}
  //void Plot(const DblVec& x, OR::EnvironmentBase& env, std::vector<OR::GraphHandlePtr>& handles);
protected:
  BeliefRobotAndDOFPtr brad_;
  VarVector theta0_vars_;
  VarVector theta1_vars_;
  VarVector u_vars_;
  ConstraintType type_;
};

typedef boost::function<VectorXd(VectorXd)> VectorOfVectorFun;
typedef boost::shared_ptr<VectorOfVectorFun> VectorOfVectorFunPtr;
Eigen::MatrixXd calcForwardNumJac(VectorOfVectorFun f, const VectorXd& x, double epsilon=1e-5);
osg::Matrix beliefToTransform(const Eigen::Vector3d& mean, const Eigen::Matrix3d& cov);

}
