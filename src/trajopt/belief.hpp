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

  Eigen::VectorXd Observe(const Eigen::VectorXd& dofs, const Eigen::VectorXd& r);
  Eigen::VectorXd Dynamics(const Eigen::VectorXd& dofs, const Eigen::VectorXd& u, const Eigen::VectorXd& q);
  Eigen::VectorXd BeliefDynamics(const Eigen::VectorXd& theta0, const Eigen::VectorXd& u0);
  Eigen::VectorXd VectorXdRand(int size);

  void composeBelief(const Eigen::VectorXd& x, const Eigen::MatrixXd& rt_S, VectorXd& theta);
  void decomposeBelief(const Eigen::VectorXd& theta, VectorXd& x, Eigen::MatrixXd& rt_S);
  void ekfUpdate(const Eigen::VectorXd& u0, const Eigen::VectorXd& xest0, const Eigen::MatrixXd& Vest0, VectorXd& xest, Eigen::MatrixXd& Vest);

  Eigen::MatrixXd EndEffectorJacobian(const Eigen::VectorXd& x0) {
  	Eigen::MatrixXd jac(3,3);
		double l1 = 0.16;
		double l2 = 0.16;
		double l3 = 0.08;
		double s1 = -l1 * sin(x0(0));
		double s2 = -l2 * sin(x0(0)+x0(1));
		double s3 = -l3 * sin(x0(0)+x0(1)+x0(2));
		double c1 = l1 * cos(x0(0));
		double c2 = l2 * cos(x0(0)+x0(1));
		double c3 = l3 * cos(x0(0)+x0(1)+x0(2));
		jac << s1+s2+s3, s2+s3, s3, c1+c2+c3, c2+c3, c3, 0, 0, 0;
		return jac;
  }

  Eigen::VectorXd sqrt_sigma;
};
typedef boost::shared_ptr<BeliefRobotAndDOF> BeliefRobotAndDOFPtr;

class BeliefDynamicsConstraint: public Constraint {
public:
	BeliefDynamicsConstraint(const VarVector& theta0_vars,	const VarVector& theta1_vars, const VarVector& u_vars, BeliefRobotAndDOFPtr brad);
  vector<double> value(const vector<double>& x);
  ConvexConstraintsPtr convex(const vector<double>& x, Model* model);
  ConstraintType type() {return type_;}
//  void Plot(const DblVec& x, OR::EnvironmentBase& env, std::vector<OR::GraphHandlePtr>& handles);
protected:
  BeliefRobotAndDOFPtr brad_;
  VarVector theta0_vars_;
  VarVector theta1_vars_;
  VarVector u_vars_;
  ConstraintType type_;
};

typedef boost::function<VectorXd(VectorXd)> VectorOfVectorFun;
typedef boost::shared_ptr<VectorOfVectorFun> VectorOfVectorFunPtr;
Eigen::MatrixXd calcNumJac(VectorOfVectorFun f, const VectorXd& x, double epsilon=0.00048828125);
osg::Matrix gaussianToTransform(const Eigen::Vector3d& mean, const Eigen::Matrix3d& cov);

}
