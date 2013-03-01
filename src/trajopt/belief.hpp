#pragma once
#include "ipi/sco/modeling.hpp"
#include "ipi/sco/modeling_utils.hpp"
#include "ipi/sco/sco_fwd.hpp"
#include <Eigen/Core>
#include "trajopt/common.hpp"

#include "osgviewer/osgviewer.hpp"
#include "robot_and_dof.hpp"
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>

namespace trajopt {

class BeliefRobotAndDOF : public RobotAndDOF {
private:
  boost::variate_generator<boost::mt19937, boost::normal_distribution<> > generator;
  int n_theta;
public:
  BeliefRobotAndDOF(OpenRAVE::RobotBasePtr _robot, const IntVec& _joint_inds, int _affinedofs=0, const OR::Vector _rotationaxis=OR::Vector());

  Eigen::MatrixXd GetDynNoise();
  inline int GetQSize() { return GetDOF(); }
//  inline int GetQSize() { return GetDynNoise().cols(); }

  Eigen::MatrixXd GetObsNoise();
  inline int GetRSize() { return GetDOF(); }
  inline int GetObsSize() { return GetDOF(); }
//  inline int GetRSize() { return GetObsNoise().cols(); }
//  inline int GetObsSize() { return GetObsNoise().rows(); }

  // UKF update functions
  MatrixXd sigmaPoints(const VectorXd& theta);
  MatrixXd sigmaPoints(const VectorXd& mean, const MatrixXd& cov);
  VectorXd sigmaPoint(const VectorXd& mean, const MatrixXd& cov, int idx);
  void ukfUpdate(const VectorXd& u0, const VectorXd& xest0, const MatrixXd& Vest0, VectorXd& xest, MatrixXd& Vest);

  Eigen::VectorXd Observe(const Eigen::VectorXd& dofs, const Eigen::VectorXd& r);
  Eigen::VectorXd Dynamics(const Eigen::VectorXd& dofs, const Eigen::VectorXd& u, const Eigen::VectorXd& q);
  Eigen::VectorXd BeliefDynamics(const Eigen::VectorXd& theta0, const Eigen::VectorXd& u0);
  Eigen::VectorXd VectorXdRand(int size);

  void composeBelief(const Eigen::VectorXd& x, const Eigen::MatrixXd& rt_S, VectorXd& theta);
  void decomposeBelief(const Eigen::VectorXd& theta, VectorXd& x, Eigen::MatrixXd& rt_S);
  void ekfUpdate(const Eigen::VectorXd& u0, const Eigen::VectorXd& x0, const Eigen::MatrixXd& rtSigma0, VectorXd& xest, Eigen::MatrixXd& Vest);

  Eigen::MatrixXd EndEffectorJacobian(const Eigen::VectorXd& x0);
  void GetEndEffectorNoiseAsGaussian(const Eigen::VectorXd& theta, Eigen::VectorXd& mean, Eigen::MatrixXd& cov);

  inline int GetNTheta() { return n_theta; }

	OR::KinBody::LinkPtr link;
	// Scaled UKF update vars
	double alpha, beta, kappa;
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

class BeliefDynamicsConstraint2 : public ConstraintFromNumDiff {
public:
	BeliefDynamicsConstraint2(const VarVector& theta0_vars,	const VarVector& theta1_vars, const VarVector& u_vars,
			BeliefRobotAndDOFPtr brad, const BoolVec& enabled=BoolVec());
};

class CovarianceCost : public Cost {
public:
	CovarianceCost(const VarVector& rtSigma_vars, const Eigen::MatrixXd& Q, BeliefRobotAndDOFPtr brad);
  virtual ConvexObjectivePtr convex(const vector<double>& x, Model* model);
  virtual double value(const vector<double>&);
private:
  VarVector rtSigma_vars_;
  Eigen::MatrixXd Q_;
  QuadExpr expr_;
  BeliefRobotAndDOFPtr brad_;
};

class ControlCost : public Cost {
public:
	ControlCost(const VarArray& traj, const VectorXd& coeffs);
  virtual ConvexObjectivePtr convex(const vector<double>& x, Model* model);
  virtual double value(const vector<double>&);
private:
  VarArray vars_;
  VectorXd coeffs_;
  QuadExpr expr_;
};

typedef boost::function<VectorXd(VectorXd)> VectorOfVectorFun;
typedef boost::shared_ptr<VectorOfVectorFun> VectorOfVectorFunPtr;
Eigen::MatrixXd calcNumJac(VectorOfVectorFun f, const VectorXd& x, double epsilon=0.00048828125);
osg::Matrix gaussianAsTransform(const Eigen::Vector3d& mean, const Eigen::Matrix3d& cov);

}
