#pragma once
#include "trajopt/belief.hpp"
#include "ipi/sco/modeling.hpp"
#include "ipi/sco/modeling_utils.hpp"
#include "ipi/sco/sco_fwd.hpp"

namespace trajopt {
using namespace ipi::sco;

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

}
