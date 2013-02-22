#include "trajopt/belief.hpp"
#include <boost/bind.hpp>
#include "ipi/sco/expr_ops.hpp"
#include "ipi/sco/modeling_utils.hpp"
#include "utils/eigen_conversions.hpp"

using namespace std;
using namespace ipi::sco;
using namespace Eigen;
using namespace OpenRAVE;
using namespace util;

namespace trajopt {

BeliefRobotAndDOF::BeliefRobotAndDOF(OpenRAVE::RobotBasePtr _robot, const IntVec& _joint_inds, int _affinedofs, const OR::Vector _rotationaxis) :
	RobotAndDOF(_robot, _joint_inds, _affinedofs, _rotationaxis),
	generator(boost::mt19937(time(NULL)+rand()), boost::normal_distribution<>(0, 1))
{}

MatrixXd BeliefRobotAndDOF::GetDynNoise() {
	VectorXd diag_noise(3);
	diag_noise << 0.08, 0.13, 0.18;
	return diag_noise.asDiagonal();
}

MatrixXd BeliefRobotAndDOF::GetObsNoise() { return VectorXd::Constant(3,0.09).asDiagonal(); }

VectorXd BeliefRobotAndDOF::Observe(VectorXd dofs, VectorXd r) {
		VectorXd z = dofs + GetObsNoise()*r;
		return z;
}

VectorXd BeliefRobotAndDOF::Dynamics(VectorXd dofs, VectorXd u, VectorXd q) {
		VectorXd dofs1 = dofs+u + GetDynNoise()*q;
		return dofs1;
}

VectorXd BeliefRobotAndDOF::beta(VectorXd theta0, VectorXd u0) {
	VectorXd theta;
	ekfUpdate(VectorXd::Zero(GetObsSize()), u0, theta0, theta);
	return theta;
}

VectorXd BeliefRobotAndDOF::VectorXdRand(int size, double sigma) {
	VectorXd v(size);
	for (int i=0; i<size; i++) v(i) = generator();
	return v;
}

void BeliefRobotAndDOF::composeBelief(VectorXd x, MatrixXd V, VectorXd& theta) {
	int n_dof = GetDOF();
	theta.resize(n_dof+n_dof*n_dof);
	theta.topRows(n_dof) = x;
	V.resize(n_dof*n_dof,1);
	theta.bottomRows(n_dof*n_dof) = V;
}

void BeliefRobotAndDOF::decomposeBelief(VectorXd theta, VectorXd& x, MatrixXd& V) {
	int n_dof = GetDOF();
	x = theta.topRows(n_dof);
	V = theta.bottomRows(n_dof*n_dof);
	V.resize(n_dof,n_dof);
}

void BeliefRobotAndDOF::ekfUpdate(VectorXd z0, VectorXd u0, VectorXd theta0, VectorXd& theta) {
	VectorXd xest0, xest;
	MatrixXd Vest0, Vest;
	decomposeBelief(theta0, xest0, Vest0);
	ekfUpdate(z0, u0, xest0, Vest0, xest, Vest);
	composeBelief(xest, Vest, theta);
}

void BeliefRobotAndDOF::ekfUpdate(VectorXd z0, VectorXd u0, VectorXd xest0, MatrixXd Vest0, VectorXd& xest, MatrixXd& Vest) {
//	cout << "EKF"<< endl;
//	cout << "z0" << endl;
//	cout << z0 << endl;
//	cout << "u0" << endl;
//	cout << u0 << endl;
//	cout << "xest0" << endl;
//	cout << xest0 << endl;
//	cout << "Vest0" << endl;
//	cout << Vest0 << endl;
	VectorXd q = VectorXd::Zero(GetQSize());

	MatrixXd A = calcForwardNumJac(boost::bind(&BeliefRobotAndDOF::Dynamics, this, _1, u0, q), xest0);
	MatrixXd Q = calcForwardNumJac(boost::bind(&BeliefRobotAndDOF::Dynamics, this, xest0, u0, _1), q);

	VectorXd xpred = Dynamics(xest0, u0, q);
	MatrixXd Vpred = A * (Vest0 * Vest0.transpose()) * A.transpose() + Q*Q;

	VectorXd r = VectorXd::Zero(GetRSize());
	MatrixXd C = calcForwardNumJac(boost::bind(&BeliefRobotAndDOF::Observe, this, _1, r), xpred);
	MatrixXd R = calcForwardNumJac(boost::bind(&BeliefRobotAndDOF::Observe, this, xpred, _1), r);

	MatrixXd C_Vpred = C*Vpred;

	MatrixXd A_K = C_Vpred*C.transpose() + R*R.transpose();
  PartialPivLU<MatrixXd> solver(A_K);
  MatrixXd L_transpose = solver.solve(C_Vpred);
  MatrixXd L = L_transpose.transpose();

	xest = xpred + L*(z0-Observe(xpred, r));

	LLT<MatrixXd> lltofVest((MatrixXd::Identity(GetDOF(),GetDOF()) - L*C) * Vpred);
	Vest = lltofVest.matrixL();

//	cout << "xest" << endl;
//	cout << xest << endl;
//	cout << "Vest" << endl;
//	cout << Vest << endl;
}



BeliefDynamicsConstraint::BeliefDynamicsConstraint(const VarVector& theta0_vars,	const VarVector& theta1_vars, const VarVector& u_vars, BeliefRobotAndDOFPtr brad) :
    Constraint("BeliefDynamics"), brad_(brad), theta0_vars_(theta0_vars), theta1_vars_(theta1_vars), u_vars_(u_vars), type_(EQ)
{}

vector<double> BeliefDynamicsConstraint::value(const vector<double>& xin) {
  VectorXd theta0_hat = getVec(xin, theta0_vars_);
  VectorXd theta1_hat = getVec(xin, theta1_vars_);
  VectorXd u_hat = getVec(xin, u_vars_);
  return toDblVec(brad_->beta(theta0_hat, u_hat) - theta1_hat);
}
ConvexConstraintsPtr BeliefDynamicsConstraint::convex(const vector<double>& xin, Model* model) {
  VectorXd theta0_hat = getVec(xin, theta0_vars_);
  VectorXd theta1_hat = getVec(xin, theta1_vars_);
  VectorXd u_hat = getVec(xin, u_vars_);

  // linearize belief dynamics around theta0_hat and u_hat
  // beta is the belief dynamics without noise
  MatrixXd A = calcForwardNumJac(boost::bind(&BeliefRobotAndDOF::beta, brad_.get(), _1, u_hat), theta0_hat);
  MatrixXd B = calcForwardNumJac(boost::bind(&BeliefRobotAndDOF::beta, brad_.get(), theta0_hat, _1), u_hat);
  VectorXd c = brad_->beta(theta0_hat, u_hat);

//  // test convexification
//  cout << "theta0" << endl;
//  VectorXd theta0 = theta0_hat + brad_->VectorXdRand(12,0.1);
//  cout << theta0 << endl;
//  cout << "u "<< endl;
//  VectorXd u = u_hat + brad_->VectorXdRand(3,0.1);
//  cout << u << endl;
//  cout << "theta1_approx" << endl;
//  VectorXd theta1_approx = A * (theta0 - theta0_hat) + B * (u - u_hat) + c;
//  cout << theta1_approx << endl;
//  cout << "theta1" << endl;
//  VectorXd theta1 = brad_->beta(theta0, u);
//  cout << theta1 << endl;

  // equality constraint
  // theta1_vars_ = A * (theta0_vars_ - theta0_hat) + B * (u_vars_ - u_hat) + c
  ConvexConstraintsPtr out(new ConvexConstraints(model));
  assert(A.rows() == B.rows());
  for (int i=0; i < A.rows(); ++i) {
      AffExpr aff_theta0;
      aff_theta0.constant = c[i] - A.row(i).dot(theta0_hat);
      aff_theta0.coeffs = toDblVec(A.row(i));
      aff_theta0.vars = theta0_vars_;
      AffExpr aff_u;
      aff_u.constant = - B.row(i).dot(u_hat);
      aff_u.coeffs = toDblVec(B.row(i));
      aff_u.vars = u_vars_;
      AffExpr aff_theta1;
      aff_u.constant = 0;
      aff_u.coeffs = vector<double>(theta1_vars_.size(),1);
      aff_u.vars = theta1_vars_;
      AffExpr aff_theta0_u = exprAdd(aff_theta0, aff_u);
      AffExpr aff = exprSub(aff_theta0_u, aff_theta1);
      aff = cleanupAff(aff);
      out->addEqCnt(aff);
  }
  return out;
}

/*
#include "osgviewer/osgviewer.hpp"
void BeliefDynamicsConstraint::Plot(const DblVec& x, OR::EnvironmentBase& env, std::vector<OR::GraphHandlePtr>& handles) {
  VectorXd theta0 = getVec(x, theta0_vars_);

  VectorXd xest;
  MatrixXd Vest;
  brad_->decomposeBelief(theta0, xest, Vest);

  OR::RobotBasePtr robot = brad_->GetRobot();
  brad_->SetDOFValues(toDblVec(xest));
  OR::KinBody::LinkPtr link = robot->GetLink("Finger");
  MatrixXd jac = brad_->PositionJacobian(3,OR::Vector(0,0,0,0));


  boost::shared_ptr<OSGViewer> viewer = OSGViewer::GetOrCreate(brad_->GetRobot()->GetEnv());

  Eigen::Matrix3d cov;
  cov << 0.025, 0.0075, 0.00175, 0.0075, 0.0070, 0.00135, 0.00175, 0.00135, 0.00043;
  Eigen::Vector3d mean(0,0,0.1);

  handles.push_back(viewer.PlotEllipsoid(beliefToTransform(mean,cov), OR::Vector(1,0,0,1)));
}
*/

MatrixXd calcForwardNumJac(VectorOfVectorFun f, const VectorXd& x, double epsilon) {
  VectorXd y = f(x);
  MatrixXd out(y.size(), x.size());
  VectorXd xpert = x;
  for (size_t i=0; i < size_t(x.size()); ++i) {
    xpert(i) = x(i) + epsilon;
    VectorXd ypert = f(xpert);
    out.col(i) = (ypert - y) / epsilon;
    xpert(i) = x(i);
  }
  return out;
}

osg::Matrix beliefToTransform(const Eigen::Vector3d& mean, const Eigen::Matrix3d& cov) {
//	Eigen::Matrix3d cov;
//	cov << 0.025, 0.0075, 0.00175, 0.0075, 0.0070, 0.00135, 0.00175, 0.00135, 0.00043;
//	Eigen::Vector3d mean(0,0,0.1);

	Eigen::EigenSolver<Eigen::Matrix3d> es(cov);
	Eigen::Matrix4d t = Eigen::Matrix4d::Identity();
	t.block(0,0,3,3) = es.eigenvectors().real() * es.eigenvalues().real().cwiseSqrt().asDiagonal();
	t.block(0,3,3,1) = mean;

	osg::Matrix osg_t;
	osg_t.set(t.data());
	return osg_t;
}

}
