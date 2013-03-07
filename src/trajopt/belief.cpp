#include "trajopt/belief.hpp"
#include <time.h>
#include <boost/bind.hpp>
#include "sco/expr_ops.hpp"
#include "sco/modeling_utils.hpp"
#include "utils/eigen_conversions.hpp"

using namespace std;
using namespace sco;
using namespace Eigen;
using namespace OpenRAVE;
using namespace util;

namespace {
	double sigmoid(double x,double mean) {
		double y = (x - mean);
		double s = (y/sqrt(1+y*y))+1.0;

		if (x < mean)
			return s*0.1;
		else
			return s*7.0;
	}
}

namespace trajopt {

BeliefRobotAndDOF::BeliefRobotAndDOF(OpenRAVE::RobotBasePtr _robot, const IntVec& _joint_inds, int _affinedofs, const OR::Vector _rotationaxis) :
			RobotAndDOF(_robot, _joint_inds, _affinedofs, _rotationaxis),
			generator(boost::mt19937(time(NULL)+rand()), boost::normal_distribution<>(0, 1)),
			n_theta(GetDOF() + GetDOF() * (GetDOF()+1) / 2)
{
	if (GetDOF() == 3) link = GetRobot()->GetLink("Finger");
	else link = GetRobot()->GetLink("Base");

	// UKF vars
	alpha = 0.5;
	beta = 2.0;
	kappa = 5.0;
}

MatrixXd BeliefRobotAndDOF::GetDynNoise() {
	int n_dof = GetDOF();
	VectorXd diag_noise(n_dof);
	if (n_dof == 3)	diag_noise << 0.08, 0.13, 0.18;
	else diag_noise << 0.05, 0.05;
	return diag_noise.asDiagonal();
}

MatrixXd BeliefRobotAndDOF::GetObsNoise() {
	int n_dof = GetDOF();
	VectorXd diag_noise(n_dof);
	if (n_dof == 3) diag_noise << 0.09, 0.09, 0.09;
	else diag_noise << 0.005, 0.005;
	assert(0);
	return diag_noise.asDiagonal();
}

VectorXd BeliefRobotAndDOF::Observe(const VectorXd& dofs, const VectorXd& r) {
	int n_dof = GetDOF();
	OR::RobotBase::RobotStateSaver saver = const_cast<BeliefRobotAndDOF*>(this)->Save();
	SetDOFValues(toDblVec(dofs));
	OR::Vector trans = link->GetTransform().trans;

	VectorXd z(n_dof);
	if (n_dof == 3) {
		z = Vector3d(trans.x, trans.y, trans.z);
		//z += sigmoid(trans.y, -0.2)*GetObsNoise()*r;
		z += (0.5*pow(trans.y+0.2,2)+1)*r; // as in the Platt paper
		//	z += ((trans.y+0.2)/0.4)*GetObsNoise()*r;
		//	if (trans.y<-0.2) z += 0.1*GetObsNoise()*r;
		//	else z += 10*GetObsNoise()*r;
	} else {
		z = Vector2d(trans.x, trans.y) + (0.5*pow(5.0 - trans.x,2)+0.001)*r; // as in the Platt paper
	}
	return z;
}

VectorXd BeliefRobotAndDOF::Dynamics(const VectorXd& dofs, const VectorXd& u, const VectorXd& q) {
	VectorXd dofs1 = dofs+u + GetDynNoise()*q;
	//for (int i=0; i<q.size(); i++) assert(q[i] == 0);
	//VectorXd dofs1 = dofs+u;
	return dofs1;
}

VectorXd BeliefRobotAndDOF::BeliefDynamics(const VectorXd& theta0, const VectorXd& u0) {
	VectorXd x0, x;
	MatrixXd rtSigma0, rtSigma;
	decomposeBelief(theta0, x0, rtSigma0);
	ekfUpdate(u0, x0, rtSigma0, x, rtSigma);
	VectorXd theta;
	composeBelief(x, rtSigma, theta);

//	cout << "rtSigma " << endl << rtSigma << endl;
//	cout << "theta " << theta.transpose() << endl << endl;
//	decomposeBelief(theta, x, rtSigma);
//	cout << "rtSigma " << endl << rtSigma << endl;
//	cout << "theta " << theta.transpose() << endl << endl;
//	cout << "---------------" << endl;

	return theta;
}

VectorXd BeliefRobotAndDOF::VectorXdRand(int size) {
	VectorXd v(size);
	for (int i=0; i<size; i++) v(i) = generator();
	return v;
}

MatrixXd BeliefRobotAndDOF::sigmaPoints(const VectorXd& theta) {
	VectorXd x;
	MatrixXd rtSigma;
	decomposeBelief(theta, x, rtSigma);
	return sigmaPoints(x, rtSigma*rtSigma.transpose());
}

MatrixXd BeliefRobotAndDOF::sigmaPoints(const VectorXd& mean, const MatrixXd& cov)
{
	int n_dof = GetDOF();
	int n_r = GetRSize();

	int L = n_dof + n_r;

	double lambda = alpha*alpha*(L + kappa) - L;
	double w = 1 / (2*(L + lambda));
	double mw = lambda / (L + lambda);
	double vw = mw + (1 - alpha*alpha + beta);

	MatrixXd sigmapts(n_dof, 2*n_dof+1);
	sigmapts.col(0) = mean;

	Eigen::JacobiSVD<MatrixXd, NoQRPreconditioner> svd((n_dof+lambda)*cov, ComputeThinU | ComputeThinV);
	MatrixXd rt_scaled_cov = svd.matrixU() * svd.singularValues().array().sqrt().matrix().asDiagonal() * svd.matrixV().transpose();

	for(int i = 0; i < n_dof; ++i){
		sigmapts.col(2*i+1) = (mean + rt_scaled_cov.col(i));
		sigmapts.col(2*i+2) = (mean - rt_scaled_cov.col(i));
	}

	return sigmapts;
}

// TODO add motion noise
void BeliefRobotAndDOF::ukfUpdate(const VectorXd& u0, const VectorXd& x0, const MatrixXd& rtSigma0, VectorXd& x, MatrixXd& rtSigma)
{
	int n_dof = GetDOF();
	int n_u = u0.rows();
	int n_r = GetRSize();
	int n_z = GetObsSize();

 	int L = n_dof + n_r;

	double lambda = alpha*alpha*(L + kappa) - L;
	double w = 1 / (2*(L + lambda));
	double mw = lambda / (L + lambda);
	double vw = mw + (1 - alpha*alpha + beta);

	MatrixXd Sigma0 = rtSigma0 * rtSigma0.transpose();

	// propagate sigma points through f
	VectorXd q = VectorXd::Zero(GetQSize());
	MatrixXd sigmapts(n_dof, 2*n_dof+1);

	sigmapts.col(0) = Dynamics(x0, u0, q);

	Eigen::JacobiSVD<MatrixXd, NoQRPreconditioner> svd_Sigma0((L+lambda)*Sigma0, ComputeThinU | ComputeThinV);
	MatrixXd rt_scaled_Sigma0 = svd_Sigma0.matrixU() * svd_Sigma0.singularValues().array().sqrt().matrix().asDiagonal() * svd_Sigma0.matrixV().transpose();

	for (int i = 0; i < n_dof; ++i) {
		sigmapts.col(2*i+1) = Dynamics(x0 + rt_scaled_Sigma0.col(i), u0, q);
		sigmapts.col(2*i+2) = Dynamics(x0 - rt_scaled_Sigma0.col(i), u0, q);
	}

	// calculate mean -- O(xDim^2)
	x = mw * sigmapts.col(0);
	for (int i = 1; i < sigmapts.cols(); ++i) {
		x += w * sigmapts.col(i);
	}

	// calculate variance -- O(xDim^3)
	MatrixXd Sigma = vw * (sigmapts.col(0) - x)*(sigmapts.col(0) - x).transpose();
	for (int i = 1; i < sigmapts.cols(); ++i) {
		Sigma += w * (sigmapts.col(i) - x)*(sigmapts.col(i) - x).transpose();
	}

	// Measurement Update
	MatrixXd Z(n_z, 2*(n_dof+n_z)+1);
	VectorXd r = VectorXd::Zero(GetRSize());

	int idx = 0;
	for (int i = 0; i < sigmapts.cols(); ++i) {
		Z.col(idx) = Observe(sigmapts.col(i), r);
		++idx;
	}

	double factor = sqrt(L + lambda);
	for (int i = sigmapts.cols(); i < 2*(n_dof+n_z)+1; ++i) {
		r(i) = factor;
		Z.col(idx) = Observe(sigmapts.col(0), r);
		++idx;
		r(i) = -factor;
		Z.col(idx) = Observe(sigmapts.col(0), r);
		++idx;
		r(i) = 0;
	}

	// calculate mean -- O(xDim*zDim + zDim^2)
	VectorXd z0(n_z);
	z0 = mw * Z.col(0);
	for (int i = 1; i < Z.cols(); ++i) {
		z0 += w * Z.col(i);
	}

	// calculate variance -- O(zDim^3 + xDim*zDim^2)
	MatrixXd Pzz = vw * (Z.col(0) - z0)*(Z.col(0) - z0).transpose();
	for (int i = 1; i < Z.cols(); ++i) {
		Pzz += w * (Z.col(i) - z0)*(Z.col(i) - z0).transpose();
	}
	// calculate cross-covariance -- O(xDim^2*zDim + xDim*zDim^2)
	MatrixXd Pxz = vw * (sigmapts.col(0) - x0)*(Z.col(0) - z0).transpose();
	for (int i = 1; i < sigmapts.cols(); ++i) {
		Pxz += w * (sigmapts.col(i) - x0)*(Z.col(i) - z0).transpose();
	}
	for (int i = sigmapts.cols(); i < Z.cols(); ++i) {
		Pxz += w * (sigmapts.col(i) - x0)*(Z.col(i) - z0).transpose();
	}

	PartialPivLU<MatrixXd> solver(Pxz);
	MatrixXd K = solver.solve(Pzz); // Pxz/Pzz? Check.
	x += K*(Z.col(0) - z0); // O(xDim*zDim)
	Sigma -= Pxz*K.transpose(); // O(xDim^2*zDim)

	Eigen::JacobiSVD<MatrixXd, NoQRPreconditioner> svd_Sigma(Sigma, ComputeThinU | ComputeThinV);
	rtSigma = svd_Sigma.matrixU() * svd_Sigma.singularValues().array().sqrt().matrix().asDiagonal() * svd_Sigma.matrixV().transpose();
}

void BeliefRobotAndDOF::ekfUpdate(const VectorXd& u0, const VectorXd& x0, const MatrixXd& rtSigma0, VectorXd& x, MatrixXd& rtSigma) {
	int n_dof = GetDOF();

	VectorXd q = VectorXd::Zero(GetQSize());
	x = Dynamics(x0, u0, q);

	MatrixXd Sigma0 = rtSigma0 * rtSigma0.transpose();

	MatrixXd A = calcNumJac(boost::bind(&BeliefRobotAndDOF::Dynamics, this, _1, u0, q), x0);
	MatrixXd Q = calcNumJac(boost::bind(&BeliefRobotAndDOF::Dynamics, this, x0, u0, _1), q);
	MatrixXd Gamma0 = A * Sigma0 * A.transpose() + Q*Q.transpose();

	VectorXd r = VectorXd::Zero(GetRSize());
	MatrixXd C = calcNumJac(boost::bind(&BeliefRobotAndDOF::Observe, this, _1, r), x0);
	MatrixXd R = calcNumJac(boost::bind(&BeliefRobotAndDOF::Observe, this, x0, _1), r);

	MatrixXd W;
	MatrixXd A_K = C*Gamma0*C.transpose() + R*R.transpose();
	PartialPivLU<MatrixXd> solver(A_K);
	MatrixXd L = solver.solve(C*Gamma0);
	MatrixXd Sigma = Gamma0 - Gamma0 * C.transpose() * L;

//	LLT<MatrixXd> lltofSigma(Sigma);
//	rtSigma = lltofSigma.matrixL();
	Eigen::JacobiSVD<MatrixXd, NoQRPreconditioner> svd(Sigma, ComputeThinU | ComputeThinV);
	rtSigma = svd.matrixU() * svd.singularValues().array().sqrt().matrix().asDiagonal() * svd.matrixV().transpose();

	VectorXd s = svd.singularValues();
	for (int i=0; i<s.size(); i++)
		assert(s(i)>=0);
}

Eigen::MatrixXd BeliefRobotAndDOF::EndEffectorJacobian(const Eigen::VectorXd& x0) {
	int n_dof = GetDOF();

	Eigen::MatrixXd jac(n_dof,n_dof);

	if (n_dof == 3) {
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
	} else {
		jac << 1,0,0,1;
	}

//		// analytical jacobian computed in openrave
//		OR::RobotBase::RobotStateSaver saver = const_cast<BeliefRobotAndDOF*>(this)->Save();
//		DblMatrix Jxyz = PositionJacobian(link->GetIndex(), link->GetTransform().trans);
//		cout << "Jxyz" << endl;
//		std::cout << Jxyz << std::endl;
//		cout << "JAC" << endl;
//		cout << jac << endl;
//		cout << "-----------" << endl;

	return jac;
}

void BeliefRobotAndDOF::GetEndEffectorNoiseAsGaussian(const VectorXd& theta, VectorXd& mean, MatrixXd& cov) {
	VectorXd x;
	MatrixXd rt_Sigma;
	decomposeBelief(theta, x, rt_Sigma);

	MatrixXd jac = EndEffectorJacobian(x);

	OR::RobotBase::RobotStateSaver saver = const_cast<BeliefRobotAndDOF*>(this)->Save();
	SetDOFValues(toDblVec(x));
	OR::Vector trans = link->GetTransform().trans;
	mean = Vector3d(trans.x, trans.y, trans.z);
	MatrixXd partial_cov = jac * rt_Sigma * rt_Sigma.transpose() * jac.transpose();
	assert(partial_cov.rows() == partial_cov.cols());
	assert(partial_cov.rows() == 2 || partial_cov.rows() == 3);
	cov = MatrixXd::Zero(3,3);
	cov.topLeftCorner(partial_cov.rows(), partial_cov.cols()) = partial_cov;
}


MatrixXd calcNumJac(VectorOfVectorFun f, const VectorXd& x, double epsilon) {
	VectorXd y = f(x);
	MatrixXd out(y.size(), x.size());
	VectorXd x_plus = x;
	VectorXd x_minus = x;
	for (size_t i=0; i < size_t(x.size()); ++i) {
		x_plus(i) = x(i) + epsilon;
		VectorXd y_plus = f(x_plus);
		x_minus(i) = x(i) - epsilon;
		VectorXd y_minus = f(x_minus);
		out.col(i) = (y_plus - y_minus) / (2*epsilon);
		x_plus(i) = x(i);
		x_minus(i) = x(i);
	}
	return out;
}

osg::Matrix gaussianAsTransform(const Eigen::Vector3d& mean, const Eigen::Matrix3d& cov) {
	//	Eigen::Matrix3d cov;
	//	cov << 0.025, 0.0075, 0.00175, 0.0075, 0.0070, 0.00135, 0.00175, 0.00135, 0.00043;
	//	Eigen::Vector3d mean(0,0,0.1);

	Eigen::EigenSolver<Eigen::Matrix3d> es(cov);
	Eigen::Matrix4d t = Eigen::Matrix4d::Identity();
	t.block(0,0,3,3) = es.eigenvectors().real() * es.eigenvalues().real().cwiseSqrt().asDiagonal();
	t.block(0,3,3,1) = mean;
	Eigen::Matrix4d t_transpose = t.transpose();

	osg::Matrix osg_t;
	osg_t.set(t_transpose.data());
	return osg_t;
}

}
