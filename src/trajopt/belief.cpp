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

	VectorXi sigma_vec_inds(GetNTheta()-GetDOF());
	for (int i=0; i<sigma_vec_inds.rows(); i++) sigma_vec_inds(i) = i;
	MatrixXi sigma_matrix_inds = toSigmaMatrix(sigma_vec_inds);
	for (int j=0; j<sigma_matrix_inds.cols(); j++)
		sigma_col_to_indices.push_back(toDblVec((VectorXi) sigma_matrix_inds.col(j)));

	// UKF vars
	alpha = 0.1;
	beta = 2.0;
	kappa = 1.0;
}

void BeliefRobotAndDOF::SetBeliefValues(const DblVec& theta) {
	assert(theta.size() == GetNTheta());
	SetDOFValues(DblVec(theta.begin(), theta.begin()+GetDOF()));
	sigma_vec = DblVec(theta.begin()+GetDOF(), theta.end());
}

DblVec BeliefRobotAndDOF::GetBeliefValues() {
	DblVec theta = GetDOFValues();
	theta.insert(theta.end(), sigma_vec.begin(), sigma_vec.end());
	assert(theta.size() == GetNTheta());
	return theta;
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
	else diag_noise << 0.01, 0.01;
	assert(0);
	return diag_noise.asDiagonal();
}

VectorXd BeliefRobotAndDOF::Observe(const VectorXd& dofs, const VectorXd& r) {
	int n_dof = GetDOF();
	OR::RobotBase::RobotStateSaver saver = const_cast<BeliefRobotAndDOF*>(this)->Save();
	SetDOFValues(toDblVec(dofs));
	OR::Vector trans = link->GetTransform().trans;

	VectorXd z;
	if (n_dof == 3) {
		z = Vector3d(trans.x, trans.y, trans.z);
		//z += sigmoid(trans.y, -0.2)*GetObsNoise()*r;
		z += (0.5*pow(trans.y+0.2,2)+1)*r; // as in the Platt paper
		//	z += ((trans.y+0.2)/0.4)*GetObsNoise()*r;
		//	if (trans.y<-0.2) z += 0.1*GetObsNoise()*r;
		//	else z += 10*GetObsNoise()*r;
	}  else if (n_dof == 7) {
		z = Vector2d(0,0);
		OR::Vector beacon;
		beacon.x = -0.75; beacon.y = 0; beacon.z = 0.75;
		double dist = (trans.x - beacon.x)*(trans.x - beacon.x) + (trans.y - beacon.y)*(trans.y - beacon.y) + (trans.z - beacon.z)*(trans.z - beacon.z);
		z[0] = 1.0/(1.0 + dist) + 0.1*dist*r;
		z[1] = dofs[0] + 0.01*r;
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
	//ekfUpdate(u0, x0, rtSigma0, x, rtSigma);
	ukfUpdate(u0, x0, rtSigma0, x, rtSigma);
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

	Eigen::JacobiSVD<MatrixXd, NoQRPreconditioner> svd((L+lambda)*cov, ComputeThinU | ComputeThinV);
	MatrixXd rt_scaled_cov = svd.matrixU() * svd.singularValues().array().sqrt().matrix().asDiagonal() * svd.matrixV().transpose();

	for(int i = 0; i < n_dof; ++i){
		sigmapts.col(2*i+1) = (mean + rt_scaled_cov.col(i));
		sigmapts.col(2*i+2) = (mean - rt_scaled_cov.col(i));
	}

	return sigmapts;
}

void BeliefRobotAndDOF::ukfUpdate(const VectorXd& u0, const VectorXd& x0, const MatrixXd& rtSigma0, VectorXd& x, MatrixXd& rtSigma)
{
	int n_dof = GetDOF();
	int n_u = u0.rows();
	int n_r = GetRSize();
	int n_z = GetObsSize();
	int n_q = GetQSize();

 	int L = n_dof + n_q + n_r;

 	double ukfalpha = 0.001; //alpha;
 	double ukfbeta = 2; //beta;
 	double ukfkappa = 0; //kappa;

	double lambda = ukfalpha*ukfalpha*(L + ukfkappa) - L;
	double w = 1 / (2*(L + lambda));
	double mw = lambda / (L + lambda);
	double vw = mw + (1 - ukfalpha*ukfalpha + ukfbeta);

	//cout << "L: " << L << endl;
	//cout << "lambda: " << lambda << endl;
	//cout << "w: " << w << endl;
	//cout << "mw: " << mw << endl;
	//cout << "vw: " << vw << endl;

	MatrixXd Sigma0 = rtSigma0 * rtSigma0.transpose();

	//cout << "Sigma0:\n" << Sigma0 << endl;

	// propagate sigma points through f
	VectorXd q = VectorXd::Zero(GetQSize());
	MatrixXd sigmapts(n_dof, 2*(n_dof+n_q)+1);

	int idx = 0;
	sigmapts.col(idx++) = Dynamics(x0, u0, q);

	Eigen::JacobiSVD<MatrixXd, NoQRPreconditioner> svd_Sigma0((L+lambda)*Sigma0, ComputeThinU | ComputeThinV);
	MatrixXd rt_scaled_Sigma0 = svd_Sigma0.matrixU() * svd_Sigma0.singularValues().array().sqrt().matrix().asDiagonal() * svd_Sigma0.matrixV().transpose();

	for (int i = 0; i < n_dof; ++i) {
		sigmapts.col(idx++) = Dynamics(x0 + rt_scaled_Sigma0.col(i), u0, q);
		sigmapts.col(idx++) = Dynamics(x0 - rt_scaled_Sigma0.col(i), u0, q);
	}

	double qstep = sqrt(L + lambda);
	for (int i = 0; i < n_q; ++i) {
		q(i) = qstep;
		sigmapts.col(idx++) = Dynamics(x0, u0, q);
		q(i) = -qstep;
		sigmapts.col(idx++) = Dynamics(x0, u0, q);
		q(i) = 0;
	}

	//cout << "sigmapts:\n" << sigmapts << endl;

	// calculate mean -- O(xDim^2)
	VectorXd xhat(n_dof);
	xhat = (mw + 2*n_r*w) * sigmapts.col(0);
	for (int i = 1; i < sigmapts.cols(); ++i) {
		xhat += w * sigmapts.col(i);
	}

	//cout << "xhat: " << xhat.transpose() << endl;

	// calculate variance -- O(xDim^3)
	MatrixXd Sigma = (vw + 2*n_r*w) * (sigmapts.col(0) - xhat)*(sigmapts.col(0) - xhat).transpose();
	for (int i = 1; i < sigmapts.cols(); ++i) {
		Sigma += w * (sigmapts.col(i) - xhat)*(sigmapts.col(i) - xhat).transpose();
	}

	//cout << "Sigma predict:\n" << Sigma << endl;

	// Measurement Update
	MatrixXd Z(n_z, 2*(n_dof+n_q+n_r)+1);
	VectorXd r = VectorXd::Zero(GetRSize());

	idx = 0;
	for (int i = 0; i < sigmapts.cols(); ++i) {
		Z.col(idx++) = Observe(sigmapts.col(i), r);
	}

	double rstep = sqrt(L + lambda);
	for (int i = 0; i < n_r; ++i) {
		r(i) = rstep;
		Z.col(idx++) = Observe(sigmapts.col(0), r);
		r(i) = -rstep;
		Z.col(idx++) = Observe(sigmapts.col(0), r);
		r(i) = 0;
	}

	//cout << "Z:\n" << Z << endl;

	// calculate mean -- O(xDim*zDim + zDim^2)
	VectorXd zhat(n_z);
	zhat = mw * Z.col(0);
	for (int i = 1; i < Z.cols(); ++i) {
		zhat += w * Z.col(i);
	}

	//cout << "zhat: " << zhat.transpose() << endl;

	// calculate variance -- O(zDim^3 + xDim*zDim^2)
	MatrixXd Pzz = vw * (Z.col(0) - zhat)*(Z.col(0) - zhat).transpose();
	for (int i = 1; i < Z.cols(); ++i) {
		Pzz += w * (Z.col(i) - zhat)*(Z.col(i) - zhat).transpose();
	}
	//cout << "Pzz:\n" << Pzz << endl;

	// calculate cross-covariance -- O(xDim^2*zDim + xDim*zDim^2)
	MatrixXd Pxz = vw * (sigmapts.col(0) - xhat)*(Z.col(0) - zhat).transpose();
	for (int i = 1; i < sigmapts.cols(); ++i) {
		Pxz += w * (sigmapts.col(i) - xhat)*(Z.col(i) - zhat).transpose();
	}
	for (int i = sigmapts.cols(); i < Z.cols(); ++i) {
		Pxz += w * (sigmapts.col(0) - xhat)*(Z.col(i) - zhat).transpose();
	}

	//cout << "Pxz:\n" << Pxz << endl;

	PartialPivLU<MatrixXd> solver(Pzz);
	MatrixXd K = solver.solve(Pxz); // Pxz\Pzz? Check.
	//x += K*(Z.col(0) - z0); // O(xDim*zDim)

	//cout << "K:\n" << K << endl;

	x = xhat;
	Sigma -= Pxz*K.transpose(); // O(xDim^2*zDim)

	//cout << "Sigma update:\n" << Sigma << endl;

	//int num;
	//cin >> num;

	//exit(-1);

	Eigen::JacobiSVD<MatrixXd, NoQRPreconditioner> svd_Sigma(Sigma, ComputeThinU | ComputeThinV);
	rtSigma = svd_Sigma.matrixU() * svd_Sigma.singularValues().array().sqrt().matrix().asDiagonal() * svd_Sigma.matrixV().transpose();
}

// A \ B
MatrixXd backlashDivide(const MatrixXd& A, const MatrixXd B) {
	PartialPivLU<MatrixXd> solver(B);
	return (MatrixXd) solver.solve(A);
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

//void BeliefRobotAndDOF::ekfUpdate(const VectorXd& u0, const VectorXd& x0, const MatrixXd& rtSigma0, VectorXd& x, MatrixXd& rtSigma) {
//	int n_dof = GetDOF();
//
//	VectorXd q = VectorXd::Zero(GetQSize());
//	VectorXd x_pred = Dynamics(x0, u0, q);
//
//	MatrixXd Sigma0 = rtSigma0 * rtSigma0.transpose();
//
//	MatrixXd A = calcNumJac(boost::bind(&BeliefRobotAndDOF::Dynamics, this, _1, u0, q), x0);
//	MatrixXd Q = calcNumJac(boost::bind(&BeliefRobotAndDOF::Dynamics, this, x0, u0, _1), q);
//	MatrixXd Sigma_pred = A * Sigma0 * A.transpose() + Q*Q.transpose();
//
//	VectorXd r = VectorXd::Zero(GetRSize());
//	MatrixXd C = calcNumJac(boost::bind(&BeliefRobotAndDOF::Observe, this, _1, r), x_pred);
//	MatrixXd R = calcNumJac(boost::bind(&BeliefRobotAndDOF::Observe, this, x_pred, _1), r);
//
//	MatrixXd L = backlashDivide((MatrixXd) Sigma_pred * C.transpose(), (MatrixXd) C*Sigma_pred*C.transpose() + R*R.transpose());
//	x = x_pred;
//	MatrixXd Sigma = Sigma_pred - L*C*Sigma_pred;
//
////	LLT<MatrixXd> lltofSigma(Sigma);
////	rtSigma = lltofSigma.matrixL();
//	Eigen::JacobiSVD<MatrixXd, NoQRPreconditioner> svd(Sigma, ComputeThinU | ComputeThinV);
//	rtSigma = svd.matrixU() * svd.singularValues().array().sqrt().matrix().asDiagonal() * svd.matrixV().transpose();
//
//	VectorXd s = svd.singularValues();
//	for (int i=0; i<s.size(); i++)
//		assert(s(i)>=0);
//}

// theta needs to be set accordingly before calling this (i.e. call SetBeliefValues)
Eigen::MatrixXd BeliefRobotAndDOF::BeliefJacobian(int link_ind, int sigma_pt_ind, const OR::Vector& pt) {
	MatrixXd pos_jac = PositionJacobian(link->GetIndex(), pt);
	MatrixXd jac = MatrixXd::Zero(3, GetNTheta());
	jac.leftCols(GetDOF()) = pos_jac;

	float lambda;
	if (sigma_pt_ind == 0) lambda = 0;
	else lambda = ((sigma_pt_ind-1)%2)==0 ? 1 : -1; // TODO non-unit magnitude for lambda

	if (lambda != 0) {
		int sigma_matrix_col = (sigma_pt_ind - 1)/2;
		vector<int> sigma_vec_inds = sigmaColToSigmaIndices(sigma_matrix_col);
		assert(sigma_vec_inds.size() == GetDOF());
		for (int i=0; i<sigma_vec_inds.size(); i++) {
			int jac_col = GetDOF() + sigma_vec_inds[i];
			jac.col(jac_col) = lambda * pos_jac.col(i);
		}
	}

//	cout << "lambda " << lambda << endl;
//	if (lambda != 0) {
//		cout << "sigma_pt_ind " << sigma_pt_ind << endl;
//		int sigma_matrix_col = (sigma_pt_ind - 1)/2;
//		vector<int> sigma_vec_inds = sigmaColToSigmaIndices(sigma_matrix_col);
//		cout << "sigma_vec_inds " << toVectorXd(sigma_vec_inds).transpose() << endl;
//	}
//	cout << jac << endl << endl;

	return jac;
}

void BeliefRobotAndDOF::GetEndEffectorNoiseAsGaussian(const VectorXd& theta, VectorXd& mean, MatrixXd& cov) {
	VectorXd x;
	MatrixXd rt_Sigma;
	decomposeBelief(theta, x, rt_Sigma);

	OR::RobotBase::RobotStateSaver saver = const_cast<BeliefRobotAndDOF*>(this)->Save();
	SetDOFValues(toDblVec(x));
	OR::Vector trans = link->GetTransform().trans;
	MatrixXd jac = PositionJacobian(link->GetIndex(), trans);
	mean = Vector3d(trans.x, trans.y, trans.z);
	cov = jac * rt_Sigma * rt_Sigma.transpose() * jac.transpose();
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
