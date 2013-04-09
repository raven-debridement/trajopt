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

namespace trajopt {

BeliefRobotAndDOF::BeliefRobotAndDOF(OpenRAVE::RobotBasePtr _robot, const IntVec& _joint_inds, int _affinedofs, const OR::Vector _rotationaxis) :
							RobotAndDOF(_robot, _joint_inds, _affinedofs, _rotationaxis),
							generator(boost::mt19937(time(NULL)+rand()), boost::normal_distribution<>(0, 1)),
							sigma_pts_scale(2)
{
	x_dim = GetDOF();
	b_dim = GetDOF() + GetDOF() * (GetDOF()+1) / 2;
	s_dim = GetDOF() * (GetDOF()+1) / 2;
	q_dim = GetDOF();
	(x_dim == 3 || x_dim == 2)? z_dim = x_dim: z_dim = 3;
	r_dim = z_dim;
	u_dim = GetDOF();

	if (x_dim == 3) endeffector = GetRobot()->GetLink("Finger");
	else if (x_dim == 7) endeffector = GetRobot()->GetLink("wam7");
	else if (x_dim == 2) endeffector = GetRobot()->GetLink("Base");

	VectorXi sigma_vec_inds(s_dim);

	for (int i=0; i<s_dim; i++) sigma_vec_inds(i) = i;

	MatrixXi sigma_matrix_inds(x_dim, x_dim);
	int idx = 0;
	for (int j = 0; j < x_dim; ++j) {
		for (int i = j; i < x_dim; ++i) {
			sigma_matrix_inds(i,j) = sigma_matrix_inds(j,i) = sigma_vec_inds(idx);
			idx++;
		}
	}

	for (int j=0; j<sigma_matrix_inds.cols(); j++)
		sigma_col_to_indices.push_back(toDblVec((VectorXi) sigma_matrix_inds.col(j)));

	// UKF vars
	alpha = 0.1;
	beta = 2.0;
	kappa = 1.0;
}

void BeliefRobotAndDOF::ForwardKinematics(const VectorXd& dofs, Vector3d& eetrans) {
	if (x_dim == 3) {
		double x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11;
		x0=cos((double)dofs[0]);
		x1=cos((double)dofs[1]);
		x2=sin((double)dofs[0]);
		x3=sin((double)dofs[1]);
		x4=sin((double)dofs[2]);
		x5=cos((double)dofs[2]);
		x6=0.160000000000000*x2;
		x7=0.0800000000000000*x2;
		x8=0.160000000000000*x0;
		x9=0.0800000000000000*x0;
		x10=((((x1)*(x9)))+((-1.00000000000000*(x3)*(x7))));
		x11=((((x1)*(x7)))+(((x3)*(x9))));
		eetrans[0]=((((x1)*(x8)))+(x8)+((-1.00000000000000*(x11)*(x4)))+((-1.00000000000000*(x3)*(x6)))+(((x10)*(x5))));
		eetrans[1]=((((x1)*(x6)))+(x6)+(((x11)*(x5)))+(((x3)*(x8)))+(((x10)*(x4))));
		eetrans[2] = 0;
	} else if (x_dim == 7) {
		Matrix<double,9,1> eerot;
		double x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47;
		x0=cos((double)dofs[0]);
		x1=cos((double)dofs[3]);
		x2=cos((double)dofs[1]);
		x3=cos((double)dofs[2]);
		x4=sin((double)dofs[1]);
		x5=sin((double)dofs[3]);
		x6=sin((double)dofs[2]);
		x7=sin((double)dofs[0]);
		x8=cos((double)dofs[4]);
		x9=sin((double)dofs[4]);
		x10=cos((double)dofs[5]);
		x11=sin((double)dofs[5]);
		x12=((double(1.00000000000000))*(x1));
		x13=((double(0.0450000000000000))*(x2));
		x14=((double(1.00000000000000))*(x5));
		x15=((double(0.0450000000000000))*(x5));
		x16=((double(0.300000000000000))*(x1));
		x17=((double(1.00000000000000))*(x2));
		x18=((x3)*(x4));
		x19=((x0)*(x4));
		x20=((x4)*(x7));
		x21=((x0)*(x2));
		x22=((x0)*(x6));
		x23=((x3)*(x7));
		x24=((x6)*(x7));
		x25=((x4)*(x6));
		x26=((x0)*(x3));
		x27=((double(1.00000000000000))*(x24));
		x28=((double(0.0450000000000000))*(x24));
		x29=((double(0.0450000000000000))*(x22));
		x30=((x12)*(x2));
		x31=((x14)*(x19));
		x32=((x12)*(x19));
		x33=((x13)*(x26));
		x34=((x14)*(x20));
		x35=((x12)*(x20));
		x36=((x13)*(x23));
		x37=((((x21)*(x3)))+(((double(-1.00000000000000))*(x27))));
		x38=((((x2)*(x23)))+(x22));
		x39=((((double(-1.00000000000000))*(x17)*(x24)))+(x26));
		x40=((((double(-1.00000000000000))*(x17)*(x26)))+(x27));
		x41=((((double(-1.00000000000000))*(x23)))+(((double(-1.00000000000000))*(x17)*(x22))));
		x42=((((double(-1.00000000000000))*(x22)))+(((double(-1.00000000000000))*(x17)*(x23))));
		x43=((((x14)*(x2)))+(((x12)*(x18))));
		x44=((x36)+(x29));
		x45=((((x25)*(x9)))+(((double(-1.00000000000000))*(x43)*(x8))));
		x46=((((x39)*(x9)))+(((x8)*(((((double(-1.00000000000000))*(x34)))+(((x1)*(x38))))))));
		x47=((((x8)*(((((double(-1.00000000000000))*(x31)))+(((x1)*(x37)))))))+(((x41)*(x9))));
		eerot[0]=((((x10)*(x47)))+(((x11)*(((((double(-1.00000000000000))*(x32)))+(((x40)*(x5))))))));
		eerot[1]=((((x9)*(((((double(-1.00000000000000))*(x12)*(x37)))+(x31)))))+(((x41)*(x8))));
		eerot[2]=((((x11)*(x47)))+(((x10)*(((((double(-1.00000000000000))*(x14)*(x40)))+(x32))))));
		eetrans[0]=((double(0.220000000000000))+(((x5)*(((((double(-0.300000000000000))*(x24)))+(((double(0.300000000000000))*(x21)*(x3)))))))+(((double(0.550000000000000))*(x19)))+(((double(-1.00000000000000))*(x28)))+(((x16)*(x19)))+(x33)+(((x1)*(((((double(-1.00000000000000))*(x33)))+(x28)))))+(((x15)*(x19))));
		eerot[3]=((((x11)*(((((double(-1.00000000000000))*(x35)))+(((x42)*(x5)))))))+(((x10)*(x46))));
		eerot[4]=((((x9)*(((((double(-1.00000000000000))*(x12)*(x38)))+(x34)))))+(((x39)*(x8))));
		eerot[5]=((((x10)*(((((double(-1.00000000000000))*(x14)*(x42)))+(x35)))))+(((x11)*(x46))));
		eetrans[1]=((double(0.140000000000000))+(((double(-1.00000000000000))*(x1)*(x44)))+(((double(0.550000000000000))*(x20)))+(((x5)*(((((double(0.300000000000000))*(x22)))+(((double(0.300000000000000))*(x2)*(x23)))))))+(((x16)*(x20)))+(x44)+(((x15)*(x20))));
		eerot[6]=((((x10)*(x45)))+(((x11)*(((((x18)*(x5)))+(((double(-1.00000000000000))*(x30))))))));
		eerot[7]=((((x43)*(x9)))+(((x25)*(x8))));
		eerot[8]=((((x11)*(x45)))+(((x10)*(((((double(-1.00000000000000))*(x14)*(x18)))+(x30))))));
		double x48=((double(0.0450000000000000))*(x18));
		eetrans[2]=((double(0.346000000000000))+(((double(-0.300000000000000))*(x18)*(x5)))+(((double(-1.00000000000000))*(x48)))+(((x16)*(x2)))+(((x1)*(x48)))+(((x13)*(x5)))+(((double(0.550000000000000))*(x2))));
	} else if (x_dim == 2) {
		eetrans[0] = dofs[0];
		eetrans[1] = dofs[1];
		eetrans[2] = 0;
	}
}

void BeliefRobotAndDOF::composeBelief(const VectorXd& x, const MatrixXd& rt_S, VectorXd& theta) {
	theta.resize(b_dim);
	theta.topRows(x_dim) = x;
	VectorXd rt_S_vec(s_dim);
	int idx = 0;
	for (int i=0; i<x_dim; i++) {
		for (int j=i; j<x_dim; j++) {
			rt_S_vec[idx] = 0.5 * (rt_S(i,j)+rt_S(j,i));
			idx++;
		}
	}

	theta.bottomRows(s_dim) = rt_S_vec;
}

void BeliefRobotAndDOF::decomposeBelief(const VectorXd& theta, VectorXd& x, MatrixXd& rt_S) {
	x = theta.topRows(x_dim);
	VectorXd rt_S_vec = theta.bottomRows(s_dim);

	int idx = 0;
	for (int j = 0; j < x_dim; ++j) {
		for (int i = j; i < x_dim; ++i) {
			rt_S(i,j) = rt_S(j,i) = rt_S_vec[idx];
			idx++;
		}
	}
}

void BeliefRobotAndDOF::SetBeliefValues(const DblVec& theta) {
	assert(theta.size() == b_dim);
	SetDOFValues(DblVec(theta.begin(), theta.begin()+x_dim));
	sigma_vec = DblVec(theta.begin()+x_dim, theta.end());
}

DblVec BeliefRobotAndDOF::GetBeliefValues() {
	DblVec theta = GetDOFValues();
	theta.insert(theta.end(), sigma_vec.begin(), sigma_vec.end());
	assert(theta.size() == b_dim);
	return theta;
}

MatrixXd BeliefRobotAndDOF::GetDynNoise() {
	MatrixXd diag_noise = MatrixXd::Identity(q_dim,q_dim);

	if (x_dim == 3)	{
		diag_noise(0,0) = 0.08; diag_noise(1,1) = 0.13; diag_noise(2,2) = 0.18;
	} else if (x_dim == 7) {
		//		double scale = 0.0025;
		double scale = 0.01;
		for(int i = 0; i < x_dim; ++i) {
			diag_noise(i,i) = scale;
		}
		//diag_noise(0,0) = 0.01; diag_noise(1,1) = 0.01; diag_noise(2,2) = 0.01;
		//diag_noise(3,3) = 0.01; diag_noise(4,4) = 0.01; diag_noise(5,5) = 0.01;
		//diag_noise(6,6) = 0.01;
	} else if (x_dim == 2) {
		diag_noise(0,0) = 0.01; diag_noise(1,1) = 0.01;
	}
	return diag_noise;
}

MatrixXd BeliefRobotAndDOF::GetObsNoise() {
	MatrixXd diag_noise = MatrixXd::Identity(z_dim, z_dim);
	if (x_dim == 3) {
		diag_noise(0,0) = 0.09; diag_noise(1,1) = 0.09; diag_noise(2,2) = 0.09;
	} else if (x_dim == 2) {
		diag_noise(0,0) = 0.01; diag_noise(1,1) = 0.01;
	} else if (x_dim == 7) {
		diag_noise(0,0) = 0.001; diag_noise(1,1) = 0.001; diag_noise(2,2) = 0.001;
	}
	return diag_noise;
}

VectorXd BeliefRobotAndDOF::Dynamics(const VectorXd& dofs, const VectorXd& u, const VectorXd& q) {
	return (dofs + u + GetDynNoise()*q);
}

MatrixXd BeliefRobotAndDOF::dfdx(const VectorXd& x, const VectorXd& u, const VectorXd& q) {
	MatrixXd J(x_dim,x_dim);
	VectorXd x_plus = x, x_minus = x;
	for (int i = 0; i < x_dim; ++i) {
		x_plus[i] += STEP; x_minus[i] -= STEP;
		J.col(i) = (Dynamics(x_plus, u, q) - Dynamics(x_minus, u, q)) / (2*STEP);
		x_plus[i] = x_minus[i] = x[i];
	}
	return J;
}

MatrixXd BeliefRobotAndDOF::dfdq(const VectorXd& x, const VectorXd& u, const VectorXd& q) {
	MatrixXd J(x_dim,q_dim);
	VectorXd q_plus = q, q_minus = q;
	for (int i = 0; i < q_dim; ++i) {
		q_plus[i] += STEP; q_minus[i] -= STEP;
		J.col(i) = (Dynamics(x, u, q_plus) - Dynamics(x, u, q_minus)) / (2*STEP);
		q_plus[i] = q_minus[i] = q[i];
	}
	return J;
}

double sigmoid(double x,double mean) {
	double y = (x - mean);
	double s = (y/sqrt(1+y*y))+1.0;

	if (x < mean)
		return s*0.1;
	else
		return s*7.0;
}

VectorXd BeliefRobotAndDOF::Observe(const VectorXd& x, const VectorXd& r) {
	//OR::RobotBase::RobotStateSaver saver = const_cast<BeliefRobotAndDOF*>(this)->Save();
	//SetDOFValues(toDblVec(x));
	//OR::Vector trans2 = endeffector->GetTransform().trans;

	Vector3d trans;
	ForwardKinematics(x, trans);

	VectorXd z(z_dim);
	if (x_dim == 3) {
		//z[0] = trans.x; z[1] = trans.y; z[2] = trans.z;
		double scale = (0.5*(trans[1]+0.2)*(trans[1]+0.2)+1);
		z[0] = trans[0] + scale*r[0]; z[1] = trans[1] + scale*r[1]; z[2] = trans[2] + scale*r[2];
	}  else if (x_dim == 7) {
		z[0] = 0; z[1] = 0;
		OR::Vector beacon;
		beacon.x = -1.0; beacon.y = 0; beacon.z = 0.5;
		//double dist = (trans.x - beacon.x)*(trans.x - beacon.x) + (trans.y - beacon.y)*(trans.y - beacon.y) + (trans.z - beacon.z)*(trans.z - beacon.z);
		//double dist = (trans[0] - beacon.x)*(trans[0] - beacon.x) + (trans[1] - beacon.y)*(trans[1] - beacon.y) + (trans[2] - beacon.z)*(trans[2] - beacon.z);
		double dist = 3.0*(trans[0] - beacon.x)*(trans[0] - beacon.x);
		z[0] = 1.0/(1.0 + dist) + 0.1*r[0];
		z[1] = x[0] + 0.01*r[1];
		z[2] = x[3] + 0.01*r[2];
	} else if (x_dim == 2) {
		// as in the Platt paper
		//z[0] = trans.x; z[1] = trans.y;
		z[0] = trans[0]; z[1] = trans[1];
		double scale = (0.5*(5.0 - trans[0])*(5.0 - trans[0])+0.001);
		z[0] += scale*r[0];
		z[1] += scale*r[1];
	}
	return z;
}

MatrixXd BeliefRobotAndDOF::dhdx(const VectorXd& x, const VectorXd& r) {
	MatrixXd J(z_dim,x_dim);
	VectorXd x_plus = x, x_minus = x;
	for (int i = 0; i < x_dim; ++i) {
		x_plus[i] += STEP; x_minus[i] -= STEP;
		J.col(i) = (Observe(x_plus, r) - Observe(x_minus, r)) / (2*STEP);
		x_plus[i] = x_minus[i] = x[i];
	}
	return J;
}

MatrixXd BeliefRobotAndDOF::dhdr(const VectorXd& x, const VectorXd& r) {
	MatrixXd J(z_dim,r_dim);
	VectorXd r_plus = r, r_minus = r;
	for (int i = 0; i < r_dim; ++i) {
		r_plus[i] += STEP; r_minus[i] -= STEP;
		J.col(i) = (Observe(x, r_plus) - Observe(x, r_minus)) / (2*STEP);
		r_plus[i] = r_minus[i] = r[i];
	}
	return J;
}


VectorXd BeliefRobotAndDOF::BeliefDynamics(const VectorXd& theta0, const VectorXd& u0) {
	VectorXd x0(x_dim), x(x_dim);
	MatrixXd rtSigma0(x_dim,x_dim), rtSigma(x_dim,x_dim);

	decomposeBelief(theta0, x0, rtSigma0);

	ekfUpdate(u0, x0, rtSigma0, x, rtSigma, false, VectorXd::Zero(z_dim));
	//ukfUpdate(u0, x0, rtSigma0, x, rtSigma);

	VectorXd theta(b_dim);
	composeBelief(x, rtSigma, theta);

	return theta;
}

MatrixXd BeliefRobotAndDOF::dgdb(const VectorXd& theta, const VectorXd& u) {
	MatrixXd J(b_dim,b_dim);
	VectorXd theta_plus = theta, theta_minus = theta;
	for (int i = 0; i < b_dim; ++i) {
		theta_plus[i] += STEP; theta_minus[i] -= STEP;
		J.col(i) = (BeliefDynamics(theta_plus, u) - BeliefDynamics(theta_minus, u)) / (2*STEP);
		theta_plus[i] = theta_minus[i] = theta[i];
	}
	return J;
}

MatrixXd BeliefRobotAndDOF::dgdu(const VectorXd& theta, const VectorXd& u) {
	MatrixXd J(b_dim,u_dim);
	VectorXd u_plus = u, u_minus = u;
	for (int i = 0; i < u_dim; ++i) {
		u_plus[i] += STEP; u_minus[i] -= STEP;
		J.col(i) = (BeliefDynamics(theta, u_plus) - BeliefDynamics(theta, u_minus)) / (2*STEP);
		u_plus[i] = u_minus[i] = u[i];
	}
	return J;
}

VectorXd BeliefRobotAndDOF::VectorXdRand(int size) {
	VectorXd v(size);
	for (int i=0; i<size; i++) {
		v(i) = generator();
		// hack!!!
		if (v(i) > 2) v(i) = 2;
		if (v(i) < -2) v(i) = -2;
	}
	return v;
}

VectorXd BeliefRobotAndDOF::mvnrnd(const VectorXd& mu, const MatrixXd& Sigma) {
	Eigen::SelfAdjointEigenSolver<MatrixXd> eigenSolver(Sigma);
	return mu + eigenSolver.eigenvectors() * eigenSolver.eigenvalues().cwiseSqrt().asDiagonal() * VectorXdRand(x_dim);
}

MatrixXd BeliefRobotAndDOF::sigmaPoints(const VectorXd& theta) {
	VectorXd x(x_dim);
	MatrixXd rtSigma(x_dim,x_dim);
	decomposeBelief(theta, x, rtSigma);
	return sigmaPoints(x, rtSigma);
}

MatrixXd BeliefRobotAndDOF::sigmaPoints(const VectorXd& mean, const MatrixXd& sqrtcov)
{
	int L = x_dim + r_dim;

	double lambda = alpha*alpha*(L + kappa) - L;
	//	double w = 1 / (2*(L + lambda));
	//	double mw = lambda / (L + lambda);
	//	double vw = mw + (1 - alpha*alpha + beta);

	MatrixXd sigmapts(x_dim, 2*x_dim+1);
	sigmapts.col(0) = mean;

	//	double scale = sqrt(L+lambda);
	for(int i = 0; i < x_dim; ++i){
		sigmapts.col(2*i+1) = (mean + sigma_pts_scale*sqrtcov.col(i));
		sigmapts.col(2*i+2) = (mean - sigma_pts_scale*sqrtcov.col(i));
	}

	return sigmapts;
}

void BeliefRobotAndDOF::ekfUpdate(const VectorXd& u0, const VectorXd& x0, const MatrixXd& rtSigma0, VectorXd& x, MatrixXd& rtSigma, bool observe, const VectorXd& z) {
	VectorXd q = VectorXd::Zero(q_dim);
	VectorXd x_pred = Dynamics(x0, u0, q);

	MatrixXd A = dfdx(x0, u0, q);
	MatrixXd Q = dfdq(x0, u0, q);
	MatrixXd Sigma_pred = A * rtSigma0 * rtSigma0.transpose() * A.transpose() + Q*Q.transpose();

	VectorXd r = VectorXd::Zero(r_dim);
	MatrixXd C = dhdx(x_pred, r);
	MatrixXd R = dhdr(x_pred, r);

	PartialPivLU<MatrixXd> solver(C*Sigma_pred*C.transpose() + R*R.transpose());
	MatrixXd K = solver.solve((Sigma_pred*C.transpose()).transpose()).transpose();

	x = x_pred;
	if (observe) { x += K*(z-Observe(x_pred, r)); }

	//Eigen::JacobiSVD<MatrixXd, QRPreconditioner> svd(Sigma_pred - K*C*Sigma_pred, ComputeThinU | ComputeThinV);
	//rtSigma = svd.matrixU() * svd.singularValues().array().sqrt().matrix().asDiagonal() * svd.matrixV().transpose();

	SelfAdjointEigenSolver<MatrixXd> es(Sigma_pred - K*C*Sigma_pred);
	rtSigma = es.eigenvectors().real() * es.eigenvalues().real().cwiseSqrt().asDiagonal() * es.eigenvectors().real().transpose();
}

void BeliefRobotAndDOF::ukfUpdate(const VectorXd& u0, const VectorXd& x0, const MatrixXd& rtSigma0, VectorXd& x, MatrixXd& rtSigma, bool observe, const VectorXd& z)
{
	int L = x_dim + q_dim + r_dim;

	double ukfalpha = 0.001; //alpha;
	double ukfbeta = 2; //beta;
	double ukfkappa = 0; //kappa;

	double lambda = ukfalpha*ukfalpha*(L + ukfkappa) - L;
	double w = 1 / (2*(L + lambda));
	double mw = lambda / (L + lambda);
	double vw = mw + (1 - ukfalpha*ukfalpha + ukfbeta);

	// propagate sigma points through f
	VectorXd q = VectorXd::Zero(q_dim);
	MatrixXd sigmapts(x_dim, 2*(x_dim+q_dim)+1);

	int idx = 0;
	sigmapts.col(idx++) = Dynamics(x0, u0, q);

	MatrixXd rt_scaled_Sigma0 = sqrt(L+lambda)*rtSigma0;

	for (int i = 0; i < x_dim; ++i) {
		sigmapts.col(idx++) = Dynamics(x0 + rt_scaled_Sigma0.col(i), u0, q);
		sigmapts.col(idx++) = Dynamics(x0 - rt_scaled_Sigma0.col(i), u0, q);
	}

	double qstep = sqrt(L + lambda);
	for (int i = 0; i < q_dim; ++i) {
		q(i) = qstep;
		sigmapts.col(idx++) = Dynamics(x0, u0, q);
		q(i) = -qstep;
		sigmapts.col(idx++) = Dynamics(x0, u0, q);
		q(i) = 0;
	}

	// calculate mean -- O(xDim^2)
	VectorXd xhat(x_dim);
	xhat = (mw + 2*r_dim*w) * sigmapts.col(0);
	for (int i = 1; i < sigmapts.cols(); ++i) {
		xhat += w * sigmapts.col(i);
	}

	// calculate variance -- O(xDim^3)
	MatrixXd Sigma = (vw + 2*r_dim*w) * (sigmapts.col(0) - xhat)*(sigmapts.col(0) - xhat).transpose();
	for (int i = 1; i < sigmapts.cols(); ++i) {
		Sigma += w * (sigmapts.col(i) - xhat)*(sigmapts.col(i) - xhat).transpose();
	}

	// Measurement Update
	MatrixXd Z(z_dim, 2*(x_dim+q_dim+r_dim)+1);
	VectorXd r = VectorXd::Zero(r_dim);

	idx = 0;
	for (int i = 0; i < sigmapts.cols(); ++i) {
		Z.col(idx++) = Observe(sigmapts.col(i), r);
	}

	double rstep = sqrt(L + lambda);
	for (int i = 0; i < r_dim; ++i) {
		r(i) = rstep;
		Z.col(idx++) = Observe(sigmapts.col(0), r);
		r(i) = -rstep;
		Z.col(idx++) = Observe(sigmapts.col(0), r);
		r(i) = 0;
	}

	// calculate mean -- O(xDim*zDim + zDim^2)
	VectorXd zhat(z_dim);
	zhat = mw * Z.col(0);
	for (int i = 1; i < Z.cols(); ++i) {
		zhat += w * Z.col(i);
	}

	// calculate variance -- O(zDim^3 + xDim*zDim^2)
	MatrixXd Pzz = vw * (Z.col(0) - zhat)*(Z.col(0) - zhat).transpose();
	for (int i = 1; i < Z.cols(); ++i) {
		Pzz += w * (Z.col(i) - zhat)*(Z.col(i) - zhat).transpose();
	}

	// calculate cross-covariance -- O(xDim^2*zDim + xDim*zDim^2)
	MatrixXd Pxz = vw * (sigmapts.col(0) - xhat)*(Z.col(0) - zhat).transpose();
	for (int i = 1; i < sigmapts.cols(); ++i) {
		Pxz += w * (sigmapts.col(i) - xhat)*(Z.col(i) - zhat).transpose();
	}
	for (int i = sigmapts.cols(); i < Z.cols(); ++i) {
		Pxz += w * (sigmapts.col(0) - xhat)*(Z.col(i) - zhat).transpose();
	}

	PartialPivLU<MatrixXd> solver(Pzz);
	MatrixXd K = solver.solve(Pxz); // Pxz/Pzz

	x = xhat;
	if (observe) {
		x += K*(Z.col(0) - z);
	}

	Sigma -= Pxz*K.transpose();

	SelfAdjointEigenSolver<MatrixXd> es(Sigma);
	rtSigma = es.eigenvectors().real() * es.eigenvalues().real().cwiseSqrt().asDiagonal() * es.eigenvectors().real().transpose();
}

// theta needs to be set accordingly before calling this (i.e. call SetBeliefValues)
MatrixXd BeliefRobotAndDOF::BeliefJacobian(int link_ind, int sigma_pt_ind, const OR::Vector& pt) {
	MatrixXd pos_jac = PositionJacobian(endeffector->GetIndex(), pt);
	MatrixXd jac = MatrixXd::Zero(3,b_dim);
	jac.leftCols(x_dim) = pos_jac;

	double lambda;
	if (sigma_pt_ind == 0) lambda = 0;
	else lambda = ((sigma_pt_ind-1)%2)==0 ? 1 : -1; // TODO non-unit magnitude for lambda

	if (lambda != 0) {
		int sigma_matrix_col = (sigma_pt_ind - 1)/2;
		vector<int> sigma_vec_inds = sigmaColToSigmaIndices(sigma_matrix_col);
		assert(sigma_vec_inds.size() == x_dim);
		for (int i=0; i<sigma_vec_inds.size(); i++) {
			int jac_col = x_dim + sigma_vec_inds[i];
			jac.col(jac_col) = lambda * pos_jac.col(i);
		}
	}
	return jac;
}

void BeliefRobotAndDOF::GetEndEffectorNoiseAsGaussian(const VectorXd& theta, Vector3d& mean, Matrix3d& cov) {
	VectorXd x(x_dim);
	MatrixXd rt_Sigma(x_dim,x_dim);
	decomposeBelief(theta, x, rt_Sigma);

	OR::RobotBase::RobotStateSaver saver = const_cast<BeliefRobotAndDOF*>(this)->Save();
	SetDOFValues(toDblVec(x));
	OR::Vector trans = endeffector->GetTransform().trans;

	MatrixXd jac = PositionJacobian(endeffector->GetIndex(), trans);
	mean = Vector3d(trans.x, trans.y, trans.z);
	cov = jac * rt_Sigma * rt_Sigma.transpose() * jac.transpose();
}

osg::Matrix gaussianAsTransform(const Vector3d& mean, const Matrix3d& cov) {
	EigenSolver<Matrix3d> es(cov);
	Matrix4d t = Matrix4d::Identity();
	t.block(0,0,3,3) = es.eigenvectors().real() * es.eigenvalues().real().cwiseSqrt().asDiagonal();
	t.block(0,3,3,1) = mean;
	Matrix4d t_transpose = t.transpose();

	osg::Matrix osg_t;
	osg_t.set(t_transpose.data());
	return osg_t;
}

}
