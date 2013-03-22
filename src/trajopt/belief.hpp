#pragma once
#include "trajopt/robot_and_dof.hpp"

#include "osgviewer/osgviewer.hpp"
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>

#include "utils/eigen_conversions.hpp"

using namespace std;
using namespace Eigen;
using namespace util;

namespace trajopt {

// hack to set parameter sizes to avoid matrix allocations at runtime
//#define X_DIM 2
//#define U_DIM 2
//#define Z_DIM 2
//#define Q_DIM 2
//#define R_DIM 2

#define X_DIM 3
#define U_DIM 3
#define Z_DIM 3
#define Q_DIM 3
#define R_DIM 3

//#define X_DIM 7
//#define U_DIM 7
//#define Z_DIM 3
//#define Q_DIM 7
//#define R_DIM 3

#define S_DIM (X_DIM*(X_DIM+1))/2
#define B_DIM (X_DIM + S_DIM)

#define N_STEPS 10

#define STEP 0.00048828125

typedef Matrix<double,X_DIM,X_DIM> MatXXd;
typedef Matrix<double,X_DIM,Q_DIM> MatXQd;
typedef Matrix<double,Z_DIM,Z_DIM> MatZZd;
typedef Matrix<double,Z_DIM,X_DIM> MatZXd;
typedef Matrix<double,X_DIM,Z_DIM> MatXZd;
typedef Matrix<double,Z_DIM,R_DIM> MatZRd;
typedef Matrix<double,B_DIM,B_DIM> MatBBd;
typedef Matrix<double,B_DIM,U_DIM> MatBUd;

typedef Matrix<double,X_DIM,1> VecXd;
typedef Matrix<double,U_DIM,1> VecUd;
typedef Matrix<double,S_DIM,1> VecSd;
typedef Matrix<double,Z_DIM,1> VecZd;
typedef Matrix<double,U_DIM,1> VecUd;
typedef Matrix<double,Q_DIM,1> VecQd;
typedef Matrix<double,R_DIM,1> VecRd;
typedef Matrix<double,B_DIM,1> VecBd;

typedef Matrix<int,X_DIM,X_DIM> MatXXi;
typedef Matrix<int,S_DIM,1> VecSi;

class TRAJOPT_API BeliefRobotAndDOF : public RobotAndDOF {
private:
	boost::variate_generator<boost::mt19937, boost::normal_distribution<> > generator;
	DblVec sigma_vec;
public:

	BeliefRobotAndDOF(OpenRAVE::RobotBasePtr _robot, const IntVec& _joint_inds, int _affinedofs=0, const OR::Vector _rotationaxis=OR::Vector());

	void ForwardKinematics(const VecXd& dofs, Vector3d& eetrans);

	void SetBeliefValues(const DblVec& theta);
	DblVec GetBeliefValues();

	MatXXd GetDynNoise();
	MatZZd GetObsNoise();

	VecZd Observe(const VecXd& x, const VecRd& r);
	VecXd Dynamics(const VecXd& x, const VecUd& u, const VecQd& q);

	MatXXd dfdx(const VecXd& x, const VecUd& u, const VecQd& q);
	MatXXd dfdq(const VecXd& x, const VecUd& u, const VecQd& q);
	MatZXd dhdx(const VecXd& x, const VecRd& r);
	MatZZd dhdr(const VecXd& x, const VecRd& r);

	VecBd BeliefDynamics(const VecBd& theta0, const VecUd& u0);
	MatBBd dgdb(const VecBd& theta, const VecUd& u);
	MatBUd dgdu(const VecBd& theta, const VecUd& u);

	VectorXd VectorXdRand(int size);

	template <typename T> Eigen::Matrix<T,Eigen::Dynamic,1> toSigmaVec(const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>& rt_S) {
		Eigen::Matrix<T,Eigen::Dynamic,1> rt_S_vec(S_DIM);
		int idx = 0;
		for (int i=0; i<X_DIM; i++) {
			for (int j=i; j<X_DIM; j++) {
				rt_S_vec[idx] = 0.5 * (rt_S(i,j)+rt_S(j,i));
				idx++;
			}
		}
		return rt_S_vec;
	}

	template <typename T> Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> toSigmaMatrix(const Eigen::Matrix<T,Eigen::Dynamic,1>& rt_S_vec) {
		Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> rt_S(X_DIM, X_DIM);
		int idx = 0;
		for (int j = 0; j < X_DIM; ++j) {
			for (int i = j; i < X_DIM; ++i) {
				rt_S(i,j) = rt_S(j,i) = rt_S_vec[idx];
				idx++;
			}
		}
		return rt_S;
	}

	void composeBelief(const VecXd& x, const MatXXd& rt_S, VecBd& theta);
	void decomposeBelief(const VecBd& theta, VecXd& x, MatXXd& rt_S);

	void ekfUpdate(const VecUd& u0, const VecXd& x0, const MatXXd& rtSigma0, VecXd& x, MatXXd& rtSigma, bool observe, const VecZd& z);

	MatrixXd sigmaPoints(const VecBd& theta);
	MatrixXd sigmaPoints(const VecXd& mean, const MatXXd& sqrtcov);
	VecXd sigmaPoint(const VecXd& mean, const MatXXd& cov, int idx);

	void ukfUpdate(const VecUd& u0, const VecXd& x0, const MatXXd& rtSigma0, VecXd& x, MatXXd& rtSigma, bool observe, const VecZd& z);

	void GetEndEffectorNoiseAsGaussian(const VecBd& theta, Vector3d& mean, Matrix3d& cov);

	// theta needs to be set accordingly before calling this (i.e. call SetBeliefValues)
	Matrix<double, 3, B_DIM> BeliefJacobian(int link_ind, int sigma_pt_ind, const OR::Vector& pt);

	OR::KinBody::LinkPtr endeffector;
	// Scaled UKF update vars
	double alpha, beta, kappa;

private:
	// returns the indices of the terms in Sigma_vec that are in the jth column of the corresponding Sigma_matrix
	inline std::vector<int> sigmaColToSigmaIndices(int j) { return sigma_col_to_indices[j]; }
	std::vector<std::vector<int> > sigma_col_to_indices;
};

typedef boost::shared_ptr<BeliefRobotAndDOF> BeliefRobotAndDOFPtr;

osg::Matrix gaussianAsTransform(const Vector3d& mean, const Matrix3d& cov);

}
