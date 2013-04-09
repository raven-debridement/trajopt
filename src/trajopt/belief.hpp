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

#define STEP 0.00048828125

class TRAJOPT_API BeliefRobotAndDOF : public RobotAndDOF {
public:
	BeliefRobotAndDOF(OpenRAVE::RobotBasePtr _robot, const IntVec& _joint_inds, int _affinedofs=0, const OR::Vector _rotationaxis=OR::Vector());

	void ForwardKinematics(const VectorXd& dofs, Vector3d& eetrans);

	void SetBeliefValues(const DblVec& theta);
	DblVec GetBeliefValues();

	MatrixXd GetDynNoise();
	MatrixXd GetObsNoise();

	VectorXd Observe(const VectorXd& x, const VectorXd& r);
	VectorXd Dynamics(const VectorXd& x, const VectorXd& u, const VectorXd& q);

	MatrixXd dfdx(const VectorXd& x, const VectorXd& u, const VectorXd& q);
	MatrixXd dfdq(const VectorXd& x, const VectorXd& u, const VectorXd& q);
	MatrixXd dhdx(const VectorXd& x, const VectorXd& r);
	MatrixXd dhdr(const VectorXd& x, const VectorXd& r);

	VectorXd BeliefDynamics(const VectorXd& theta0, const VectorXd& u0);
	MatrixXd dgdb(const VectorXd& theta, const VectorXd& u);
	MatrixXd dgdu(const VectorXd& theta, const VectorXd& u);

	VectorXd VectorXdRand(int size);
	VectorXd mvnrnd(const VectorXd& mu, const MatrixXd& Sigma);

	template <typename T> Matrix<T,Dynamic,1> toSigmaVec(const Matrix<T,Dynamic,Dynamic>& rt_S) {
		Matrix<T,Dynamic,1> rt_S_vec(s_dim);
		int idx = 0;
		for (int i=0; i < x_dim; i++) {
			for (int j=i; j < x_dim; j++) {
				rt_S_vec[idx] = 0.5 * (rt_S(i,j)+rt_S(j,i));
				idx++;
			}
		}
		return rt_S_vec;
	}

	template <typename T> Matrix<T,Dynamic,Dynamic> toSigmaMatrix(const Matrix<T,Dynamic,1>& rt_S_vec) {
		Matrix<T,Dynamic,Dynamic> rt_S(x_dim, x_dim);
		int idx = 0;
		for (int j = 0; j < x_dim; ++j) {
			for (int i = j; i < x_dim; ++i) {
				rt_S(i,j) = rt_S(j,i) = rt_S_vec[idx];
				idx++;
			}
		}
		return rt_S;
	}

	void composeBelief(const VectorXd& x, const MatrixXd& rt_S, VectorXd& theta);
	void decomposeBelief(const VectorXd& theta, VectorXd& x, MatrixXd& rt_S);

	void ekfUpdate(const VectorXd& u0, const VectorXd& x0, const MatrixXd& rtSigma0, VectorXd& x, MatrixXd& rtSigma, bool observe, const VectorXd& z);

	MatrixXd sigmaPoints(const VectorXd& theta);
	MatrixXd sigmaPoints(const VectorXd& mean, const MatrixXd& sqrtcov);
	VectorXd sigmaPoint(const VectorXd& mean, const MatrixXd& cov, int idx);

	void ukfUpdate(const VectorXd& u0, const VectorXd& x0, const MatrixXd& rtSigma0, VectorXd& x, MatrixXd& rtSigma, bool observe, const VectorXd& z);

	void GetEndEffectorNoiseAsGaussian(const VectorXd& theta, Vector3d& mean, Matrix3d& cov);

	// theta needs to be set accordingly before calling this (i.e. call SetBeliefValues)
	MatrixXd BeliefJacobian(int link_ind, int sigma_pt_ind, const OR::Vector& pt);

	void SetSigmaPointsScale(double scale) { sigma_pts_scale = scale; }

	inline int GetXDim() { return x_dim; }
	inline int GetBDim() { return b_dim; }
	inline int GetSDim() { return s_dim; }
	inline int GetQDim() { return q_dim; }
	inline int GetZDim() { return z_dim; }
	inline int GetRDim() { return r_dim; }
	inline int GetUDim() { return u_dim; }

	OR::KinBody::LinkPtr endeffector;
	// Scaled UKF update vars
	double alpha, beta, kappa;

private:
	double sigma_pts_scale;
	// returns the indices of the terms in Sigma_vec that are in the jth column of the corresponding Sigma_matrix
	inline vector<int> sigmaColToSigmaIndices(int j) { return sigma_col_to_indices[j]; }
	vector<vector<int> > sigma_col_to_indices;
	boost::variate_generator<boost::mt19937, boost::normal_distribution<> > generator;
	DblVec sigma_vec;
	int x_dim, b_dim, s_dim, z_dim, q_dim, r_dim, u_dim;
};

typedef boost::shared_ptr<BeliefRobotAndDOF> BeliefRobotAndDOFPtr;

osg::Matrix gaussianAsTransform(const Vector3d& mean, const Matrix3d& cov);

}
