#pragma once
#include "trajopt/robot_and_dof.hpp"

#include "osgviewer/osgviewer.hpp"
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>

#include "utils/eigen_conversions.hpp"

namespace trajopt {

class BeliefRobotAndDOF : public RobotAndDOF {
private:
  boost::variate_generator<boost::mt19937, boost::normal_distribution<> > generator;
  int n_theta;
  DblVec sigma_vec;
public:
  BeliefRobotAndDOF(OpenRAVE::RobotBasePtr _robot, const IntVec& _joint_inds, int _affinedofs=0, const OR::Vector _rotationaxis=OR::Vector());
  Eigen::Vector3d EndEffectorPosition(const Eigen::VectorXd& dofs) {
		OR::RobotBase::RobotStateSaver saver = const_cast<BeliefRobotAndDOF*>(this)->Save();
		SetDOFValues(util::toDblVec(dofs));
		OR::Vector trans = link->GetTransform().trans;
		return Eigen::Vector3d(trans.x, trans.y, trans.z);
  }

  void SetBeliefValues(const DblVec& theta);
  DblVec GetBeliefValues();

  Eigen::MatrixXd GetDynNoise();
  inline int GetQSize() { return GetDOF(); }
//  inline int GetQSize() { return GetDynNoise().cols(); }

  Eigen::MatrixXd GetObsNoise();
  inline int GetRSize() { return GetDOF(); }
  inline int GetObsSize() { return GetDOF(); }
//  inline int GetRSize() { return GetObsNoise().cols(); }
//  inline int GetObsSize() { return GetObsNoise().rows(); }

  // UKF update functions
  Eigen::MatrixXd sigmaPoints(const VectorXd& theta);
  Eigen::MatrixXd sigmaPoints(const VectorXd& mean, const Eigen::MatrixXd& cov);
  VectorXd sigmaPoint(const VectorXd& mean, const Eigen::MatrixXd& cov, int idx);
  void ukfUpdate(const VectorXd& u0, const VectorXd& xest0, const Eigen::MatrixXd& Vest0, VectorXd& xest, Eigen::MatrixXd& Vest);

  Eigen::VectorXd Observe(const Eigen::VectorXd& dofs, const Eigen::VectorXd& r);
  Eigen::VectorXd Dynamics(const Eigen::VectorXd& dofs, const Eigen::VectorXd& u, const Eigen::VectorXd& q);
  Eigen::VectorXd BeliefDynamics(const Eigen::VectorXd& theta0, const Eigen::VectorXd& u0);
  Eigen::VectorXd VectorXdRand(int size);

  template <typename T> Eigen::Matrix<T,Eigen::Dynamic,1> toSigmaVec(const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>& rt_S) {
  	Eigen::Matrix<T,Eigen::Dynamic,1> rt_S_vec(GetNTheta()-GetDOF());
		int idx = 0;
		for (int i=0; i<GetDOF(); i++) {
			for (int j=i; j<GetDOF(); j++) {
				rt_S_vec(idx) = 0.5 * (rt_S(i,j)+rt_S(j,i));
				idx++;
			}
		}
		return rt_S_vec;
	}
  template <typename T> Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> toSigmaMatrix(const Eigen::Matrix<T,Eigen::Dynamic,1>& rt_S_vec) {
  	Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> rt_S(GetDOF(), GetDOF());
		int idx = 0;
		for (int j = 0; j < GetDOF(); ++j) {
			for (int i = j; i < GetDOF(); ++i) {
				rt_S(i,j) = rt_S(j,i) = rt_S_vec(idx);
				idx++;
			}
		}
		return rt_S;
	}
  template <typename T> void composeBelief(const Eigen::Matrix<T,Eigen::Dynamic,1>& x, const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>& rt_S,
  		Eigen::Matrix<T,Eigen::Dynamic,1>& theta) {
  	theta.resize(GetNTheta());
  	theta.topRows(GetDOF()) = x;
  	theta.bottomRows(GetNTheta()-GetDOF()) = toSigmaVec(rt_S);
  }
  template <typename T> void decomposeBelief(const Eigen::Matrix<T,Eigen::Dynamic,1>& theta,
  		Eigen::Matrix<T,Eigen::Dynamic,1>& x, Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>& rt_S) {
  	x = theta.topRows(GetDOF());
  	rt_S = toSigmaMatrix((Eigen::Matrix<T,Eigen::Dynamic,1>) theta.bottomRows(GetNTheta()-GetDOF()));
  }

  void ekfUpdate(const Eigen::VectorXd& u0, const Eigen::VectorXd& x0, const Eigen::MatrixXd& rtSigma0, VectorXd& xest, Eigen::MatrixXd& Vest);

  void GetEndEffectorNoiseAsGaussian(const Eigen::VectorXd& theta, Eigen::VectorXd& mean, Eigen::MatrixXd& cov);

  // theta needs to be set accordingly before calling this (i.e. call SetBeliefValues)
  Eigen::MatrixXd BeliefJacobian(int link_ind, int sigma_pt_ind, const OR::Vector& pt);

  inline int GetNTheta() { return n_theta; }

	OR::KinBody::LinkPtr link;
	// Scaled UKF update vars
	double alpha, beta, kappa;

private:
	// returns the indices of the terms in Sigma_vec that are in the jth column of the corresponding Sigma_matrix
	inline std::vector<int> sigmaColToSigmaIndices(int j) { return sigma_col_to_indices[j]; }
	std::vector<std::vector<int> > sigma_col_to_indices;
};
typedef boost::shared_ptr<BeliefRobotAndDOF> BeliefRobotAndDOFPtr;

typedef boost::function<VectorXd(VectorXd)> VectorOfVectorFun;
typedef boost::shared_ptr<VectorOfVectorFun> VectorOfVectorFunPtr;
Eigen::MatrixXd calcNumJac(VectorOfVectorFun f, const VectorXd& x, double epsilon=0.00048828125);
osg::Matrix gaussianAsTransform(const Eigen::Vector3d& mean, const Eigen::Matrix3d& cov);

}
