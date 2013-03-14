#include "trajopt/common.hpp"
#include <openrave-core.h>
#include "trajopt/collision_checker.hpp"
#include "osgviewer/osgviewer.hpp"
#include "utils/eigen_conversions.hpp"
#include "trajopt/utils.hpp"
#include <boost/algorithm/string.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/foreach.hpp>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Geometry>
using namespace OpenRAVE;
using namespace std;
using namespace trajopt;
using namespace util;
using namespace Eigen;

typedef Matrix<double,6,1> Vector6d;
typedef Matrix<double,6,Eigen::Dynamic> Matrix6Xd;
typedef Matrix<double,6,6> Matrix6d;
typedef Matrix<double,12,1> Vector12d;
typedef Matrix<double,12,Eigen::Dynamic> Matrix12Xd;
typedef Matrix<double,12,12> Matrix12d;

CollisionCheckerPtr cc;
vector<GraphHandlePtr> handles;
EnvironmentBasePtr env;
OSGViewerPtr viewer;

// copied from rave_utils
RobotBase::ManipulatorPtr GetManipulatorByName(RobotBase& robot, const std::string& name) {
  vector<RobotBase::ManipulatorPtr> manips = robot.GetManipulators();
  BOOST_FOREACH(RobotBase::ManipulatorPtr& manip, manips) {
    if (manip->GetName()==name) return manip;
  }
  return RobotBase::ManipulatorPtr();
}

// copied from problem_description
BeliefRobotAndDOFPtr RADFromName(const string& name, RobotBasePtr robot) {
  if (name == "active") {
    return BeliefRobotAndDOFPtr(new BeliefRobotAndDOF(robot, robot->GetActiveDOFIndices(), robot->GetAffineDOF(), robot->GetAffineRotationAxis()));
  }
  vector<int> dof_inds;
  int affinedofs = 0;
  Vector rotationaxis(0,0,1);
  vector<string> components;
  boost::split(components, name, boost::is_any_of("+"));
  for (int i=0; i < components.size(); ++i) {
    std::string& component = components[i];
    if (RobotBase::ManipulatorPtr manip = GetManipulatorByName(*robot, component)) {
      vector<int> inds = manip->GetArmIndices();
      dof_inds.insert(dof_inds.end(), inds.begin(), inds.end());
    }
    else if (component == "base3d") {
    	affinedofs |= DOF_XYZ | DOF_Rotation3D;
    }
    else if (component == "base") {
      affinedofs |= DOF_X | DOF_Y | DOF_RotationAxis;
    }
    else if (component == "base_point") {
      affinedofs |= DOF_X | DOF_Y;
    }
    else if (KinBody::JointPtr joint = robot->GetJoint(component)) {
      dof_inds.push_back(joint->GetDOFIndex());
    }
    else PRINT_AND_THROW( boost::format("error in reading manip description: %s must be a manipulator, link, or 'base'")%component );
  }
  return BeliefRobotAndDOFPtr(new BeliefRobotAndDOF(robot, dof_inds, affinedofs, rotationaxis));
}

float alpha = 1;
void AdjustTransparency(float da) {
  alpha += da;
  alpha = fmin(alpha, 1);
  alpha = fmax(alpha, 0);
  viewer->SetAllTransparency(alpha);
}

bool isTrajectoryInCollision(CollisionCheckerPtr cc, TrajArray traj, BeliefRobotAndDOFPtr rad) {
	vector<Collision> collisions;
	cc->DiscreteCheckTrajectory(traj, rad, collisions);
	for (int i=0; i<collisions.size(); i++) {
		if (collisions[i].distance < 0) return true;
	}
	return false;
}

OR::Vector toQuat(const Matrix3d& m) {
	Quaterniond rq(m);
	return OR::Vector(rq.w(), rq.x(), rq.y(), rq.z());
}
OR::Vector toQuat(const AngleAxisd& aa) {
	return toQuat(aa.toRotationMatrix());
}
Eigen::Matrix3d toMatrix3d(OR::RaveVector<float> rq) {
	return Quaterniond(rq[0], rq[1], rq[2], rq[3]).toRotationMatrix();
}

void renderSigmaPts(BeliefRobotAndDOFPtr rad, const MatrixXd& sigma_pts, const osg::Vec4f& colorvec) {
  vector<KinBody::LinkPtr> links;
	vector<int> joint_inds;
	rad->GetAffectedLinks(links, true, joint_inds);

	// render sigma points
	for (int j=0; j<sigma_pts.cols(); j++) {
		cout << sigma_pts.col(j).transpose() << endl;
		rad->SetDOFValues(toDblVec(sigma_pts.col(j)));
		handles.push_back(viewer->PlotKinBody(rad->GetRobot()));
		if (j==0) SetColor(handles.back(), osg::Vec4f(0,0,1,1));
		else //SetColor(handles.back(), osg::Vec4f(0,0,1,1));
			SetColor(handles.back(), colorvec);
	}

	// render convex hulls of sigma points
	vector<DblVec> dofvals(sigma_pts.cols());
	for (int i=0; i<sigma_pts.cols(); i++)
		dofvals[i] = toDblVec(sigma_pts.col(i));
	cc->SetContactDistance(100);
//	cc->PlotCastHull(*rad, links, dofvals, handles);
}

Matrix3d skewSymmetric(const Vector3d& v) {
	Matrix3d m;
	m << 0, -v(2), v(1), v(2), 0, -v(0), -v(1), v(0), 0;
	return m;
}

//VectorXd psi0, VectorXd psi1,
// psi0 is robot
// psi1 is obstacle
Vector6d delta(const Vector6d& psi0, const Vector6d& psi1) {
	Vector3d pos0 = psi0.topRows(3);
	Vector3d pos1 = psi1.topRows(3);
	Matrix3d rot0, rot1;
	if (psi0.bottomRows(3).norm() > 1e-5)
		rot0 = (Matrix3d) AngleAxisd(psi0.bottomRows(3).norm(), psi0.bottomRows(3).normalized());
	else
		rot0 = Matrix3d::Identity();
	if (psi1.bottomRows(3).norm() > 1e-5)
		rot1 = (Matrix3d) AngleAxisd(psi1.bottomRows(3).norm(), psi1.bottomRows(3).normalized());
	else
		rot1 = Matrix3d::Identity();

	Vector6d rel;
	rel.topRows(3) = rot1 * (pos0 - pos1);
	Vector3d rel_rot = Vector3d::Zero();
	for (int j=0; j<3; j++) {
		rel_rot += skewSymmetric(rot0.col(j)) * rot1.col(j);
	}
	rel.bottomRows(3) = 0.5 * rel_rot;

	return rel;
}


// psi1 robot
// psi2 obstacle
Matrix6Xd relativeSigmaPts(const Matrix6Xd& psi1, const Matrix6Xd& psi2) {
	Vector6d mean = delta(psi1.col(0), psi2.col(0));

	Matrix6Xd sigmapts(6,13);
	sigmapts.col(0) = mean;

	for(int i = 0; i < 6; ++i) {
		sigmapts.col(2*i+1) = (mean + (delta(psi1.col(2*i+1), psi1.col(0)) + delta(psi2.col(2*i+1), psi2.col(0))));
		sigmapts.col(2*i+2) = (mean + (delta(psi1.col(2*i+2), psi1.col(0)) + delta(psi2.col(2*i+2), psi2.col(0))));
	}

	return sigmapts;
}

Matrix6Xd jointSigmaPts(const Vector12d& mean, const Matrix12d& rt_Sigma) {
	int n_dof = 12;
//	int n_r = 6;
//	int L = n_dof + n_r;
	int L = 12;
	double alpha = 0.5;
	double kappa = 5.0;
	double lambda = alpha*alpha*(L + kappa) - L;
	cout << "lambda" << endl;
	cout << lambda << endl;

	MatrixXd rt_scaled_cov = sqrt(L + lambda) * rt_Sigma;

	Matrix6d cov = Matrix6d::Zero();
	Matrix6d var1 = Matrix6d::Zero();
	Matrix6d var2 = Matrix6d::Zero();

	Matrix6Xd sigmapts(6,25);
	sigmapts.col(0) = delta(mean.topRows(6), mean.bottomRows(6));
	for(int i = 0; i < 12; ++i) {
		VectorXd sp = mean + rt_scaled_cov.col(i);
		VectorXd sm = mean - rt_scaled_cov.col(i);
		sigmapts.col(2*i+1) = (sigmapts.col(0) + (delta(sp.topRows(6), mean.topRows(6)) - delta(sp.bottomRows(6), mean.bottomRows(6))));
		sigmapts.col(2*i+2) = (sigmapts.col(0) + (delta(sm.topRows(6), mean.topRows(6)) - delta(sm.bottomRows(6), mean.bottomRows(6))));

//		VectorXd spo = mean + rt_scaled_cov.col((i+6)%12);
//		VectorXd smo = mean - rt_scaled_cov.col((i+6)%12);
//		sigmapts.col(2*i+1) = (sigmapts.col(0) + (delta(sp.topRows(6), mean.topRows(6)) - delta(smo.bottomRows(6), mean.bottomRows(6))));
//		sigmapts.col(2*i+2) = (sigmapts.col(0) + (delta(sm.topRows(6), mean.topRows(6)) - delta(spo.bottomRows(6), mean.bottomRows(6))));

		cov += delta(sp.topRows(6), mean.topRows(6)) * delta(sp.bottomRows(6), mean.bottomRows(6)).transpose();
		cov += delta(sm.topRows(6), mean.topRows(6)) * delta(sm.bottomRows(6), mean.bottomRows(6)).transpose();

		var1 += delta(sp.topRows(6), mean.topRows(6)) * delta(sp.topRows(6), mean.topRows(6)).transpose();
		var1 += delta(sm.topRows(6), mean.topRows(6)) * delta(sm.topRows(6), mean.topRows(6)).transpose();

		var2 += delta(sp.bottomRows(6), mean.bottomRows(6)) * delta(sp.bottomRows(6), mean.bottomRows(6)).transpose();
		var2 += delta(sm.bottomRows(6), mean.bottomRows(6)) * delta(sm.bottomRows(6), mean.bottomRows(6)).transpose();
	}
	Matrix6d var_total = Matrix6d::Zero();
	for (int i=1; i<sigmapts.cols(); i++) {
		var_total += delta(sigmapts.col(i),sigmapts.col(0)) * delta(sigmapts.col(i),sigmapts.col(0)).transpose();
	}
	cout << "var_total" << endl;
	cout << var_total << endl;

	cov /= 12.0;
//	var1 /= 12.0;
//	var2 /= 12.0;
	cout << "cov" << endl;
	cout << cov << endl;
	cout << "var1" << endl;
	cout << var1 << endl;
	cout << "var2" << endl;
	cout << var2 << endl;
	cout << "sigmapts" << endl;
	cout << sigmapts << endl;
	return sigmapts;
}

int main() {
  RaveInitialize(false, OpenRAVE::Level_Debug);
  env = RaveCreateEnvironment();
  env->StopSimulation();
  env->Load(DATA_DIR"/boxes.env.xml");


  vector<RobotBasePtr> robots;
  env->GetRobots(robots);
  RobotBasePtr robot = robots[0];
  vector<RobotBase::ManipulatorPtr> manips = robot->GetManipulators();

  cc = CollisionChecker::GetOrCreate(*env);
  viewer.reset(new OSGViewer(env));
  env->AddViewer(viewer);
  dynamic_cast<OSGViewer::EventHandler*>(viewer->GetViewer().getCameraManipulator())->setTransformation(osg::Vec3d(0,0,2), osg::Vec3d(0,0,0), osg::Vec3d(0,1,0));

  viewer->AddKeyCallback('=', boost::bind(&AdjustTransparency, .05));
  viewer->AddKeyCallback('-', boost::bind(&AdjustTransparency, -.05));

  BeliefRobotAndDOFPtr rad = RADFromName("base3d", robot);
  const int n_dof = rad->GetDOF();
  const int n_theta = rad->GetNTheta();
  const int n_steps = 10;

  BeliefRobotAndDOFPtr obstacle = RADFromName("base3d", robots[1]);

	//////////////////////////////////////////////

  VectorXd x(6);
//  x << -0.3,0.4,-0.4,0,M_PI_4,M_PI_4/2.0;
  x << -0.3,0.4,-0.4,0,0,0;
	MatrixXd rt_Sigma = Eigen::MatrixXd::Identity(6,6)*0.1;
	rt_Sigma.topLeftCorner(3,3) = MatrixXd::Identity(3,3)*0.1;
//	rt_Sigma.bottomRightCorner(3,3) = MatrixXd::Identity(3,3);
//	rt_Sigma.topRightCorner(3,3) = MatrixXd::Identity(3,3)*0.02;
//	rt_Sigma.bottomLeftCorner(3,3) = MatrixXd::Identity(3,3)*0.02;
	VectorXd theta;
	rad->composeBelief(x, rt_Sigma, theta);

	VectorXd obstacle_x = toVectorXd(obstacle->GetDOFValues());
	MatrixXd obstacle_rt_Sigma = Eigen::MatrixXd::Identity(6,6)*0.1;
	obstacle_rt_Sigma.topLeftCorner(3,3) = MatrixXd::Identity(3,3)*0.05;
//	obstacle_rt_Sigma.bottomRightCorner(3,3) = MatrixXd::Identity(3,3)*0.05;
//	obstacle_rt_Sigma.topRightCorner(3,3) = MatrixXd::Identity(3,3)*0.01;
//	obstacle_rt_Sigma.bottomLeftCorner(3,3) = MatrixXd::Identity(3,3)*0.01;
	VectorXd obstacle_theta;
	obstacle->composeBelief(obstacle_x, obstacle_rt_Sigma, obstacle_theta);

	renderSigmaPts(rad, rad->sigmaPoints(theta), osg::Vec4f(1,0,0,0.5));
//	SetTransparency(handles.back(), 1);
	renderSigmaPts(obstacle, obstacle->sigmaPoints(obstacle_theta),osg::Vec4f(0,1,0,0.5));
//	SetTransparency(handles.back(), 1);

	obstacle->SetDOFValues(toDblVec(obstacle_x));

	// plot collisions
//	vector<Collision> collisions;
//	cc->DiscreteCheckSigma(rad, rad->sigmaPoints(theta), collisions);
//	cout << "collision " << collisions.size() << " " << collisions.back().distance << endl;

	VectorXd joint_x(12);
	joint_x.topRows(6) = x;
	joint_x.bottomRows(6) = obstacle_x;
	MatrixXd joint_rt_sigma = MatrixXd::Identity(12,12);
	joint_rt_sigma.topLeftCorner(6,6) = rt_Sigma;
	joint_rt_sigma.bottomRightCorner(6,6) = obstacle_rt_Sigma;

//	MatrixXd  rel_sigmapts = relativeSigmaPts(rad->sigmaPoints(theta), obstacle->sigmaPoints(obstacle_theta));
	MatrixXd  rel_sigmapts = jointSigmaPts(joint_x, joint_rt_sigma);
	renderSigmaPts(rad, rel_sigmapts, osg::Vec4f(0,0,1,0.5));
	SetColor(handles.back(), osg::Vec4f(0,0.6,0.6,0.2));
//	SetTransparency(handles.back(), 1);

	OR::Transform T_identity;
	T_identity.identity();
	KinBodyPtr obstacle_body = env->GetKinBody("obstacle");
	cout << obstacle_body->GetTransform().trans << endl;
	obstacle_body->SetTransform(T_identity);
	cout << obstacle_body->GetTransform().trans << endl;

	// plot collisions
//	cc->DiscreteCheckSigma(rad, rel_sigmapts, collisions);
//	cout << "collision " << collisions.size() << " " << collisions.back().distance << endl;
////	cc->DiscreteCheckTrajectory(theta.transpose(), rad, collisions);
//	PlotCollisions(collisions, *env, handles, 0);

//	handles.push_back(viewer->PlotSphere(OR::Vector(x[0],x[1],x[2]), 0.1, OR::RaveVector<float>(1,0,0,0.8)));

//	OR::Transform T = robot->GetTransform();
//	T.rot = toQuat(AngleAxisd(M_PI_4, Vector3d(1,0,0)));
//	cout << toMatrix3d(toQuat(AngleAxisd(M_PI_4, Vector3d(1,0,0)))) << endl << endl;
//	cout << (Matrix3d) AngleAxisd(M_PI_4, Vector3d(1,0,0)) << endl << endl;
//	robot->SetTransform(T);
	rad->SetDOFValues(toDblVec(VectorXd::Ones(6)*1000));
	obstacle->SetDOFValues(toDblVec(VectorXd::Ones(6)*1000));
  viewer->Idle();

  env.reset();
  viewer.reset();
  RaveDestroy();
}
