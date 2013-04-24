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
#include <Eigen/Dense>

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

void DiscreteCheckSigmaSigma(BeliefRobotAndDOFPtr rad0, const MatrixXd& sigma_pts0, BeliefRobotAndDOFPtr rad1, const MatrixXd& sigma_pts1, vector<Collision>& collisions) {
	OR::RobotBase::RobotStateSaver saver0 = rad0->Save();
	OR::RobotBase::RobotStateSaver saver1 = rad1->Save();

	vector<int> inds;

	KinBody::LinkPtr link0, link1;
	vector<KinBody::LinkPtr> links;
	rad0->GetAffectedLinks(links, true, inds);
	assert(links.size() == 1);
	link0 = links[0];
	rad1->GetAffectedLinks(links, true, inds);
	assert(links.size() == 1);
	link1 = links[0];

	vector<OR::Transform> tf0(sigma_pts0.cols()), tf1(sigma_pts1.cols());

	for (int i=0; i<sigma_pts0.cols(); i++) {
		rad0->SetDOFValues(toDblVec(sigma_pts0.col(i)));
		tf0[i] = link0->GetTransform();
	}
	for (int i=0; i<sigma_pts1.cols(); i++) {
		rad1->SetDOFValues(toDblVec(sigma_pts1.col(i)));
		tf1[i] = link1->GetTransform();
	}
	cc->MultiCastVsMultiCast(link0, tf0, link1, tf1, collisions);
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
	OR::RobotBase::RobotStateSaver saver = rad->Save();

	vector<KinBody::LinkPtr> links;
	vector<int> joint_inds;
	rad->GetAffectedLinks(links, true, joint_inds);

	// render sigma points
	for (int j=0; j<sigma_pts.cols(); j++) {
		rad->SetDOFValues(toDblVec(sigma_pts.col(j)));
		handles.push_back(viewer->PlotKinBody(rad->GetRobot()));
		SetColor(handles.back(), colorvec);
	}

	// render convex hulls of sigma points
	vector<DblVec> dofvals(sigma_pts.cols());
	for (int i=0; i<sigma_pts.cols(); i++)
		dofvals[i] = toDblVec(sigma_pts.col(i));
	cc->SetContactDistance(100);
	cc->PlotCastHull(*rad, links, dofvals, handles);
	SetTransparency(handles.back(), 0.2);
}

Matrix3d skewSymmetric(const Vector3d& v) {
	Matrix3d m;
	m << 0, -v(2), v(1), v(2), 0, -v(0), -v(1), v(0), 0;
	return m;
}

Matrix3d poseToRot(VectorXd x) {
	Vector3d aa;
	if (x.rows() == 3) aa = x;
	else if (x.rows() == 6) aa = x.bottomRows(3);
	Matrix3d rot;
	if (aa.norm() > 1e-5)
		rot = (Matrix3d) AngleAxisd(aa.norm(), aa.normalized());
	else
		rot = Matrix3d::Identity();
	return rot;
}
Vector3d toEigenVector(OR::Vector v) {
	return Vector3d(v[0], v[1], v[2]);
}
OR::Vector toORVector(Vector3d v) {
	return OR::Vector(v[0], v[1], v[2]);
}

Collision transformedCollision(const Collision& col_rhs, const Vector6d& pose) {
	Matrix3d rot = poseToRot(pose);
	Vector3d pos = pose.topRows(3);
	Collision col = col_rhs;
	col.ptA = toORVector(rot * toEigenVector(col.ptA) + pos);
	col.ptB = toORVector(rot * toEigenVector(col.ptB) + pos);
	col.normalB2A = toORVector(rot * toEigenVector(col.normalB2A));
	return col;
}

//VectorXd psi0, VectorXd psi1,
// psi0 is robot
// psi1 is obstacle
Vector6d delta(const Vector6d& psi0, const Vector6d& psi1) {
	Vector3d pos0 = psi0.topRows(3);
	Vector3d pos1 = psi1.topRows(3);
	Matrix3d rot0 = poseToRot(psi0);
	Matrix3d rot1 = poseToRot(psi1);
//	if (psi0.bottomRows(3).norm() > 1e-5)
//		rot0 = (Matrix3d) AngleAxisd(psi0.bottomRows(3).norm(), psi0.bottomRows(3).normalized());
//	else
//		rot0 = Matrix3d::Identity();
//	if (psi1.bottomRows(3).norm() > 1e-5)
//		rot1 = (Matrix3d) AngleAxisd(psi1.bottomRows(3).norm(), psi1.bottomRows(3).normalized());
//	else
//		rot1 = Matrix3d::Identity();

//	Vector6d rel_old;
//	rel_old.topRows(3) = rot1 * (pos0 - pos1);
//	Vector3d rel_rot_old = Vector3d::Zero();
//	for (int j=0; j<3; j++) {
//		rel_rot_old += skewSymmetric(rot0.col(j)) * rot1.col(j);
//	}
//	rel_old.bottomRows(3) = 0.5 * rel_rot_old;

	Matrix4d tf0 = Matrix4d::Identity();
	tf0.topLeftCorner(3,3) = rot0;
	tf0.topRightCorner(3,1) = pos0;
	Matrix4d tf1 = Matrix4d::Identity();
	tf1.topLeftCorner(3,3) = rot1;
	tf1.topRightCorner(3,1) = pos1;
	Matrix4d rel_tf = tf1.inverse() * tf0;

	Vector6d rel;
	Matrix3d rel_rot = rot1.transpose() * rot0;
//	Matrix3d rel_rot = rel_tf.topLeftCorner(3,3);
	Vector3d rel_pos = rot1.transpose() * (pos0 - pos1);
//	Vector3d rel_pos = rel_tf.topRightCorner(3,1);
	AngleAxisd rel_rot_aa(rel_rot);
	rel.topRows(3) = rel_pos;
	rel.bottomRows(3) = rel_rot_aa.axis() * rel_rot_aa.angle();

	return rel;
}

Matrix6Xd jointSigmaPts(const Vector12d& mean, const Matrix12d& rt_Sigma, BeliefRobotAndDOFPtr rad) {
	double lambda = 2;
	Matrix12d rt_scaled_cov = lambda * rt_Sigma;

	Matrix6d cov = Matrix6d::Zero();
	Matrix6d var1 = Matrix6d::Zero();
	Matrix6d var2 = Matrix6d::Zero();

	double w = (1/(2*lambda*lambda));

	Matrix6Xd sigmapts(6,13);
	for(int i = 0; i < 12; ++i) {
		Vector12d sp = mean + rt_scaled_cov.col(i);
		Vector12d sm = mean - rt_scaled_cov.col(i);

		cov += w * delta(sp.topRows(6), mean.topRows(6)) * delta(sp.bottomRows(6), mean.bottomRows(6)).transpose();
		cov += w * delta(sm.topRows(6), mean.topRows(6)) * delta(sm.bottomRows(6), mean.bottomRows(6)).transpose();

		var1 += w * delta(sp.topRows(6), mean.topRows(6)) * delta(sp.topRows(6), mean.topRows(6)).transpose();
		var1 += w * delta(sm.topRows(6), mean.topRows(6)) * delta(sm.topRows(6), mean.topRows(6)).transpose();

		var2 += w * delta(sp.bottomRows(6), mean.bottomRows(6)) * delta(sp.bottomRows(6), mean.bottomRows(6)).transpose();
		var2 += w * delta(sm.bottomRows(6), mean.bottomRows(6)) * delta(sm.bottomRows(6), mean.bottomRows(6)).transpose();
	}

	Vector6d rel_mean = delta(mean.topRows(6), mean.bottomRows(6));
	Matrix6d rel_var = var1+var2-2*cov;

	cout << "rel_var" << endl;
	cout << rel_var << endl;

	Eigen::SelfAdjointEigenSolver<MatrixXd> es(rel_var);
	Matrix6d sqrt_rel_var = es.eigenvectors().real() * es.eigenvalues().real().cwiseSqrt().asDiagonal() * es.eigenvectors().real().transpose();
	sigmapts = rad->sigmaPoints(rel_mean, sqrt_rel_var);

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

Quaterniond quatAdd(const Quaterniond& q0, const Quaterniond& q1) {
	return Quaterniond(q0.w()+q1.w(), q0.x()+q1.x(), q0.y()+q1.y(), q0.z()+q1.z());
}

class DualQuat {
	Quaterniond m_qr;
	Quaterniond m_qt;
public:
	DualQuat(const Quaterniond& qr, const Quaterniond& qt) : m_qr(qr), m_qt(qt) {}
	DualQuat operator* (const DualQuat& dq) {
		return DualQuat(m_qr*dq.m_qr, quatAdd(m_qr*dq.m_qt, m_qt*dq.m_qr));
	}
	VectorXd asVectorXd() {
		VectorXd v(8);
		v << m_qr.w(), m_qr.x(), m_qr.y(), m_qr.z(), m_qt.w(), m_qt.x(), m_qt.y(), m_qt.z();
		return v;
	}
};

DualQuat DOFToDualQuat(const VectorXd& x) {
	DualQuat dqt(Quaterniond::Identity(), Quaterniond(0, x(0)/2.0, x(1)/2.0, x(2)/2.0));
	Quaterniond qr = Quaterniond::Identity();
	if (x.bottomRows(3).norm() > 1e-5)
		qr = Quaterniond(AngleAxisd(x.bottomRows(3).norm(), x.bottomRows(3).normalized()));
	DualQuat dqr(qr, Quaterniond(0,0,0,0));
	cout << "---" << endl;
	cout << dqt.asVectorXd().transpose() << endl;
	cout << dqr.asVectorXd().transpose() << endl;
	cout << "---" << endl;
	return dqt*dqr;
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

	BeliefRobotAndDOFPtr rad = RADFromName("base3d", robot);
	const int n_steps = 10;

	BeliefRobotAndDOFPtr obstacle = RADFromName("base3d", robots[1]);

	cc = CollisionChecker::GetOrCreate(*env);
	viewer.reset(new OSGViewer(env));
	env->AddViewer(viewer);
	dynamic_cast<OSGViewer::EventHandler*>(viewer->GetViewer().getCameraManipulator())->setTransformation(osg::Vec3d(0,0,2), osg::Vec3d(0,0,0), osg::Vec3d(0,1,0));

	viewer->AddKeyCallback('=', boost::bind(&AdjustTransparency, .05));
	viewer->AddKeyCallback('-', boost::bind(&AdjustTransparency, -.05));

	//////////////////////////////////////////////

	VectorXd x = toVectorXd(rad->GetDOFValues());
	//  x << -0.3,0.4,-0.4,0,M_PI_4,M_PI_4/2.0;
	//x << -0.3,0.5+0.4,0,0,0,0;
	rad->SetDOFValues(toDblVec(x));
	VectorXd rt_Sigma_diag(6);
	rt_Sigma_diag << 0.1, 0.1, 0.1, 0.1, 0.1, 0.1;
	MatrixXd rt_Sigma = rt_Sigma_diag.asDiagonal();
	VectorXd theta;
	rad->composeBelief(x, rt_Sigma, theta);

	VectorXd obstacle_x = toVectorXd(obstacle->GetDOFValues());
	//obstacle_x << 0.2,0.5+0.25,0,0,0,0;
	obstacle->SetDOFValues(toDblVec(obstacle_x));
	VectorXd obstacle_rt_Sigma_diag(6);
	obstacle_rt_Sigma_diag << 0.1, 0.1, 0.1, 0.1, 0.1, 0.1;
	MatrixXd obstacle_rt_Sigma = obstacle_rt_Sigma_diag.asDiagonal();
	VectorXd obstacle_theta;
	obstacle->composeBelief(obstacle_x, obstacle_rt_Sigma, obstacle_theta);

	renderSigmaPts(rad, rad->sigmaPoints(theta), osg::Vec4f(0,0,1,0.5));
	renderSigmaPts(obstacle, obstacle->sigmaPoints(obstacle_theta), osg::Vec4f(0,0,1,0.5));

	// plot collisions
	vector<Collision> collisions;
	DiscreteCheckSigmaSigma(rad, rad->sigmaPoints(theta), obstacle, obstacle->sigmaPoints(obstacle_theta), collisions);
	cout << "collisions " << collisions.size() << ": ";
	for (int i=0; i<collisions.size(); i++) cout << collisions[i].distance << " "; cout << endl;

	{
	Vector12d joint_x;
	joint_x.topRows(6) = x;
	joint_x.bottomRows(6) = obstacle_x;
	Matrix12d joint_rt_sigma = Matrix12d::Identity();
	joint_rt_sigma.topLeftCorner(6,6) = rt_Sigma;
	joint_rt_sigma.bottomRightCorner(6,6) = obstacle_rt_Sigma;
	//joint_rt_sigma.bottomLeftCorner(6,6) = 0.099*Matrix6d::Identity();
	//joint_rt_sigma.topRightCorner(6,6) = 0.099*Matrix6d::Identity();

	Matrix6Xd  rel_sigmapts = jointSigmaPts(joint_x, joint_rt_sigma, rad);

	renderSigmaPts(rad, rel_sigmapts, osg::Vec4f(0,0,1,0.5));
	renderSigmaPts(obstacle, Vector6d::Zero(), osg::Vec4f(0,0,1,0.5));

	DiscreteCheckSigmaSigma(rad, rel_sigmapts, obstacle, Vector6d::Zero(), collisions);
	collisions.push_back(transformedCollision(collisions.back(), obstacle_x));
	cout << "collisions " << collisions.size() << ": ";
	for (int i=0; i<collisions.size(); i++) cout << collisions[i].distance << " "; cout << endl;
	}

	{
	Vector12d joint_x;
	joint_x.topRows(6) = obstacle_x;
	joint_x.bottomRows(6) = x;
	Matrix12d joint_rt_sigma = Matrix12d::Identity();
	joint_rt_sigma.topLeftCorner(6,6) = obstacle_rt_Sigma;
	joint_rt_sigma.bottomRightCorner(6,6) = rt_Sigma;

	Matrix6Xd  rel_sigmapts = jointSigmaPts(joint_x, joint_rt_sigma, obstacle);

	renderSigmaPts(rad, Vector6d::Zero(), osg::Vec4f(0,0,1,0.5));
	renderSigmaPts(obstacle, rel_sigmapts, osg::Vec4f(0,0,1,0.5));

	DiscreteCheckSigmaSigma(obstacle, rel_sigmapts, rad, Vector6d::Zero(), collisions);
	collisions.push_back(transformedCollision(collisions.back(), x));
	cout << "collisions " << collisions.size() << ": ";
	for (int i=0; i<collisions.size(); i++) cout << collisions[i].distance << " "; cout << endl;
	}

	PlotCollisions(collisions, *env, handles, 0);

	viewer->Idle();

	env.reset();
	viewer.reset();
	RaveDestroy();
}
