#include "trajopt/common.hpp"
#include <openrave-core.h>
#include "trajopt/collision_checker.hpp"
#include "osgviewer/osgviewer.hpp"
#include "utils/eigen_conversions.hpp"
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

Matrix3d skewSymmetric(const Vector3d& v) {
	Matrix3d m;
	m << 0, -v(2), v(1), v(2), 0, -v(0), -v(1), v(0), 0;
	return m;
}

Eigen::Matrix3d toMatrix3d(OR::RaveVector<float> rq) {
	return Quaterniond(rq[0], rq[1], rq[2], rq[3]).toRotationMatrix();
}

int main() {
  RaveInitialize(false, OpenRAVE::Level_Debug);
  env = RaveCreateEnvironment();
  env->StopSimulation();
  env->Load(DATA_DIR"/three_links_obstacles.env.xml");

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

  BeliefRobotAndDOFPtr rad = RADFromName("arm", robot);

  ////////////// INTERESTING STUFF STARTS HERE //////////////

  const int n_dof = rad->GetDOF();
  const int n_steps = 10;

  // consider this trajectory
  TrajArray traj(n_steps, n_dof); // aka Eigen::MatrixXd
  VectorXd start = toVectorXd(rad->GetDOFValues());
  VectorXd end = VectorXd::Ones(n_dof);
  for (int idof = 0; idof < n_dof; ++idof) {
		traj.col(idof) = VectorXd::LinSpaced(n_steps, start(idof), end(idof));
	}

  // render trajectory
  for (int i=0; i<n_steps; i++) {
		rad->SetDOFValues(toDblVec(traj.block(i,0,1,n_dof).transpose()));
		handles.push_back(viewer->PlotKinBody(rad->GetRobot()));
		double color_param = (((double)i)/((double)n_steps-1.0));
		SetColor(handles.back(), osg::Vec4f(0,color_param,1.0-color_param,1));
	}

  // check for collisions in the trajectory
  cout << "is trajectory in collisions " << isTrajectoryInCollision(cc, traj, rad) << endl;

  viewer->Idle(); // pause until 'p' is pressed
  handles.clear(); // clear the previous renderings

  // consider this dof_values
  Vector3d dof_values(0.6,0.4,0.2);

  // render a particular dof value
  rad->SetDOFValues(toDblVec(dof_values));
  handles.push_back(viewer->PlotKinBody(rad->GetRobot()));

  // check for collisions in for a particular dof
  cout << "is dof in collisions " << isTrajectoryInCollision(cc, dof_values.transpose(), rad) << endl;

	OR::Vector pt(0.2,0.36,0);
	handles.push_back(viewer->PlotSphere(pt, 0.01, OR::RaveVector<float>(1,1,0,1)));

	MatrixXd dpdq(3,6);
	dpdq.leftCols(3) = MatrixXd::Identity(3,3);
	dpdq.rightCols(3) = skewSymmetric(-toVector3d(pt));
	cout << dpdq << endl;
	cout << endl;

	MatrixXd jac_pos = rad->PositionJacobian(3, Vector(0,0,0));
	MatrixXd jac_rot = rad->RotationJacobian(3);
	MatrixXd dqdtheta(6,rad->GetDOF());
	dqdtheta.topRows(3) = jac_pos;
	dqdtheta.bottomRows(3) = jac_rot;
	cout << dqdtheta << endl;
	cout << endl;

	MatrixXd dpdtheta_rave = rad->PositionJacobian(3, pt);
	MatrixXd dpdtheta = dpdq * dqdtheta;
	cout << dpdtheta_rave << endl;
	cout << endl;
	cout << dpdtheta << endl;

  viewer->Idle();

  ////////////// INTERESTING STUFF ENDS HERE //////////////

  env.reset();
  viewer.reset();
  RaveDestroy();
}
