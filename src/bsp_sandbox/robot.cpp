#include "bsp/bsp.hpp"
//#include "robot.hpp"
#include <openrave-core.h>
#include <openrave/openrave.h>

using namespace std;
using namespace OpenRAVE;
using namespace trajopt;

int main() {
  RaveInitialize();
  EnvironmentBasePtr env = RaveCreateEnvironment();
  env->StopSimulation();
  env->Load(string("/home/rocky/trajopt/data") + "/barrett.env.xml");
  OSGViewerPtr viewer = OSGViewer::GetOrCreate(env);
  RobotBasePtr robot = GetRobot(*env);
  //viewer.reset(new OSGViewer(env));
  //viewer->UpdateSceneData();

  {
    double robot_dofs_array[] = { 1.3, 1.3, 1.3, 0.5 };
    int robot_dof_idx_array[] = { 7, 8, 9, 10 };
    robot->SetDOFValues(vector<double>(robot_dofs_array, end(robot_dofs_array)), false, vector<int>(robot_dof_idx_array, end(robot_dof_idx_array)));
  }

  int robot_active_dof_idx_array[] = { 0, 1, 2, 3, 4, 5, 6 };
  vector<int> robot_active_dof_idx(robot_active_dof_idx_array, end(robot_active_dof_idx_array));

  for (int i= 0; i < robot_active_dof_idx.size(); ++i) cout << robot_active_dof_idx[i]<<" ";

  robot->SetActiveDOFs(robot_active_dof_idx);

  { // set camera transformation
    osg::Vec3d osg_eye(0, 0, 4);
    osg::Vec3d osg_center(0, 0, 0);
    osg::Vec3d osg_up(0, 1, 0);
    viewer->m_handler->setTransformation(osg_eye, osg_center, osg_up);
  }

  {
    double active_dof_array[] = { -0.81103357, 0.5414401, 0., 1.88999729, -4.08813702, 1.10511549, 2.1282108 };
      //-0.81103357, 1.58598118, 0., -1.42894471, 0.81721394, 1.46287215, -0.11429623 };
    robot->SetActiveDOFValues(vector<double>(active_dof_array, end(active_dof_array)));//, false, robot_active_dof_idx);
  }

  //env->AddViewer(viewer);
  viewer->Idle();

  RaveDestroy();

  return 0;
}
