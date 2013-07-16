#include "bsp/bsp.hpp"
//#include "robot.hpp"
#include <openrave-core.h>
#include <openrave/openrave.h>
#include <array>

using namespace std;
using namespace OpenRAVE;
using namespace trajopt;

template<typename T, size_t n>
vector<T> vec(const std::array<T, n>& arr) {
  return vector<T>(arr.begin(), arr.end());
}

int main() {
  RaveInitialize();
  EnvironmentBasePtr env = RaveCreateEnvironment();
  env->StopSimulation();
  env->Load(string("/home/rocky/trajopt/data") + "/barrett.env.xml");
  OSGViewerPtr viewer = OSGViewer::GetOrCreate(env);
  RobotBasePtr robot = GetRobot(*env);

  robot->SetDOFValues(vec(std::array<double, 4>{{1.3, 1.3, 1.3, 0.5}}), false, vec(std::array<int, 4>{{7, 8, 9, 10}}));

  robot->SetActiveDOFs(vec(std::array<int, 7>{{0, 1, 2, 3, 4, 5, 6}}));

  { // set camera transformation
    osg::Vec3d osg_eye(0, 0, 4);
    osg::Vec3d osg_center(0, 0, 0);
    osg::Vec3d osg_up(0, 1, 0);
    viewer->m_handler->setTransformation(osg_eye, osg_center, osg_up);
  }

  robot->SetActiveDOFValues(vec(std::array<double, 7>{{-0.81103357, 0.5414401, 0., 1.88999729, -4.08813702, 1.10511549, 2.1282108}}));

  viewer->Idle();

  RaveDestroy();

  return 0;
}
