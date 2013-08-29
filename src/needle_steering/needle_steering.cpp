#include "needle_steering.hpp"
#include "strtk.hpp"
#include <stdio.h>  /* defines FILENAME_MAX */
#ifdef WINDOWS
  #include <direct.h>
  #define GetCurrentDir _getcwd
#else
  #include <unistd.h>
  #define GetCurrentDir getcwd
#endif

using namespace Needle;

string get_current_directory() {
  char cCurrentPath[FILENAME_MAX];
  if (!GetCurrentDir(cCurrentPath, sizeof(cCurrentPath))) {
    throw std::runtime_error("cannot get current path");
  }
  return string(cCurrentPath);
}

int main(int argc, char** argv) {

  NeedleProblemHelperPtr helper(new NeedleProblemHelper());

  bool plotting = false;
  bool verbose = false;
  bool plot_final_result = false;
  double env_transparency = 0.1;
  boost::shared_ptr<Needle::TrajPlotter> plotter;

  string data_dir = get_current_directory() + "/../data";


  vector<string> start_string_vec;
  start_string_vec.push_back("0,0,0,0,0,0");
  start_string_vec.push_back("1,2,0,0,0,0");
  vector<string> goal_string_vec;
  goal_string_vec.push_back("1,2,7,0.78,0,0");
  goal_string_vec.push_back("0,0,9,0.78,0,0");

  {
    Config config;
    config.add(new Parameter<bool>("plotting", &plotting, "plotting"));
    config.add(new Parameter<bool>("plot_final_result", &plot_final_result, "plot_final_result"));
    config.add(new Parameter<bool>("verbose", &verbose, "verbose"));
    config.add(new Parameter<double>("env_transparency", &env_transparency, "env_transparency"));
    config.add(new Parameter<string>("data_dir", &data_dir, "data_dir"));
    config.add(new Parameter< vector<string> >("start_vec", &start_string_vec, "s"));
    config.add(new Parameter< vector<string> >("goal_vec", &goal_string_vec, "g"));
    CommandParser parser(config);
    parser.read(argc, argv, true);
  }

  if (start_string_vec.size() != goal_string_vec.size()) {
    throw std::runtime_error("The number of start and goal vectors must be the same!");
  }

  if (start_string_vec.size() == 0) {
    throw std::runtime_error("You must provide at least 1 start and 1 goal vector.");
  }

  RaveInitialize(false, verbose ? Level_Debug : Level_Info);
  EnvironmentBasePtr env = RaveCreateEnvironment();
  env->StopSimulation();

  OSGViewerPtr viewer;
  if (plotting || plot_final_result) {
    viewer = OSGViewer::GetOrCreate(env);
    assert(viewer);
  }

  env->Load(data_dir + "/channel.env.xml");

  if (plotting || plot_final_result) {
    viewer->SetAllTransparency(env_transparency);
  }


  int n_needles = start_string_vec.size();

  vector<Vector6d> starts;
  vector<Vector6d> goals;

  for (int i = 0; i < n_needles; ++i) {
    DblVec start_vec;
    DblVec goal_vec;
    strtk::parse(start_string_vec[i], ",", start_vec);
    strtk::parse(goal_string_vec[i], ",", goal_vec);
    Vector6d start = toVectorXd(start_vec), goal = toVectorXd(goal_vec);
    Matrix4d start_pose, goal_pose;
    start_pose.topLeftCorner<3, 3>() = rotMat(start.tail<3>());
    start_pose.topRightCorner<3, 1>() = start.head<3>();
    start_pose.bottomLeftCorner<1, 3>() = Vector3d::Zero();
    start_pose(3, 3) = 1;
    goal_pose.topLeftCorner<3, 3>() = rotMat(goal.tail<3>());
    goal_pose.topRightCorner<3, 1>() = goal.head<3>();
    goal_pose.bottomLeftCorner<1, 3>() = Vector3d::Zero();
    goal_pose(3, 3) = 1;
    starts.push_back(logDown(start_pose));
    goals.push_back(logDown(goal_pose));
  }


  vector<KinBodyPtr> robots;

  vector<VectorXd> sol;

  for (int i = 0; i < n_needles; ++i) {
    // one needle at a time
    helper->InitParametersFromConsole(argc, argv);
    helper->n_needles = 1;
    helper->starts.push_back(starts[i]);
    helper->goals.push_back(goals[i]);
    helper->robots.push_back(env->ReadRobotURI(RobotBasePtr(), data_dir + "/needlebot.xml"));
    env->Add(helper->robots.back(), true);
    OptProbPtr prob(new OptProb());
    helper->ConfigureProblem(*prob);
    OptimizerT opt(prob);
    helper->ConfigureOptimizer(opt);
    if (plotting || plot_final_result) {
      plotter.reset(new Needle::TrajPlotter(helper->pis));
    }
    if (plotting) {
      opt.addCallback(boost::bind(&Needle::TrajPlotter::OptimizerCallback, boost::ref(plotter), _1, _2, helper));
    } 
    opt.optimize();
    for (int j = 0; j < helper->robots.size(); ++j) {
      robots.push_back(helper->robots[j]);
    }

    sol.push_back(helper->GetSolutions(opt).front());
    
    cout << "adding needles" << endl;
    helper->AddNeedlesToBullet(opt);
    cout << "finished adding needles" << endl;
  }

  trajopt::SetUserData(*env, "trajopt_cc", OpenRAVE::UserDataPtr());
  helper->InitParametersFromConsole(argc, argv);
  helper->n_needles = n_needles;
  helper->starts = starts;
  helper->goals = goals;
  for (int i = 0; i < n_needles; ++i) {
    helper->robots.push_back(env->ReadRobotURI(RobotBasePtr(), data_dir + "/needlebot.xml"));
    env->Add(helper->robots.back(), true);
  }

  OptProbPtr prob(new OptProb());
  helper->ConfigureProblem(*prob);
  OptimizerT opt(prob);
  helper->ConfigureOptimizer(opt);
  helper->SetSolutions(sol, opt);

  if (plotting || plot_final_result) {
    plotter.reset(new Needle::TrajPlotter(helper->pis));
  }
  if (plotting) {
    opt.addCallback(boost::bind(&Needle::TrajPlotter::OptimizerCallback, boost::ref(plotter), _1, _2, helper));
  }

  opt.optimize();
  RaveDestroy();

  return 0;
  
}
