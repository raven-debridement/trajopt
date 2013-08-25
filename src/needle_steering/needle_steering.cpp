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

template<typename T, size_t N>
T* end(T (&ra)[N]) {
  return ra + N;
}

void printVector(VectorXd x) {
  for (int i = 0; i < x.size(); ++i) {
    cout << x(i);
    if (i < x.size() - 1) cout << ", ";
  }
}

int main(int argc, char** argv)
{

  NeedleProblemHelperPtr helper(new NeedleProblemHelper());

  bool plotting = false;
  bool verbose = false;
  bool plot_final_result = false;
  double env_transparency = 0.1;
  boost::shared_ptr<Needle::TrajPlotter> plotter;

  string data_dir = get_current_directory() + "/../data";

  helper->T = 25;
  helper->r_min = 2.98119536;
  helper->n_dof = 6;

  helper->formulation = NeedleProblemHelper::Form1;
  helper->curvature_constraint = NeedleProblemHelper::ConstantRadius;
  helper->speed_formulation = NeedleProblemHelper::ConstantSpeed;
  helper->method = NeedleProblemHelper::Colocation;
  helper->curvature_formulation = NeedleProblemHelper::UseRadius;
  helper->rotation_cost = NeedleProblemHelper::UseRotationQuadraticCost;
  helper->use_speed_deviation_constraint = false;
  helper->use_speed_deviation_cost = false;
  helper->continuous_collision = true;
  helper->explicit_controls = true;
  helper->control_constraints = true;

  // parameters for the optimizer
  helper->improve_ratio_threshold = 0.1;
  helper->trust_shrink_ratio = 0.9;
  helper->trust_expand_ratio = 1.3;
  helper->record_trust_region_history = false;
  helper->merit_error_coeff = 10;
  helper->max_merit_coeff_increases = 10;

  vector<string> start_string_vec;
  start_string_vec.push_back("-6.60848,12.6176,-8.06651,2.53666,-0.868663,1.31701");
  start_string_vec.push_back("-6.60848,12.6176,-8.06651,2.53666,-0.868663,1.31701");
  vector<string> goal_string_vec;
  goal_string_vec.push_back("-3.21932,6.87362,-1.21877,0,0,0");
  goal_string_vec.push_back("-2.71912,8.00334,-1.12736,0,0,0");
  

  
  helper->coeff_rotation = 1.;
  helper->coeff_speed = 1.;
  helper->coeff_rotation_regularization = 0.1;
  helper->coeff_orientation_error = 1;
  helper->collision_dist_pen = 0.05;
  helper->collision_coeff = 10;

  const char *ignored_kinbody_c_strs[] = { "KinBodyProstate", "KinBodyDermis", "KinBodyEpidermis", "KinBodyHypodermis" };
  helper->ignored_kinbody_names = vector<string>(ignored_kinbody_c_strs, end(ignored_kinbody_c_strs));

  {
    Config config;
    config.add(new Parameter<bool>("plotting", &plotting, "plotting"));
    config.add(new Parameter<bool>("plot_final_result", &plot_final_result, "plot_final_result"));
    config.add(new Parameter<bool>("verbose", &verbose, "verbose"));
    config.add(new Parameter<double>("env_transparency", &env_transparency, "env_transparency"));
    config.add(new Parameter<int>("T", &helper->T, "T"));
    config.add(new Parameter<int>("formulation", &helper->formulation, "formulation"));
    config.add(new Parameter<int>("curvature_constraint", &helper->curvature_constraint, "curvature_constraint"));
    config.add(new Parameter<int>("method", &helper->method, "method"));
    config.add(new Parameter<int>("curvature_formulation", &helper->curvature_formulation, "curvature_formulation"));
    config.add(new Parameter<int>("speed_formulation", &helper->speed_formulation, "speed_formulation"));
    config.add(new Parameter<int>("rotation_cost", &helper->rotation_cost, "rotation_cost"));
    config.add(new Parameter<double>("coeff_rotation_regularization", &helper->coeff_rotation_regularization, "coeff_rotation_regularization"));
    config.add(new Parameter<double>("coeff_rotation", &helper->coeff_rotation, "coeff_rotation"));
    config.add(new Parameter<double>("coeff_speed", &helper->coeff_speed, "coeff_speed"));
    config.add(new Parameter<double>("coeff_orientation_error", &helper->coeff_orientation_error, "coeff_orientation_error"));
    config.add(new Parameter<double>("r_min", &helper->r_min, "r_min"));
    config.add(new Parameter<double>("improve_ratio_threshold", &helper->improve_ratio_threshold, "improve_ratio_threshold"));
    config.add(new Parameter<double>("trust_shrink_ratio", &helper->trust_shrink_ratio, "trust_shrink_ratio"));
    config.add(new Parameter<double>("trust_expand_ratio", &helper->trust_expand_ratio, "trust_expand_ratio"));
    config.add(new Parameter<double>("collision_dist_pen", &helper->collision_dist_pen, "collision_dist_pen"));
    config.add(new Parameter<bool>("use_speed_deviation_constraint", &helper->use_speed_deviation_constraint, "use_speed_deviation_constraint"));
    config.add(new Parameter<bool>("use_speed_deviation_cost", &helper->use_speed_deviation_cost, "use_speed_deviation_cost"));
    config.add(new Parameter<bool>("record_trust_region_history", &helper->record_trust_region_history, "record_trust_region_history"));
    config.add(new Parameter<bool>("explicit_controls", &helper->explicit_controls, "explicit_controls"));
    config.add(new Parameter<bool>("continuous_collision", &helper->continuous_collision, "continuous_collision"));
    config.add(new Parameter<bool>("control_constraints", &helper->control_constraints, "control_constraints"));
    config.add(new Parameter<string>("data_dir", &data_dir, "data_dir"));
    config.add(new Parameter< vector<string> >("start", &start_string_vec, "s"));
    config.add(new Parameter< vector<string> >("goal", &goal_string_vec, "g"));
    CommandParser parser(config);
    parser.read(argc, argv);
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

  env->Load(data_dir + "/prostate.env.xml");

  if (plotting || plot_final_result) {
    viewer->SetAllTransparency(env_transparency);
  }

  helper->robot = GetRobot(*env);

  helper->n_needles = start_string_vec.size();

  helper->starts.clear();
  helper->goals.clear();

  for (int i = 0; i < helper->n_needles; ++i) {
    DblVec start;
    DblVec goal;
    strtk::parse(start_string_vec[i], ",", start);
    strtk::parse(goal_string_vec[i], ",", goal);
    helper->starts.push_back(toVectorXd(start));
    helper->goals.push_back(toVectorXd(goal));
    cout << "start: " << toVectorXd(start).transpose() << endl;
    cout << "goal: " << toVectorXd(goal).transpose() << endl;
  }

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

  RaveDestroy();

  return 0;
  
}
