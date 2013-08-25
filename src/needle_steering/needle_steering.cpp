#include "needle_steering.hpp"

using namespace Needle;

void printVector(VectorXd x) {
  for (int i = 0; i < x.size(); ++i) {
    cout << x(i);
    if (i < x.size() - 1) cout << ", ";
  }
}

int main(int argc, char** argv)
{

  NeedleProblemHelperPtr helper(new NeedleProblemHelper());

  double plotting = false;
  double plot_final_result = false;

  helper->plotting=false;
  helper->verbose=false;
  helper->plot_final_result=false;
  helper->env_transparency = 0.5;

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
  helper->improve_ratio_threshold = 0.1;//0.25;
  helper->trust_shrink_ratio = 0.9;//0.7;
  helper->trust_expand_ratio = 1.3;//1.2;
  helper->record_trust_region_history = false;
  
  double start_vec_array[] = {-6.60848, 12.6176, -8.06651, 2.53666, -0.868663, 1.31701};//{-12.82092, 6.80976, 0.06844, 0, 0, 0};//{0, 0, 0, 0, 0, 0};
  double goal_vec_array[] = {-3.21932, 6.87362, -1.21877, 0, 0, 0};//{0, 0.896343312427, 7.49334469032, 0, 0, 0};

  helper->coeff_rotation = 1.;
  helper->coeff_speed = 1.;
  helper->coeff_rotation_regularization = 0.1;
  helper->coeff_orientation_error = 1;
  helper->collision_dist_pen = 0.05;
  helper->collision_coeff = 10;

  vector<double> start_vec(start_vec_array, start_vec_array + n_dof);
  vector<double> goal_vec(goal_vec_array, goal_vec_array + n_dof);
  {
    Config config;
    config.add(new Parameter<bool>("plotting", &helper->plotting, "plotting"));
    config.add(new Parameter<bool>("plot_final_result", &helper->plot_final_result, "plot_final_result"));
    config.add(new Parameter<bool>("verbose", &helper->verbose, "verbose"));
    config.add(new Parameter<double>("env_transparency", &helper->env_transparency, "env_transparency"));
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
    config.add(new Parameter< vector<double> >("s", &start_vec, "s"));
    config.add(new Parameter< vector<double> >("g", &goal_vec, "g"));
    CommandParser parser(config);
    parser.read(argc, argv);
  }

  RaveInitialize(false, helper->verbose ? Level_Debug : Level_Info);
  EnvironmentBasePtr env = RaveCreateEnvironment();
  env->StopSimulation();

  OSGViewerPtr viewer;
  if (helper->plotting || helper->plot_final_result) {
    viewer = OSGViewer::GetOrCreate(env);
    assert(viewer);
  }

  env->Load(string(DATA_DIR) + "/prostate.env.xml");

  if (helper->plotting || helper->plot_final_result) {
    viewer->SetAllTransparency(helper->env_transparency);
  }

  RobotBasePtr robot = GetRobot(*env);

  for (int i = 0; i < n_dof; ++i) helper->start[i] = start_vec[i];
  for (int i = 0; i < n_dof; ++i) helper->goal[i] = goal_vec[i];

  const char *ignored_kinbody_c_strs[] = { "KinBodyProstate", "KinBodyDermis", "KinBodyEpidermis", "KinBodyHypodermis" };
  helper->ignored_kinbody_names = vector<string>(ignored_kinbody_c_strs, end(ignored_kinbody_c_strs));

  
  helper->robot = robot;

  OptProbPtr prob;
  OptimizerT opt;

  helper->Initialize(argc, argv);
  helper->max_merit_coeff_increases = 10;
  prob.reset(new OptProb());
  helper->ConfigureProblem(*prob);
  opt = OptimizerT(prob);
  helper->ConfigureOptimizer(opt);
  opt.optimize();
  

  RaveDestroy();

  

  return 0;
  
}
