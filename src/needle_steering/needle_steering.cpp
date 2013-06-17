#include "needle_steering.hpp"

using namespace Needle;

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

  bool plotting=false, verbose=false;
  bool plot_final_result=false;
  double env_transparency = 0.5;

  int T = 25;
  double r_min = 2.98119536;
  int n_dof = 6;

  int formulation = NeedleProblemHelper::Form2;
  int curvature_constraint = NeedleProblemHelper::BoundedRadius;
  int speed_formulation = NeedleProblemHelper::ConstantSpeed;
  int method = NeedleProblemHelper::Shooting;
  int curvature_formulation = NeedleProblemHelper::UseRadius;
  bool use_speed_deviation_constraint = false;
  bool use_speed_deviation_cost = false;

  double improve_ratio_threshold = 0.1;//0.25;
  double trust_shrink_ratio = 0.9;//0.7;
  double trust_expand_ratio = 1.3;//1.2;
  bool record_trust_region_history = false;
  
  double start_vec_array[] = {0, 0, 0, 0, 0, 0};
  double goal_vec_array[] = {0, 0.896343312427, 7.49334469032, 0, 0, 0};

  vector<double> start_vec(start_vec_array, start_vec_array + n_dof);
  vector<double> goal_vec(goal_vec_array, goal_vec_array + n_dof);

  {
    Config config;
    config.add(new Parameter<bool>("plotting", &plotting, "plotting"));
    config.add(new Parameter<bool>("plot_final_result", &plot_final_result, "plot_final_result"));
    config.add(new Parameter<bool>("verbose", &verbose, "verbose"));
    config.add(new Parameter<double>("env_transparency", &env_transparency, "env_transparency"));
    config.add(new Parameter<int>("T", &T, "T"));
    config.add(new Parameter<int>("formulation", &formulation, "formulation"));
    config.add(new Parameter<int>("curvature_constraint", &curvature_constraint, "curvature_constraint"));
    config.add(new Parameter<int>("method", &method, "method"));
    config.add(new Parameter<int>("curvature_formulation", &curvature_formulation, "curvature_formulation"));
    config.add(new Parameter<int>("speed_formulation", &speed_formulation, "speed_formulation"));
    config.add(new Parameter<double>("r_min", &r_min, "r_min"));
    config.add(new Parameter<double>("improve_ratio_threshold", &improve_ratio_threshold, "improve_ratio_threshold"));
    config.add(new Parameter<double>("trust_shrink_ratio", &trust_shrink_ratio, "trust_shrink_ratio"));
    config.add(new Parameter<double>("trust_expand_ratio", &trust_expand_ratio, "trust_expand_ratio"));
    config.add(new Parameter<bool>("use_speed_deviation_constraint", &use_speed_deviation_constraint, "use_speed_deviation_constraint"));
    config.add(new Parameter<bool>("use_speed_deviation_cost", &use_speed_deviation_cost, "use_speed_deviation_cost"));
    config.add(new Parameter<bool>("record_trust_region_history", &record_trust_region_history, "record_trust_region_history"));
    config.add(new Parameter< vector<double> >("s", &start_vec, "s"));
    config.add(new Parameter< vector<double> >("g", &goal_vec, "g"));
    CommandParser parser(config);
    parser.read(argc, argv);
  }

  double coeff_rotation = 1.;
  double coeff_speed = 1.;

  RaveInitialize(false, verbose ? Level_Debug : Level_Info);
  EnvironmentBasePtr env = RaveCreateEnvironment();
  env->StopSimulation();
  OSGViewerPtr viewer = OSGViewer::GetOrCreate(env);
  assert(viewer);

  env->Load(string(DATA_DIR) + "/prostate.env.xml");//needleprob.env.xml");
  viewer->SetAllTransparency(env_transparency);
  RobotBasePtr robot = GetRobot(*env);

  Vector6d start; for (int i = 0; i < n_dof; ++i) start[i] = start_vec[i];
  Vector6d goal; for (int i = 0; i < n_dof; ++i) goal[i] = goal_vec[i];

  const char *ignored_kinbody_c_strs[] = { "KinBodyProstate", "KinBodyDermis", "KinBodyEpidermis", "KinBodyHypodermis" };
  vector<string> ignored_kinbody_names(ignored_kinbody_c_strs, end(ignored_kinbody_c_strs));

  OptProbPtr prob(new OptProb());

  NeedleProblemHelperPtr helper(new NeedleProblemHelper());
  helper->start = start;
  helper->goal = goal;
  helper->coeff_rotation = coeff_rotation;
  helper->coeff_speed = coeff_speed;
  helper->T = T;
  helper->r_min = r_min;
  helper->n_dof = n_dof;
  helper->ignored_kinbody_names = ignored_kinbody_names;
  helper->collision_dist_pen = 0.025;
  helper->collision_coeff = 20;
  helper->formulation = formulation;
  helper->curvature_constraint = curvature_constraint;
  helper->speed_formulation = speed_formulation;
  helper->method = method;
  helper->curvature_formulation = curvature_formulation;
  helper->use_speed_deviation_constraint = use_speed_deviation_constraint;
  helper->use_speed_deviation_cost = use_speed_deviation_cost;
  helper->robot = robot;
  helper->ConfigureProblem(*prob);

  BasicTrustRegionSQP opt(prob);
  opt.max_iter_ = 500;    
  opt.improve_ratio_threshold_ = improve_ratio_threshold;
  opt.trust_shrink_ratio_ = trust_shrink_ratio;
  opt.trust_expand_ratio_ = trust_expand_ratio;
  opt.record_trust_region_history_ = record_trust_region_history;

  helper->ConfigureOptimizer(opt);

  boost::shared_ptr<TrajPlotter> plotter;
  plotter.reset(new TrajPlotter(helper->local_configs, helper->twistvars));
  if (plotting) {
    opt.addCallback(boost::bind(&TrajPlotter::OptimizerCallback, boost::ref(plotter), _1, _2, helper));
  }

  opt.optimize();
  
  if (plot_final_result) plotter->PlotBothTrajectories(prob, opt, helper);

  if (opt.record_trust_region_history_) {
    for (int i = 0; i < opt.trust_region_size_history.size(); ++i) {
      cout << "Trust region size history for iteration #" << i << " ";
      for (int j = 0; j < opt.trust_region_size_history[i].size(); ++j) {
        cout << opt.trust_region_size_history[i][j] << " ";
      }
      cout << endl;
    }
    for (int i = 0; i < opt.log_trust_region_size_history.size(); ++i) {
      cout << "Trust region logarithm size history for iteration #" << i << " ";
      for (int j = 0; j < opt.log_trust_region_size_history[i].size(); ++j) {
        cout << opt.log_trust_region_size_history[i][j] << " ";
      }
      cout << endl;
    }
  }

  RaveDestroy();

}
