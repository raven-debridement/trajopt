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

  int formulation = NeedleProblemHelper::Form1;
  int curvature_constraint = NeedleProblemHelper::ConstantRadius;
  int speed_formulation = NeedleProblemHelper::ConstantSpeed;
  int method = NeedleProblemHelper::Colocation;
  int curvature_formulation = NeedleProblemHelper::UseRadius;
  int rotation_cost = NeedleProblemHelper::UseRotationQuadraticCost;
  bool use_speed_deviation_constraint = false;
  bool use_speed_deviation_cost = false;

  double improve_ratio_threshold = 0.1;//0.25;
  double trust_shrink_ratio = 0.9;//0.7;
  double trust_expand_ratio = 1.3;//1.2;
  bool record_trust_region_history = false;
  
  double start_vec_array[] = {-6.60848, 12.6176, -8.06651, 2.53666, -0.868663, 1.31701};//{-12.82092, 6.80976, 0.06844, 0, 0, 0};//{0, 0, 0, 0, 0, 0};
  double goal_vec_array[] = {-3.21932, 6.87362, -1.21877, 0, 0, 0};//{0, 0.896343312427, 7.49334469032, 0, 0, 0};

  double coeff_rotation = 1.;
  double coeff_speed = 1.;
  double coeff_rotation_regularization = 0.1;

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
    config.add(new Parameter<int>("rotation_cost", &rotation_cost, "rotation_cost"));
    config.add(new Parameter<double>("coeff_rotation_regularization", &coeff_rotation_regularization, "coeff_rotation_regularization"));
    config.add(new Parameter<double>("coeff_rotation", &coeff_rotation, "coeff_rotation"));
    config.add(new Parameter<double>("coeff_speed", &coeff_speed, "coeff_speed"));
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

  

  RaveInitialize(false, verbose ? Level_Debug : Level_Info);
  EnvironmentBasePtr env = RaveCreateEnvironment();
  env->StopSimulation();
  OSGViewerPtr viewer;
  if (plotting || plot_final_result) {
    viewer = OSGViewer::GetOrCreate(env);
    assert(viewer);
  }

  env->Load(string(DATA_DIR) + "/prostate.env.xml");//needleprob.env.xml");

  if (plotting || plot_final_result) {
    viewer->SetAllTransparency(env_transparency);
  }
  RobotBasePtr robot = GetRobot(*env);

  Vector6d start; for (int i = 0; i < n_dof; ++i) start[i] = start_vec[i];
  Vector6d goal; for (int i = 0; i < n_dof; ++i) goal[i] = goal_vec[i];

  const char *ignored_kinbody_c_strs[] = { "KinBodyProstate", "KinBodyDermis", "KinBodyEpidermis", "KinBodyHypodermis" };
  vector<string> ignored_kinbody_names(ignored_kinbody_c_strs, end(ignored_kinbody_c_strs));


  NeedleProblemHelperPtr helper(new NeedleProblemHelper());
  helper->start = start;
  helper->goal = goal;
  helper->coeff_rotation = coeff_rotation;
  helper->coeff_rotation_regularization = coeff_rotation_regularization;
  helper->coeff_speed = coeff_speed;
  helper->T = T;
  helper->r_min = r_min;
  helper->n_dof = n_dof;
  helper->ignored_kinbody_names = ignored_kinbody_names;
  helper->collision_dist_pen = 0.05;
  helper->collision_coeff = 10;
  helper->formulation = formulation;
  helper->curvature_constraint = curvature_constraint;
  helper->speed_formulation = speed_formulation;
  helper->method = method;
  helper->curvature_formulation = curvature_formulation;
  helper->rotation_cost = rotation_cost;
  helper->use_speed_deviation_constraint = use_speed_deviation_constraint;
  helper->use_speed_deviation_cost = use_speed_deviation_cost;
  helper->robot = robot;

  OptProbPtr prob(new OptProb());
  helper->ConfigureProblem(*prob, true);

  OptimizerT opt(prob);
  opt.max_iter_ = 500;    
  opt.improve_ratio_threshold_ = improve_ratio_threshold;
  opt.trust_shrink_ratio_ = trust_shrink_ratio;
  opt.trust_expand_ratio_ = trust_expand_ratio;
  opt.record_trust_region_history_ = record_trust_region_history;

  helper->ConfigureOptimizer(opt);

  boost::shared_ptr<Needle::TrajPlotter> plotter;
  if (plotting || plot_final_result) {
    plotter.reset(new Needle::TrajPlotter(helper->local_configs, helper->twistvars));
  }
  if (plotting) {
    opt.addCallback(boost::bind(&Needle::TrajPlotter::OptimizerCallback, boost::ref(plotter), _1, _2, helper));
  }

  OptStatus status = opt.optimize();

  //DblVec x = opt.x();
  //vector<KinBodyPtr> bodies = helper->local_configs[0]->GetBodies();
  //MatrixXd vals = getTraj(x, helper->twistvars);
  //for (int i=0; i < vals.rows(); ++i) {
  //  helper->local_configs[i]->pose = helper->local_configs[i]->pose * expUp(vals.row(i));
  //  //helper->local_configs[i]->SetDOFValues(toDblVec(Vector6d::Zero()));
  //}

  //setVec(x, helper->twistvars.m_data, DblVec(helper->twistvars.size(), 0));

  //double total_collision_cost = 0.0;
  //for (int i = 0; i < helper->collision_costs.size(); ++i) {
  //  total_collision_cost += helper->collision_costs[i]->value(x, opt.getModel().get());
  //}

  //cout << "total collision cost: " << total_collision_cost << endl;

  //if (status != OPT_CONVERGED || total_collision_cost > opt.cnt_tolerance_) {
  //  NeedleProblemHelperPtr helper2(new NeedleProblemHelper());
  //  helper2->start = start;
  //  helper2->goal = goal;
  //  helper2->coeff_rotation = coeff_rotation;
  //  helper2->coeff_rotation_regularization = coeff_rotation_regularization;
  //  helper2->coeff_speed = coeff_speed;
  //  helper2->T = T;
  //  helper2->r_min = r_min;
  //  helper2->n_dof = n_dof;
  //  helper2->ignored_kinbody_names = ignored_kinbody_names;
  //  helper2->collision_dist_pen = 0.05;
  //  helper2->collision_coeff = 10;
  //  helper2->formulation = formulation;
  //  helper2->curvature_constraint = curvature_constraint;
  //  helper2->speed_formulation = speed_formulation;
  //  helper2->method = method;
  //  helper2->curvature_formulation = curvature_formulation;
  //  helper2->rotation_cost = rotation_cost;
  //  helper2->use_speed_deviation_constraint = use_speed_deviation_constraint;
  //  helper2->use_speed_deviation_cost = use_speed_deviation_cost;
  //  helper2->robot = robot;
  //  OptProbPtr prob2(new OptProb());
  //  helper2->ConfigureProblem(*prob2);
  //  OptimizerT opt2(prob2);
  //  opt2.max_iter_ = 500;    
  //  opt2.improve_ratio_threshold_ = improve_ratio_threshold;
  //  opt2.trust_shrink_ratio_ = trust_shrink_ratio;
  //  opt2.trust_expand_ratio_ = trust_expand_ratio;
  //  opt2.record_trust_region_history_ = record_trust_region_history;

  //  opt2.merit_error_coeff_ *= opt2.merit_coeff_increase_ratio_;
  //  opt2.max_merit_coeff_increases_ -= 1;

  //  helper2->ConfigureOptimizer(opt2);

  //  opt2.initialize(x);

  //  for (int i = 0; i < helper2->local_configs.size(); ++i) {
  //    helper2->local_configs[i]->pose = helper->local_configs[i]->pose;
  //  }

  //  boost::shared_ptr<Needle::TrajPlotter> plotter2;
  //  if (plotting || plot_final_result) {
  //    plotter2.reset(new Needle::TrajPlotter(helper2->local_configs, helper2->twistvars));
  //  }
  //  if (plotting) {
  //    opt2.addCallback(boost::bind(&Needle::TrajPlotter::OptimizerCallback, boost::ref(plotter2), _1, _2, helper2));
  //  }
  //  opt2.optimize();


  //  cout << logDown(helper2->local_configs[0]->pose).transpose() << endl;
  //  
  //  if (plot_final_result) plotter->PlotBothTrajectories(prob2, opt2, helper2);

  //  if (opt2.record_trust_region_history_) {
  //    for (int i = 0; i < opt2.trust_region_size_history.size(); ++i) {
  //      cout << "Trust region size history for iteration #" << i << " ";
  //      for (int j = 0; j < opt2.trust_region_size_history[i].size(); ++j) {
  //        cout << opt2.trust_region_size_history[i][j] << " ";
  //      }
  //      cout << endl;
  //    }
  //    for (int i = 0; i < opt2.log_trust_region_size_history.size(); ++i) {
  //      cout << "Trust region logarithm size history for iteration #" << i << " ";
  //      for (int j = 0; j < opt2.log_trust_region_size_history[i].size(); ++j) {
  //        cout << opt2.log_trust_region_size_history[i][j] << " ";
  //      }
  //      cout << endl;
  //    }
  //  }
  //} else {
  //  cout << "Constraints are satisfied already in stage 1!" << endl;
  //}

  //RaveDestroy();

}
