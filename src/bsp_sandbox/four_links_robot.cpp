#include "bsp/bsp.hpp"
#include "four_links_robot.hpp"
#include "sco/sco_fwd.hpp"

using namespace std;
using namespace OpenRAVE;
using namespace trajopt;

namespace FourLinksRobotBSP {

  template<typename T, size_t n>
  vector<T> vec(const std::array<T, n>& arr) {
    return vector<T>(arr.begin(), arr.end());
  }

  void initialize_robot(RobotBasePtr robot, const Vector4d& start) {
    robot->SetDOFValues(toDblVec(start));
  }

  void initialize_viewer(OSGViewerPtr viewer) {
    osg::Vec3d osg_eye(0, 0, 20);
    osg::Vec3d osg_center(0, 0, -3);
    osg::Vec3d osg_up(0, 1, 0);
    viewer->m_handler->setTransformation(osg_eye, osg_center, osg_up);
  }

  FourLinksRobotBSPPlanner::FourLinksRobotBSPPlanner() : BSPPlanner<FourLinksRobotBSPProblemHelper>() {}

  void FourLinksRobotBSPPlanner::initialize_optimizer_parameters(BSPTrustRegionSQP& opt, bool is_first_time) {
    opt.max_iter_                   = 350;
    opt.merit_error_coeff_          = 10;
    opt.merit_coeff_increase_ratio_ = 10;
    opt.max_merit_coeff_increases_  = 5;
    //opt.trust_shrink_ratio_         = .1;
    opt.trust_shrink_ratio_         = .8;
    opt.trust_expand_ratio_         = 1.5;
    //opt.trust_expand_ratio_         = 1.2;
    opt.min_trust_box_size_         = 1e-4;
    opt.min_approx_improve_         = 1e-4;
    opt.min_approx_improve_frac_    = -INFINITY;
    opt.improve_ratio_threshold_    = 0.15;
    //opt.improve_ratio_threshold_    = 0.15;
    opt.trust_box_size_             = 1e-1;
    opt.cnt_tolerance_              = 1e-4;
  }

  void FourLinksRobotBSPPlanner::initialize() {
    BSPPlanner<FourLinksRobotBSPProblemHelper>::initialize(); 
    helper->robot = robot;
    helper->rad = rad;
    helper->link = link;
    helper->goal_pos = goal_pos;
    helper->base_config = base_config;
    helper->link_lengths = link_lengths;
    helper->sigma_pts_scale = sigma_pts_scale;
    vector<double> lbs, ubs;
    rad->GetDOFLimits(lbs, ubs);
    cout << "state bounds: " << toVectorXd(lbs).transpose() << endl;
    cout << "state bounds: " << toVectorXd(ubs).transpose() << endl;
    helper->set_state_bounds(lbs, ubs);
  }

  FourLinksRobotGoalError::FourLinksRobotGoalError(FourLinksRobotBSPProblemHelperPtr helper) : helper(helper) {}

  VectorXd FourLinksRobotGoalError::operator()(const VectorXd& a) const {
    return helper->goal_pos - helper->angle_to_endpoint_position((Vector4d) a);
  }

  TransT FourLinksRobotBSPProblemHelper::angle_to_transform(const Vector4d& angle) const {
    TransT mat(2, 4);
    double t0 = base_config(2) + angle(0);
    double t1 = t0 + angle(1);
    double t2 = t1 + angle(2);
    double t3 = t2 + angle(3);
    mat << cos(t0), cos(t1), cos(t2), cos(t3),
           sin(t0), sin(t1), sin(t2), sin(t3);
    return mat;
  }

  Vector10d FourLinksRobotBSPProblemHelper::angle_to_joint_positions(const Vector4d& angle) const {
    Vector4d l0; l0 << link_lengths(0),               0,               0,               0;
    Vector4d l1; l1 << link_lengths(0), link_lengths(1),               0,               0;
    Vector4d l2; l2 << link_lengths(0), link_lengths(1), link_lengths(2),               0;
    Vector4d l3; l3 << link_lengths(0), link_lengths(1), link_lengths(2), link_lengths(3);
    TransT mat = angle_to_transform(angle);
    Vector10d res;
    res.segment<2>(0) = base_config.head<2>();
    res.segment<2>(2) = base_config.head<2>() + mat * l0;
    res.segment<2>(4) = base_config.head<2>() + mat * l1;
    res.segment<2>(6) = base_config.head<2>() + mat * l2;
    res.segment<2>(8) = base_config.head<2>() + mat * l3;
    return res;
  }

  Vector2d FourLinksRobotBSPProblemHelper::angle_to_endpoint_position(const Vector4d& angle) const {
    Vector4d l3; l3 << link_lengths(0), link_lengths(1), link_lengths(2), link_lengths(3);
    TransT mat = angle_to_transform(angle);
    return base_config.head<2>() + mat * l3;
  }

  FourLinksRobotBSPProblemHelper::FourLinksRobotBSPProblemHelper() : BSPProblemHelper<FourLinksRobotBeliefFunc>() {

    set_state_dim(4);
    set_sigma_dof(10);
    set_observe_dim(2);
    set_control_dim(4);

    //set_control_bounds( vector<double>(4, -1), vector<double>(4, 1) );
    set_control_bounds( vec(std::array<double, 4>{{-0.1, -0.3, -0.6, -1}}), vec(std::array<double, 4>{{0.1, 0.3, 0.6, 1}}));

    set_variance_cost(VarianceT::Identity(state_dim, state_dim) * 100);
    set_final_variance_cost(VarianceT::Identity(state_dim, state_dim) * 100);
    set_control_cost(ControlCostT::Identity(control_dim, control_dim));
  }

  ControlL1CostError::ControlL1CostError(const MatrixXd& R) : R(R) {}

  VectorXd ControlL1CostError::operator()(const VectorXd& a) const {
    assert(a.size() == R.cols());   
    assert(a.size() == R.rows());   
    return R * a;
  }

  void FourLinksRobotBSPProblemHelper::add_goal_constraint(OptProb& prob) {
    VectorOfVectorPtr f(new FourLinksRobotGoalError(boost::static_pointer_cast<FourLinksRobotBSPProblemHelper>(this->shared_from_this())));
    Vector2d coeffs = Vector2d::Ones();
    prob.addConstraint(ConstraintPtr(new ConstraintFromFunc(f, state_vars.row(T), coeffs, EQ, "goal")));
  }

  void FourLinksRobotBSPProblemHelper::add_control_cost(OptProb& prob) {
    for (int i = 0; i < T; ++i) {
      VectorOfVectorPtr f(new ControlL1CostError(R));
      VectorXd coeffs = VectorXd::Ones(control_dim);
      prob.addCost(CostPtr(new CostFromErrFunc(f, control_vars.row(i), coeffs, ABS, "control_cost")));
    }
  }

  void FourLinksRobotBSPProblemHelper::add_collision_term(OptProb& prob) {
    for (int i = 0; i < T; ++i) {
    //for (int i = 0; i < T; ++i) {
      //prob.addIneqConstraint(ConstraintPtr(new BeliefCollisionConstraint<FourLinksRobotBeliefFunc>(0.05, 30, rad, belief_vars.row(i), belief_func, link)));
      //prob.addIneqConstraint(ConstraintPtr(new BeliefCollisionConstraint<FourLinksRobotBeliefFunc>(0.05, 30, rad, belief_vars.row(i), belief_vars.row(i+1), belief_func, link)));
      prob.addIneqConstraint(ConstraintPtr(new BeliefCollisionConstraint<FourLinksRobotBeliefFunc>(0.05, 1, rad, belief_vars.row(i), belief_vars.row(i+1), belief_func, link)));
      //prob.addIneqConstraint(ConstraintPtr(new CollisionConstraint(0.05, 30, rad, state_vars.row(i), state_vars.row(i+1))));
      //prob.addCost(CostPtr(new BeliefCollisionCost<FourLinksRobotBeliefFunc>(0.05, 30, rad, belief_vars.row(i), belief_vars.row(i+1), belief_func, link)));
    }
    BeliefCollisionCheckerPtr cc = BeliefCollisionChecker::GetOrCreate(*(rad->GetEnv()));
    //CollisionChecker::GetOrCreate(*(rad->GetEnv()))->SetContactDistance(0.09);
    cc->SetContactDistance(0.09);
  }

  void FourLinksRobotBSPProblemHelper::configure_problem(OptProb& prob) {
    BSPProblemHelper<FourLinksRobotBeliefFunc>::configure_problem(prob);
    add_collision_term(prob);
  }

  void FourLinksRobotBSPProblemHelper::initialize() {
    BSPProblemHelper<FourLinksRobotBeliefFunc>::initialize();
    this->belief_func->sigma_pts_scale = sigma_pts_scale;
  }

  FourLinksRobotStateFunc::FourLinksRobotStateFunc() : StateFunc<StateT, ControlT, StateNoiseT>() {}

  FourLinksRobotStateFunc::FourLinksRobotStateFunc(BSPProblemHelperBasePtr helper) :
                            StateFunc<StateT, ControlT, StateNoiseT>(helper), four_links_robot_helper(boost::static_pointer_cast<FourLinksRobotBSPProblemHelper>(helper)) {}

  StateT FourLinksRobotStateFunc::operator()(const StateT& x, const ControlT& u, const StateNoiseT& m) const {
    Vector4d noise_scale(0.01, 0.03, 0.03, 0.04);
    return x + u + noise_scale.asDiagonal() * m;
  }

  FourLinksRobotObserveFunc::FourLinksRobotObserveFunc() : ObserveFunc<StateT, ObserveT, ObserveNoiseT>() {}

  FourLinksRobotObserveFunc::FourLinksRobotObserveFunc(BSPProblemHelperBasePtr helper) :
              ObserveFunc<StateT, ObserveT, ObserveNoiseT>(helper), four_links_robot_helper(boost::static_pointer_cast<FourLinksRobotBSPProblemHelper>(helper)) {}

  ObserveT FourLinksRobotObserveFunc::operator()(const StateT& x, const ObserveNoiseT& n) const {
    Vector2d pos = four_links_robot_helper->angle_to_endpoint_position(x);
    //double scale = 0.05 * (pos - four_links_robot_helper->goal_pos).norm() + 0.01;
    double scale = 0.01 * fabs(pos.x() + 3) + 0.001;
    //double scale = 0.2 * exp((pos.x() - 5) * 0.35) + 0.001;// 0.5*(pos.y()-3)*(pos.y()-3)+0.01;
    //double scale = 0.5*(pos.y()-3)*(pos.y()-3)+0.01;
    return pos + scale * n;
  }

  FourLinksRobotBeliefFunc::FourLinksRobotBeliefFunc() : EkfBeliefFunc<FourLinksRobotStateFunc, FourLinksRobotObserveFunc, BeliefT>() {}

  FourLinksRobotBeliefFunc::FourLinksRobotBeliefFunc(BSPProblemHelperBasePtr helper, StateFuncPtr f, ObserveFuncPtr h) :
             EkfBeliefFunc<FourLinksRobotStateFunc, FourLinksRobotObserveFunc, BeliefT>(helper, f, h),
             four_links_robot_helper(boost::static_pointer_cast<FourLinksRobotBSPProblemHelper>(helper)) {}

}

using namespace FourLinksRobotBSP;

int main(int argc, char *argv[]) {
  int T = 10;
  bool sim_plotting = false;
  bool sim_result_plotting = false;
  bool stage_plotting = false;
  bool stage_result_plotting = false;
  bool first_step_only = false;
  double noise_level = 1.;
  double sigma_pts_scale = 1.5;

  string data_dir = get_current_directory() + "/../data";

  {
    Config config;
    config.add(new Parameter<bool>("sim_plotting", &sim_plotting, "sim_plotting"));
    config.add(new Parameter<bool>("sim_result_plotting", &sim_result_plotting, "sim_result_plotting"));
    config.add(new Parameter<bool>("stage_plotting", &stage_plotting, "stage_plotting"));
    config.add(new Parameter<bool>("stage_result_plotting", &stage_result_plotting, "stage_result_plotting"));
    config.add(new Parameter<bool>("first_step_only", &first_step_only, "first_step_only"));
    config.add(new Parameter<double>("noise_level", &noise_level, "noise_level"));
    config.add(new Parameter<double>("sigma_pts_scale", &sigma_pts_scale, "sigma_pts_scale"));
    config.add(new Parameter<string>("data_dir", &data_dir, "data_dir"));
    CommandParser parser(config);
    parser.read(argc, argv, true);
  }
  
  string manip_name("arm");
  string link_name("Finger");

  RaveInitialize();
  EnvironmentBasePtr env = RaveCreateEnvironment();
  env->StopSimulation();
  env->Load(data_dir + "/four_links.env.xml");
  //env->Load(data_dir + "/four_links_2.env.xml");
  OSGViewerPtr viewer;
  RobotBasePtr robot = GetRobot(*env);


  double initial_skew_angle = PI/8;
  //Vector4d start(PI/2, 0, -PI/2, 0);//-initial_skew_angle, -PI + 2*initial_skew_angle, -initial_skew_angle);
  //Vector4d start(PI/2, initial_skew_angle, PI - 2*initial_skew_angle, initial_skew_angle);
  Vector4d start(-5*PI/6, PI/6, PI/3, PI/6);//initial_skew_angle, PI - 2*initial_skew_angle, initial_skew_angle);

  //Vector4d start(-PI/4, 0, 0, 0);
  Matrix4d start_sigma = Vector4d(0.001, 0.02, 0.03, 0.05).asDiagonal();//Matrix4d::Identity() * 0.01;
  deque<Vector4d> initial_controls;
  for (int i = 0; i < T; ++i) {
    initial_controls.push_back(Vector4d::Zero());
  }

  Vector2d goal_pos(3.3, -3);
  //Vector2d goal_pos(-2.5, 4);

  initialize_robot(robot, start);

  FourLinksRobotBSPPlannerPtr planner(new FourLinksRobotBSPPlanner());

  planner->start = start;
  planner->start_sigma = start_sigma;
  planner->goal_pos = goal_pos;
  planner->link_lengths = Vector4d(2, 2, 2, 2);
  planner->sigma_pts_scale = sigma_pts_scale;
  //planner->link_lengths = Vector4d(2, 2, 2, 2);
  planner->T = T;
  planner->env = env;
  planner->controls = initial_controls;
  planner->robot = robot;
  planner->base_config = Vector3d(0, 0, 0);
  planner->noise_level = noise_level;
  planner->rad = RADFromName(manip_name, robot);
  planner->link = planner->rad->GetRobot()->GetLink(link_name);
  planner->method = BSP::DiscontinuousBeliefSpace;
  planner->initialize();


  boost::function<void(OptProb*, DblVec&)> opt_callback;
  if (stage_plotting || sim_plotting || stage_result_plotting || sim_result_plotting) {
    viewer = OSGViewer::GetOrCreate(env);
    initialize_viewer(viewer);
  }

  if (stage_plotting) {
    opt_callback = boost::bind(&OpenRAVEPlotterMixin<FourLinksRobotBSPPlanner>::stage_plot_callback, 
                               planner, planner->helper->rad, viewer, _1, _2);
  }

  while (!planner->finished()) {
    planner->solve(opt_callback, 1, 1);
    if (stage_result_plotting) {
      OpenRAVEPlotterMixin<FourLinksRobotBSPPlanner>::stage_plot_callback(planner, planner->helper->rad, viewer, &(*(planner->prob)), planner->result);
    }
    planner->simulate_execution();
    if (first_step_only) break;
    if (sim_plotting) {
      OpenRAVEPlotterMixin<FourLinksRobotBSPPlanner>::sim_plot_callback(planner, planner->rad, viewer);
    }
  }

  if (planner->finished()) {
    if (is_sim_trajectory_in_collision(planner, planner->rad, env) {
      cout << "IN COLLISION" << endl;
    } else {
      cout << "NOT IN COLLISION" << endl;
    }
  }
  RaveDestroy();
  return 0;
}
