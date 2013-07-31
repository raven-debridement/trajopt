#include "bsp/bsp.hpp"
#include "barrett_robot.hpp"
#include "barrett_robot_forward_kinematics.hpp"
#include "sco/sco_fwd.hpp"

using namespace std;
using namespace OpenRAVE;
using namespace trajopt;

namespace BarrettRobotBSP {

  template<typename T, size_t n>
  vector<T> vec(const std::array<T, n>& arr) {
    return vector<T>(arr.begin(), arr.end());
  }

  vector<Vector7d> get_initial_trajectory(const string& data_dir, int T) {
    vector<Vector7d> ret;
    ifstream traj_file(data_dir + "/barrett_traj.txt", ifstream::in);
    if (!traj_file.is_open()) {
      cout << "error while loading initial trajectory!" << endl;
      RaveDestroy();
      exit(1);
    }
    for (int i = 0; i <= T; ++i) {
      Vector7d cur_traj;
      for (int j = 0; j < 7; ++j) {
        traj_file >> cur_traj(j);
      }
      ret.push_back(cur_traj);
    }
    return ret;
  }

  deque<Vector7d> get_initial_controls(const vector<Vector7d>& initial_trajectory) {
    deque<Vector7d> ret;
    for (int i = 0; i < (int)initial_trajectory.size() - 1; ++i) {
      ret.push_back(initial_trajectory[i+1] - initial_trajectory[i]);
    }
    return ret;
  }

  void initialize_robot(RobotBasePtr robot, const Vector7d& start) {
    robot->SetDOFValues(vec(std::array<double, 4>{{1.3, 1.3, 1.3, 0.5}}), false, vec(std::array<int, 4>{{7, 8, 9, 10}}));
    robot->SetActiveDOFs(vec(std::array<int, 7>{{0, 1, 2, 3, 4, 5, 6}}));
    robot->SetActiveDOFValues(toDblVec(start));
  }

  void initialize_viewer(OSGViewerPtr viewer) {
    osg::Vec3d osg_eye(0, 0, 4);
    osg::Vec3d osg_center(0, 0, 0);
    osg::Vec3d osg_up(0, 1, 0);
    viewer->m_handler->setTransformation(osg_eye, osg_center, osg_up);
  }

  BarrettRobotBSPPlanner::BarrettRobotBSPPlanner() : BSPPlanner<BarrettRobotBSPProblemHelper>() {}

  void BarrettRobotBSPPlanner::initialize_optimizer_parameters(BSPTrustRegionSQP& opt, bool is_first_time) {
    opt.max_iter_                   = 350;
    opt.merit_error_coeff_          = 10;
    opt.merit_coeff_increase_ratio_ = 10;
    opt.max_merit_coeff_increases_  = 5;
    opt.trust_shrink_ratio_         = .1;
    opt.trust_expand_ratio_         = 1.5;
    opt.min_trust_box_size_         = 1e-4;
    opt.min_approx_improve_         = 1e-4;
    opt.min_approx_improve_frac_    = -INFINITY;
    opt.improve_ratio_threshold_    = 0.25;
    opt.trust_box_size_             = 1e-1;
    opt.cnt_tolerance_              = 1e-4;
  }

  void BarrettRobotBSPPlanner::initialize() {
    BSPPlanner<BarrettRobotBSPProblemHelper>::initialize(); 
    helper->robot = robot;
    helper->rad = rad;
    helper->link = link;
    helper->goal_trans = goal_trans;
    vector<double> lbs, ubs;
    rad->GetDOFLimits(lbs, ubs);
    helper->set_state_bounds(lbs, ubs);
  }

  BarrettRobotBSPProblemHelper::BarrettRobotBSPProblemHelper() : BSPProblemHelper<BarrettRobotBeliefFunc>() {

    set_state_dim(7);
    set_sigma_dof(28);
    set_observe_dim(3);
    set_control_dim(7);

    set_control_bounds( vector<double>(7, -0.3), vector<double>(7, 0.3) );

    set_variance_cost(VarianceT::Identity(state_dim, state_dim) * sqrt(10));
    set_final_variance_cost(VarianceT::Identity(state_dim, state_dim) * sqrt(10));
    set_control_cost(ControlCostT::Identity(control_dim, control_dim)*0.1);
  }

  void BarrettRobotBSPProblemHelper::add_goal_constraint(OptProb& prob) {
    VectorXd coeffs(6); coeffs << 1, 1, 1, 1, 1, 1;
    VectorOfVectorPtr f(new CartPoseErrCalculator(matrixToTransform(goal_trans), rad, link));
    prob.addConstraint(ConstraintPtr(new ConstraintFromFunc(f, state_vars.row(T), coeffs, EQ, "goal")));
  }

  void BarrettRobotBSPProblemHelper::add_collision_term(OptProb& prob) {
    for (int i = 0; i <= T; ++i) {
    //for (int i = 0; i < T; ++i) {
      //prob.addIneqConstraint(ConstraintPtr(new BeliefCollisionConstraint<BarrettRobotBeliefFunc>(0.025, 1, rad, belief_vars.row(i), belief_func, link)));
      prob.addCost(CostPtr(new BeliefCollisionCost<BarrettRobotBeliefFunc>(0.025, 1, rad, belief_vars.row(i), belief_func, link)));
      //prob.addCost(CostPtr(new BeliefCollisionCost<BarrettRobotBeliefFunc>(0.025, 1, rad, belief_vars.row(i), belief_vars.row(i+1), belief_func, link)));
    }
    BeliefCollisionCheckerPtr cc = BeliefCollisionChecker::GetOrCreate(*(rad->GetEnv()));
    cc->SetContactDistance(0.065);
  }

  void BarrettRobotBSPProblemHelper::configure_problem(OptProb& prob) {
    BSPProblemHelper<BarrettRobotBeliefFunc>::configure_problem(prob);
    add_collision_term(prob);
  }

  BarrettRobotStateFunc::BarrettRobotStateFunc() : StateFunc<StateT, ControlT, StateNoiseT>() {}

  BarrettRobotStateFunc::BarrettRobotStateFunc(BSPProblemHelperBasePtr helper) :
                            StateFunc<StateT, ControlT, StateNoiseT>(helper), barrett_robot_helper(boost::static_pointer_cast<BarrettRobotBSPProblemHelper>(helper)) {}

  StateT BarrettRobotStateFunc::operator()(const StateT& x, const ControlT& u, const StateNoiseT& m) const {
    return x + u + 0.01 * m;
  }

  BarrettRobotObserveFunc::BarrettRobotObserveFunc() : ObserveFunc<StateT, ObserveT, ObserveNoiseT>() {}

  BarrettRobotObserveFunc::BarrettRobotObserveFunc(BSPProblemHelperBasePtr helper) :
              ObserveFunc<StateT, ObserveT, ObserveNoiseT>(helper), barrett_robot_helper(boost::static_pointer_cast<BarrettRobotBSPProblemHelper>(helper)) {}

  ObserveT BarrettRobotObserveFunc::operator()(const StateT& x, const ObserveNoiseT& n) const {
    ObserveT ret(observe_dim);
    Vector3d trans = forward_kinematics(x);
    Vector3d beacon(-1.0, 0.0, 0.5);
    double dist = 3.0 * (trans(0) - beacon(0)) * (trans(0) - beacon(0));
    ret(0) = 1.0 / (1.0 + dist) + 0.1 * n(0);
    ret(1) = x(0) + 0.01 * n(1);
    ret(2) = x(3) + 0.01 * n(2);
    return ret;
  }

  BarrettRobotBeliefFunc::BarrettRobotBeliefFunc() : EkfBeliefFunc<BarrettRobotStateFunc, BarrettRobotObserveFunc, BeliefT>() {}

  BarrettRobotBeliefFunc::BarrettRobotBeliefFunc(BSPProblemHelperBasePtr helper, StateFuncPtr f, ObserveFuncPtr h) :
             EkfBeliefFunc<BarrettRobotStateFunc, BarrettRobotObserveFunc, BeliefT>(helper, f, h),
             barrett_robot_helper(boost::static_pointer_cast<BarrettRobotBSPProblemHelper>(helper)) {}

}

using namespace BarrettRobotBSP;

int main(int argc, char *argv[]) {
  int T = 26;
  bool sim_plotting = false;
  bool stage_plotting = false;
  bool first_step_only = false;

  string data_dir = get_current_directory() + "/../data";

  {
    Config config;
    config.add(new Parameter<bool>("sim_plotting", &sim_plotting, "sim_plotting"));
    config.add(new Parameter<bool>("stage_plotting", &stage_plotting, "stage_plotting"));
    config.add(new Parameter<bool>("first_step_only", &first_step_only, "first_step_only"));
    config.add(new Parameter<string>("data_dir", &data_dir, "data_dir"));
    CommandParser parser(config);
    parser.read(argc, argv, true);
  }

  string manip_name("arm");
  string link_name("wam7");

  RaveInitialize();
  EnvironmentBasePtr env = RaveCreateEnvironment();
  env->StopSimulation();
  env->Load(data_dir + "/barrett.env.xml");
  OSGViewerPtr viewer;
  RobotBasePtr robot = GetRobot(*env);

  auto initial_trajectory = get_initial_trajectory(data_dir, T);
  auto initial_controls = get_initial_controls(initial_trajectory);

  Vector7d start = initial_trajectory[0];
  Matrix7d start_sigma = Matrix7d::Identity() * 0.22 * 0.22;

  Matrix4d goal_trans;
  goal_trans <<  0,  0, 1, 0.6,
                 0,  1, 0, 0.4,
                -1,  0, 0, 0.6,
                 0,  0, 0,   1;

  initialize_robot(robot, start);

  BarrettRobotBSPPlannerPtr planner(new BarrettRobotBSPPlanner());

  planner->start = start;
  planner->start_sigma = start_sigma;
  planner->goal_trans = goal_trans;
  planner->T = T;
  planner->controls = initial_controls;
  planner->robot = robot;
  planner->rad = RADFromName(manip_name, robot);
  planner->link = planner->rad->GetRobot()->GetLink(link_name);
  planner->method = BSP::DiscontinuousBeliefSpace;
  planner->initialize();


  boost::function<void(OptProb*, DblVec&)> opt_callback;
  if (stage_plotting || sim_plotting) {
    viewer = OSGViewer::GetOrCreate(env);
    initialize_viewer(viewer);
  }
  if (stage_plotting) {
    opt_callback = boost::bind(&OpenRAVEPlotterMixin<BarrettRobotBSPPlanner>::stage_plot_callback, 
                               planner, planner->helper->rad, viewer, _1, _2);
  }

  while (!planner->finished()) {
    planner->solve(opt_callback, 1, 1);
    planner->simulate_execution();
    if (first_step_only) break;
    if (sim_plotting) {
      OpenRAVEPlotterMixin<BarrettRobotBSPPlanner>::sim_plot_callback(planner, planner->rad, viewer);
    }
  }

  RaveDestroy();
  return 0;
}
