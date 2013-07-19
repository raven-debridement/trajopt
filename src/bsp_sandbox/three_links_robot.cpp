#include "bsp/bsp.hpp"


#include "three_links_robot.hpp"
#include "three_links_robot_forward_kinematics.hpp"
#include "sco/sco_fwd.hpp"

using namespace std;
using namespace OpenRAVE;
using namespace trajopt;

namespace ThreeLinksRobotBSP {

  template<typename T, size_t n>
  vector<T> vec(const std::array<T, n>& arr) {
    return vector<T>(arr.begin(), arr.end());
  }

  void initialize_robot(RobotBasePtr robot, const Vector3d& start) {
    robot->SetDOFValues(toDblVec(start));
  }

  void initialize_viewer(OSGViewerPtr viewer) {
    //osg::Vec3d osg_eye(0, 0, 4);
    //osg::Vec3d osg_center(0, 0, 0);
    //osg::Vec3d osg_up(0, 1, 0);
    //viewer->m_handler->setTransformation(osg_eye, osg_center, osg_up);
  }

  ThreeLinksRobotBSPPlanner::ThreeLinksRobotBSPPlanner() : BSPPlanner<ThreeLinksRobotBSPProblemHelper>() {}

  void ThreeLinksRobotBSPPlanner::initialize_optimizer_parameters(BSPTrustRegionSQP& opt, bool is_first_time) {
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

  void ThreeLinksRobotBSPPlanner::initialize() {
    BSPPlanner<ThreeLinksRobotBSPProblemHelper>::initialize(); 
    helper->robot = robot;
    helper->rad = rad;
    helper->link = link;
    helper->goal_trans = goal_trans;
    vector<double> lbs, ubs;
    rad->GetDOFLimits(lbs, ubs);
    helper->set_state_bounds(lbs, ubs);
  }

  ThreeLinksRobotBSPProblemHelper::ThreeLinksRobotBSPProblemHelper() : BSPProblemHelper<ThreeLinksRobotBeliefFunc>() {

    set_state_dim(3);
    set_sigma_dof(6);
    set_observe_dim(3);
    set_control_dim(3);

    set_control_bounds( vector<double>(3, -0.4), vector<double>(3, 0.4) );

    set_variance_cost(VarianceT::Identity(state_dim, state_dim));
    set_final_variance_cost(VarianceT::Identity(state_dim, state_dim));
    set_control_cost(ControlCostT::Identity(control_dim, control_dim));
  }

  void ThreeLinksRobotBSPProblemHelper::add_goal_constraint(OptProb& prob) {
    VectorXd coeffs(6); coeffs << 1, 1, 1, 1, 1, 1;
    VectorOfVectorPtr f(new CartPoseErrCalculator(matrixToTransform(goal_trans), rad, link));
    prob.addConstraint(ConstraintPtr(new ConstraintFromFunc(f, state_vars.row(T), coeffs, EQ, "goal")));
  }

  void ThreeLinksRobotBSPProblemHelper::add_collision_term(OptProb& prob) {
    for (int i = 0; i <= T; ++i) {
      //prob.addIneqConstraint(ConstraintPtr(new BeliefCollisionConstraint<ThreeLinksRobotBeliefFunc>(0.025, 1, rad, belief_vars.row(i), belief_func, link)));
      prob.addCost(CostPtr(new BeliefCollisionCost<ThreeLinksRobotBeliefFunc>(0.05, 30, rad, belief_vars.row(i), belief_func, link)));
    }
    BeliefCollisionCheckerPtr cc = BeliefCollisionChecker::GetOrCreate(*(rad->GetEnv()));
    cc->SetContactDistance(0.09);
  }

  void ThreeLinksRobotBSPProblemHelper::configure_problem(OptProb& prob) {
    BSPProblemHelper<ThreeLinksRobotBeliefFunc>::configure_problem(prob);
    add_collision_term(prob);
  }

  ThreeLinksRobotStateFunc::ThreeLinksRobotStateFunc() : StateFunc<StateT, ControlT, StateNoiseT>() {}

  ThreeLinksRobotStateFunc::ThreeLinksRobotStateFunc(BSPProblemHelperBasePtr helper) :
                            StateFunc<StateT, ControlT, StateNoiseT>(helper), three_links_robot_helper(boost::static_pointer_cast<ThreeLinksRobotBSPProblemHelper>(helper)) {}

  StateT ThreeLinksRobotStateFunc::operator()(const StateT& x, const ControlT& u, const StateNoiseT& m) const {
    Vector3d noise_scale(0.08, 0.13, 0.18);
    return x + u + noise_scale.asDiagonal() * m;
  }

  ThreeLinksRobotObserveFunc::ThreeLinksRobotObserveFunc() : ObserveFunc<StateT, ObserveT, ObserveNoiseT>() {}

  ThreeLinksRobotObserveFunc::ThreeLinksRobotObserveFunc(BSPProblemHelperBasePtr helper) :
              ObserveFunc<StateT, ObserveT, ObserveNoiseT>(helper), three_links_robot_helper(boost::static_pointer_cast<ThreeLinksRobotBSPProblemHelper>(helper)) {}

  ObserveT ThreeLinksRobotObserveFunc::operator()(const StateT& x, const ObserveNoiseT& n) const {
    Vector3d trans = forward_kinematics(x);
    double scale = 0.5*(trans(1)+0.2)*(trans(1)+0.2)+1;
    return trans + scale * n;
  }

  ThreeLinksRobotBeliefFunc::ThreeLinksRobotBeliefFunc() : EkfBeliefFunc<ThreeLinksRobotStateFunc, ThreeLinksRobotObserveFunc, BeliefT>() {}

  ThreeLinksRobotBeliefFunc::ThreeLinksRobotBeliefFunc(BSPProblemHelperBasePtr helper, StateFuncPtr f, ObserveFuncPtr h) :
             EkfBeliefFunc<ThreeLinksRobotStateFunc, ThreeLinksRobotObserveFunc, BeliefT>(helper, f, h),
             three_links_robot_helper(boost::static_pointer_cast<ThreeLinksRobotBSPProblemHelper>(helper)) {}

  ThreeLinksRobotOptimizerTask::ThreeLinksRobotOptimizerTask(QObject* parent) : BSPOptimizerTask(parent) {}

  ThreeLinksRobotOptimizerTask::ThreeLinksRobotOptimizerTask(int argc, char **argv, QObject* parent) : BSPOptimizerTask(argc, argv, parent) {}

  void ThreeLinksRobotOptimizerTask::run() {
    int T = 25;
    bool sim_plotting = false;
    bool stage_plotting = false;
    bool first_step_only = false;
    double noise_level = 1.;

    {
      Config config;
      config.add(new Parameter<bool>("sim_plotting", &sim_plotting, "sim_plotting"));
      config.add(new Parameter<bool>("stage_plotting", &stage_plotting, "stage_plotting"));
      config.add(new Parameter<bool>("first_step_only", &first_step_only, "first_step_only"));
      config.add(new Parameter<double>("noise_level", &noise_level, "noise_level"));
      CommandParser parser(config);
      parser.read(argc, argv, true);
    }

    string manip_name("arm");
    string link_name("Finger");

    RaveInitialize();
    EnvironmentBasePtr env = RaveCreateEnvironment();
    env->StopSimulation();
    env->Load(string(DATA_DIR) + "/three_links.env.xml");
    OSGViewerPtr viewer;
    RobotBasePtr robot = GetRobot(*env);

    Vector3d start = Vector3d::Zero();
    Matrix3d start_sigma = Matrix3d::Identity();
    deque<Vector3d> initial_controls;
    for (int i = 0; i < T; ++i) {
      initial_controls.push_back(Vector3d::Zero());
    }

    Matrix4d goal_trans;
    goal_trans <<  0, -1, 0, -0.1,
                   1,  0, 0,  0.2,
                   0,  0, 1,    0,
                   0,  0, 0,    1;

    initialize_robot(robot, start);

    ThreeLinksRobotBSPPlannerPtr planner(new ThreeLinksRobotBSPPlanner());

    planner->start = start;
    planner->start_sigma = start_sigma;
    planner->goal_trans = goal_trans;
    planner->T = T;
    planner->controls = initial_controls;
    planner->robot = robot;
    planner->noise_level = noise_level;
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
      opt_callback = boost::bind(&OpenRAVEPlotterMixin<ThreeLinksRobotBSPPlanner>::stage_plot_callback, 
                                 planner, planner->helper->rad, viewer, _1, _2);
    }

    while (!planner->finished()) {
      planner->solve(opt_callback, 1, 1);
      planner->simulate_execution();
      if (first_step_only) break;
      if (sim_plotting) {
        OpenRAVEPlotterMixin<ThreeLinksRobotBSPPlanner>::sim_plot_callback(planner, planner->rad, viewer);
      }
    }

    RaveDestroy();
  }

}

using namespace ThreeLinksRobotBSP;

int main(int argc, char *argv[]) {
	ThreeLinksRobotOptimizerTask task(argc, argv);
  task.run();
	return 0;
}
