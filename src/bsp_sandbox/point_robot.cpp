#include "bsp/bsp.hpp"
#include "point_robot.hpp"
#include "point_robot_forward_kinematics.hpp"
#include "sco/sco_fwd.hpp"

using namespace std;
using namespace OpenRAVE;
using namespace trajopt;

namespace PointRobotBSP {

  template<typename T, size_t n>
  vector<T> vec(const std::array<T, n>& arr) {
    return vector<T>(arr.begin(), arr.end());
  }

  void initialize_robot(RobotBasePtr robot, const Vector2d& start) {
    robot->SetDOFValues(toDblVec(start));
  }

  void initialize_viewer(OSGViewerPtr viewer) {
    //osg::Vec3d osg_eye(0, 0, 4);
    //osg::Vec3d osg_center(0, 0, 0);
    //osg::Vec3d osg_up(0, 1, 0);
    //viewer->m_handler->setTransformation(osg_eye, osg_center, osg_up);
  }

  PointRobotBSPPlanner::PointRobotBSPPlanner() : BSPPlanner<PointRobotBSPProblemHelper>() {}

  void PointRobotBSPPlanner::initialize_optimizer_parameters(BSPTrustRegionSQP& opt, bool is_first_time) {
    opt.max_iter_                   = 350;
    opt.merit_error_coeff_          = 10;
    opt.merit_coeff_increase_ratio_ = 10;
    opt.max_merit_coeff_increases_  = 5;
    opt.trust_shrink_ratio_         = .1;
    opt.trust_expand_ratio_         = 1.5;
    opt.min_trust_box_size_         = 1e-4;
    opt.min_approx_improve_         = 1e-4;
    opt.min_approx_improve_frac_    = -INFINITY;
    opt.improve_ratio_threshold_    = 0.15;
    opt.trust_box_size_             = 1e-1;
    opt.cnt_tolerance_              = 1e-4;
  }

  void PointRobotBSPPlanner::initialize() {
    BSPPlanner<PointRobotBSPProblemHelper>::initialize(); 
    helper->robot = robot;
    helper->rad = rad;
    helper->link = link;
    helper->goal_trans = goal_trans;
  }

  PointRobotBSPProblemHelper::PointRobotBSPProblemHelper() : BSPProblemHelper<PointRobotBeliefFunc>() {

    set_state_dim(2);
    set_sigma_dof(3);
    set_observe_dim(2);
    set_control_dim(2);

    set_state_bounds( vector<double>(2, -INFINITY), vector<double>(2, INFINITY) );
    set_control_bounds( vector<double>(2, -10), vector<double>(2, 10) );

    set_variance_cost(VarianceT::Identity(state_dim, state_dim) * 100);
    set_final_variance_cost(VarianceT::Identity(state_dim, state_dim) * 100);
    set_control_cost(ControlCostT::Identity(control_dim, control_dim)*0.05);
  }

  void PointRobotBSPProblemHelper::add_goal_constraint(OptProb& prob) {
    VectorXd coeffs(6); coeffs << 1, 1, 1, 1, 1, 1;
    VectorOfVectorPtr f(new CartPoseErrCalculator(matrixToTransform(goal_trans), rad, link));
    prob.addConstraint(ConstraintPtr(new ConstraintFromFunc(f, state_vars.row(T), coeffs, EQ, "goal")));
  }

  void PointRobotBSPProblemHelper::add_collision_term(OptProb& prob) {
    for (int i = 0; i < T; ++i) {
    //for (int i = 0; i <= T; ++i) {
      prob.addIneqConstraint(ConstraintPtr(new BeliefCollisionConstraint<PointRobotBeliefFunc>(0.1, 1, rad, belief_vars.row(i), belief_vars.row(i+1), belief_func, link)));
      //prob.addIneqConstraint(ConstraintPtr(new BeliefCollisionConstraint<PointRobotBeliefFunc>(0.1, 1, rad, belief_vars.row(i), belief_func, link)));
      //prob.addIneqConstraint(ConstraintPtr(new CollisionConstraint(0.1, 1, rad, state_vars.row(i), state_vars.row(i+1))));
      //prob.addCost(CostPtr(new BeliefCollisionCost<PointRobotBeliefFunc>(0.1, 1, rad, belief_vars.row(i), belief_vars.row(i+1), belief_func, link)));
    }
    BeliefCollisionCheckerPtr cc = BeliefCollisionChecker::GetOrCreate(*(rad->GetEnv()));
    cc->SetContactDistance(0.14);
  }

  void PointRobotBSPProblemHelper::configure_problem(OptProb& prob) {
    BSPProblemHelper<PointRobotBeliefFunc>::configure_problem(prob);
    add_collision_term(prob);
  }

  PointRobotStateFunc::PointRobotStateFunc() : StateFunc<StateT, ControlT, StateNoiseT>() {}

  PointRobotStateFunc::PointRobotStateFunc(BSPProblemHelperBasePtr helper) :
                            StateFunc<StateT, ControlT, StateNoiseT>(helper), point_robot_helper(boost::static_pointer_cast<PointRobotBSPProblemHelper>(helper)) {}

  StateT PointRobotStateFunc::operator()(const StateT& x, const ControlT& u, const StateNoiseT& m) const {
    return x + u + 0.01 * m;
  }

  PointRobotObserveFunc::PointRobotObserveFunc() : ObserveFunc<StateT, ObserveT, ObserveNoiseT>() {}

  PointRobotObserveFunc::PointRobotObserveFunc(BSPProblemHelperBasePtr helper) :
              ObserveFunc<StateT, ObserveT, ObserveNoiseT>(helper), point_robot_helper(boost::static_pointer_cast<PointRobotBSPProblemHelper>(helper)) {}

  ObserveT PointRobotObserveFunc::operator()(const StateT& x, const ObserveNoiseT& n) const {
    Vector3d trans = forward_kinematics(x);
		double scale = (0.5*(5.0 - trans[0])*(5.0 - trans[0])+0.01);
    return x.head<2>() + scale * n;
  }

  PointRobotBeliefFunc::PointRobotBeliefFunc() : EkfBeliefFunc<PointRobotStateFunc, PointRobotObserveFunc, BeliefT>() {}

  PointRobotBeliefFunc::PointRobotBeliefFunc(BSPProblemHelperBasePtr helper, StateFuncPtr f, ObserveFuncPtr h) :
             EkfBeliefFunc<PointRobotStateFunc, PointRobotObserveFunc, BeliefT>(helper, f, h),
             point_robot_helper(boost::static_pointer_cast<PointRobotBSPProblemHelper>(helper)) {}

  void stage_plot_callback(boost::shared_ptr<PointRobotBSPPlanner> planner, OSGViewerPtr viewer, OptProb* prob, DblVec& x) {
    vector<GraphHandlePtr> handles;
    handles.clear();
    OpenRAVEPlotterMixin<PointRobotBSPPlanner>::plot_opt_trajectory(planner, planner->rad, viewer, prob, x, &handles);
    vector<double> color_params;
    for (int i = 0; i <= planner->T; ++i) {
      color_params.push_back(((double)i)/((double)planner->T-1.0));
    }
    StateT state;
    VarianceT sigma;
    Vector3d mean;
    Matrix3d cov;
    for (int i = 0; i <= planner->T; ++i) {
      BeliefT b = getVec(x, planner->helper->belief_vars.row(i));
      planner->helper->belief_func->extract_state_and_sigma(b, &state, &sigma);
      belief_to_endeffector_noise(planner->rad, planner->link, state, sigma, &mean, &cov);
      handles.push_back(viewer->PlotEllipseXYContour(gaussian_as_transform(mean, cov), OR::Vector(0,color_params[i],1.0-color_params[i],1)));
    }
    BeliefT b = getVec(x, planner->helper->belief_vars.row(0));
    for (int i = 0; i <= planner->T; ++i) {
      planner->helper->belief_func->extract_state_and_sigma(b, &state, &sigma);
      belief_to_endeffector_noise(planner->rad, planner->link, state, sigma, &mean, &cov);
      handles.push_back(viewer->PlotEllipseXYContour(gaussian_as_transform(mean, cov), OR::Vector(0,color_params[i],1.0-color_params[i],1), true)); 
      if (i < planner->T) b = planner->helper->belief_func->call(b, getVec(x, planner->helper->control_vars.row(i)));
    }
    viewer->Idle();
  }

}

using namespace PointRobotBSP;

int main(int argc, char *argv[]) {
  int T = 4;
  bool sim_plotting = false;
  bool stage_plotting = false;
  bool first_step_only = false;

  {
    Config config;
    config.add(new Parameter<bool>("sim_plotting", &sim_plotting, "sim_plotting"));
    config.add(new Parameter<bool>("stage_plotting", &stage_plotting, "stage_plotting"));
    config.add(new Parameter<bool>("first_step_only", &first_step_only, "first_step_only"));
    CommandParser parser(config);
    parser.read(argc, argv, true);
  }

  string manip_name("base_point");
  string link_name("Base");

  RaveInitialize();
  EnvironmentBasePtr env = RaveCreateEnvironment();
  env->StopSimulation();
  env->Load(string(DATA_DIR) + "/point.env.xml");
  OSGViewerPtr viewer;
  RobotBasePtr robot = GetRobot(*env);


  Vector2d start = Vector2d::Zero();
  Matrix2d start_sigma = Matrix2d::Identity() * 0.2 * 0.2;

  deque<Vector2d> initial_controls;
  for (int i = 0; i < T; ++i) {
    //initial_controls.push_back((Vector2d(0, 5) - start) / T);//Vector2d::Zero());
    initial_controls.push_back(Vector2d::Zero());
  }

  Matrix4d goal_trans;
  goal_trans <<  1, 0, 0, 0,
                 0, 1, 0, 5,
                 0, 0, 1, 0,
                 0, 0, 0, 1;

  initialize_robot(robot, start);

  PointRobotBSPPlannerPtr planner(new PointRobotBSPPlanner());

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

  cout << "dof: " << planner->rad->GetDOF() << endl;

  boost::function<void(OptProb*, DblVec&)> opt_callback;
  if (stage_plotting || sim_plotting) {
    viewer = OSGViewer::GetOrCreate(env);
    initialize_viewer(viewer);
  }

  if (stage_plotting) {
    opt_callback = boost::bind(&stage_plot_callback, 
                               planner, viewer, _1, _2);
  }

  while (!planner->finished()) {
    planner->solve(opt_callback, 1, 1);
    planner->simulate_execution();
    if (first_step_only) break;
    if (sim_plotting) {
      OpenRAVEPlotterMixin<PointRobotBSPPlanner>::sim_plot_callback(planner, planner->rad, viewer);
    }
  }

  RaveDestroy();
  return 0;
}
