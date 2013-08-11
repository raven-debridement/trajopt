#include "bsp/bsp.hpp"
#include "raven.hpp"
#include "raven_forward_kinematics.hpp"
#include "sco/sco_fwd.hpp"

using namespace std;
using namespace OpenRAVE;
using namespace trajopt;

namespace RavenBSP {

  template<typename T, size_t n>
  vector<T> vec(const std::array<T, n>& arr) {
    return vector<T>(arr.begin(), arr.end());
  }

  vector<Vector6d> get_initial_trajectory(const string& data_dir, int T) {
    vector<Vector6d> ret;
    ifstream traj_file(data_dir + "/raven_traj.txt", ifstream::in);
    if (!traj_file.is_open()) {
      cout << "error while loading initial trajectory!" << endl;
      RaveDestroy();
      exit(1);
    }
    for (int i = 0; i <= T; ++i) {
      Vector6d cur_traj;
      for (int j = 0; j < 6; ++j) {
        traj_file >> cur_traj(j);
      }
      ret.push_back(cur_traj);
    }
    return ret;
  }

  deque<Vector6d> get_initial_controls(const vector<Vector6d>& initial_trajectory) {
    deque<Vector6d> ret;
    for (int i = 0; i < (int)initial_trajectory.size() - 1; ++i) {
      ret.push_back(initial_trajectory[i+1] - initial_trajectory[i]);
    }
    return ret;
  }

  void initialize_robot(RobotBasePtr robot, std::string manip_name, const Vector6d& start) {
	  robot->SetActiveManipulator(manip_name);
	  robot->SetDOFValues(toDblVec(start),0,robot->GetActiveManipulator()->GetArmIndices());
  }

  void initialize_viewer(OSGViewerPtr viewer) {
    osg::Vec3d osg_eye(0, 0, .4)// 0 0 4
    osg::Vec3d osg_center(0, 0, 0); // 0 0 0
    osg::Vec3d osg_up(0, 0, 0); // 0 1 0
    viewer->m_handler->setTransformation(osg_eye, osg_center, osg_up);
    /*const osg::Matrixd camera_mat;
    camera_mat << 0, 0, 1,// -.89,
    			  -1, 0, 0,// -.3,
			   	   0, -1, 0;// .4,
				   //0,  0, 0,  1;*/
    //osg::Vec3d eye(0,0,.4);
    //osg::Quat quat(.5, -.5, -.5, .5);
    //viewer->m_handler->setTransformation(eye, quat);
  }

  RavenBSPPlanner::RavenBSPPlanner() : BSPPlanner<RavenBSPProblemHelper>() {}

  void RavenBSPPlanner::initialize_optimizer_parameters(BSPTrustRegionSQP& opt, bool is_first_time) {
	  //not robot-specific, but need to be tweaked
    opt.max_iter_                   = 100;
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
    opt.cnt_tolerance_              = 1e-5;
  }

  void RavenBSPPlanner::initialize() {
    BSPPlanner<RavenBSPProblemHelper>::initialize();
    helper->robot = robot;
    helper->rad = rad;
    helper->link = link;
    helper->goal_trans = goal_trans;
    vector<double> lbs, ubs;
    rad->GetDOFLimits(lbs, ubs);
    helper->set_state_bounds(lbs, ubs);

    helper->sigma_pts_scale = sigma_pts_scale;
    helper->sigma_pts_scale_vec = sigma_pts_scale_vec;
  }

  RavenBSPProblemHelper::RavenBSPProblemHelper() : BSPProblemHelper<RavenBeliefFunc>() {

    set_state_dim(6); //TODO: 6
    set_sigma_dof(21); //TODO: 21
    set_observe_dim(6); //TODO: 6 + 6
    set_control_dim(6); //TODO: 6

    set_control_bounds( vector<double>(6, -0.3), vector<double>(6, 0.3) );

    set_variance_cost(VarianceT::Identity(state_dim, state_dim) * sqrt(10));
    set_final_variance_cost(VarianceT::Identity(state_dim, state_dim) * sqrt(10));
    set_control_cost(ControlCostT::Identity(control_dim, control_dim)*0.1);
  }

  void RavenBSPProblemHelper::add_goal_constraint(OptProb& prob) {
    VectorXd coeffs(6); coeffs << 1, 1, 1, 1, 1, 1;
    VectorOfVectorPtr f(new CartPoseErrCalculator(matrixToTransform(goal_trans), rad, link));
    prob.addConstraint(ConstraintPtr(new ConstraintFromFunc(f, state_vars.row(T), coeffs, EQ, "goal")));
  }

  void RavenBSPProblemHelper::add_collision_term(OptProb& prob) {
    for (int i = 0; i <= T; ++i) {
    //for (int i = 0; i < T; ++i) {
      //prob.addIneqConstraint(ConstraintPtr(new BeliefCollisionConstraint<RavenBeliefFunc>(0.025, 1, rad, belief_vars.row(i), belief_func, link)));
      prob.addCost(CostPtr(new BeliefCollisionCost<RavenBeliefFunc>(0.001, 1, rad, belief_vars.row(i), belief_func, link)));
      //prob.addCost(CostPtr(new BeliefCollisionCost<RavenBeliefFunc>(0.025, 1, rad, belief_vars.row(i), belief_vars.row(i+1), belief_func, link)));
    }
    BeliefCollisionCheckerPtr cc = BeliefCollisionChecker::GetOrCreate(*(rad->GetEnv()));
    cc->SetContactDistance(0.065);
  }

  void RavenBSPProblemHelper::configure_problem(OptProb& prob) {
    BSPProblemHelper<RavenBeliefFunc>::configure_problem(prob);
    add_collision_term(prob);
  }

  void RavenBSPProblemHelper::initialize() {
    BSPProblemHelper<RavenBeliefFunc>::initialize();
    this->belief_func->sigma_pts_scale = sigma_pts_scale;
    this->belief_func->sigma_pts_scale_vec = sigma_pts_scale_vec;
  }

  RavenStateFunc::RavenStateFunc() : StateFunc<StateT, ControlT, StateNoiseT>() {}

  RavenStateFunc::RavenStateFunc(BSPProblemHelperBasePtr helper) :
                            StateFunc<StateT, ControlT, StateNoiseT>(helper), barrett_robot_helper(boost::static_pointer_cast<RavenBSPProblemHelper>(helper)) {}

  StateT RavenStateFunc::operator()(const StateT& x, const ControlT& u, const StateNoiseT& m) const {
    return x + u + 0.01 * m;
  }

  RavenObserveFunc::RavenObserveFunc() : ObserveFunc<StateT, ObserveT, ObserveNoiseT>() {}

  RavenObserveFunc::RavenObserveFunc(BSPProblemHelperBasePtr helper) :
              ObserveFunc<StateT, ObserveT, ObserveNoiseT>(helper), barrett_robot_helper(boost::static_pointer_cast<RavenBSPProblemHelper>(helper)) {}

  ObserveT RavenObserveFunc::operator()(const StateT& x, const ObserveNoiseT& n) const {
	  //TODO: important update
	  /*
    ObserveT ret(observe_dim);
    Vector3d trans(0,0,1); //forward_kinematics(x);
    Vector3d beacon(-1.0, 0.0, 0.5);
    double dist = 3.0 * (trans(0) - beacon(0)) * (trans(0) - beacon(0));
    ret(0) = 1.0 / (1.0 + dist) + 0.1 * n(0);
    ret(1) = x(0) + 0.01 * n(1);
    ret(2) = x(3) + 0.01 * n(2);*/
	  ObserveT ret(observe_dim);
	  ret = x + 0.01*n;
    return ret;
  }
  //TODO: observation_masks, signed distance

  RavenBeliefFunc::RavenBeliefFunc() : EkfBeliefFunc<RavenStateFunc, RavenObserveFunc, BeliefT>() {}

  RavenBeliefFunc::RavenBeliefFunc(BSPProblemHelperBasePtr helper, StateFuncPtr f, ObserveFuncPtr h) :
             EkfBeliefFunc<RavenStateFunc, RavenObserveFunc, BeliefT>(helper, f, h),
             barrett_robot_helper(boost::static_pointer_cast<RavenBSPProblemHelper>(helper)) {}

}

using namespace RavenBSP;

RavenBSPWrapper::RavenBSPWrapper() : manip_name("arm"), link_name("wam7"), sigma_pts_scale(2), insertion_factor(0.1) {
}

void RavenBSPWrapper::setEnvironment(const EnvironmentBasePtr& env) {
	this->env = env;
}
void RavenBSPWrapper::setViewer(const OSGViewerPtr& viewer, bool sim_plotting, bool stage_plotting) {
	this->viewer = viewer;
	this->sim_plotting = sim_plotting;
	this->stage_plotting = stage_plotting;
}

void RavenBSPWrapper::initialize() {
	env->StopSimulation();
	RobotBasePtr robot = GetRobot(*env);

	planner = RavenBSPPlannerPtr(new RavenBSPPlanner());

	planner->start = start;
	planner->start_sigma = start_sigma;
	planner->goal_trans = goal_trans;
	planner->T = T;
	planner->sigma_pts_scale = 0;
	planner->sigma_pts_scale_vec = sigma_pts_scale * StateT::Ones();
	planner->sigma_pts_scale_vec(2) *= insertion_factor;
	planner->controls = controls;
	planner->robot = robot;
	planner->rad = RADFromName(manip_name, robot);
	planner->link = planner->rad->GetRobot()->GetLink(link_name);
	planner->method = BSP::DiscontinuousBeliefSpace;
	planner->initialize();


	if (stage_plotting || sim_plotting) {
		viewer = OSGViewer::GetOrCreate(env);
		initialize_viewer(viewer);
	}
	if (stage_plotting) {
		opt_callback = boost::bind(&OpenRAVEPlotterMixin<RavenBSPPlanner>::stage_plot_callback,
				planner, planner->helper->rad, viewer, _1, _2);
	}
}
bool RavenBSPWrapper::finished() {
	return planner->finished();
}
void RavenBSPWrapper::solve() {
	planner->solve(opt_callback, 1, 1);
}
void RavenBSPWrapper::simulate_execution() {
	planner->simulate_execution();
}

void RavenBSPWrapper::run() {
	initialize();

	while (!finished()) {
		solve();
		simulate_execution();
		if (sim_plotting) {
			OpenRAVEPlotterMixin<RavenBSPPlanner>::sim_plot_callback(planner, planner->rad, viewer);
		}
	}
}



int main(int argc, char *argv[]) {
  int T = 10;
  bool sim_plotting = false;
  bool stage_plotting = false;
  bool first_step_only = false;
  double insertion_factor = 0.1;

  string data_dir = get_current_directory(argv) + "/data";

  {
    Config config;
    config.add(new Parameter<bool>("sim_plotting", &sim_plotting, "sim_plotting"));
    config.add(new Parameter<bool>("stage_plotting", &stage_plotting, "stage_plotting"));
    config.add(new Parameter<bool>("first_step_only", &first_step_only, "first_step_only"));
    config.add(new Parameter<string>("data_dir", &data_dir, "data_dir"));
    config.add(new Parameter<double>("insertion_factor", &insertion_factor, "insertion_factor"));
    CommandParser parser(config);
    parser.read(argc, argv, true);
  }

  string manip_name("right_arm");
  string link_name("tool_R");
  // note: tool_R is 1cm off in z direction

  RaveInitialize();
  EnvironmentBasePtr env = RaveCreateEnvironment();
  env->StopSimulation();
  env->Load(data_dir + "/raven.env.xml");
  OSGViewerPtr viewer;
  RobotBasePtr robot = GetRobot(*env);

  auto initial_trajectory = get_initial_trajectory(data_dir, T);
  auto initial_controls = get_initial_controls(initial_trajectory);

  Vector6d start = initial_trajectory[0];
  Vector6d end = initial_trajectory.back();
  Matrix6d start_sigma = Matrix6d::Identity() *pow(0.00001,2);

  initialize_robot(robot, manip_name, start);

  OpenRAVE::geometry::RaveTransform<double> rave_start_trans = robot->GetActiveManipulator()->GetEndEffectorTransform();

  // set end EE transform
  OpenRAVE::geometry::RaveTransform<double> rave_goal_trans(rave_start_trans);
  rave_goal_trans.trans.y += .04;
  Matrix4d goal_trans = transformToMatrix(rave_goal_trans);

//#define RAVEN_CREATE_OBSTACLES
#ifdef RAVEN_CREATE_OBSTACLES
  // create box obstacle
  OpenRAVE::geometry::RaveTransform<double> box_trans;
  box_trans.identity();
  box_trans.trans.x = (rave_start_trans.trans.x + rave_goal_trans.trans.x) / 2;
  box_trans.trans.y = (rave_start_trans.trans.y + rave_goal_trans.trans.y) / 2;
  box_trans.trans.z = (rave_start_trans.trans.z + rave_goal_trans.trans.z) / 2;
  Vector3d box_extents(.01,.0025,.02);
  addOpenraveBox(env,"obstacle_box",box_trans,box_extents);

  // create workspace floor
  OpenRAVE::geometry::RaveTransform<double> floor_trans(box_trans);
  floor_trans.trans.z -= box_extents[2];
  Vector3d floor_extents(.1, .1, .0025);
  addOpenraveBox(env,"workspace_floor",floor_trans,floor_extents);
#endif

  RavenBSPPlannerPtr planner(new RavenBSPPlanner());


  planner->start = start;
  planner->start_sigma = start_sigma;
  planner->goal_trans = goal_trans;
  planner->T = T;

  double sigma_scale = 2;
  planner->sigma_pts_scale = 0;
  planner->sigma_pts_scale_vec = sigma_scale * StateT::Ones();
  planner->sigma_pts_scale_vec(2) *= insertion_factor;

  planner->controls = initial_controls;
  planner->robot = robot;
  planner->rad = RADFromName(manip_name, robot);
  planner->link = planner->rad->GetRobot()->GetLink(link_name);
  planner->method = BSP::DiscontinuousBeliefSpace;
  planner->initialize();

  vector<GraphHandlePtr> handles;
  boost::function<void(OptProb*, DblVec&)> opt_callback;
  if (stage_plotting || sim_plotting) {
    viewer = OSGViewer::GetOrCreate(env);
    initialize_viewer(viewer);

    handles.push_back(viewer->PlotAxes(robot->GetActiveManipulator()->GetEndEffectorTransform(),.05));
    handles.push_back(viewer->PlotAxes(matrixToTransform(goal_trans),.05));

    cout << "Initial setup. Press p to continue to BSP" << endl;
    viewer->Idle();
  }
  if (stage_plotting) {
    opt_callback = boost::bind(&OpenRAVEPlotterMixin<RavenBSPPlanner>::stage_plot_callback,
                               planner, planner->helper->rad, viewer, _1, _2);
  }


  while (!planner->finished()) {
    planner->solve(opt_callback, 1, 1);
    planner->simulate_execution();
    if (first_step_only) break;
    if (sim_plotting) {
      OpenRAVEPlotterMixin<RavenBSPPlanner>::sim_plot_callback(planner, planner->rad, viewer);
      handles.push_back(viewer->PlotAxes(robot->GetActiveManipulator()->GetEndEffectorTransform(),.025));
    }
  }

  RaveDestroy();
  return 0;
}

