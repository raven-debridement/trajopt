#include "bsp/bsp.hpp"
#include "raven.hpp"
#include "raven_forward_kinematics.hpp"
#include "sco/sco_fwd.hpp"
#include <array>

using namespace std;
using namespace OpenRAVE;
using namespace trajopt;

namespace RavenBSP {

  template<typename T, size_t n>
  vector<T> vec(const std::array<T, n>& arr) {
    return vector<T>(arr.begin(), arr.end());
  }

  vector<Vector12d> get_initial_trajectory(const string& data_dir, int T) {
    vector<Vector12d> ret;
    ifstream traj_file(data_dir + "/raven_traj.txt", ifstream::in);
    if (!traj_file.is_open()) {
      cout << "error while loading initial trajectory!" << endl;
      RaveDestroy();
      exit(1);
    }
    for (int i = 0; i <= T; ++i) {
      Vector12d cur_traj;
      for (int j = 0; j < 6; ++j) {
        traj_file >> cur_traj(j);
      }
      ret.push_back(cur_traj);
    }
    return ret;
  }

  deque<Vector12d> get_initial_controls(const vector<Vector12d>& initial_trajectory) {
    deque<Vector12d> ret;
    for (int i = 0; i < (int)initial_trajectory.size() - 1; ++i) {
      Vector12d vec;
      //vec << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
      ret.push_back(vec);
      ret.push_back(initial_trajectory[i+1] - initial_trajectory[i]);
    }
    return ret;
  }

  void initialize_robot(RobotBasePtr robot, const Vector12d& start) {
	  Matrix<int, 12, 1> indv;
	  indv << 2,  3,  4,  5,  6,  7, 15, 16, 17, 18, 19, 20;
	  std::vector<int> inds;
	  for (int i=0;i<12;i++) {
		  inds.push_back(indv(i));
	  }
	  robot->SetActiveDOFs(inds);
	  robot->SetActiveDOFValues(toDblVec(start));
  }

  void initialize_viewer(OSGViewerPtr viewer) {
    osg::Vec3d osg_eye(1,0, .8);// 0 0 4
    osg::Vec3d osg_center(0, 0, 0); // 0 0 0
    osg::Vec3d osg_up(0, 1, 0); // 0 1 0
    viewer->m_handler->setTransformation(osg_eye, osg_center, osg_up);
    viewer->m_handler->setRotation(osg::Quat(.5,-.5,-.5,.5));
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
    helper->link_L = link_L;
    helper->link_R = link_R;
    helper->goal_trans_L = goal_trans_L;
    helper->goal_trans_R = goal_trans_R;
    vector<double> lbs, ubs;
    rad->GetDOFLimits(lbs, ubs);
    helper->set_state_bounds(lbs, ubs);
  }

  RavenBSPProblemHelper::RavenBSPProblemHelper() : BSPProblemHelper<RavenBeliefFunc>() {

    set_state_dim(12); //TODO: 6
    set_sigma_dof(78); //TODO: 21
    set_observe_dim(12); //TODO: 6 + 6
    set_control_dim(12); //TODO: 6

    set_control_bounds( vector<double>(12, -0.3), vector<double>(12, 0.3) );

    set_variance_cost(VarianceT::Identity(state_dim, state_dim) * sqrt(10));
    set_final_variance_cost(VarianceT::Identity(state_dim, state_dim) * sqrt(10));
    set_control_cost(ControlCostT::Identity(control_dim, control_dim)*0.1);
  }

  void RavenBSPProblemHelper::add_goal_constraint(OptProb& prob) {
    VectorXd coeffs(6); coeffs << 1, 1, 1, 1, 1, 1;
    {
    	VectorOfVectorPtr f(new CartPoseErrCalculator(matrixToTransform(goal_trans_L), rad, link_L));
    	prob.addConstraint(ConstraintPtr(new ConstraintFromFunc(f, state_vars.row(T), coeffs, EQ, "goal_L")));
    }
    {
    	VectorOfVectorPtr f(new CartPoseErrCalculator(matrixToTransform(goal_trans_R), rad, link_R));
    	prob.addConstraint(ConstraintPtr(new ConstraintFromFunc(f, state_vars.row(T), coeffs, EQ, "goal_R")));
    }
  }

  void RavenBSPProblemHelper::add_collision_term(OptProb& prob) {
//	  //discrete State space:
//	  for (int i = 0; i <= T; ++i) {
//		  prob.addIneqConstraint(ConstraintPtr(new CollisionConstraint(0.001, 1, rad, state_vars.row(i))));
//	  }
//	  CollisionCheckerPtr cc = CollisionChecker::GetOrCreate(*(rad->GetEnv()));
//
//	  //cont. State space:
//	  for (int i = 0; i < T; ++i) {
//		  prob.addIneqConstraint(ConstraintPtr(new CollisionConstraint(0.001, 1, rad, state_vars.row(i), state_vars.row(i+1))));
//	  }
//	  CollisionCheckerPtr cc = CollisionChecker::GetOrCreate(*(rad->GetEnv()));
//
//	  //discrete bsp
//	  for (int i = 0; i <= T; ++i) {
//		  //for (int i = 0; i < T; ++i) {
//		  prob.addIneqConstraint(ConstraintPtr(new BeliefCollisionConstraint<RavenBeliefFunc>(0.001, 1, rad, belief_vars.row(i), belief_func, link)));
//	  }
//	  BeliefCollisionCheckerPtr cc = BeliefCollisionChecker::GetOrCreate(*(rad->GetEnv()));

	  //cont. bsp
	  for (int i = 0; i < T; ++i) {
		  prob.addIneqConstraint(ConstraintPtr(new BeliefCollisionConstraint<RavenBeliefFunc>(0.001, 1, rad, belief_vars.row(i), belief_vars.row(i+1), belief_func, link_R)));
		  prob.addIneqConstraint(ConstraintPtr(new BeliefCollisionConstraint<RavenBeliefFunc>(0.001, 1, rad, belief_vars.row(i), belief_vars.row(i+1), belief_func, link_L)));
	  }
	  BeliefCollisionCheckerPtr cc = BeliefCollisionChecker::GetOrCreate(*(rad->GetEnv()));

	  cc->ExcludeCollisionPair(*(robot->GetLink("grasper1_R")),*(robot->GetLink("grasper2_R")));
	  cc->ExcludeCollisionPair(*(robot->GetLink("grasper1_L")),*(robot->GetLink("grasper2_L")));
	  cc->SetContactDistance(0.001);
  }

  void RavenBSPProblemHelper::configure_problem(OptProb& prob) {
    BSPProblemHelper<RavenBeliefFunc>::configure_problem(prob);
    add_collision_term(prob);
  }

  void RavenBSPProblemHelper::initialize() {
    BSPProblemHelper<RavenBeliefFunc>::initialize();
  }

  RavenStateFunc::RavenStateFunc() : StateFunc<StateT, ControlT, StateNoiseT>() {}

  RavenStateFunc::RavenStateFunc(BSPProblemHelperBasePtr helper) :
                            StateFunc<StateT, ControlT, StateNoiseT>(helper), barrett_robot_helper(boost::static_pointer_cast<RavenBSPProblemHelper>(helper)) {}

  StateT RavenStateFunc::operator()(const StateT& x, const ControlT& u, const StateNoiseT& m) const {
	  // StateNoiseT adj_noise = 0.01 * m;
	  // adj_noise(2) *= .01;
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

RavenBSPWrapper::RavenBSPWrapper() : manip_name_L("left_arm"), manip_name_R("right_arm"), link_name_L("tool_L"), link_name_R("tool_R"), insertion_factor(0.1) {
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
	planner->goal_trans_R = goal_trans_L;
	planner->goal_trans_R = goal_trans_R;
	planner->T = T;
	planner->controls = controls;
	planner->robot = robot;
	planner->rad = RADFromName("active", robot);
	planner->link_L = planner->rad->GetRobot()->GetLink(link_name_L);
	planner->link_R = planner->rad->GetRobot()->GetLink(link_name_R);
	planner->method = BSP::ContinuousBeliefSpace;
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
	controls = planner->controls;
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
  int T = 50;
  bool sim_plotting = false;
  bool stage_plotting = false;
  bool first_step_only = false;
  double insertion_factor = 0.1;

  double box_x = 0.005;
  double box_y = .005;
  double box_z = .005;

  double box_trans_x = 0.0;
  double box_trans_y = 0.0;
  double box_trans_z = 0.0;

  bool add_box = false;
  bool add_floor = false;

  string data_dir = get_current_directory(argv) + "/data";

  {
    Config config;
    config.add(new Parameter<bool>("sim_plotting", &sim_plotting, "sim_plotting"));
    config.add(new Parameter<bool>("stage_plotting", &stage_plotting, "stage_plotting"));
    config.add(new Parameter<bool>("first_step_only", &first_step_only, "first_step_only"));
    config.add(new Parameter<string>("data_dir", &data_dir, "data_dir"));
    config.add(new Parameter<double>("insertion_factor", &insertion_factor, "insertion_factor"));
    config.add(new Parameter<int>("T", &T, "T"));
    config.add(new Parameter<double>("box_x", &box_x, "box_x"));
    config.add(new Parameter<double>("box_y", &box_y, "box_y"));
    config.add(new Parameter<double>("box_z", &box_z, "box_z"));
    config.add(new Parameter<double>("box_trans_x", &box_trans_x, "box_trans_x"));
    config.add(new Parameter<double>("box_trans_y", &box_trans_y, "box_trans_y"));
    config.add(new Parameter<double>("box_trans_z", &box_trans_z, "box_trans_z"));
    config.add(new Parameter<bool>("add_box", &add_box, "add_box"));
    config.add(new Parameter<bool>("add_floor", &add_floor, "add_floor"));
    CommandParser parser(config);
    parser.read(argc, argv, true);
  }

  string manip_name_R("right_arm");
  string link_name_R("tool_R");
  string manip_name_L("left_arm");
  string link_name_L("tool_L");

  RaveInitialize();
  EnvironmentBasePtr env = RaveCreateEnvironment();
  env->StopSimulation();
  env->Load(data_dir + "/raven.env.xml");
  OSGViewerPtr viewer;
  RobotBasePtr robot = GetRobot(*env);


  Vector12d startJoints;
  startJoints <<
		  0.35669398307800293,
		  1.9775397777557373,
		  -0.14986588060855865,
		  0.6145784854888916,
		  -0.14158375561237335,
		  0.0,
		  0.7972189784049988,
		  2.0593690872192383,
		  -0.14696653187274933,
		  -3.0865819454193115,
		  0.6097292304039001,
		  0.0;

  Vector12d endJoints;
  endJoints <<
		  0.7972189784049988,
		  2.0593690872192383,
		  -0.14696653187274933,
		  -3.0865819454193115,
		  0.6097292304039001,
		  0.0,
		  0.35669398307800293,
		  1.9775397777557373,
		  -0.14986588060855865,
		  0.6145784854888916,
		  -0.14158375561237335,
		  0.0;


  vector<Vector12d> initial_trajectory;
  for (int i=0;i< T;i++) {
	  float factor = ((float) i)/T;
	  Vector12d value = startJoints + factor * (endJoints-startJoints);
	  initial_trajectory.push_back(value);
  }

  auto initial_controls = get_initial_controls(initial_trajectory);

  Vector12d start = initial_trajectory[0];
  Vector12d end = initial_trajectory.back();
  Matrix12d start_sigma = Matrix12d::Identity() *pow(0.00001,2);
  start_sigma(2,2) *= insertion_factor;

  initialize_robot(robot, start);

  robot->SetActiveManipulator(manip_name_L);
  OpenRAVE::geometry::RaveTransform<double> rave_start_trans_L = robot->GetActiveManipulator()->GetEndEffectorTransform();

  robot->SetActiveManipulator(manip_name_R);
  OpenRAVE::geometry::RaveTransform<double> rave_start_trans_R = robot->GetActiveManipulator()->GetEndEffectorTransform();

  // set end EE transform
  OpenRAVE::geometry::RaveTransform<double> rave_goal_trans_L(rave_start_trans_R);
  Matrix4d goal_trans_L = transformToMatrix(rave_goal_trans_L);

  OpenRAVE::geometry::RaveTransform<double> rave_goal_trans_R(rave_start_trans_L);
  Matrix4d goal_trans_R = transformToMatrix(rave_goal_trans_R);

#define RAVEN_CREATE_OBSTACLES
#ifdef RAVEN_CREATE_OBSTACLES
  // create box obstacle
  OpenRAVE::geometry::RaveTransform<double> box_trans;
  box_trans.identity();
  box_trans.trans.x = ((rave_start_trans_R.trans.x + rave_goal_trans_R.trans.x) / 2) + box_trans_x;
  box_trans.trans.y = ((rave_start_trans_R.trans.y + rave_goal_trans_R.trans.y) / 2) + box_trans_y;
  box_trans.trans.z = ((rave_start_trans_R.trans.z + rave_goal_trans_R.trans.z) / 2) + box_trans_z;
  Vector3d box_extents(box_x,box_y,box_z);//Vector3d box_extents(.01,.0025,.015);
  if (add_box) {
	  addOpenraveBox(env,"obstacle_box",box_trans,box_extents);
  }

  // create workspace floor
  OpenRAVE::geometry::RaveTransform<double> floor_trans(box_trans);
  floor_trans.trans.z -= box_extents[2] + box_z/2.0 - .0025;
  Vector3d floor_extents(.1, .1, .0025);

  if (add_floor) {
	  addOpenraveBox(env,"workspace_floor",floor_trans,floor_extents);
  }
#endif

  RavenBSPPlannerPtr planner(new RavenBSPPlanner());


  planner->start = start;
  planner->start_sigma = start_sigma;
  planner->goal_trans_L = goal_trans_L;
  planner->goal_trans_R = goal_trans_R;
  planner->T = T;

  planner->controls = initial_controls;
  planner->robot = robot;
  planner->rad = RADFromName("active", robot);
  planner->link_L = planner->rad->GetRobot()->GetLink(link_name_L);
  planner->link_R = planner->rad->GetRobot()->GetLink(link_name_R);
  planner->method = BSP::ContinuousBeliefSpace;
  planner->initialize();

  // so not in self-collision. might not need anymore
  int index = robot->GetJointIndex("grasper_joint_1_R");
  cout << index << endl;
  std::vector<int> indices;
  indices.push_back(index);
  std::vector<double> values;
  values.push_back(.5);
  robot->SetDOFValues(values,0,indices);

  vector<GraphHandlePtr> handles;
  boost::function<void(OptProb*, DblVec&)> opt_callback;
  if (stage_plotting || sim_plotting) {
    viewer = OSGViewer::GetOrCreate(env);
    initialize_viewer(viewer);

    //robot->SetActiveManipulator(manip_name_L);
    //handles.push_back(viewer->PlotAxes(robot->GetActiveManipulator()->GetEndEffectorTransform(),.05));
    //handles.push_back(viewer->PlotAxes(matrixToTransform(goal_trans_L),.05));

//    robot->SetActiveManipulator(manip_name_R);
//	handles.push_back(viewer->PlotAxes(robot->GetActiveManipulator()->GetEndEffectorTransform(),.05));
//    handles.push_back(viewer->PlotAxes(matrixToTransform(goal_trans_R),.05));

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

    for (int i=0; i < planner->controls.size(); i++) {
    	cout << planner->controls[i] << endl;
    }

    if (first_step_only) break;
    if (sim_plotting) {
      OpenRAVEPlotterMixin<RavenBSPPlanner>::sim_plot_callback(planner, planner->rad, viewer);
      //handles.push_back(viewer->PlotAxes(robot->GetActiveManipulator()->GetEndEffectorTransform(),.015));
    }
  }

  robot->SetActiveDOFValues(toDblVec(startJoints));
  deque<ControlT> controls = planner->controls;
  std::vector<double> currentJoints;
  robot->GetActiveDOFValues(currentJoints);
  viewer = OSGViewer::GetOrCreate(env);
  initialize_viewer(viewer);
  for(int i=0; i < controls.size(); i++) {
	  ControlT control = controls[i];
	  cout << "control " << i << ": " << control << endl;
	  cout << "currentJoints: ";
	  for(int j=0; j < control.size(); j++) {
		  currentJoints[j] += control[j];
		  cout << currentJoints[j] << " ";
	  }
	  cout << endl;
	  robot->SetActiveDOFValues(currentJoints);
	  viewer->Idle();
  }

  RaveDestroy();
  return 0;
}

