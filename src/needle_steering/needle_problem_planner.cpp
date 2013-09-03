#include "needle_steering.hpp"

namespace Needle {

  inline double rndnum() {
	  return ((double) random()) / RAND_MAX;
  }

  inline double normal() {
    double u_1 = 0;
    while (u_1 == 0) {
      u_1 = rndnum();
    }
    double u_2 = 0;
    while (u_2 == 0) {
      u_2 = rndnum();
    }
    return sqrt(-2*log(u_1)) * sin(2*PI*u_2);
  }

  NeedleProblemPlanner::NeedleProblemPlanner(int argc, char **argv) :
    argc(argc),
    argv(argv),
    helper(new NeedleProblemHelper()),
    plotting(false),
    verbose(false),
    plot_final_result(false),
    env_transparency(0.1),
    data_dir(get_current_directory() + "/../data") {

    vector<string> start_string_vec;
    start_string_vec.push_back("-11.67067,5.54934,0,0,0.78,0");
    start_string_vec.push_back("-11.17067,5.04934,0,0,0.78,0");
    vector<string> goal_string_vec;
    goal_string_vec.push_back("-2.71912,8.00334,-1.12736,0,0.78,0");
    goal_string_vec.push_back("-3.926,7.78291,-1.08402,0,0.75,0");

    int T = 25;
    
    Config config;
    config.add(new Parameter<bool>("plotting", &this->plotting, "plotting"));
    config.add(new Parameter<bool>("plot_final_result", &this->plot_final_result, "plot_final_result"));
    config.add(new Parameter<bool>("verbose", &this->verbose, "verbose"));
    config.add(new Parameter<double>("env_transparency", &this->env_transparency, "env_transparency"));
    config.add(new Parameter<string>("data_dir", &this->data_dir, "data_dir"));
    config.add(new Parameter<string>("env_file_path", &this->env_file_path, "env_file_path"));
    config.add(new Parameter<string>("robot_file_path", &this->robot_file_path, "robot_file_path"));
    config.add(new Parameter< vector<string> >("start_vec", &start_string_vec, "s"));
    config.add(new Parameter< vector<string> >("goal_vec", &goal_string_vec, "g"));
    config.add(new Parameter<int>("T", &T, "T"));
    CommandParser parser(config);
    parser.read(argc, argv, true);

    if (this->env_file_path.length() == 0) {
      this->env_file_path = data_dir + "/prostate.env.xml";
    }
    if (this->robot_file_path.length() == 0) {
      this->robot_file_path = data_dir + "/needlebot.xml";
    }

    if (start_string_vec.size() != goal_string_vec.size()) {
      throw std::runtime_error("The number of start and goal vectors must be the same!");
    }

    if (start_string_vec.size() == 0) {
      throw std::runtime_error("You must provide at least 1 start and 1 goal vector.");
    }

    RaveInitialize(false, verbose ? Level_Debug : Level_Info);
    this->env = RaveCreateEnvironment();
    this->env->StopSimulation();

    OSGViewerPtr viewer;
    if (this->plotting || this->plot_final_result) {
      viewer = OSGViewer::GetOrCreate(env);
      assert(viewer);
    }

    this->env->Load(this->env_file_path);

    if (this->plotting || this->plot_final_result) {
      viewer->SetAllTransparency(this->env_transparency);
    }


    this->n_needles = start_string_vec.size();
    this->starts.clear();
    this->goals.clear();

    for (int i = 0; i < n_needles; ++i) {
      DblVec start_vec;
      DblVec goal_vec;
      strtk::parse(start_string_vec[i], ",", start_vec);
      strtk::parse(goal_string_vec[i], ",", goal_vec);
      Vector6d start = toVectorXd(start_vec), goal = toVectorXd(goal_vec);
      Matrix4d start_pose, goal_pose;
      start_pose.topLeftCorner<3, 3>() = rotMat(start.tail<3>());
      start_pose.topRightCorner<3, 1>() = start.head<3>();
      start_pose.bottomLeftCorner<1, 3>() = Vector3d::Zero();
      start_pose(3, 3) = 1;
      goal_pose.topLeftCorner<3, 3>() = rotMat(goal.tail<3>());
      goal_pose.topRightCorner<3, 1>() = goal.head<3>();
      goal_pose.bottomLeftCorner<1, 3>() = Vector3d::Zero();
      goal_pose(3, 3) = 1;
      starts.push_back(logDown(start_pose));
      goals.push_back(logDown(goal_pose));
    }

    for (int i = 0; i < n_needles; ++i) {
      this->Ts.push_back(T);
    }

  }

  vector<VectorXd> NeedleProblemPlanner::Solve(const vector<VectorXd>& initial) {
    trajopt::SetUserData(*this->env, "trajopt_cc", OpenRAVE::UserDataPtr());
    helper->InitParametersFromConsole(this->argc, this->argv);
    helper->n_needles = this->n_needles;
    helper->starts = this->starts;
    helper->goals = this->goals;
    helper->Ts = this->Ts;
    for (int i = 0; i < n_needles; ++i) {
      helper->robots.push_back(this->env->ReadRobotURI(RobotBasePtr(), this->robot_file_path));
      this->env->Add(helper->robots.back(), true);
    }

    OptProbPtr prob(new OptProb());
    helper->ConfigureProblem(*prob);
    OptimizerT opt(prob);
    helper->ConfigureOptimizer(opt);

    if (initial.size() == helper->n_needles) {
      cout << "set solution initialization" << endl;
      helper->SetSolutions(initial, opt);
    }

    if (this->plotting || this->plot_final_result) {
      this->plotter.reset(new Needle::TrajPlotter(helper->pis));
    }
    if (this->plotting) {
      opt.addCallback(boost::bind(&Needle::TrajPlotter::OptimizerCallback, boost::ref(this->plotter), _1, _2, helper));
    }

    opt.optimize();

    this->x = opt.x();
    
    return helper->GetSolutions(opt);
  }

  Vector6d NeedleProblemPlanner::PerturbState(const Vector6d& state) {
    Vector6d ret;
    for (int i = 0; i < 3; ++i) {
      ret(i) += normal() * 0.01;
    }
    return ret;
  }

  vector<Vector6d> NeedleProblemPlanner::SimulateExecution(const vector<Vector6d>& current_states) {
    vector<Vector6d> ret;
    for (int i = 0; i < current_states.size() - n_needles; ++i) { // leave as it is
      ret.push_back(current_states[i]);
    }

    double phi = helper->GetPhi(this->x, 0, helper->pis.front());
    double Delta = helper->GetDelta(this->x, 0, helper->pis.front());
    double curvature_or_radius = helper->GetCurvatureOrRadius(this->x, 0, helper->pis.front());
    Vector6d state_to_change = current_states[current_states.size() - n_needles];
    ret.push_back(PerturbState(logDown(helper->TransformPose(expUp(state_to_change), phi, Delta, curvature_or_radius))));

    for (int i = current_states.size() - n_needles + 1, j = 0; i < current_states.size(); ++i, ++j) { // leave as it is
      ret.push_back(current_states[i]);
    }

    vector<Vector6d> prev_starts = starts;
    vector<Vector6d> prev_goals = goals;

    starts.clear();
    goals.clear();

    if (Ts.front() > 1) {
      starts.push_back(ret[current_states.size() - n_needles]);
      goals.push_back(prev_goals.front());
      --Ts.front();
    } else {
      // get rid of first time step
      vector<int> new_Ts;
      for (int i = 1; i < Ts.size(); ++i) {
        new_Ts.push_back(Ts[i]);
      }
      this->Ts = new_Ts;
      --n_needles;
    }

    for (int i = current_states.size() - n_needles, j = 1; i < current_states.size(); ++i, ++j) {
      this->starts.push_back(ret[i]);
      this->goals.push_back(prev_goals[j]);
    }

    return ret;
  }

  vector<VectorXd> NeedleProblemPlanner::GetSolutionsWithoutFirstTimestep(const vector<VectorXd>& sol) {
    if (sol.size() > 0) {
      return helper->GetSolutionsWithoutFirstTimestep(sol);
    } else {
      return sol;
    }
  }

  bool NeedleProblemPlanner::Finished() const {
    return this->Ts.size() == 0;
  }

  NeedleProblemPlanner::~NeedleProblemPlanner() {
    RaveDestroy();
  }
}
