#include "needle_steering.hpp"

namespace Needle {
  template<typename T, size_t N>
  T* end(T (&ra)[N]) {
    return ra + N;
  }

  void NeedleProblemHelper::AddRotationCost(OptProb& prob) {
    switch (rotation_cost) {
      case UseRotationQuadraticCost: {
        prob.addCost(CostPtr(new RotationQuadraticCost(phivars.col(0), coeff_rotation, shared_from_this())));
        break;
      }
      case UseRotationL1Cost: {
        prob.addCost(CostPtr(new RotationL1Cost(phivars.col(0), coeff_rotation_regularization, shared_from_this())));
        break;
      }
      SWITCH_DEFAULT;
    }
  }

  void NeedleProblemHelper::AddSpeedCost(OptProb& prob) {
    switch (speed_formulation) {
      case ConstantSpeed:
        prob.addCost(CostPtr(new ConstantSpeedCost(Deltavar, coeff_speed, shared_from_this())));
        break;
      case VariableSpeed:
        prob.addCost(CostPtr(new VariableSpeedCost(Deltavars.col(0), coeff_speed, shared_from_this())));
        if (use_speed_deviation_cost) {
          prob.addCost(CostPtr(new SpeedDeviationCost(Deltavars.col(0), Delta_lb, coeff_speed, shared_from_this())));
        }
        break;
      SWITCH_DEFAULT;
    }
  }

  void NeedleProblemHelper::ConfigureProblem(OptProb& prob) {
    CreateVariables(prob);
    InitLocalConfigurations(robot, prob);
    InitTrajectory(prob);
    AddRotationCost(prob);
    AddSpeedCost(prob);
    AddStartConstraint(prob);
    AddGoalConstraint(prob);
    AddSpeedConstraint(prob);
    if (control_constraints) {
      if (this->explicit_controls) {
        AddControlConstraint(prob);
      } else {
        AddPoseConstraint(prob);
      }
    }
    AddCollisionConstraint(prob);
  }

  void NeedleProblemHelper::InitOptimizeVariables(OptimizerT& opt) {
    if (this->initVec.size() == 0) {
      // Initialize twistvars
      for (int i = 0; i <= T; ++i) {
        for (int j = 0; j < n_dof; ++j) {
          this->initVec.push_back(0.);
        }
      }
      // Initialize phivars
      for (int i = 0; i < T; ++i) {
        this->initVec.push_back(0.);
      }
      switch (speed_formulation) {
        case ConstantSpeed:
          // Initialize Delta
          initVec.push_back(Delta_lb);
          break;
        case VariableSpeed:
          for (int i = 0; i < T; ++i) {
            this->initVec.push_back(Delta_lb);
          }
          break;
        SWITCH_DEFAULT;
      }
      // Initialize time frame radii
      if (curvature_constraint == BoundedRadius) {
        for (int i = 0; i < T; ++i) {
          switch (curvature_formulation) {
            case UseCurvature:
              this->initVec.push_back(1.0 / r_min);
              break;
            case UseRadius:
              this->initVec.push_back(r_min);
              break;
            SWITCH_DEFAULT;
          }
        }
      }
    }
    opt.initialize(this->initVec);
  }

  Matrix4d NeedleProblemHelper::TransformPose(const Matrix4d& pose, double phi, double Delta, double curvature_or_radius) const {
    double theta;
    switch (curvature_formulation) {
      case NeedleProblemHelper::UseCurvature:
        theta = Delta * curvature_or_radius;
        break;
      case NeedleProblemHelper::UseRadius:
        theta = Delta / curvature_or_radius;
        break;
      SWITCH_DEFAULT;
    }
    switch (formulation) {
      case NeedleProblemHelper::Form1: {
        Vector6d w; w << 0, 0, 0, 0, 0, phi;
        Vector6d v; v << 0, 0, Delta, theta, 0, 0;
        return pose * expUp(w) * expUp(v);
      }
      case NeedleProblemHelper::Form2: {
        Vector6d w; w << 0, 0, Delta, theta, 0, phi;
        return pose * expUp(w);
      }
      SWITCH_DEFAULT;
    }
  }

  double NeedleProblemHelper::GetPhi(const DblVec& x, int i) const {
    return x[phivars.row(i)[0].var_rep->index];
  }

  double NeedleProblemHelper::GetDelta(const DblVec& x, int i) const {
    switch (speed_formulation) {
      case ConstantSpeed:
        return x[Deltavar.var_rep->index];
      case VariableSpeed:
        return x[Deltavars.row(i)[0].var_rep->index];
      SWITCH_DEFAULT;
    }
  }

  double NeedleProblemHelper::GetCurvatureOrRadius(const DblVec& x, int i) const {
    switch (curvature_formulation) {
      case UseCurvature: 
        return GetCurvature(x, i);
      case UseRadius:
        return GetRadius(x, i);
      SWITCH_DEFAULT;
    }
  }
  
  double NeedleProblemHelper::GetCurvature(const DblVec& x, int i) const {
    assert (curvature_formulation == UseCurvature);
    switch (curvature_constraint) {
      case ConstantRadius:
        return 1.0 / r_min;
      case BoundedRadius:
        return x[curvature_or_radius_vars.row(i)[0].var_rep->index];
      SWITCH_DEFAULT;
    }
  }

  double NeedleProblemHelper::GetRadius(const DblVec& x, int i) const {
    assert (curvature_formulation == UseRadius);
    switch (curvature_constraint) {
      case ConstantRadius:
        return r_min;
      case BoundedRadius:
        return x[curvature_or_radius_vars.row(i)[0].var_rep->index];
      SWITCH_DEFAULT;
    }
  }

  #ifdef NEEDLE_TEST
  void NeedleProblemHelper::checkAlignment(DblVec& x) {
    double diff = 0.;
    for (int i = 0; i < T; ++i) {
      double phi = GetPhi(x, i);
      double Delta = GetDelta(x, i);
      double curvature_or_radius = GetCurvatureOrRadius(x, i);
      diff += (local_configs[i+1]->pose_ - TransformPose(local_configs[i]->pose_, phi, Delta, curvature_or_radius)).norm();
    }
  }
  #endif

  void NeedleProblemHelper::OptimizerCallback(OptProb*, DblVec& x) {
    
    switch (method) {
      case Colocation: {
        MatrixXd twistvals = getTraj(x, twistvars);
        for (int i = 0; i < local_configs.size(); ++i) {
          local_configs[i]->pose = local_configs[i]->pose * expUp(twistvals.row(i));
        }
        #ifdef NEEDLE_TEST
        checkAlignment(x);
        #endif
        setVec(x, twistvars.m_data, DblVec(twistvars.size(), 0));
        break;
      }
      case Shooting: {
        // execute the control input to set local configuration poses
        local_configs[0]->pose = expUp(start);
        for (int i = 0; i < T; ++i) {
          double phi = GetPhi(x, i);
          double Delta = GetDelta(x, i);
          double curvature_or_radius = GetCurvatureOrRadius(x, i);
          local_configs[i+1]->pose = TransformPose(local_configs[i]->pose, phi, Delta, curvature_or_radius);
        }
        setVec(x, twistvars.m_data, DblVec(twistvars.size(), 0));
        break;
      }
      SWITCH_DEFAULT;
    }
  }

  void NeedleProblemHelper::ConfigureOptimizer(OptimizerT& opt) {
    opt.max_iter_ = 500;    
    opt.improve_ratio_threshold_ = this->improve_ratio_threshold;
    opt.trust_shrink_ratio_ = this->trust_shrink_ratio;
    opt.trust_expand_ratio_ = this->trust_expand_ratio;
    opt.record_trust_region_history_ = this->record_trust_region_history;
    opt.max_merit_coeff_increases_ = this->max_merit_coeff_increases;
    opt.merit_error_coeff_ = this->dynamics_coeff;

    InitOptimizeVariables(opt);

    opt.addCallback(boost::bind(&Needle::NeedleProblemHelper::OptimizerCallback, this, _1, _2));

    if (this->plotting || this->plot_final_result) {
      this->plotter.reset(new Needle::TrajPlotter(this->local_configs, this->twistvars));
    }
    if (this->plotting) {
      opt.addCallback(boost::bind(&Needle::TrajPlotter::OptimizerCallback, boost::ref(this->plotter), _1, _2, shared_from_this()));
    }
  }

  void NeedleProblemHelper::CreateVariables(OptProb& prob) {
    // Time frame varies from 0 to T instead of from 0 to T-1
    AddVarArray(prob, T+1, n_dof, "twist", twistvars);
    AddVarArray(prob, T, 1, -PI, PI, "phi", phivars);
    Delta_lb = (goal.topRows(3) - start.topRows(3)).norm() / T / r_min;
    switch (speed_formulation) {
      case ConstantSpeed:
        //Deltavar = prob.createVariables(singleton<string>("Delta"), singleton<double>(Delta_lb),singleton<double>(INFINITY))[0];
        Deltavar = prob.createVariables(singleton<string>("Delta"), singleton<double>(Delta_lb),singleton<double>(INFINITY))[0];
        break;
      case VariableSpeed:
        //AddVarArray(prob, T, 1, Delta_lb * 0.5, INFINITY, "speed", Deltavars); // TODO should we have a resonable upper bound for this?
        AddVarArray(prob, T, 1, Delta_lb*0.5, INFINITY, "speed", Deltavars); // TODO should we have a resonable upper bound for this?
        //AddVarArray(prob, T, 1, Delta_lb*2, INFINITY, "speed", Deltavars); // TODO should we have a resonable upper bound for this?
        break;
      SWITCH_DEFAULT;
    }
    // Only the twist variables are incremental (i.e. their trust regions should be around zero)
    //prob.setIncremental(twistvars.flatten());
    if (curvature_constraint == BoundedRadius) {
      switch (curvature_formulation) {
        case UseCurvature:
          AddVarArray(prob, T, 1, 0.01, 1. / r_min, "curvature", curvature_or_radius_vars);
          break;
        case UseRadius:
          AddVarArray(prob, T, 1, r_min, 100.0, "radius", curvature_or_radius_vars);
          break;
        SWITCH_DEFAULT;
      }
    }
  }

  void NeedleProblemHelper::InitLocalConfigurations(const KinBodyPtr robot, OptProb& prob) {
    for (int i = 0; i <= T; ++i) {
      local_configs.push_back(LocalConfigurationPtr(new LocalConfiguration(robot)));
    }
  }

  void NeedleProblemHelper::InitTrajectory(OptProb& prob) {
    MatrixXd initTraj(T+1, n_dof);
    for (int idof = 0; idof < n_dof; ++idof) {
      initTraj.col(idof) = VectorXd::LinSpaced(T+1, start[idof], goal[idof]);
    }
    for (int i = 0; i <= T; ++i) {
      local_configs[i]->pose = expUp(initTraj.row(i));
    }
  }

  void NeedleProblemHelper::AddStartConstraint(OptProb& prob) {
    VarVector vars = twistvars.row(0);
    VectorOfVectorPtr f(new Needle::PositionError(local_configs[0], start, shared_from_this()));
    Vector6d coeffs; coeffs << 1., 1., 1., this->coeff_orientation_error, this->coeff_orientation_error, this->coeff_orientation_error;
    prob.addConstraint(ConstraintPtr(new ConstraintFromFunc(f, vars, coeffs, EQ, "entry")));
    dynamics_constraints.push_back(prob.getConstraints().back());
  }

  void NeedleProblemHelper::AddGoalConstraint(OptProb& prob) {
    VarVector vars = twistvars.row(T);
    VectorOfVectorPtr f(new Needle::PositionError(local_configs[T], goal, shared_from_this()));
    Vector6d coeffs; coeffs << 1., 1., 1., 0., 0., 0.;
    prob.addConstraint(ConstraintPtr(new ConstraintFromFunc(f, vars, coeffs, EQ, "goal")));
    dynamics_constraints.push_back(prob.getConstraints().back());
  }

  void NeedleProblemHelper::AddSpeedConstraint(OptProb& prob) {
    if (speed_formulation == VariableSpeed && use_speed_deviation_constraint) {
      VectorOfVectorPtr f(new Needle::SpeedDeviationError(Delta_lb, shared_from_this()));
      Vector1d coeffs = Vector1d::Ones();
      prob.addConstraint(ConstraintPtr(new ConstraintFromFunc(f, Deltavars.col(0), coeffs, INEQ, "speed deviation")));
    }
  }

  void NeedleProblemHelper::AddControlConstraint(OptProb& prob) {
    for (int i = 0; i < T; ++i) {
      VarVector vars = concat(concat(twistvars.row(i), twistvars.row(i+1)), phivars.row(i));
      switch (speed_formulation) {
        case ConstantSpeed:
          vars.push_back(Deltavar);
          break;
        case VariableSpeed:
          vars = concat(vars, Deltavars.row(i));
          break;
        SWITCH_DEFAULT;
      }
      if (curvature_constraint == BoundedRadius) {
        vars = concat(vars, curvature_or_radius_vars.row(i));
      }
      VectorOfVectorPtr f(new Needle::ControlError(local_configs[i], local_configs[i+1], shared_from_this()));
      VectorXd coeffs = VectorXd::Ones(boost::static_pointer_cast<Needle::ControlError>(f)->outputSize());
      coeffs.tail<3>() *= this->coeff_orientation_error;
      prob.addConstraint(ConstraintPtr(new ConstraintFromFunc(f, vars, coeffs, EQ, (boost::format("control%i")%i).str())));
      dynamics_constraints.push_back(prob.getConstraints().back());
    }
  }

  void NeedleProblemHelper::AddPoseConstraint(OptProb& prob) {
    for (int i = 0; i < T; ++i) {
      VarVector vars = concat(twistvars.row(i), twistvars.row(i+1));
      switch (speed_formulation) {
        case ConstantSpeed:
          vars.push_back(Deltavar);
          break;
        case VariableSpeed:
          vars = concat(vars, Deltavars.row(i));
          break;
        SWITCH_DEFAULT;
      }
      if (curvature_constraint == BoundedRadius) {
        vars = concat(vars, curvature_or_radius_vars.row(i));
      }
      VectorOfVectorPtr f(new Needle::PoseError(local_configs[i], local_configs[i+1], shared_from_this()));
      VectorXd coeffs = VectorXd::Ones(6);
      coeffs(3) = coeffs(4) = 0;
      coeffs.tail<3>() *= this->coeff_orientation_error;
      prob.addConstraint(ConstraintPtr(new ConstraintFromFunc(f, vars, coeffs, EQ, (boost::format("control%i")%i).str())));
      dynamics_constraints.push_back(prob.getConstraints().back());
    }
  }

  void NeedleProblemHelper::AddCollisionConstraint(OptProb& prob) {
    Str2Dbl tag2dist_pen(collision_dist_pen), tag2coeff(10);//collision_coeff); // the coefficient is essentially multiplied by opt.merit_error_coeff_
    cout << "dist pen: " << collision_dist_pen << endl;
    cout << "collision coeff: " << 10 << endl;//collision_coeff / dynamics_coeff * 10 << endl;
    for (int i = 0; i < ignored_kinbody_names.size(); ++i) {
      tag2coeff.insert( std::pair<string, double>(ignored_kinbody_names[i], 0.0) );
    }
    if (continuous_collision) {
      for (int i = 0; i < T; ++i) {
        prob.addConstraint(ConstraintPtr(new CollisionTaggedConstraint(tag2dist_pen, tag2coeff, local_configs[i], local_configs[i+1], twistvars.row(i), twistvars.row(i+1))));
        collision_constraints.push_back(prob.getConstraints().back());
      }
    } else {
      for (int i = 0; i <= T; ++i) {
        prob.addConstraint(ConstraintPtr(new CollisionTaggedConstraint(tag2dist_pen, tag2coeff, local_configs[i], twistvars.row(i))));
        collision_constraints.push_back(prob.getConstraints().back());
      }
    }

    EnvironmentBasePtr env = local_configs[0]->GetEnv();
    vector<KinBodyPtr> bodies; env->GetBodies(bodies);
    CollisionChecker::GetOrCreate(*env)->SetContactDistance(collision_dist_pen + 0.05);

    for (int i=0; i < bodies.size(); ++i) {
      if (bodies[i]->GetName() == "KinBodyProstate" ||
          bodies[i]->GetName() == "KinBodyDermis" ||
          bodies[i]->GetName() == "KinBodyEpidermis" ||
          bodies[i]->GetName() == "KinBodyHypodermis") {
        CollisionChecker::GetOrCreate(*env)->ExcludeCollisionPair(*bodies[i]->GetLinks()[0], *robot->GetLinks()[0]);
      }
    }
  }

  DblVec NeedleProblemHelper::EvaluateDynamicsViolations(const DblVec& x, Model* model) const {
    DblVec out(this->dynamics_constraints.size());
    for (size_t i = 0; i < this->dynamics_constraints.size(); ++i) {
      out[i] = this->dynamics_constraints[i]->violation(x, model);
    }
    return out;
  }

  DblVec NeedleProblemHelper::EvaluateCollisionViolations(const DblVec& x, Model* model) const {
    DblVec out(this->collision_constraints.size());
    for (size_t i = 0; i < this->collision_constraints.size(); ++i) {
      out[i] = this->collision_constraints[i]->violation(x, model);
    }
    return out;
  }

  bool NeedleProblemHelper::HasDynamicsViolations(const DblVec& x, Model* model) const {
    DblVec viols = EvaluateDynamicsViolations(x, model);
    return viols.size() > 0 && vecMax(viols) > 1e-4;
  }

  bool NeedleProblemHelper::HasCollisionViolations(const DblVec& x, Model* model) const {
    DblVec viols = EvaluateCollisionViolations(x, model);
    cout << "total violations: " << vecMax(viols) << endl;
    for (int i = 0; i < x.size(); ++i) {
      cout << x[i] << " ";
    }cout << endl;
    return viols.size() > 0 && vecMax(viols) > 1e-4 * (collision_coeff / dynamics_coeff * 10);
  }

  void NeedleProblemHelper::Initialize(int argc, char** argv) {

    this->plotting=false;
    this->verbose=false;
    this->plot_final_result=false;
    this->env_transparency = 0.5;

    this->T = 25;
    this->r_min = 2.98119536;
    this->n_dof = 6;

    this->formulation = NeedleProblemHelper::Form1;
    this->curvature_constraint = NeedleProblemHelper::ConstantRadius;
    this->speed_formulation = NeedleProblemHelper::ConstantSpeed;
    this->method = NeedleProblemHelper::Colocation;
    this->curvature_formulation = NeedleProblemHelper::UseRadius;
    this->rotation_cost = NeedleProblemHelper::UseRotationQuadraticCost;
    this->use_speed_deviation_constraint = false;
    this->use_speed_deviation_cost = false;
    this->continuous_collision = true;
    this->explicit_controls = true;
    this->control_constraints = true;

    // parameters for the optimizer
    this->improve_ratio_threshold = 0.1;//0.25;
    this->trust_shrink_ratio = 0.9;//0.7;
    this->trust_expand_ratio = 1.3;//1.2;
    this->record_trust_region_history = false;
    
    double start_vec_array[] = {-6.60848, 12.6176, -8.06651, 2.53666, -0.868663, 1.31701};//{-12.82092, 6.80976, 0.06844, 0, 0, 0};//{0, 0, 0, 0, 0, 0};
    double goal_vec_array[] = {-3.21932, 6.87362, -1.21877, 0, 0, 0};//{0, 0.896343312427, 7.49334469032, 0, 0, 0};

    this->coeff_rotation = 1.;
    this->coeff_speed = 1.;
    this->coeff_rotation_regularization = 0.1;
    this->coeff_orientation_error = 1;
    this->collision_dist_pen = 0.05;
    this->collision_coeff = 10;

    vector<double> start_vec(start_vec_array, start_vec_array + n_dof);
    vector<double> goal_vec(goal_vec_array, goal_vec_array + n_dof);
    {
      Config config;
      config.add(new Parameter<bool>("plotting", &this->plotting, "plotting"));
      config.add(new Parameter<bool>("plot_final_result", &this->plot_final_result, "plot_final_result"));
      config.add(new Parameter<bool>("verbose", &this->verbose, "verbose"));
      config.add(new Parameter<double>("env_transparency", &this->env_transparency, "env_transparency"));
      config.add(new Parameter<int>("T", &this->T, "T"));
      config.add(new Parameter<int>("formulation", &this->formulation, "formulation"));
      config.add(new Parameter<int>("curvature_constraint", &this->curvature_constraint, "curvature_constraint"));
      config.add(new Parameter<int>("method", &this->method, "method"));
      config.add(new Parameter<int>("curvature_formulation", &this->curvature_formulation, "curvature_formulation"));
      config.add(new Parameter<int>("speed_formulation", &this->speed_formulation, "speed_formulation"));
      config.add(new Parameter<int>("rotation_cost", &this->rotation_cost, "rotation_cost"));
      config.add(new Parameter<double>("coeff_rotation_regularization", &this->coeff_rotation_regularization, "coeff_rotation_regularization"));
      config.add(new Parameter<double>("coeff_rotation", &this->coeff_rotation, "coeff_rotation"));
      config.add(new Parameter<double>("coeff_speed", &this->coeff_speed, "coeff_speed"));
      config.add(new Parameter<double>("coeff_orientation_error", &this->coeff_orientation_error, "coeff_orientation_error"));
      config.add(new Parameter<double>("r_min", &this->r_min, "r_min"));
      config.add(new Parameter<double>("improve_ratio_threshold", &this->improve_ratio_threshold, "improve_ratio_threshold"));
      config.add(new Parameter<double>("trust_shrink_ratio", &this->trust_shrink_ratio, "trust_shrink_ratio"));
      config.add(new Parameter<double>("trust_expand_ratio", &this->trust_expand_ratio, "trust_expand_ratio"));
      config.add(new Parameter<double>("collision_dist_pen", &this->collision_dist_pen, "collision_dist_pen"));
      config.add(new Parameter<bool>("use_speed_deviation_constraint", &this->use_speed_deviation_constraint, "use_speed_deviation_constraint"));
      config.add(new Parameter<bool>("use_speed_deviation_cost", &this->use_speed_deviation_cost, "use_speed_deviation_cost"));
      config.add(new Parameter<bool>("record_trust_region_history", &this->record_trust_region_history, "record_trust_region_history"));
      config.add(new Parameter<bool>("explicit_controls", &this->explicit_controls, "explicit_controls"));
      config.add(new Parameter<bool>("continuous_collision", &this->continuous_collision, "continuous_collision"));
      config.add(new Parameter<bool>("control_constraints", &this->control_constraints, "control_constraints"));
      config.add(new Parameter< vector<double> >("s", &start_vec, "s"));
      config.add(new Parameter< vector<double> >("g", &goal_vec, "g"));
      CommandParser parser(config);
      parser.read(argc, argv);
    }

    RaveInitialize(false, this->verbose ? Level_Debug : Level_Info);
    EnvironmentBasePtr env = RaveCreateEnvironment();
    env->StopSimulation();
    OSGViewerPtr viewer;
    if (this->plotting || this->plot_final_result) {
      viewer = OSGViewer::GetOrCreate(env);
      assert(viewer);
    }

    env->Load(string(DATA_DIR) + "/prostate.env.xml");

    if (this->plotting || this->plot_final_result) {
      viewer->SetAllTransparency(this->env_transparency);
    }

    RobotBasePtr robot = GetRobot(*env);

    for (int i = 0; i < n_dof; ++i) this->start[i] = start_vec[i];
    for (int i = 0; i < n_dof; ++i) this->goal[i] = goal_vec[i];

    const char *ignored_kinbody_c_strs[] = { "KinBodyProstate", "KinBodyDermis", "KinBodyEpidermis", "KinBodyHypodermis" };
    this->ignored_kinbody_names = vector<string>(ignored_kinbody_c_strs, end(ignored_kinbody_c_strs));

    
    this->robot = robot;

  }
}
