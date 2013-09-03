#include "needle_steering.hpp"

namespace Needle {
  template<typename T, size_t N>
  T* end(T (&ra)[N]) {
    return ra + N;
  } 

  void NeedleProblemHelper::AddRotationCost(OptProb& prob, NeedleProblemInstancePtr pi) {
    switch (rotation_cost) {
      case UseRotationQuadraticCost: {
        prob.addCost(CostPtr(new RotationQuadraticCost(pi->phivars.col(0), coeff_rotation, shared_from_this())));
        break;
      }
      case UseRotationL1Cost: {
        prob.addCost(CostPtr(new RotationL1Cost(pi->phivars.col(0), coeff_rotation_regularization, shared_from_this())));
        break;
      }
      SWITCH_DEFAULT;
    }
  }

  void NeedleProblemHelper::AddSpeedCost(OptProb& prob, NeedleProblemInstancePtr pi) {
    switch (speed_formulation) {
      case ConstantSpeed:
        prob.addCost(CostPtr(new ConstantSpeedCost(pi->Deltavar, coeff_speed, shared_from_this(), pi)));
        break;
      case VariableSpeed:
        prob.addCost(CostPtr(new VariableSpeedCost(pi->Deltavars.col(0), coeff_speed, shared_from_this())));
        if (use_speed_deviation_cost) {
          prob.addCost(CostPtr(new SpeedDeviationCost(pi->Deltavars.col(0), pi->Delta_lb, coeff_speed, shared_from_this())));
        }
        break;
      SWITCH_DEFAULT;
    }
  }

  void NeedleProblemHelper::ConfigureProblem(OptProb& prob) {
    for (int i = 0; i < n_needles; ++i) {
      NeedleProblemInstancePtr pi(new NeedleProblemInstance());
      pi->start = starts[i];
      pi->goal = goals[i];
      pi->T = Ts[i];
      CreateVariables(prob, pi);
      InitLocalConfigurations(this->robots[i], prob, pi);
      InitTrajectory(prob, pi);
      AddRotationCost(prob, pi);
      AddSpeedCost(prob, pi);
      AddStartConstraint(prob, pi);
      AddGoalConstraint(prob, pi);
      AddSpeedConstraint(prob, pi);
      if (control_constraints) {
        if (this->explicit_controls) {
          AddControlConstraint(prob, pi);
        } else {
          AddPoseConstraint(prob, pi);
        }
      }
      AddCollisionConstraint(prob, pi);
      this->pis.push_back(pi);
    }

    for (int i = 0; i < n_needles; ++i) {
      for (int j = i+1; j < n_needles; ++j) {
        AddSelfCollisionConstraint(prob, pis[i], pis[j]);
      }
    }

    InitializeCollisionEnvironment();
  }

  void NeedleProblemHelper::InitOptimizeVariables(OptimizerT& opt) {
    DblVec initVec;
    for (int k = 0; k < n_needles; ++k) {
      NeedleProblemInstancePtr pi = pis[k];
      if (pi->initVec.size() == 0) {
        // Initialize twistvars
        for (int i = 0; i <= pi->T; ++i) {
          for (int j = 0; j < n_dof; ++j) {
            pi->initVec.push_back(0.);
          }
        }
        // Initialize phivars
        for (int i = 0; i < pi->T; ++i) {
          pi->initVec.push_back(0.);
        }
        switch (speed_formulation) {
          case ConstantSpeed:
            // Initialize Delta
            pi->initVec.push_back(pi->Delta_lb);
            break;
          case VariableSpeed:
            for (int i = 0; i < pi->T; ++i) {
              pi->initVec.push_back(pi->Delta_lb);
            }
            break;
          SWITCH_DEFAULT;
        }
        // Initialize time frame radii
        if (curvature_constraint == BoundedRadius) {
          for (int i = 0; i < pi->T; ++i) {
            switch (curvature_formulation) {
              case UseCurvature:
                pi->initVec.push_back(1.0 / r_min);
                break;
              case UseRadius:
                pi->initVec.push_back(r_min);
                break;
              SWITCH_DEFAULT;
            }
          }
        }

        for (int i = 0; i < pi->initVec.size(); ++i) {
          initVec.push_back(pi->initVec[i]);
        }
      }
    }
    opt.initialize(initVec);
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

  double NeedleProblemHelper::GetPhi(const DblVec& x, int i, NeedleProblemInstancePtr pi) const {
    return x[pi->phivars.row(i)[0].var_rep->index];
  }

  double NeedleProblemHelper::GetDelta(const DblVec& x, int i, NeedleProblemInstancePtr pi) const {
    switch (speed_formulation) {
      case ConstantSpeed:
        return x[pi->Deltavar.var_rep->index];
      case VariableSpeed:
        return x[pi->Deltavars.row(i)[0].var_rep->index];
      SWITCH_DEFAULT;
    }
  }

  double NeedleProblemHelper::GetCurvatureOrRadius(const DblVec& x, int i, NeedleProblemInstancePtr pi) const {
    switch (curvature_formulation) {
      case UseCurvature: 
        return GetCurvature(x, i, pi);
      case UseRadius:
        return GetRadius(x, i, pi);
      SWITCH_DEFAULT;
    }
  }
  
  double NeedleProblemHelper::GetCurvature(const DblVec& x, int i, NeedleProblemInstancePtr pi) const {
    assert (curvature_formulation == UseCurvature);
    switch (curvature_constraint) {
      case ConstantRadius:
        return 1.0 / r_min;
      case BoundedRadius:
        return x[pi->curvature_or_radius_vars.row(i)[0].var_rep->index];
      SWITCH_DEFAULT;
    }
  }

  double NeedleProblemHelper::GetRadius(const DblVec& x, int i, NeedleProblemInstancePtr pi) const {
    assert (curvature_formulation == UseRadius);
    switch (curvature_constraint) {
      case ConstantRadius:
        return r_min;
      case BoundedRadius:
        return x[pi->curvature_or_radius_vars.row(i)[0].var_rep->index];
      SWITCH_DEFAULT;
    }
  }

  void NeedleProblemHelper::OptimizerCallback(OptProb*, DblVec& x) {
    
    switch (method) {
      case Colocation: {
        for (int i = 0; i < n_needles; ++i) {
          MatrixXd twistvals = getTraj(x, pis[i]->twistvars);
          for (int j = 0; j < pis[i]->local_configs.size(); ++j) {
            pis[i]->local_configs[j]->pose = pis[i]->local_configs[j]->pose * expUp(twistvals.row(j));
          }
          setVec(x, pis[i]->twistvars.m_data, DblVec(pis[i]->twistvars.size(), 0));
        }
        break;
      }
      case Shooting: {
        // execute the control input to set local configuration poses
        for (int i = 0; i < n_needles; ++i) {
          pis[i]->local_configs[0]->pose = expUp(pis[i]->start);
          for (int j = 0; j < pis[i]->T; ++j) {
            double phi = GetPhi(x, j, pis[i]);
            double Delta = GetDelta(x, j, pis[i]);
            double curvature_or_radius = GetCurvatureOrRadius(x, j, pis[i]);
            pis[i]->local_configs[j+1]->pose = TransformPose(pis[i]->local_configs[j]->pose, phi, Delta, curvature_or_radius);
            setVec(x, pis[i]->twistvars.m_data, DblVec(pis[i]->twistvars.size(), 0));
          }
        }
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
    opt.merit_error_coeff_ = this->merit_error_coeff;

    InitOptimizeVariables(opt);

    opt.addCallback(boost::bind(&Needle::NeedleProblemHelper::OptimizerCallback, this, _1, _2));
  }

  void NeedleProblemHelper::CreateVariables(OptProb& prob, NeedleProblemInstancePtr pi) {
    // Time frame varies from 0 to T instead of from 0 to T-1
    AddVarArray(prob, pi->T+1, n_dof, "twist", pi->twistvars);
    AddVarArray(prob, pi->T, 1, -PI, PI, "phi", pi->phivars);
    pi->Delta_lb = (pi->goal.topRows(3) - pi->start.topRows(3)).norm() / pi->T / r_min;
    switch (speed_formulation) {
      case ConstantSpeed:
        pi->Deltavar = prob.createVariables(singleton<string>("Delta"), singleton<double>(pi->Delta_lb),singleton<double>(INFINITY))[0];
        break;
      case VariableSpeed:
        //AddVarArray(prob, T, 1, Delta_lb * 0.5, INFINITY, "speed", Deltavars); // TODO should we have a resonable upper bound for this?
        AddVarArray(prob, pi->T, 1, pi->Delta_lb*0.5, INFINITY, "speed", pi->Deltavars); // TODO should we have a resonable upper bound for this?
        //AddVarArray(prob, T, 1, Delta_lb*2, INFINITY, "speed", Deltavars); // TODO should we have a resonable upper bound for this?
        break;
      SWITCH_DEFAULT;
    }
    // Only the twist variables are incremental (i.e. their trust regions should be around zero)
    //prob.setIncremental(twistvars.flatten());
    if (curvature_constraint == BoundedRadius) {
      switch (curvature_formulation) {
        case UseCurvature:
          AddVarArray(prob, pi->T, 1, 0.01, 1. / r_min, "curvature", pi->curvature_or_radius_vars);
          break;
        case UseRadius:
          AddVarArray(prob, pi->T, 1, r_min, 100.0, "radius", pi->curvature_or_radius_vars);
          break;
        SWITCH_DEFAULT;
      }
    }
  }

  void NeedleProblemHelper::InitLocalConfigurations(const KinBodyPtr robot, OptProb& prob, NeedleProblemInstancePtr pi) {
    for (int i = 0; i <= pi->T; ++i) {
      pi->local_configs.push_back(LocalConfigurationPtr(new LocalConfiguration(robot)));
    }
  }

  void NeedleProblemHelper::InitTrajectory(OptProb& prob, NeedleProblemInstancePtr pi) {
    MatrixXd initTraj(pi->T+1, n_dof);
    for (int idof = 0; idof < n_dof; ++idof) {
      initTraj.col(idof) = VectorXd::LinSpaced(pi->T+1, pi->start[idof], pi->goal[idof]);
    }
    for (int i = 0; i <= pi->T; ++i) {
      pi->local_configs[i]->pose = expUp(initTraj.row(i));
    }
  }

  void NeedleProblemHelper::AddStartConstraint(OptProb& prob, NeedleProblemInstancePtr pi) {
    VarVector vars = pi->twistvars.row(0);
    VectorOfVectorPtr f(new Needle::PositionError(pi->local_configs[0], pi->start, this->start_position_error_relax, this->start_orientation_deviation_relax, shared_from_this()));
    Vector6d coeffs; coeffs << 1., 1., 1., this->coeff_orientation_error, this->coeff_orientation_error, this->coeff_orientation_error;
    prob.addConstraint(ConstraintPtr(new ConstraintFromFunc(f, vars, coeffs, EQ, "entry")));
    pi->dynamics_constraints.push_back(prob.getConstraints().back());
  }

  void NeedleProblemHelper::AddGoalConstraint(OptProb& prob, NeedleProblemInstancePtr pi) {
    VarVector vars = pi->twistvars.row(pi->T);
    VectorOfVectorPtr f(new Needle::PositionError(pi->local_configs[pi->T], pi->goal, Vector3d::Zero(), shared_from_this()));
    Vector6d coeffs; 
    if (goal_orientation_constraint) {
      coeffs << 1., 1., 1., this->coeff_orientation_error, this->coeff_orientation_error, this->coeff_orientation_error;
    } else {
      coeffs << 1., 1., 1., 0, 0, 0;
    }
    prob.addConstraint(ConstraintPtr(new ConstraintFromFunc(f, vars, coeffs, EQ, "goal")));
    pi->dynamics_constraints.push_back(prob.getConstraints().back());
  }

  void NeedleProblemHelper::AddSpeedConstraint(OptProb& prob, NeedleProblemInstancePtr pi) {
    if (speed_formulation == VariableSpeed && use_speed_deviation_constraint) {
      VectorOfVectorPtr f(new Needle::SpeedDeviationError(pi->Delta_lb, shared_from_this()));
      Vector1d coeffs = Vector1d::Ones();
      prob.addConstraint(ConstraintPtr(new ConstraintFromFunc(f, pi->Deltavars.col(0), coeffs, INEQ, "speed deviation")));
    }
  }

  void NeedleProblemHelper::AddControlConstraint(OptProb& prob, NeedleProblemInstancePtr pi) {
    for (int i = 0; i < pi->T; ++i) {
      VarVector vars = concat(concat(pi->twistvars.row(i), pi->twistvars.row(i+1)), pi->phivars.row(i));
      switch (speed_formulation) {
        case ConstantSpeed:
          vars.push_back(pi->Deltavar);
          break;
        case VariableSpeed:
          vars = concat(vars, pi->Deltavars.row(i));
          break;
        SWITCH_DEFAULT;
      }
      if (curvature_constraint == BoundedRadius) {
        vars = concat(vars, pi->curvature_or_radius_vars.row(i));
      }
      VectorOfVectorPtr f(new Needle::ControlError(pi->local_configs[i], pi->local_configs[i+1], shared_from_this()));
      VectorXd coeffs = VectorXd::Ones(boost::static_pointer_cast<Needle::ControlError>(f)->outputSize());
      coeffs.tail<3>() *= this->coeff_orientation_error;
      prob.addConstraint(ConstraintPtr(new ConstraintFromFunc(f, vars, coeffs, EQ, (boost::format("control%i")%i).str())));
      pi->dynamics_constraints.push_back(prob.getConstraints().back());
    }
  }

  void NeedleProblemHelper::AddPoseConstraint(OptProb& prob, NeedleProblemInstancePtr pi) {
    for (int i = 0; i < pi->T; ++i) {
      VarVector vars = concat(pi->twistvars.row(i), pi->twistvars.row(i+1));
      switch (speed_formulation) {
        case ConstantSpeed:
          vars.push_back(pi->Deltavar);
          break;
        case VariableSpeed:
          vars = concat(vars, pi->Deltavars.row(i));
          break;
        SWITCH_DEFAULT;
      }
      if (curvature_constraint == BoundedRadius) {
        vars = concat(vars, pi->curvature_or_radius_vars.row(i));
      }
      VectorOfVectorPtr f(new Needle::PoseError(pi->local_configs[i], pi->local_configs[i+1], shared_from_this()));
      VectorXd coeffs = VectorXd::Ones(6);
      coeffs(3) = coeffs(4) = 0;
      coeffs.tail<3>() *= this->coeff_orientation_error;
      prob.addConstraint(ConstraintPtr(new ConstraintFromFunc(f, vars, coeffs, EQ, (boost::format("control%i")%i).str())));
      pi->dynamics_constraints.push_back(prob.getConstraints().back());
    }
  }

  void NeedleProblemHelper::AddCollisionConstraint(OptProb& prob, NeedleProblemInstancePtr pi) {
    if (continuous_collision) {
      for (int i = 0; i < pi->T; ++i) {
        prob.addConstraint(ConstraintPtr(new CollisionConstraint(collision_dist_pen, collision_coeff, pi->local_configs[i], pi->local_configs[i+1], pi->twistvars.row(i), pi->twistvars.row(i+1))));
        pi->collision_constraints.push_back(prob.getConstraints().back());
      }
    } else {
      for (int i = 0; i <= pi->T; ++i) {
        prob.addConstraint(ConstraintPtr(new CollisionConstraint(collision_dist_pen, collision_coeff, pi->local_configs[i], pi->twistvars.row(i))));
        pi->collision_constraints.push_back(prob.getConstraints().back());
      }
    }
  }

  void NeedleProblemHelper::AddSelfCollisionConstraint(OptProb& prob, NeedleProblemInstancePtr piA, NeedleProblemInstancePtr piB) {
    if (continuous_collision) {
      for (int i = 0; i < piA->T; ++i) {
        for (int j = 0; j < piB->T; ++j) {
          prob.addConstraint(ConstraintPtr(new CollisionConstraint(collision_dist_pen, collision_coeff, piA->local_configs[i], piA->local_configs[i+1], piB->local_configs[j], piB->local_configs[j+1], piA->twistvars.row(i), piA->twistvars.row(i+1), piB->twistvars.row(j), piB->twistvars.row(j+1))));
          self_collision_constraints.push_back(prob.getConstraints().back());
        }
      }
    } else {
      throw std::runtime_error("Not implemented");
    }
  }

  void NeedleProblemHelper::InitializeCollisionEnvironment() {
    EnvironmentBasePtr env = pis[0]->local_configs[0]->GetEnv();
    vector<KinBodyPtr> bodies; env->GetBodies(bodies);
    CollisionChecker::GetOrCreate(*env)->SetContactDistance(collision_dist_pen + 0.05);

    for (int k = 0; k < n_needles; ++k) {
      for (int i=0; i < bodies.size(); ++i) {
        if (std::find(ignored_kinbody_names.begin(), ignored_kinbody_names.end(), bodies[i]->GetName()) != ignored_kinbody_names.end()) {
          BOOST_FOREACH(const KinBody::LinkPtr& robot_link, robots[k]->GetLinks()) {
            BOOST_FOREACH(const KinBody::LinkPtr& body_link, bodies[i]->GetLinks()) {
              CollisionChecker::GetOrCreate(*env)->ExcludeCollisionPair(*body_link, *robot_link);//*bodies[i]->GetLinks()[0], *robots[k]->GetLinks()[0]);
            }
          }
        }
      }
    }
  }

  void NeedleProblemHelper::InitParameters() {
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
    this->goal_orientation_constraint = true;

    // parameters for the optimizer
    this->improve_ratio_threshold = 0.1;
    this->trust_shrink_ratio = 0.9;
    this->trust_expand_ratio = 1.3;
    this->record_trust_region_history = false;
    this->merit_error_coeff = 10;
    this->max_merit_coeff_increases = 10;

    this->coeff_rotation = 1.;
    this->coeff_speed = 1.;
    this->coeff_rotation_regularization = 0.1;
    this->coeff_orientation_error = 1;
    this->collision_dist_pen = 0.05;
    this->collision_coeff = 10;

    this->start_position_error_relax = Vector3d(1.25, 1.25, 0);
    this->start_orientation_deviation_relax = 0.0873;

    const char *ignored_kinbody_c_strs[] = { "KinBodyProstate", "KinBodyDermis", "KinBodyEpidermis", "KinBodyHypodermis" };
    this->ignored_kinbody_names = vector<string>(ignored_kinbody_c_strs, end(ignored_kinbody_c_strs));
  }

  void NeedleProblemHelper::InitParametersFromConsole(int argc, char** argv) {
    Clear();
    InitParameters();
    double start_position_error_relax_x = this->start_a;
    Config config;
    //config.add(new Parameter<int>("T", &this->T, "T"));
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
    config.add(new Parameter<double>("merit_error_coeff", &this->merit_error_coeff, "merit_error_coeff"));
    config.add(new Parameter<bool>("use_speed_deviation_constraint", &this->use_speed_deviation_constraint, "use_speed_deviation_constraint"));
    config.add(new Parameter<bool>("use_speed_deviation_cost", &this->use_speed_deviation_cost, "use_speed_deviation_cost"));
    config.add(new Parameter<bool>("record_trust_region_history", &this->record_trust_region_history, "record_trust_region_history"));
    config.add(new Parameter<bool>("explicit_controls", &this->explicit_controls, "explicit_controls"));
    config.add(new Parameter<bool>("continuous_collision", &this->continuous_collision, "continuous_collision"));
    config.add(new Parameter<bool>("control_constraints", &this->control_constraints, "control_constraints"));
    config.add(new Parameter<bool>("goal_orientation_constraint", &this->goal_orientation_constraint, "goal_orientation_constraint"));
    
    CommandParser parser(config);
    parser.read(argc, argv, true);
  }

  void NeedleProblemHelper::Clear() {
    starts.clear();
    goals.clear();
    if (n_needles > 0) {
      EnvironmentBasePtr env = pis[0]->local_configs[0]->GetEnv();
      for (int i = 0; i < n_needles; ++i) {
        env->Remove(robots[i]);
        trajopt::RemoveUserData(*robots[i], "bt");
      }
      robots.clear();
      self_collision_constraints.clear();
      pis.clear();
      n_needles = 0;
    }
  }

  void NeedleProblemHelper::AddNeedlesToBullet(OptimizerT& opt) {
    for (int i = 0; i < n_needles; ++i) {
      AddNeedleToBullet(pis[i], opt);
    }
  }

  void NeedleProblemHelper::AddNeedleToBullet(NeedleProblemInstancePtr pi, OptimizerT& opt) {
    
    EnvironmentBasePtr env = pi->local_configs[0]->GetEnv();
    DblVec x = opt.x();
    MatrixXd twistvals = getTraj(x, pi->twistvars);
    boost::shared_ptr<BulletCollisionChecker> cc = boost::dynamic_pointer_cast<BulletCollisionChecker>(CollisionChecker::GetOrCreate(*env));
    for (int i = 0; i < pi->T; ++i) {
      vector<KinBody::LinkPtr> links;
      vector<int> inds;
      pi->local_configs[i]->GetAffectedLinks(links, true, inds);
      cc->AddCastHullShape(*pi->local_configs[i], *pi->local_configs[i+1], links, toDblVec(twistvals.row(i)), toDblVec(twistvals.row(i+1)));
    }
  }

  vector<VectorXd> NeedleProblemHelper::GetSolutions(OptimizerT& opt) {
    vector<VectorXd> ret;
    for (int i = 0; i < n_needles; ++i) {
      ret.push_back(pis[i]->GetSolution(opt));
    }
    return ret;
  }

  void NeedleProblemHelper::SetSolutions(const vector<VectorXd>& sol, OptimizerT& opt) {
    assert (sol.size() == n_needles);
    for (int i = 0; i < n_needles; ++i) {
      pis[i]->SetSolution(sol[i], opt);
    }
  }

  vector<VectorXd> NeedleProblemHelper::GetSolutionsWithoutFirstTimestep(const vector<VectorXd>& sol) {
    vector<VectorXd> ret;
    assert (sol.size() == n_needles);
    if (this->Ts.front() > 1) {
      ret.push_back(pis.front()->GetSolutionWithoutFirstTimestep(sol.front()));
    }
    for (int i = 1; i < n_needles; ++i) {
      ret.push_back(sol[i]);
    }
    return ret;
  }
}
