#include "needle_steering.hpp"

namespace Needle {
  
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
        prob.addCost(CostPtr(new ConstantSpeedCost(pi->Deltavar, coeff_speed, shared_from_this())));
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
      CreateVariables(prob, pi);
      InitLocalConfigurations(robot, prob, pi);
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

    InitializeCollisionEnvironment();
  }

  void NeedleProblemHelper::InitOptimizeVariables(OptimizerT& opt) {
    DblVec initVec;
    for (int k = 0; k < n_needles; ++k) {
      NeedleProblemInstancePtr pi = pis[k];
      if (pi->initVec.size() == 0) {
        // Initialize twistvars
        for (int i = 0; i <= T; ++i) {
          for (int j = 0; j < n_dof; ++j) {
            pi->initVec.push_back(0.);
          }
        }
        // Initialize phivars
        for (int i = 0; i < T; ++i) {
          pi->initVec.push_back(0.);
        }
        switch (speed_formulation) {
          case ConstantSpeed:
            // Initialize Delta
            pi->initVec.push_back(pi->Delta_lb);
            break;
          case VariableSpeed:
            for (int i = 0; i < T; ++i) {
              pi->initVec.push_back(pi->Delta_lb);
            }
            break;
          SWITCH_DEFAULT;
        }
        // Initialize time frame radii
        if (curvature_constraint == BoundedRadius) {
          for (int i = 0; i < T; ++i) {
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
          for (int j = 0; j < T; ++j) {
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
    AddVarArray(prob, T+1, n_dof, "twist", pi->twistvars);
    AddVarArray(prob, T, 1, -PI, PI, "phi", pi->phivars);
    pi->Delta_lb = (pi->goal.topRows(3) - pi->start.topRows(3)).norm() / T / r_min;
    switch (speed_formulation) {
      case ConstantSpeed:
        pi->Deltavar = prob.createVariables(singleton<string>("Delta"), singleton<double>(pi->Delta_lb),singleton<double>(INFINITY))[0];
        break;
      case VariableSpeed:
        //AddVarArray(prob, T, 1, Delta_lb * 0.5, INFINITY, "speed", Deltavars); // TODO should we have a resonable upper bound for this?
        AddVarArray(prob, T, 1, pi->Delta_lb*0.5, INFINITY, "speed", pi->Deltavars); // TODO should we have a resonable upper bound for this?
        //AddVarArray(prob, T, 1, Delta_lb*2, INFINITY, "speed", Deltavars); // TODO should we have a resonable upper bound for this?
        break;
      SWITCH_DEFAULT;
    }
    // Only the twist variables are incremental (i.e. their trust regions should be around zero)
    //prob.setIncremental(twistvars.flatten());
    if (curvature_constraint == BoundedRadius) {
      switch (curvature_formulation) {
        case UseCurvature:
          AddVarArray(prob, T, 1, 0.01, 1. / r_min, "curvature", pi->curvature_or_radius_vars);
          break;
        case UseRadius:
          AddVarArray(prob, T, 1, r_min, 100.0, "radius", pi->curvature_or_radius_vars);
          break;
        SWITCH_DEFAULT;
      }
    }
  }

  void NeedleProblemHelper::InitLocalConfigurations(const KinBodyPtr robot, OptProb& prob, NeedleProblemInstancePtr pi) {
    for (int i = 0; i <= T; ++i) {
      pi->local_configs.push_back(LocalConfigurationPtr(new LocalConfiguration(robot)));
    }
  }

  void NeedleProblemHelper::InitTrajectory(OptProb& prob, NeedleProblemInstancePtr pi) {
    MatrixXd initTraj(T+1, n_dof);
    for (int idof = 0; idof < n_dof; ++idof) {
      initTraj.col(idof) = VectorXd::LinSpaced(T+1, pi->start[idof], pi->goal[idof]);
    }
    for (int i = 0; i <= T; ++i) {
      pi->local_configs[i]->pose = expUp(initTraj.row(i));
    }
  }

  void NeedleProblemHelper::AddStartConstraint(OptProb& prob, NeedleProblemInstancePtr pi) {
    VarVector vars = pi->twistvars.row(0);
    VectorOfVectorPtr f(new Needle::PositionError(pi->local_configs[0], pi->start, shared_from_this()));
    Vector6d coeffs; coeffs << 1., 1., 1., this->coeff_orientation_error, this->coeff_orientation_error, this->coeff_orientation_error;
    prob.addConstraint(ConstraintPtr(new ConstraintFromFunc(f, vars, coeffs, EQ, "entry")));
    pi->dynamics_constraints.push_back(prob.getConstraints().back());
  }

  void NeedleProblemHelper::AddGoalConstraint(OptProb& prob, NeedleProblemInstancePtr pi) {
    VarVector vars = pi->twistvars.row(T);
    VectorOfVectorPtr f(new Needle::PositionError(pi->local_configs[T], pi->goal, shared_from_this()));
    Vector6d coeffs; coeffs << 1., 1., 1., 0., 0., 0.;
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
    for (int i = 0; i < T; ++i) {
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
    for (int i = 0; i < T; ++i) {
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
    Str2Dbl tag2dist_pen(collision_dist_pen), tag2coeff(collision_coeff);
    for (int i = 0; i < ignored_kinbody_names.size(); ++i) {
      tag2coeff.insert( std::pair<string, double>(ignored_kinbody_names[i], 0.0) );
    }
    if (continuous_collision) {
      for (int i = 0; i < T; ++i) {
        prob.addConstraint(ConstraintPtr(new CollisionTaggedConstraint(tag2dist_pen, tag2coeff, pi->local_configs[i], pi->local_configs[i+1], pi->twistvars.row(i), pi->twistvars.row(i+1))));
        pi->collision_constraints.push_back(prob.getConstraints().back());
      }
    } else {
      for (int i = 0; i <= T; ++i) {
        prob.addConstraint(ConstraintPtr(new CollisionTaggedConstraint(tag2dist_pen, tag2coeff, pi->local_configs[i], pi->twistvars.row(i))));
        pi->collision_constraints.push_back(prob.getConstraints().back());
      }
    }
  }

  void NeedleProblemHelper::InitializeCollisionEnvironment() {
    EnvironmentBasePtr env = pis[0]->local_configs[0]->GetEnv();
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
}
