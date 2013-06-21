#include "needle_steering.hpp"

namespace Needle {

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
    AddControlConstraint(prob);
    AddCollisionConstraint(prob);
  }

  void NeedleProblemHelper::InitOptimizeVariables(BasicTrustRegionSQP& opt) {
    DblVec initVec;
    // Initialize twistvars
    for (int i = 0; i <= T; ++i) {
      for (int j = 0; j < n_dof; ++j) {
        initVec.push_back(0.);
      }
    }
    // Initialize phivars
    for (int i = 0; i < T; ++i) {
      initVec.push_back(0.);
    }
    switch (speed_formulation) {
      case ConstantSpeed:
        // Initialize Delta
        initVec.push_back(Delta_lb);
        break;
      case VariableSpeed:
        for (int i = 0; i < T; ++i) {
          initVec.push_back(Delta_lb);
        }
        break;
      SWITCH_DEFAULT;
    }
    // Initialize time frame radii
    if (curvature_constraint == BoundedRadius) {
      for (int i = 0; i < T; ++i) {
        switch (curvature_formulation) {
          case UseCurvature:
            initVec.push_back(1.0 / r_min);
            break;
          case UseRadius:
            initVec.push_back(r_min);
            break;
          SWITCH_DEFAULT;
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

  void NeedleProblemHelper::ConfigureOptimizer(BasicTrustRegionSQP& opt) {
    InitOptimizeVariables(opt);
    opt.addCallback(boost::bind(&Needle::NeedleProblemHelper::OptimizerCallback, this, _1, _2));
  }

  void NeedleProblemHelper::CreateVariables(OptProb& prob) {
    // Time frame varies from 0 to T instead of from 0 to T-1
    AddVarArray(prob, T+1, n_dof, "twist", twistvars);
    AddVarArray(prob, T, 1, -PI, PI, "phi", phivars);
    Delta_lb = (goal.topRows(3) - start.topRows(3)).norm() / T / r_min;
    cout << Delta_lb << endl;
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
    prob.setIncremental(twistvars.flatten());
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
    VectorXd coeffs = Vector6d::Ones();
    prob.addConstraint(ConstraintPtr(new ConstraintFromFunc(f, vars, coeffs, EQ, "entry")));
  }

  void NeedleProblemHelper::AddGoalConstraint(OptProb& prob) {
    VarVector vars = twistvars.row(T);
    VectorOfVectorPtr f(new Needle::PositionError(local_configs[T], goal, shared_from_this()));
    Vector6d coeffs; coeffs << 1., 1., 1., 0., 0., 0.;
    prob.addConstraint(ConstraintPtr(new ConstraintFromFunc(f, vars, coeffs, EQ, "goal")));
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
      prob.addConstraint(ConstraintPtr(new ConstraintFromFunc(f, vars, coeffs, EQ, (boost::format("control%i")%i).str())));
    }
  }

  void NeedleProblemHelper::AddCollisionConstraint(OptProb& prob) {
    //Str2Dbl tag2dist_pen(collision_dist_pen), tag2coeff(collision_coeff);
    //for (int i = 0; i < ignored_kinbody_names.size(); ++i) {
    //  tag2coeff.insert( std::pair<string, double>(ignored_kinbody_names[i], 0.0) );
    //}
    //for (int i = 0; i <= T; ++i) {
    //  prob.addConstraint(ConstraintPtr(new CollisionTaggedConstraint(tag2dist_pen, tag2coeff, local_configs[i], twistvars.row(i))));
    //}
  }

}
