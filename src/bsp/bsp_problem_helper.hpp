#pragma once

#include "common.hpp"
#include "belief_func.hpp"
#include "state_func.hpp"
#include "observe_func.hpp"
#include "bsp_problem_helper_base.hpp"
#include "constraints.hpp"
#include "costs.hpp"
#include "optimizers.hpp"

namespace BSP {

  template< class _BeliefFuncT >
  class BSPProblemHelper : public boost::enable_shared_from_this< BSPProblemHelper<_BeliefFuncT> >, public BSPProblemHelperBase {
  public:
    /** begin typedefs */
    typedef _BeliefFuncT BeliefFuncT;
    typedef typename BeliefConstraint<BeliefFuncT>::Ptr BeliefConstraintPtr;
    typedef Matrix<typename BeliefFuncT::scalar_type, BeliefFuncT::_control_dim, BeliefFuncT::_control_dim> ControlCostT;
    typedef typename BeliefFuncT::StateT StateT;
    typedef typename BeliefFuncT::ControlT ControlT;
    typedef typename BeliefFuncT::StateNoiseT StateNoiseT;
    typedef typename BeliefFuncT::StateGradT StateGradT;
    typedef typename BeliefFuncT::ControlGradT ControlGradT;
    typedef typename BeliefFuncT::StateNoiseGradT StateNoiseGradT;

    typedef typename BeliefFuncT::StateT ObserveStateT;
    typedef typename BeliefFuncT::ObserveT ObserveT;
    typedef typename BeliefFuncT::ObserveNoiseT ObserveNoiseT;
    typedef typename BeliefFuncT::ObserveStateGradT ObserveStateGradT;
    typedef typename BeliefFuncT::ObserveNoiseGradT ObserveNoiseGradT;

    typedef typename BeliefFuncT::BeliefT BeliefT;
    typedef typename BeliefFuncT::BeliefGradT BeliefGradT;
    typedef typename BeliefFuncT::BeliefControlGradT BeliefControlGradT;
    typedef typename BeliefFuncT::VarianceT VarianceT;
    typedef typename BeliefFuncT::StateFuncT StateFuncT;
    typedef typename BeliefFuncT::ObserveFuncT ObserveFuncT;
    typedef typename BeliefFuncT::StateFuncPtr StateFuncPtr;
    typedef typename BeliefFuncT::ObserveFuncPtr ObserveFuncPtr;
    typedef typename BeliefFuncT::Ptr BeliefFuncPtr;
    /** end typedefs */

    double noise_level;
    int state_dim;
    int control_dim;
    int observe_dim;
    int state_noise_dim;
    int observe_noise_dim;
    int belief_dim;
    int sigma_dof;
    int T; // Note: this will make the time frame from 0 to T (inclusive)!
    StateT start;
    StateT goal;
    VarianceT start_sigma;
    deque<ControlT> initial_controls;
    VarianceT Q; // variance cost
    VarianceT QF; // final variance cost
    ControlCostT R; // control cost

    DblVec state_lbs;
    DblVec state_ubs;
    DblVec control_lbs;
    DblVec control_ubs;

    VarArray state_vars;
    VarArray sqrt_sigma_vars;
    VarArray control_vars;
    VarArray belief_vars;

    StateFuncPtr state_func;
    ObserveFuncPtr observe_func;
    BeliefFuncPtr belief_func;

    vector< BeliefConstraintPtr > belief_constraints;

    BSPProblemHelper() {} 
    
    virtual void initialize() {
      state_func.reset(new StateFuncT(this->shared_from_this()));
      observe_func.reset(new ObserveFuncT(this->shared_from_this()));
      belief_func.reset(new BeliefFuncT(this->shared_from_this(), state_func, observe_func));
    }

    virtual int get_state_dim() const { return state_dim; }
    virtual int get_control_dim() const { return control_dim; }
    virtual int get_observe_dim() const { return observe_dim; }
    virtual int get_state_noise_dim() const { return state_noise_dim; }
    virtual int get_observe_noise_dim() const { return observe_noise_dim; }
    virtual int get_belief_dim() const { return belief_dim; }
    virtual int get_sigma_dof() const { return sigma_dof; }
    virtual int get_T() const { return T; }

    void set_state_dim(int new_state_dim) {
      state_dim = new_state_dim;
      state_noise_dim = new_state_dim;
    }

    void set_sigma_dof(int new_sigma_dof) {
      sigma_dof = new_sigma_dof;
      belief_dim = state_dim + sigma_dof;
    }

    void set_observe_dim(int new_observe_dim) {
      observe_dim = new_observe_dim;
      observe_noise_dim = new_observe_dim;
    }

    void set_control_dim(int new_control_dim) {
      control_dim = new_control_dim;
    }

    void set_state_bounds(const DblVec& lbs, const DblVec& ubs) {
      assert (lbs.size() == state_dim);
      assert (ubs.size() == state_dim);
      state_lbs = lbs;
      state_ubs = ubs;
    }

    void set_control_bounds(const DblVec& lbs, const DblVec& ubs) {
      assert (lbs.size() == control_dim);
      assert (ubs.size() == control_dim);
      control_lbs = lbs;
      control_ubs = ubs;
    }

    virtual void configure_problem(OptProb& prob) {
      create_state_variables(prob);
      create_sigma_variables(prob);
      create_control_variables(prob);
      create_belief_variables(prob);
      add_variance_cost(prob);
      add_control_cost(prob);
      add_start_state_constraint(prob);
      add_start_sigma_constraint(prob);
      add_goal_constraint(prob);
      add_belief_constraint(prob);
    }

    virtual void create_state_variables(OptProb& prob) {
      BSP::AddVarArray(prob, T+1, state_dim, "state", state_vars);
      for (int i = 0; i <= T; ++i) {
        for (int j = 0; j < state_dim; ++j) {
          prob.setLowerBounds(vector<double>(1, state_lbs[j]), vector<Var>(1, state_vars.at(i, j)));
          prob.setUpperBounds(vector<double>(1, state_ubs[j]), vector<Var>(1, state_vars.at(i, j)));
        }
      }
    }

    virtual void create_sigma_variables(OptProb& prob) {
      BSP::AddVarArray(prob, T+1, sigma_dof, "sigma", sqrt_sigma_vars);
    }

    virtual void create_control_variables(OptProb& prob) {
      BSP::AddVarArray(prob, T, control_dim, "control", control_vars);
      for (int i = 0; i < T; ++i) {
        for (int j = 0; j < control_dim; ++j) {
          prob.setLowerBounds(vector<double>(1, control_lbs[j]), vector<Var>(1, control_vars.at(i, j)));
          prob.setUpperBounds(vector<double>(1, control_ubs[j]), vector<Var>(1, control_vars.at(i, j)));
        }
      }
    }

    virtual void create_belief_variables(OptProb& prob) {
      belief_vars.resize(T+1, belief_dim);
      for (int i = 0; i <= T; ++i) {
        for (int j = 0; j < state_dim; ++j) {
          belief_vars(i, j) = state_vars(i, j);
        }
      }
      for (int i = 0; i <= T; ++i) {
        for (int j = 0; j < sigma_dof; ++j) {
          belief_vars(i, state_dim+j) = sqrt_sigma_vars(i, j);
        }
      }
    }

    virtual void add_variance_cost(OptProb& prob) {
      for (int i = 0; i < T; ++i) {
        prob.addCost(CostPtr(new VarianceCost(sqrt_sigma_vars.row(i), Q)));
      }
      prob.addCost(CostPtr(new VarianceCost(sqrt_sigma_vars.row(T), QF)));
    }

    virtual void add_control_cost(OptProb& prob) {
      for (int i = 0; i < T; ++i) {
        prob.addCost(CostPtr(new ControlCost(control_vars.row(i), R)));
      }
    }

    virtual void add_start_state_constraint(OptProb& prob) {
      for (int i = 0; i < state_dim; ++i) {
        prob.addLinearConstraint(exprSub(AffExpr(state_vars.at(0, i)), start(i)), EQ);
      }
      
    }

    virtual void add_start_sigma_constraint(OptProb& prob) {
      VarianceT sqrt_start_sigma = matrix_sqrt(start_sigma);
      for (int index = 0, i = 0; i < state_dim; ++i) {
        for (int j = i; j < state_dim; ++j) {
          prob.addLinearConstraint(exprSub(AffExpr(sqrt_sigma_vars.at(0, index++)), sqrt_start_sigma(i, j)), EQ);
        }
      }
    }

    virtual void add_start_constraint(OptProb& prob) {
      add_start_state_constraint(prob);
      add_start_sigma_constraint(prob);
    }

    virtual void add_goal_constraint(OptProb& prob) {
      for (int i = 0; i < state_dim; ++i) {
        prob.addLinearConstraint(exprSub(AffExpr(state_vars.at(T, i)), goal(i)), EQ);
      }
    }

    virtual void add_belief_constraint(OptProb& prob) {
      //state_func.reset(new StateFuncT(this->shared_from_this()));
      //observe_func.reset(new ObserveFuncT(this->shared_from_this()));
      //belief_func.reset(new BeliefFuncT(this->shared_from_this(), state_func, observe_func));

      for (int i = 0; i < T; ++i) {
        belief_constraints.push_back(BeliefConstraintPtr(new BeliefConstraint<BeliefFuncT>(belief_vars.row(i), control_vars.row(i), belief_vars.row(i+1), belief_func)));
        belief_constraints.back()->setName((boost::format("belief_%i")%i).str());
        prob.addConstraint(ConstraintPtr(belief_constraints.back()));
      }
    }

    virtual void init_optimize_variables(OptProb& prob, BSPTrustRegionSQP& opt) {
      DblVec x(prob.getNumVars()); 

      for (int i = 0; i < T; ++i) {
        for (int j = 0; j < control_dim; ++j) {
          x[control_vars.at(i, j).var_rep->index] = initial_controls[i](j);
        }
      }


      assert(belief_func->helper);
      vector<BeliefT> init_beliefs;
      BeliefT cur_belief(belief_dim);
      belief_func->compose_belief(start, matrix_sqrt(start_sigma), &cur_belief);
      init_beliefs.push_back(cur_belief);
      for (int i = 0; i < T; ++i) {
        cur_belief = belief_func->call(cur_belief, initial_controls[i]);
        init_beliefs.push_back(cur_belief);
      }
      for (int i = 0; i <= T; ++i) {
        for (int j = 0; j < belief_dim; ++j) {
          x[belief_vars.at(i, j).var_rep->index] = init_beliefs[i](j);
        }
      }
      opt.initialize(x);
    }

    virtual void merit_done_callback(OptProb*, DblVec& x) {
      //belief_func->scale_approx_factor(2.);

      BeliefT cur_belief;
      belief_func->compose_belief(start, matrix_sqrt(start_sigma), &cur_belief);
      for (int i = 0; i < belief_dim; ++i) {
        x[belief_vars.at(0, i).var_rep->index] = cur_belief(i);
      }
      for (int i = 0; i < T; ++i) {
        cur_belief = belief_func->call(cur_belief, (ControlT) getVec(x, control_vars.row(i)));
        for (int j = 0; j < belief_dim; ++j) {
          x[belief_vars.at(i+1, j).var_rep->index] = cur_belief(j);
        }
      }
    }

    virtual void add_optimizer_callback(OptProb& prob, BSPTrustRegionSQP& opt) {
      opt.addMeritDoneCallback(boost::bind(&BSPProblemHelper::merit_done_callback, this, _1, _2));
    }

    void configure_optimizer(OptProb& prob, BSPTrustRegionSQP& opt) {
      init_optimize_variables(prob, opt);  
      add_optimizer_callback(prob, opt);
    }

    void set_variance_cost(const VarianceT& cost) {
      Q = cost;
    }

    void set_final_variance_cost(const VarianceT& cost) {
      QF = cost;
    }

    void set_control_cost(const ControlCostT& cost) {
      R = cost;
    }

    virtual void add_state_constraint(OptProb& prob) {
      for (int i = 0; i < T; ++i) {
        VarVector vars = concat(concat(state_vars.row(i), control_vars.row(i)), state_vars.row(i+1));
        VectorOfVectorPtr f(new StateError<StateFuncT>(state_func, state_dim, control_dim, state_noise_dim));
        VectorXd coeffs = VectorXd::Ones(state_dim);
        prob.addConstraint(ConstraintPtr(new ConstraintFromFunc(f, vars, coeffs, EQ, (boost::format("state_%i")%i).str())));
      }
    }

    void initialize_controls_in_state_space() {
      OptProbPtr prob(new OptProb());
      create_state_variables(*prob);
      create_control_variables(*prob);
      add_start_state_constraint(*prob);
      add_goal_constraint(*prob);
      add_state_constraint(*prob);
      BasicTrustRegionSQP opt(prob);
      DblVec x(prob->getNumVars(), 0);
      opt.initialize(x);
      opt.optimize();
      cout << "optimized" << endl;
      DblVec result = opt.x();
      initial_controls.clear();
      for (int i = 0; i < get_T(); ++i) {
        ControlT uvec = (ControlT) getVec(result, control_vars.row(i));
        initial_controls.push_back(uvec);
      }
    }
  };
}
