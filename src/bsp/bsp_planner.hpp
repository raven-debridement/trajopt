#pragma once

#include "common.hpp"

namespace BSP {
template< class BSPProblemHelperT >
class BSPPlanner {
public:
  /** begin annoying typedefs */
  typedef boost::shared_ptr<BSPProblemHelperT> BSPProblemHelperPtr;
  typedef typename BSPProblemHelperT::StateT StateT;
  typedef typename BSPProblemHelperT::ControlT ControlT;
  typedef typename BSPProblemHelperT::StateNoiseT StateNoiseT;
  typedef typename BSPProblemHelperT::StateGradT StateGradT;
  typedef typename BSPProblemHelperT::ControlGradT ControlGradT;
  typedef typename BSPProblemHelperT::StateNoiseGradT StateNoiseGradT;
  typedef typename BSPProblemHelperT::ObserveStateT ObserveStateT;
  typedef typename BSPProblemHelperT::ObserveT ObserveT;
  typedef typename BSPProblemHelperT::ObserveNoiseT ObserveNoiseT;
  typedef typename BSPProblemHelperT::ObserveStateGradT ObserveStateGradT;
  typedef typename BSPProblemHelperT::ObserveNoiseGradT ObserveNoiseGradT;
  typedef typename BSPProblemHelperT::BeliefT BeliefT;
  typedef typename BSPProblemHelperT::BeliefGradT BeliefGradT;
  typedef typename BSPProblemHelperT::BeliefControlGradT BeliefControlGradT;
  typedef typename BSPProblemHelperT::VarianceT VarianceT;
  typedef typename BSPProblemHelperT::StateFuncT StateFuncT;
  typedef typename BSPProblemHelperT::ObserveFuncT ObserveFuncT;
  typedef typename BSPProblemHelperT::StateFuncPtr StateFuncPtr;
  typedef typename BSPProblemHelperT::ObserveFuncPtr ObserveFuncPtr;
  typedef typename BSPProblemHelperT::BeliefFuncPtr BeliefFuncPtr;

  const static int _state_noise_dim = StateFuncT::_state_noise_dim;
  const static int _observe_noise_dim = ObserveFuncT::_observe_noise_dim;

  typedef typename BSPProblemHelperT::BeliefFuncT::scalar_type scalar_type;

  typedef Matrix<scalar_type, _state_noise_dim, _state_noise_dim> StateNoiseCovT;
  typedef Matrix<scalar_type, _observe_noise_dim, _observe_noise_dim> ObserveNoiseCovT;
  /** end annoying typedefs */

  StateT start;
  VarianceT start_sigma;
  StateT goal;
  int T;
  bool initialized;
  int n_alpha_iterations;
  BSPProblemHelperPtr helper;

  enum Method { StateSpace = 0, ContinuousBeliefSpace = 1, DiscontinuousBeliefSpace = 2};
  int method;
  double noise_level;

  deque<ControlT> controls;
  StateNoiseT state_noise_mean;
  StateNoiseCovT state_noise_cov;
  ObserveNoiseT observe_noise_mean;
  ObserveNoiseCovT observe_noise_cov;
  StateT current_position; // current position initially sampled from the Gaussian distribution,
                           // and later simulated with sampled noise
  vector<StateT> simulated_positions; 
  DblVec result;

  BSPPlanner() : initialized(false) {
    n_alpha_iterations = 3;
  }

  bool finished() {
    return (helper->T <= 0);
  }

  virtual void initialize_optimizer_parameters(BSPTrustRegionSQP& opt) {
    opt.max_iter_                   = 100;
    opt.merit_error_coeff_          = 100;
    opt.merit_coeff_increase_ratio_ = 10;
    opt.max_merit_coeff_increases_  = 2;
    opt.trust_shrink_ratio_         = 0.1;
    opt.trust_expand_ratio_         = 1.5;
    opt.min_trust_box_size_         = 0.001;
    opt.min_approx_improve_         = 1e-2;
    opt.min_approx_improve_frac_    = 1e-4;
    opt.improve_ratio_threshold_    = 0.25;
    opt.trust_box_size_             = 1;
    opt.cnt_tolerance_              = 1e-06;
  }

  virtual void initialize() {
    assert(!initialized);
    helper.reset(new BSPProblemHelperT());
    helper->start = start;
    helper->goal = goal;
    helper->start_sigma = start_sigma;
    helper->T = T;
    helper->initial_controls = controls;
    helper->noise_level = noise_level;
    helper->initialize();
    state_noise_mean = StateNoiseT::Zero(helper->state_noise_dim);
    state_noise_cov = StateNoiseCovT::Identity(helper->state_noise_dim, helper->state_noise_dim);
    observe_noise_mean = ObserveNoiseT::Zero(helper->observe_noise_dim);
    observe_noise_cov = ObserveNoiseCovT::Identity(helper->observe_noise_dim, helper->observe_noise_dim);
    current_position = start;//sample_gaussian(start, start_sigma, noise_level);
    simulated_positions.push_back(current_position);
    initialized = true;
  }

  void solve(boost::function<void(OptProb*, DblVec&)>& opt_callback) {
    assert (initialized);

    if (method == StateSpace) {
      return;
    }

    OptProbPtr prob(new OptProb());
    helper->start = start;
    helper->start_sigma = start_sigma;
    helper->initialize();
    helper->configure_problem(*prob);

    if (method == ContinuousBeliefSpace) {
      n_alpha_iterations = 1;
      helper->belief_func->approx_factor = 1;
    }

    for (int i = 0; i < n_alpha_iterations; ++i) {
      BSPTrustRegionSQP opt(prob);
      initialize_optimizer_parameters(opt);
      helper->configure_optimizer(*prob, opt);
      if (opt_callback) {
        opt.addCallback(opt_callback);
      }
      opt.optimize();
      this->result = opt.x();
      helper->initial_controls.clear();
      for (int j = 0; j < helper->get_T(); ++j) {
        ControlT uvec = (ControlT) getVec(result, helper->control_vars.row(j));
        helper->initial_controls.push_back(uvec);
      }
      helper->belief_func->approx_factor *= 3;
    }

    controls = helper->initial_controls;
    
  }

  virtual void custom_simulation_update(StateT* state, VarianceT* sigma, const StateT& actual_state) {}
  
  void simulate_executions(int nsteps) {
    assert (initialized);
    if (nsteps <= 0 || nsteps > helper->T) {
      return;
    }
    StateFuncPtr state_func = helper->state_func;
    ObserveFuncPtr observe_func = helper->observe_func;
    BeliefFuncPtr belief_func = helper->belief_func;
    for (int i = 0; i < nsteps; ++i) {
      StateNoiseT state_noise = sample_gaussian(state_noise_mean, state_noise_cov, noise_level);
      ObserveNoiseT observe_noise = sample_gaussian(observe_noise_mean, observe_noise_cov, noise_level);
      // update actual position
      current_position = state_func->call(current_position, controls.front(), state_noise);
      simulated_positions.push_back(current_position);
      // update observation
      ObserveT observe = observe_func->real_observation(current_position, observe_noise);
      // update Kalman filter
      BeliefT belief(helper->belief_dim);
      belief_func->compose_belief(start, matrix_sqrt(start_sigma), &belief);
      {
        double current_approx_factor = belief_func->approx_factor;
        belief_func->approx_factor = -1;
        belief = belief_func->call(belief, controls.front(), observe, true);
        belief_func->approx_factor = current_approx_factor;
      }
      belief_func->extract_state(belief, &start);
      belief_func->extract_sigma(belief, &start_sigma);
      custom_simulation_update(&start, &start_sigma, current_position);
      controls.pop_front();
    }
    helper->initial_controls = controls;
    helper->T = controls.size();
    //cout << "Finished simulating execution, remaining horizon: " << helper->T << endl;
  }

  void simulate_execution() {
    simulate_executions(1);
  }
  
};

}
