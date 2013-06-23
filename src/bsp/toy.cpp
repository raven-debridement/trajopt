#include "toy.hpp"

using namespace BSP;
using namespace ToyBSP;

ToyBSPProblemHelper::ToyBSPProblemHelper() {
  input_dt = 1;
  set_state_dim(2);
  set_observe_dim(2);
  set_control_dim(2);
  double x_min = -10, x_max = 10, u_min = -0.9, u_max = 0.9;
  double goal_eps = 0.1;
}

void ToyBSPProblemHelper::set_state_dim(int new_state_dim) {
  state_dim = new_state_dim;
  state_noise_dim = new_state_dim;
  sigma_dof = (new_state_dim * (new_state_dim + 1)) / 2;
  belief_dim = state_dim + sigma_dof;
}

void ToyBSPProblemHelper::set_observe_dim(int new_observe_dim) {
  observe_dim = new_observe_dim;
  observe_noise_dim = new_observe_dim;
}

void ToyBSPProblemHelper::set_control_dim(int new_control_dim) {
  control_dim = new_control_dim;
}

void ToyBSPProblemHelper::configure_problem(OptProb& prob) {
  create_variables(prob);
  add_variance_cost(prob);
  add_control_cost(prob);
  add_start_constraint(prob);
  add_goal_constraint(prob);
  add_belief_constraint(prob);
}

void ToyBSPProblemHelper::create_variables(OptProb& prob) {
  AddVarArray(prob, T+1, state_dim, "state", state_vars);
  AddVarArray(prob, T+1, sigma_dof, "sigma", sigma_vars);
  AddVarArray(prob, T, control_dim, "control", control_vars);
  belief_vars.resize(T+1, belief_dim);
  for (int i = 0; i <= T; ++i) {
    for (int j = 0; j < state_dim; ++j) {
      belief_vars(i, j) = state_vars(i, j);
    }
  }
  for (int i = 0; i <= T; ++i) {
    for (int j = 0; j < sigma_dof; ++j) {
      belief_vars(i, state_dim+j) = sigma_vars(i, j);
    }
  }
}

void ToyBSPProblemHelper::add_variance_cost(OptProb& prob) {
  VarianceCostT Q = VarianceCostT::Identity(state_dim, state_dim);
  VarianceCostT QF = VarianceCostT::Identity(state_dim, state_dim) * 10;
  for (int i = 0; i < T; ++i) {
    prob.addCost(CostPtr(new VarianceCost(sigma_vars.row(i), Q)));
  }
  prob.addCost(CostPtr(new VarianceCost(sigma_vars.row(T), QF)));
}

void ToyBSPProblemHelper::add_control_cost(OptProb& prob) {
  ControlCostT R = ControlCostT::Identity(control_dim, control_dim);
  for (int i = 0; i < T; ++i) {
    prob.addCost(CostPtr(new ControlCost(control_vars.row(i), R)));
  }
}

void ToyBSPProblemHelper::add_start_constraint(OptProb& prob) {
  for (int i = 0; i < n_dof; ++i) {
    prob.addLinearConstraint(exprSub(AffExpr(state_vars.at(0, i)), start(i)), EQ);
  }
}

void ToyBSPProblemHelper::add_goal_constraint(OptProb& prob) {
  for (int i = 0; i < n_dof; ++i) {
    prob.addLinearConstraint(exprSub(AffExpr(state_vars.at(T, i)), goal(i)), EQ);
  }
}

void ToyBSPProblemHelper::add_belief_constraint(OptProb& prob) {
  state_func.reset(new ToyStateFunc(shared_from_this()));
  observe_func.reset(new ToyObserveFunc(shared_from_this()));
  belief_func.reset(new ToyBeliefFunc(shared_from_this(), state_func, observe_func));
  belief_func->set_alpha(0.5);
  belief_func->set_tol(0.1);

  for (int i = 0; i < T; ++i) {
    belief_constraints.push_back(BeliefConstraintPtr(new BeliefConstraint<ToyBeliefFunc>(belief_vars.row(i), control_vars.row(i), belief_vars.row(i+1), belief_func)));
    belief_constraints.back()->setName((boost::format("belief_%i")%i).str());
    prob.addConstraint(ConstraintPtr(belief_constraints.back()));
  }
}

void ToyBSPProblemHelper::init_optimize_variables(OptProb& prob, BSPTrustRegionSQP& opt) {
  DblVec x(prob.getNumVars()); 

  vector<ControlT> init_controls;
  ControlT control_step;
  control_step.resize(control_dim);
  control_step(0) = (goal(0) - start(0)) / T;
  control_step(1) = (goal(1) - start(1)) / T;
  for (int i = 0; i < T; ++i) {
    init_controls.push_back(control_step);
  }
  for (int i = 0; i < T; ++i) {
    for (int j = 0; j < control_dim; ++j) {
      x[control_vars.at(i, j).var_rep->index] = init_controls[i](j);
    }
  }

  assert(belief_func->helper);
  vector<BeliefT> init_beliefs;
  BeliefT cur_belief(belief_dim);
  belief_func->compose_belief(start, VarianceT::Identity(state_dim, state_dim), &cur_belief);
  init_beliefs.push_back(cur_belief);
  for (int i = 0; i < T; ++i) {
    cur_belief = belief_func->call(cur_belief, init_controls[i]);
    init_beliefs.push_back(cur_belief);
  }
  for (int i = 0; i <= T; ++i) {
    for (int j = 0; j < belief_dim; ++j) {
      x[belief_vars.at(i, j).var_rep->index] = init_beliefs[i](j);
    }
  }
  opt.initialize(x);
}

void ToyBSPProblemHelper::merit_done_callback(OptProb*, DblVec& x) {
  belief_func->scale_alpha(4.);
}

void ToyBSPProblemHelper::add_optimizer_callback(OptProb& prob, BSPTrustRegionSQP& opt) {
  opt.addMeritDoneCallback(boost::bind(&ToyBSPProblemHelper::merit_done_callback, this, _1, _2));
}

void ToyBSPProblemHelper::configure_optimizer(OptProb& prob, BSPTrustRegionSQP& opt) {
  init_optimize_variables(prob, opt);  
  add_optimizer_callback(prob, opt);
}

ToyStateFunc::ToyStateFunc() : StateFunc<StateT, ControlT, StateNoiseT>() {}

ToyStateFunc::ToyStateFunc(ToyBSPProblemHelperPtr helper) : StateFunc<StateT, ControlT, StateNoiseT>(helper), toy_helper(helper) {}

StateT ToyStateFunc::operator()(const StateT& x, const ControlT& u, const StateNoiseT& m) const {
  StateT new_x(helper->state_dim);
  new_x(0) = x(0) + u(0)*toy_helper->input_dt + sqrt(0.01*u(0)*u(0)*toy_helper->input_dt + 0.001) * m(0);
  new_x(1) = x(1) + u(1)*toy_helper->input_dt + sqrt(0.01*u(1)*u(1)*toy_helper->input_dt + 0.001) * m(1);
  return new_x;
}

ToyObserveFunc::ToyObserveFunc() : ObserveFunc<StateT, ObserveT, ObserveNoiseT>() {}

ToyObserveFunc::ToyObserveFunc(ToyBSPProblemHelperPtr helper) : ObserveFunc<StateT, ObserveT, ObserveNoiseT>(helper), toy_helper(helper) {}

ObserveT ToyObserveFunc::operator()(const StateT& x, const ObserveNoiseT& n) const {
  return x + 0.1 * n;
}

ToyBeliefFunc::ToyBeliefFunc() : BeliefFunc<ToyStateFunc, ToyObserveFunc>() {}

ToyBeliefFunc::ToyBeliefFunc(ToyBSPProblemHelperPtr helper, StateFuncPtr f, ObserveFuncPtr h) : 
  BeliefFunc<ToyStateFunc, ToyObserveFunc>(helper, f, h), toy_helper(helper) {}

BeliefT ToyBeliefFunc::operator()(const BeliefT& b, const ControlT& u) const {
  StateT x(helper->state_dim), new_x(helper->state_dim);
  BeliefT new_b(helper->belief_dim);
  VarianceT sigma(helper->state_dim, helper->state_dim);
  StateGradT A(helper->state_dim, helper->state_dim);
  StateNoiseGradT M(helper->state_dim, helper->state_noise_dim);
  GammaT gamma(helper->state_dim, helper->state_dim);
  ObserveStateGradT H(helper->observe_dim, helper->state_dim);
  ObserveNoiseGradT N(helper->observe_dim, helper->observe_noise_dim);
  KalmanT K(helper->state_dim, helper->state_dim);

  StateNoiseT zero_state_noise = StateNoiseT::Zero(helper->state_noise_dim);
  ObserveNoiseT zero_observe_noise = ObserveNoiseT::Zero(helper->observe_noise_dim);

  extract_state(b, &x);
  extract_sigma(b, &sigma);

  // linearize the state propagation function around current point
  f->linearize(x, u, zero_state_noise, &A, (ToyStateFunc::ControlGradT*) NULL, &M);

  sigma = A*sigma*A.transpose() + M*M.transpose();

  new_x = f->call(x, u, zero_state_noise);

  gamma = get_gamma(new_x); // TODO make this more generalizable

  // linearize the observation function around current point
  h->linearize(new_x, zero_observe_noise, &H, &N); 

  H = gamma * H;

  PartialPivLU<KalmanT> solver();
  K = matrix_div((KalmanT) (sigma*H.transpose()), (KalmanT) (H*sigma* H.transpose() + N*N.transpose()));

  sigma = sigma - gamma*K*(H*sigma);

  compose_belief(new_x, matrix_sqrt(sigma), &new_b);

  return new_b;
}

bool ToyBeliefFunc::sgndist(const StateT& x, StateT* dists) const {
  StateT p1(helper->state_dim); p1 << 0, 2;
  StateT p2(helper->state_dim); p2 << 0, 0;
  (*dists)(0) = (x - p1).norm() - 0.5;
  (*dists)(1) = (x - p2).norm() - 0.5;
  return (*dists)(0) < 0 || (*dists)(1) < 0;
}

GammaT ToyBeliefFunc::get_gamma(const StateT& x) const {
  StateT dists(helper->state_dim);
  sgndist(x, &dists);
  double gamma1 = 1. - (1./(1.+exp(-alpha*(dists(0)+tol))));
  double gamma2 = 1. - (1./(1.+exp(-alpha*(dists(1)+tol))));
  GammaT gamma(helper->state_dim, helper->state_dim);
  gamma << gamma1, 0,
           0, gamma2;
  return gamma;
}

void ToyBeliefFunc::set_alpha(double new_alpha) {
  alpha = new_alpha;
}

void ToyBeliefFunc::scale_alpha(double scale) {
  alpha *= scale;
}

void ToyBeliefFunc::set_tol(double new_tol) {
  tol = new_tol;
}
