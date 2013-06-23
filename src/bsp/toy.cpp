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
  Matrix2d Q = Matrix2d::Identity();
  Matrix2d QF = Matrix2d::Identity() * 10;
  for (int i = 0; i < T; ++i) {
    prob.addCost(CostPtr(new VarianceCost(sigma_vars.row(i), Q)));
  }
  prob.addCost(CostPtr(new VarianceCost(sigma_vars.row(T), QF)));
}

void ToyBSPProblemHelper::add_control_cost(OptProb& prob) {
  Matrix2d R = Matrix2d::Identity();
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

  for (int i = 0; i < T; ++i) {
    belief_constraints.push_back(BeliefConstraintPtr(new BeliefConstraint<ToyBeliefFunc>(belief_vars.row(i), control_vars.row(i), belief_vars.row(i+1), belief_func)));
    prob.addConstraint(ConstraintPtr(belief_constraints.back()));
  }
}

void ToyBSPProblemHelper::configure_optimizer(OptProb& prob, BasicTrustRegionSQP& opt) {
  
}

ToyStateFunc::ToyStateFunc() {}

ToyStateFunc::ToyStateFunc(ToyBSPProblemHelperPtr helper) : helper(helper) {}

StateT ToyStateFunc::operator()(const StateT& x, const ControlT& u, const StateNoiseT& m) const {
  StateT new_x;
  new_x(0) = x(0) + u(0)*helper->input_dt + sqrt(0.01*u(0)*u(0)*helper->input_dt + 0.001) * m(0);
  new_x(1) = x(1) + u(1)*helper->input_dt + sqrt(0.01*u(1)*u(1)*helper->input_dt + 0.001) * m(1);
  return new_x;
}

ToyObserveFunc::ToyObserveFunc() {}

ToyObserveFunc::ToyObserveFunc(ToyBSPProblemHelperPtr helper) : helper(helper) {}

ObserveT ToyObserveFunc::operator()(const StateT& x, const ObserveNoiseT& n) const {
  return x + 0.1 * n;
}

ToyBeliefFunc::ToyBeliefFunc() {}

ToyBeliefFunc::ToyBeliefFunc(ToyBSPProblemHelperPtr helper, StateFuncPtr f, ObserveFuncPtr h) : helper(helper), f(f), h(h) {}

BeliefT ToyBeliefFunc::operator()(const BeliefT& b, const ControlT& u) const {
  Vector2d x, new_x;
  Vector5d new_b;
  Matrix2d sigma, A, M, gamma, H, N, K;
  Vector2d zero = Vector2d::Zero();

  extract_state(b, &x);
  extract_sigma(b, &sigma);

  // linearize the state propagation function around current point
  f->linearize(x, u, zero, &A, (ToyStateFunc::ControlGradT*) NULL, &M);

  sigma = A*sigma*A.transpose() + M*M.transpose();

  new_x = f->call(x, u, zero);

  gamma = get_gamma(new_x); // TODO make this more generalizable

  // linearize the observation function around current point
  h->linearize(new_x, zero, &H, &N); 

  H = gamma * H;

  PartialPivLU<Matrix2d> solver();
  K = matrix_div((Matrix2d) (sigma*H.transpose()), (Matrix2d) (H*sigma* H.transpose() + N*N.transpose()));

  sigma = sigma - gamma*K*(H*sigma);

  compose_belief(new_x, matrix_sqrt(sigma), &new_b);

  return new_b;
}

bool ToyBeliefFunc::sgndist(const Vector2d& x, Vector2d* dists) const {
  Vector2d p1; p1 << 0, 2;
  Vector2d p2; p2 << 0, 0;
  (*dists)(0) = (x - p1).norm() - 0.5;
  (*dists)(1) = (x - p2).norm() - 0.5;
  return (*dists)(0) < 0 || (*dists)(1) < 0;
}

Matrix2d ToyBeliefFunc::get_gamma(const Vector2d& x) const {
  Vector2d dists;
  sgndist(x, &dists);
  double gamma1 = 1. - (1./(1.+exp(-alpha*(dists(0)+tol))));
  double gamma2 = 1. - (1./(1.+exp(-alpha*(dists(1)+tol))));
  Matrix2d gamma;
  gamma << gamma1, 0,
           0, gamma2;
  return gamma;
}
