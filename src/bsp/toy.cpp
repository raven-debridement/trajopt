#include "toy.hpp"

using namespace BSP;
using namespace ToyBSP;

ToyBSPProblemHelper::ToyBSPProblemHelper() {
  input_dt = 1;
  set_state_dim(2);
  set_observe_dim(2);
  set_control_dim(2);
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
  
}

void ToyBSPProblemHelper::configure_optimizer(OptProb& prob, BasicTrustRegionSQP& opt) {

}

Vector2d ToyStateFunc::operator()(const Vector2d& x // state
                                , const Vector2d& u // control
                                , const Vector2d& m // state noise
                                 ) const {
  Vector2d new_x;
  new_x(0) = x(0) + u(0)*helper->input_dt + sqrt(0.01*u(0)*u(0)*helper->input_dt + 0.001) * m(0);
  new_x(1) = x(1) + u(1)*helper->input_dt + sqrt(0.01*u(1)*u(1)*helper->input_dt + 0.001) * m(1);
  return new_x;
}

Vector2d ToyObserveFunc::operator()(const Vector2d& x // state
                                  , const Vector2d& n // observe noise
                                   ) const {
  return x + 0.1 * n;
}



Vector5d ToyBeliefFunc::operator()(const Vector5d& b    // current belief
                                 , const Vector2d& u    // control
                                  ) const {

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
