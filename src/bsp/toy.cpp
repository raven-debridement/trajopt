#include "toy.hpp"

using namespace BSP;

namespace ToyBSP {

  ToyBSPProblemHelper::ToyBSPProblemHelper() : BSPProblemHelper<ToyBeliefFunc>() {
    input_dt = 1.0;
    set_state_dim(2);
    set_sigma_dof(3);
    set_observe_dim(2);
    set_control_dim(2);
    set_state_bounds(DblVec(2, -10), DblVec(2, 10));
    set_control_bounds(DblVec(2, -0.9), DblVec(2, 0.9));
    set_variance_cost(VarianceT::Identity(state_dim, state_dim));
    set_final_variance_cost(VarianceT::Identity(state_dim, state_dim) * 10);
    set_control_cost(ControlCostT::Identity(control_dim, control_dim));
  }

  void ToyBSPProblemHelper::init_control_values(vector<ControlT>* output_init_controls) const {
    assert (output_init_controls != NULL);
    ControlT control_step;
    control_step.resize(control_dim);
    control_step(0) = (goal(0) - start(0)) / T;
    control_step(1) = (goal(1) - start(1)) / T;
    output_init_controls->clear();
    for (int i = 0; i < T; ++i) {
      output_init_controls->push_back(control_step);
    }
  }

  ToyStateFunc::ToyStateFunc() : StateFunc<StateT, ControlT, StateNoiseT>() {}

  ToyStateFunc::ToyStateFunc(BSPProblemHelperBasePtr helper) :
    StateFunc<StateT, ControlT, StateNoiseT>(helper), toy_helper(boost::static_pointer_cast<ToyBSPProblemHelper>(helper)) {}

  StateT ToyStateFunc::operator()(const StateT& x, const ControlT& u, const StateNoiseT& m) const {
    StateT new_x(state_dim);
    new_x(0) = x(0) + u(0)*toy_helper->input_dt + sqrt(0.01*u(0)*u(0)*toy_helper->input_dt + 0.001) * m(0);
    new_x(1) = x(1) + u(1)*toy_helper->input_dt + sqrt(0.01*u(1)*u(1)*toy_helper->input_dt + 0.001) * m(1);
    return new_x;
  }

  ToyObserveFunc::ToyObserveFunc() : ObserveFunc<StateT, ObserveT, ObserveNoiseT>() {}

  ToyObserveFunc::ToyObserveFunc(BSPProblemHelperBasePtr helper) :
    ObserveFunc<StateT, ObserveT, ObserveNoiseT>(helper), toy_helper(boost::static_pointer_cast<ToyBSPProblemHelper>(helper)) {}

  ObserveT ToyObserveFunc::operator()(const StateT& x, const ObserveNoiseT& n) const {
    return x + 0.1 * n;
  }

  ToyBeliefFunc::ToyBeliefFunc() : BeliefFunc<ToyStateFunc, ToyObserveFunc, BeliefT>(), tol(0.1), alpha(0.5) {}

  ToyBeliefFunc::ToyBeliefFunc(BSPProblemHelperBasePtr helper, StateFuncPtr f, ObserveFuncPtr h) :
    tol(0.1), alpha(0.5), BeliefFunc<ToyStateFunc, ToyObserveFunc, BeliefT>(helper, f, h), toy_helper(boost::static_pointer_cast<ToyBSPProblemHelper>(helper)) {}

  bool ToyBeliefFunc::sgndist(const StateT& x, StateT* dists) const {
    StateT p1(state_dim); p1 << 0, 2;
    StateT p2(state_dim); p2 << 0, 0;
    (*dists)(0) = (x - p1).norm() - 0.5;
    (*dists)(1) = (x - p2).norm() - 0.5;
    return (*dists)(0) < 0 || (*dists)(1) < 0;
  }

  VarianceT ToyBeliefFunc::compute_gamma(const StateT& x) const {
    StateT dists(state_dim);
    sgndist(x, &dists);
    double gamma1 = 1. - (1./(1.+exp(-alpha*(dists(0)+tol))));
    double gamma2 = 1. - (1./(1.+exp(-alpha*(dists(1)+tol))));
    VarianceT gamma(state_dim, state_dim);
    gamma << gamma1, 0,
             0, gamma2;
    return gamma;
  }
}
