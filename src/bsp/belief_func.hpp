#pragma once

#include "common.hpp"
#include "traits.hpp"
#include "utils.hpp"
#include "state_func.hpp"
#include "observe_func.hpp"
#include "bsp_problem_helper_base.hpp"
#include "problem_state.hpp"
#include "types/belief_func.hpp"

namespace BSP {

  template< class _StateFuncT=StateFunc<>, class _ObserveFuncT=ObserveFunc<>, class _BeliefT=VectorXd >
  class BeliefFunc : public ProblemState {
  public:
    define_belief_func_types();
    typedef boost::shared_ptr< BeliefFunc<_StateFuncT, _ObserveFuncT, _BeliefT> > Ptr;

    double epsilon;
    double approx_factor;
    BSPProblemHelperBasePtr helper;
    StateFuncPtr f;
    ObserveFuncPtr h;

    BeliefFunc() : epsilon(DefaultEpsilon), approx_factor(0.5) {}

    BeliefFunc(BSPProblemHelperBasePtr helper, StateFuncPtr f, ObserveFuncPtr h) :
      ProblemState(helper), helper(helper), f(f), h(h), epsilon(DefaultEpsilon), approx_factor(0.5) {}

    void set_approx_factor(double new_approx_factor) {
      approx_factor = new_approx_factor;
    }
    
    void scale_approx_factor(double scale_factor) {
      approx_factor *= scale_factor;
    }

    void get_approx_factor() {
      return approx_factor;
    }

    virtual BSPProblemHelperBasePtr get_helper() const {
      return helper;
    }

    virtual BeliefT operator()(const BeliefT& b, const ControlT& u, ObserveT* z_ptr=nullptr, ObserveT* observation_masks_ptr=nullptr, StateT* state_error_out=nullptr) const = 0;

    BeliefT call(const BeliefT& b, const ControlT& u, ObserveT* z_ptr=nullptr, ObserveT* observation_masks_ptr=nullptr, StateT* state_error_out=nullptr) const {
      return operator()(b, u, z_ptr, observation_masks_ptr, state_error_out);
    }

    virtual void extract_state(const BeliefT& belief, StateT* output_state) const {
      assert (belief.size() == belief_dim);
      assert (output_state != nullptr);
      *output_state = belief.head(state_dim);
    }

    virtual void extract_sqrt_sigma(const BeliefT& belief, VarianceT* output_sqrt_sigma) const {
      assert (belief.size() == belief_dim);
      assert (output_sqrt_sigma != nullptr);
      sqrt_sigma_vec_to_sqrt_sigma(belief.tail(sigma_dof), output_sqrt_sigma, state_dim);
    }

    virtual void compose_belief(const StateT& state, const VarianceT& sqrt_sigma, BeliefT* output_belief) const {
      assert (state.size() == state_dim);
      assert (output_belief != nullptr);
      output_belief->resize(belief_dim);
      output_belief->head(state_dim) = state;
      for (int index = state_dim, i = 0; i < state_dim; ++i) {
        for (int j = i; j < state_dim; ++j) {
          (*output_belief)(index++) = 0.5 * (sqrt_sigma(i, j) + sqrt_sigma(j, i));
        }
      }
    }

    virtual void extract_sigma(const BeliefT& belief, VarianceT* output_sigma) const {
      assert (belief.size() == belief_dim);
      assert (output_sigma != nullptr);
      sqrt_sigma_vec_to_sigma(belief.tail(sigma_dof), output_sigma, state_dim);
    }
  };

}
