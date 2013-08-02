#pragma once

#include "common.hpp"
#include "traits.hpp"
#include "utils.hpp"
#include "bsp_problem_helper_base.hpp"
#include "problem_state.hpp"
#include "types/observe_func.hpp"

namespace BSP {
  template<class _StateT=VectorXd, class _ObserveT=VectorXd, class _ObserveNoiseT=VectorXd>
  class ObserveFunc : public ProblemState {
  public:
    define_observe_func_types();
    typedef boost::shared_ptr< ObserveFunc<StateT, ObserveT, ObserveNoiseT> > Ptr;

    double epsilon;
    BSPProblemHelperBasePtr helper;

    ObserveFunc() : epsilon(DefaultEpsilon) {}
    ObserveFunc(BSPProblemHelperBasePtr helper) : ProblemState(helper), helper(helper), epsilon(DefaultEpsilon) {}

    virtual ObserveT operator()(const StateT& x, const ObserveNoiseT& n) const = 0;

    // Gives the observation in real case
    virtual ObserveT real_observation(const StateT& x, const ObserveNoiseT& n) const {
      return operator()(x, n);
    }

    virtual ObserveT observation_masks(const StateT& x, double approx_factor=-1) const {
      ObserveT masks(observe_dim);
      for (int i = 0; i < observe_dim; ++i) {
        masks(i) = 1;
      }
      return masks;
    }

    virtual ObserveT real_observation_masks(const StateT& x) const {
      return observation_masks(x);
    }

    virtual BSPProblemHelperBasePtr get_helper() const {
      return helper;
    }

    virtual void linearize(const StateT& x
                         , const ObserveNoiseT& n
                         , ObserveStateGradT* output_H
                         , ObserveNoiseGradT* output_N
                          ) const {
      if (output_H) {
        boost::function<ObserveT (const StateT& )> f_x;
        f_x = boost::bind(&ObserveFunc::operator(), this, _1, n);
        num_diff(f_x, x, observe_dim, this->epsilon, output_H);
      }
      if (output_N) {
        boost::function<ObserveT (const ObserveNoiseT& )> f_n;
        f_n = boost::bind(&ObserveFunc::operator(), this, x, _1);
        num_diff(f_n, n, observe_dim, this->epsilon, output_N);
      }
    }

    ObserveT call(const StateT& x, const ObserveNoiseT& n) const {
      return operator()(x, n);
    }

  };
}
