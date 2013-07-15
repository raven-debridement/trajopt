#pragma once

#include "common.hpp"
#include "traits.hpp"
#include "utils.hpp"
#include "bsp_problem_helper_base.hpp"
#include "problem_state.hpp"
#include "types/state_func.hpp"

namespace BSP {
  template<class _StateT=VectorXd, class _ControlT=VectorXd, class _StateNoiseT=VectorXd>
  class StateFunc : public ProblemState {
  public:
    define_state_func_types();
    typedef boost::shared_ptr< StateFunc<StateT, ControlT, StateNoiseT> > Ptr;

    double epsilon;
    BSPProblemHelperBasePtr helper;

    StateFunc() : epsilon(DefaultEpsilon) {}
    StateFunc(BSPProblemHelperBasePtr helper) : ProblemState(helper), helper(helper), epsilon(DefaultEpsilon) {}

    virtual StateT operator()(const StateT& x, const ControlT& u, const StateNoiseT& m) const = 0;

    virtual BSPProblemHelperBasePtr get_helper() const {
      return helper;
    }

    virtual void linearize(const StateT& x // state
                         , const ControlT& u // control
                         , const StateNoiseT& m // state noise
                         , StateGradT* output_A // df/dx
                         , ControlGradT* output_B // df/du
                         , StateNoiseGradT* output_M // df/dm
                          ) const {
      if (output_A) {
        boost::function<StateT (const StateT& )> f_x;
        f_x = boost::bind(&StateFunc::operator(), this, _1, u, m);
        num_diff(f_x, x, state_dim, this->epsilon, output_A);
       }
      if (output_B) {
        boost::function<StateT (const ControlT& )> f_u;
        f_u = boost::bind(&StateFunc::operator(), this, x, _1, m);
        num_diff(f_u, u, state_dim, this->epsilon, output_B);
      }
      if (output_M) {
        boost::function<StateT (const StateNoiseT& )> f_m;
        f_m = boost::bind(&StateFunc::operator(), this, x, u, _1);
        num_diff(f_m, m, state_dim, this->epsilon, output_M);
      }
    }

    StateT call(const StateT& x, const ControlT& u, const StateNoiseT& m) const {
      return operator()(x, u, m);
    }
  };
}
