#pragma once

#include "common.hpp"
#include "traits.hpp"
#include "utils.hpp"
#include "bsp_problem_helper_base.hpp"
#include "problem_state.hpp"

namespace BSP {
  template<class _StateT=VectorXd, class _ControlT=VectorXd, class _StateNoiseT=VectorXd>
  class StateFunc : public ProblemState {
  public:
    /** begin typedefs */
    typedef _StateT StateT;
    typedef _ControlT ControlT;
    typedef _StateNoiseT StateNoiseT;
    typedef boost::shared_ptr< StateFunc<StateT, ControlT, StateNoiseT> > Ptr;

    ENSURE_VECTOR(StateT);
    ENSURE_VECTOR(ControlT);
    ENSURE_VECTOR(StateNoiseT);
    ENSURE_SAME_SCALAR_TYPE(StateT, ControlT);
    ENSURE_SAME_SCALAR_TYPE(StateT, StateNoiseT);

    typedef typename MatrixTraits<StateT>::scalar_type scalar_type;

    const static int _state_dim = MatrixTraits<StateT>::rows;
    const static int _control_dim = MatrixTraits<ControlT>::rows;
    const static int _state_noise_dim = MatrixTraits<StateNoiseT>::rows;

    // types for the gradient matrices
    typedef Matrix<scalar_type, _state_dim, _state_dim> StateGradT;
    typedef Matrix<scalar_type, _state_dim, _control_dim> ControlGradT;
    typedef Matrix<scalar_type, _state_dim, _state_noise_dim> StateNoiseGradT;
    /** end typedefs */

    double epsilon;
    BSPProblemHelperBasePtr helper;

    StateFunc() : epsilon(BSP_DEFAULT_EPSILON) {}
    StateFunc(BSPProblemHelperBasePtr helper) : ProblemState(helper), helper(helper), epsilon(BSP_DEFAULT_EPSILON) {}

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
      if (output_A) num_diff((boost::function<StateT (const StateT& )>) boost::bind(&StateFunc::operator(), this, _1, u, m), x, state_dim, this->epsilon, output_A);
      if (output_B) num_diff((boost::function<StateT (const ControlT& )>) boost::bind(&StateFunc::operator(), this, x, _1, m), u, state_dim, this->epsilon, output_B);
      if (output_M) num_diff((boost::function<StateT (const StateNoiseT& )>) boost::bind(&StateFunc::operator(), this, x, u, _1), m, state_dim, this->epsilon, output_M);
    }

    StateT call(const StateT& x, const ControlT& u, const StateNoiseT& m) const {
      return operator()(x, u, m);
    }
  };
}
