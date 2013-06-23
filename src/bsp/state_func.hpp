#pragma once

#include "common.hpp"
#include "traits.hpp"
#include "bsp_problem_helper.hpp"

namespace BSP {
  template<class _StateT=VectorXd, class _ControlT=VectorXd, class _StateNoiseT=VectorXd>
  class StateFunc {
  public:
    /** begin typedefs */
    typedef _StateT StateT;
    typedef _ControlT ControlT;
    typedef _StateNoiseT StateNoiseT;
    typedef boost::shared_ptr< StateFunc<StateT, ControlT, StateNoiseT> > Ptr;

    // make sure all three types are indeed vectors
    BOOST_STATIC_ASSERT( (MatrixTraits<StateT>::cols == 1) );
    BOOST_STATIC_ASSERT( (MatrixTraits<ControlT>::cols == 1) );
    BOOST_STATIC_ASSERT( (MatrixTraits<StateNoiseT>::cols == 1) );

    typedef typename MatrixTraits<StateT>::scalar_type state_scalar_type;
    typedef typename MatrixTraits<ControlT>::scalar_type control_scalar_type;
    typedef typename MatrixTraits<StateNoiseT>::scalar_type state_noise_scalar_type;

    // make sure all three vectors are of the same scalar type
    BOOST_STATIC_ASSERT( (boost::is_same< state_scalar_type, control_scalar_type>::value) );
    BOOST_STATIC_ASSERT( (boost::is_same< state_scalar_type, state_noise_scalar_type>::value) );

    typedef state_scalar_type scalar_type;
    const static int _state_dim = MatrixTraits<StateT>::rows;
    const static int _control_dim = MatrixTraits<ControlT>::rows;
    const static int _state_noise_dim = MatrixTraits<StateNoiseT>::rows;

    // types for the gradient matrices
    typedef Matrix<scalar_type, _state_dim, _state_dim> StateGradT;
    typedef Matrix<scalar_type, _control_dim, _state_dim> ControlGradT;
    typedef Matrix<scalar_type, _state_noise_dim, _state_dim> StateNoiseGradT;
    /** end typedefs */

    double epsilon;

    StateFunc(BSPProblemHelperPtr helper);
    virtual StateT operator()(const StateT& x, const ControlT& u, const StateNoiseT& m) const = 0;

    virtual void linearize(const StateT& x // state
                         , const ControlT& u // control
                         , const StateNoiseT& m // state noise
                         , StateGradT* output_A // df/dx
                         , ControlGradT* output_B // df/du
                         , StateNoiseGradT* output_M // df/dm
                          ) const {
      if (output_A) num_diff((boost::function<StateT (const StateT& )>) boost::bind(&StateFunc::operator(), this, _1, u, m), x, helper->state_dim, this->epsilon, output_A);
      if (output_B) num_diff((boost::function<StateT (const ControlT& )>) boost::bind(&StateFunc::operator(), this, x, _1, m), u, helper->state_dim, this->epsilon, output_B);
      if (output_M) num_diff((boost::function<StateT (const StateNoiseT& )>) boost::bind(&StateFunc::operator(), this, x, u, _1), m, helper->state_dim, this->epsilon, output_M);
    }

    StateT call(const StateT& x, const ControlT& u, const StateNoiseT& m) const {
      return operator()(x, u, m);
    }
  protected:
    BSPProblemHelperPtr helper;
  };
}
