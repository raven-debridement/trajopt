#pragma once

#include "common.hpp"
#include "traits.hpp"
#include "utils.hpp"
#include "bsp_problem_helper_base.hpp"
#include "problem_state.hpp"

namespace BSP {
  template<class _StateT=VectorXd, class _ObserveT=VectorXd, class _ObserveNoiseT=VectorXd>
  class ObserveFunc : public ProblemState {
  public:
    /** begin typedefs */
    typedef _StateT StateT;
    typedef _ObserveT ObserveT;
    typedef _ObserveNoiseT ObserveNoiseT;

    typedef boost::shared_ptr< ObserveFunc<StateT, ObserveT, ObserveNoiseT> > Ptr;

    ENSURE_VECTOR(StateT);
    ENSURE_VECTOR(ObserveT);
    ENSURE_VECTOR(ObserveNoiseT);
    ENSURE_SAME_SCALAR_TYPE(StateT, ObserveT);
    ENSURE_SAME_SCALAR_TYPE(StateT, ObserveNoiseT);

    typedef typename MatrixTraits<StateT>::scalar_type scalar_type;

    const static int _state_dim = MatrixTraits<StateT>::rows;
    const static int _observe_dim = MatrixTraits<ObserveT>::rows;
    const static int _observe_noise_dim = MatrixTraits<ObserveNoiseT>::rows;

    typedef Matrix<scalar_type, _observe_dim, _state_dim> ObserveStateGradT;
    typedef Matrix<scalar_type, _observe_dim, _observe_noise_dim> ObserveNoiseGradT;
    /** end typedefs */

    double epsilon;
    BSPProblemHelperBasePtr helper;

    ObserveFunc() : epsilon(DefaultEpsilon) {}
    ObserveFunc(BSPProblemHelperBasePtr helper) : ProblemState(helper), helper(helper), epsilon(DefaultEpsilon) {}

    // This should give the exact measurement
    virtual ObserveT operator()(const StateT& x, const ObserveNoiseT& n) const = 0;

    // This should give the smooth approximation of the measurement.
    // the larger the approx_factor (> 0), the better the smooth approximation should be
    virtual ObserveT operator()(const StateT& x, const ObserveNoiseT& n, double approx_factor) const = 0;

    // Gives the observation in real case
    virtual ObserveT real_observation(const StateT& x, const ObserveNoiseT& n) const {
      return operator()(x, n);
    }

    virtual BSPProblemHelperBasePtr get_helper() const {
      return helper;
    }

    virtual void linearize(const StateT& x
                         , const ObserveNoiseT& n
                         , ObserveStateGradT* output_H
                         , ObserveNoiseGradT* output_N
                          ) const {
      if (output_H) num_diff((boost::function<ObserveT (const StateT& )>) boost::bind(&ObserveFunc::operator(), this, _1, n), x, observe_dim, this->epsilon, output_H);
      if (output_N) num_diff((boost::function<ObserveT (const ObserveNoiseT& )>) boost::bind(&ObserveFunc::operator(), this, x, _1), n, observe_dim, this->epsilon, output_N);
    }

    virtual void linearize(const StateT& x
                         , const ObserveNoiseT& n
                         , ObserveStateGradT* output_H
                         , ObserveNoiseGradT* output_N
                         , double approx_factor
                          ) const {
      if (output_H) num_diff((boost::function<ObserveT (const StateT& )>) boost::bind(&ObserveFunc::operator(), this, _1, n, approx_factor), x, observe_dim, this->epsilon, output_H);
      if (output_N) num_diff((boost::function<ObserveT (const ObserveNoiseT& )>) boost::bind(&ObserveFunc::operator(), this, x, _1, approx_factor), n, observe_dim, this->epsilon, output_N);
    }

    ObserveT call(const StateT& x, const ObserveNoiseT& n) const {
      return operator()(x, n);
    }

    ObserveT call(const StateT& x, const ObserveNoiseT& n, double approx_factor) const {
      return operator()(x, n, approx_factor);
    }
  };
}
