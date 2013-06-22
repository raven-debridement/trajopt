#pragma once

#include "common.hpp"
#include "traits.hpp"
#include "bsp_problem_helper.hpp"

namespace BSP {
  template<class _StateT=VectorXd, class _ObserveT=VectorXd, class _ObserveNoiseT=VectorXd>
  class ObserveFunc {
  public:
    /** begin typedefs */
    typedef _StateT StateT;
    typedef _ObserveT ObserveT;
    typedef _ObserveNoiseT ObserveNoiseT;

    typedef boost::shared_ptr< ObserveFunc<StateT, ObserveT, ObserveNoiseT> > Ptr;

    BOOST_STATIC_ASSERT( (MatrixTraits<StateT>::cols == 1) );
    BOOST_STATIC_ASSERT( (MatrixTraits<ObserveT>::cols == 1) );
    BOOST_STATIC_ASSERT( (MatrixTraits<ObserveNoiseT>::cols == 1) );

    typedef typename MatrixTraits<StateT>::scalar_type state_scalar_type;
    typedef typename MatrixTraits<ObserveT>::scalar_type observe_scalar_type;
    typedef typename MatrixTraits<ObserveNoiseT>::scalar_type observe_noise_scalar_type;

    BOOST_STATIC_ASSERT( (boost::is_same< state_scalar_type, observe_scalar_type>::value) );
    BOOST_STATIC_ASSERT( (boost::is_same< state_scalar_type, observe_noise_scalar_type>::value) );

    typedef state_scalar_type scalar_type;

    ObserveFunc(BSPProblemHelperPtr helper);
    const static int state_dim = MatrixTraits<StateT>::rows;
    const static int observe_dim = MatrixTraits<ObserveT>::rows;
    const static int observe_noise_dim = MatrixTraits<ObserveNoiseT>::rows;

    typedef Matrix<scalar_type, state_dim, state_dim> StateGradT;
    typedef Matrix<scalar_type, observe_noise_dim, state_dim> ObserveNoiseGradT;
    /** end typedefs */

    virtual int output_size() const = 0;

    virtual ObserveT operator()(const StateT& x, const ObserveNoiseT& n) const = 0;

    virtual void linearize(const StateT& x
                         , const ObserveT& n
                         , StateGradT* output_H
                         , ObserveNoiseGradT* output_N
                          ) const {
      if (output_H) num_diff((boost::function<StateT (const StateT& )>) boost::bind(&ObserveFunc::operator(), this, _1, n), x, this->output_size(), this->epsilon, output_H);
      if (output_N) num_diff((boost::function<StateT (const ObserveT& )>) boost::bind(&ObserveFunc::operator(), this, x, _1), n, this->output_size(), this->epsilon, output_N);
    }

    ObserveT call(const StateT& x, const ObserveNoiseT& n) const {
      return operator()(x, n);
    }

    double epsilon;
  protected:
    BSPProblemHelperPtr helper;
  };
}
