#pragma once

#include "common.hpp"
#include "traits.hpp"
#include "bsp_problem_helper.hpp"

namespace BSP {
  template< class StateFuncT=StateFunc<>, class ObserveFuncT=ObserveFunc<> >
  class BeliefFunc {
  public:
    /** begin typedefs */
    typedef typename StateFuncT::StateT StateT;
    typedef typename StateFuncT::ControlT ControlT;
    typedef typename StateFuncT::StateNoiseT StateNoiseT;

    typedef typename ObserveFuncT::StateT ObserveStateT;
    typedef typename ObserveFuncT::ObserveT ObserveT;
    typedef typename ObserveFuncT::ObserveNoiseT ObserveNoiseT;

    typedef typename StateFuncT::Ptr StateFuncPtr;
    typedef typename ObserveFuncT::Ptr ObserveFuncPtr;

    typedef boost::shared_ptr< BeliefFunc<StateFuncT, ObserveFuncT> > Ptr;

    BOOST_STATIC_ASSERT( (boost::is_same< ObserveStateT, StateT >::value) );

    const static int _state_dim = StateFuncT::_state_dim;
    const static int _control_dim = StateFuncT::_control_dim;
    const static int _state_noise_dim = StateFuncT::_state_noise_dim;
    const static int _observe_dim = ObserveFuncT::_observe_dim;
    const static int _observe_noise_dim = ObserveNoiseT::_observe_noise_dim;
    const static int _belief_dim = _state_dim < 0 ? _state_dim : (_state_dim + (_state_dim * (_state_dim + 1) / 2));

    typedef typename StateFuncT::scalar_type scalar_type;

    typedef Matrix<scalar_type, _belief_dim, 1> BeliefT;
    typedef Matrix<scalar_type, _belief_dim, _belief_dim> BeliefGradT;
    typedef Matrix<scalar_type, _control_dim, _belief_dim> ControlGradT;
    typedef Matrix<scalar_type, _state_dim, _state_dim> VarianceT;
    /** end typedefs */

    double epsilon;

    BeliefFunc(BSPProblemHelperPtr helper, StateFuncPtr f, ObserveFuncPtr h);
    virtual BeliefT operator()(const BeliefT& b, const ControlT& u) const = 0;
    void linearize(const BeliefT& b
                 , const ControlT& u
                 , BeliefGradT* output_A
                 , ControlGradT* output_B
                 , BeliefT output_c
                  ) const {
      if (output_A) num_diff((boost::function<StateT (const BeliefT& )>) boost::bind(&BeliefFunc::operator(), this, _1, u), b, helper->belief_dim, this->epsilon, output_A);
      if (output_B) num_diff((boost::function<StateT (const ControlT& )>) boost::bind(&BeliefFunc::operator(), this, b, _1), u, helper->belief_dim, this->epsilon, output_B);
      if (output_c) *output_c = this->call(b, u);
    }

    BeliefT call(const BeliefT& b, const ControlT& u) const {
      return operator()(b, u);
    }
  protected:
    BSPProblemHelperPtr helper;
    StateFuncPtr f;
    ObserveFuncPtr h;

    void extract_state(const BeliefT& belief, StateT* output_state) const {
      assert (belief.size() == helper->belief_dim);
      *output_state = belief.head(helper->state_dim);
    }

    void extract_sqrt_sigma(const BeliefT& belief, VarianceT* output_sqrt_sigma) const {
      assert (belief.size() == helper->belief_dim);
      sqrt_sigma_vec_to_sqrt_sigma(belief.tail(helper->sigma_dof), output_sqrt_sigma);
    }

    void compose_belief(const StateT& state, const VarianceT& sqrt_sigma, BeliefT* output_belief) const {
      assert (state.size() == helper->state_dim);
      assert (output_belief->resize(helper->belief_dim));
      output_belief->head(helper->state_dim) = state;
      for (int index = helper->state_dim, i = 0; i < helper->state_dim; ++i) {
        for (int j = i; j < helper->state_dim; ++j) {
          (*output_belief)(index++) = 0.5 * (sqrt_sigma(i, j) + sqrt_sigma(j, i));
        }
      }
    }

    void extract_sigma(const BeliefT& belief, VarianceT* output_sigma) const {
      assert (belief.size() == helper->belief_dim);
      sqrt_sigma_vec_to_sigma(belief.tail(helper->sigma_dof), output_sigma);
    }
  };

  
}
