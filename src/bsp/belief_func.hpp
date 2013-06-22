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

    const static int state_dim = StateFuncT::state_dim;
    const static int control_dim = StateFuncT::control_dim;
    const static int state_noise_dim = StateFuncT::state_noise_dim;
    const static int observe_dim = ObserveFuncT::observe_dim;
    const static int observe_noise_dim = ObserveNoiseT::observe_noise_dim;
    const static int belief_dim = state_dim < 0 ? state_dim : (state_dim + (state_dim * (state_dim + 1) / 2));

    typedef typename StateFuncT::scalar_type scalar_type;

    typedef Matrix<scalar_type, belief_dim, 1> BeliefT;
    typedef Matrix<scalar_type, belief_dim, belief_dim> BeliefGradT;
    typedef Matrix<scalar_type, control_dim, belief_dim> ControlGradT;
    typedef Matrix<scalar_type, state_dim, state_dim> VarianceT;
    /** end typedefs */

    double epsilon;
    BeliefFunc(BSPProblemHelperPtr helper, StateFuncPtr f, ObserveFuncPtr h);
    virtual int output_size() const = 0;
    virtual BeliefT operator()(const BeliefT& b, const ControlT& u) const = 0;
    void linearize(const BeliefT& b
                 , const ControlT& u
                 , BeliefGradT* output_A
                 , ControlGradT* output_B
                 , BeliefT output_c
                  ) const {
      if (output_A) num_diff((boost::function<StateT (const BeliefT& )>) boost::bind(&BeliefFunc::operator(), this, _1, u), b, this->output_size(), this->epsilon, output_A);
      if (output_B) num_diff((boost::function<StateT (const ControlT& )>) boost::bind(&BeliefFunc::operator(), this, b, _1), u, this->output_size(), this->epsilon, output_B);
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
      output_sqrt_sigma->resize(helper->state_dim, helper->state_dim);
      for (int index = helper->state_dim, i = 0; i < helper->state_dim; ++i) {
        for (int j = i; j < helper->state_dim; ++j) {
          // the upper diagonal entries of sqrt_sigma are stored in row order in the belief vector
          (*output_sqrt_sigma)(i, j) = (*output_sqrt_sigma)(j, i) = belief(index++);
        }
      }
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
      extract_sqrt_sigma(belief, output_sigma);
      (*output_sigma) *= output_sigma->transpose();
    }
  };

  
}
