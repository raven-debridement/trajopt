#pragma once

#include "common.hpp"
#include "traits.hpp"
#include "utils.hpp"
#include "state_func.hpp"
#include "observe_func.hpp"
#include "bsp_problem_helper_base.hpp"
#include "problem_state.hpp"

namespace BSP {

  template< class _StateFuncT=StateFunc<>, class _ObserveFuncT=ObserveFunc<>, class _BeliefT=VectorXd >
  class BeliefFunc : public ProblemState {
  public:
    /** begin typedefs */
    typedef _StateFuncT StateFuncT;
    typedef _ObserveFuncT ObserveFuncT;
    typedef _BeliefT BeliefT;
    typedef typename StateFuncT::StateT StateT;
    typedef typename StateFuncT::ControlT ControlT;
    typedef typename StateFuncT::StateNoiseT StateNoiseT;
    typedef typename StateFuncT::StateGradT StateGradT;
    typedef typename StateFuncT::ControlGradT ControlGradT;
    typedef typename StateFuncT::StateNoiseGradT StateNoiseGradT;

    typedef typename ObserveFuncT::StateT ObserveStateT;
    typedef typename ObserveFuncT::ObserveT ObserveT;
    typedef typename ObserveFuncT::ObserveNoiseT ObserveNoiseT;
    typedef typename ObserveFuncT::ObserveStateGradT ObserveStateGradT;
    typedef typename ObserveFuncT::ObserveNoiseGradT ObserveNoiseGradT;

    typedef typename StateFuncT::Ptr StateFuncPtr;
    typedef typename ObserveFuncT::Ptr ObserveFuncPtr;

    typedef boost::shared_ptr< BeliefFunc<StateFuncT, ObserveFuncT, BeliefT> > Ptr;

    BOOST_STATIC_ASSERT( (boost::is_same< ObserveStateT, StateT >::value) );

    ENSURE_VECTOR(BeliefT);
    ENSURE_SAME_SCALAR_TYPE(StateT, ObserveT);
    ENSURE_SAME_SCALAR_TYPE(StateT, BeliefT);

    const static int _state_dim = StateFuncT::_state_dim;
    const static int _control_dim = StateFuncT::_control_dim;
    const static int _state_noise_dim = StateFuncT::_state_noise_dim;
    const static int _observe_dim = ObserveFuncT::_observe_dim;
    const static int _observe_noise_dim = ObserveFuncT::_observe_noise_dim;
    const static int _belief_dim = MatrixTraits<BeliefT>::rows;
    
    BOOST_STATIC_ASSERT( (_belief_dim == -1 || _belief_dim >= _state_dim) );

    typedef typename StateFuncT::scalar_type scalar_type;

    typedef Matrix<scalar_type, _belief_dim, _belief_dim> BeliefGradT;
    typedef Matrix<scalar_type, _belief_dim, _control_dim> BeliefControlGradT;
    typedef Matrix<scalar_type, _state_dim, _state_dim> VarianceT;
    typedef Matrix<scalar_type, _observe_dim, _observe_dim> ObserveMatT;
    typedef Matrix<scalar_type, _state_dim, _observe_dim> KalmanT;
    /** end typedefs */

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

    virtual BeliefT operator()(const BeliefT& b, const ControlT& u, ObserveT* z_ptr=NULL, ObserveT* observation_masks_ptr=NULL) const {
      StateT             x(state_dim), new_x(state_dim);
      BeliefT            new_b(belief_dim);
      VarianceT          sigma(state_dim, state_dim);
      StateGradT         A(state_dim, state_dim);
      StateNoiseGradT    M(state_dim, state_noise_dim);
      ObserveStateGradT  H(observe_dim, state_dim);
      ObserveNoiseGradT  N(observe_dim, observe_noise_dim);
      KalmanT            K(state_dim, observe_dim);

      StateNoiseT        zero_state_noise = StateNoiseT::Zero(state_noise_dim);
      ObserveNoiseT      zero_observe_noise = ObserveNoiseT::Zero(observe_noise_dim);

      extract_state(b, &x);
      extract_sigma(b, &sigma);

      f->linearize(x, u, zero_state_noise, &A, (ControlGradT*) NULL, &M);

      sigma = A*sigma*A.transpose() + M*M.transpose();

      new_x = f->call(x, u, zero_state_noise);

      ObserveT observation_masks;
      if (observation_masks_ptr != NULL) {
        observation_masks = *observation_masks_ptr;
      } else {
        observation_masks = h->observation_masks(new_x, approx_factor);
      }

      h->linearize(new_x, zero_observe_noise, &H, &N);

      K = matrix_div((KalmanT) (sigma*H.transpose()), (ObserveMatT) (H*sigma*H.transpose() + N*N.transpose())) * observation_masks.asDiagonal();

      if (z_ptr != NULL) {
        new_x = new_x + K * ((*z_ptr) - h->call(new_x, zero_observe_noise));
      }

      sigma = ensure_precision((VarianceT) (sigma - K*(H*sigma)));

      compose_belief(new_x, matrix_sqrt(sigma), &new_b);

      return new_b;
    }

    void linearize(const BeliefT& b
                 , const ControlT& u
                 , BeliefGradT* output_A
                 , BeliefControlGradT* output_B
                 , BeliefT* output_c
                  ) const {
      if (output_A) num_diff((boost::function<BeliefT (const BeliefT& )>) boost::bind(&BeliefFunc<StateFuncT, ObserveFuncT, BeliefT>::operator(), this, _1, u, (ObserveT*) NULL, (ObserveT*) NULL), b, belief_dim, this->epsilon, output_A);
      if (output_B) num_diff((boost::function<BeliefT (const ControlT& )>) boost::bind(&BeliefFunc<StateFuncT, ObserveFuncT, BeliefT>::operator(), this, b, _1, (ObserveT*) NULL, (ObserveT*) NULL), u, belief_dim, this->epsilon, output_B);
      if (output_c) *output_c = this->call(b, u);
    }


    //BeliefT call(const BeliefT& b, const ControlT& u) const {
    //  return operator()(b, u);
    //}

    BeliefT call(const BeliefT& b, const ControlT& u, ObserveT* z_ptr=NULL, ObserveT* observation_masks_ptr=NULL) const {
      return operator()(b, u, z_ptr, observation_masks_ptr);
    }

    virtual void extract_state(const BeliefT& belief, StateT* output_state) const {
      assert (belief.size() == belief_dim);
      assert (output_state != NULL);
      *output_state = belief.head(state_dim);
    }

    virtual void extract_sqrt_sigma(const BeliefT& belief, VarianceT* output_sqrt_sigma) const {
      assert (belief.size() == belief_dim);
      assert (output_sqrt_sigma != NULL);
      sqrt_sigma_vec_to_sqrt_sigma(belief.tail(sigma_dof), output_sqrt_sigma, state_dim);
    }

    virtual void compose_belief(const StateT& state, const VarianceT& sqrt_sigma, BeliefT* output_belief) const {
      assert (state.size() == state_dim);
      assert (output_belief != NULL);
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
      assert (output_sigma != NULL);
      sqrt_sigma_vec_to_sigma(belief.tail(sigma_dof), output_sigma, state_dim);
    }
  };

}
