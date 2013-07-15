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
  class EkfBeliefFunc : public BeliefFunc<_StateFuncT, _ObserveFuncT, _BeliefT> {
  public:
    define_belief_func_types();
    typedef boost::shared_ptr< EkfBeliefFunc<_StateFuncT, _ObserveFuncT, _BeliefT> > Ptr;

    EkfBeliefFunc() : BeliefFunc<_StateFuncT, _ObserveFuncT, _BeliefT>() {}

    EkfBeliefFunc(BSPProblemHelperBasePtr helper, StateFuncPtr f, ObserveFuncPtr h) :
      BeliefFunc<_StateFuncT, _ObserveFuncT, _BeliefT>() {}

    virtual BeliefT operator()(const BeliefT& b, const ControlT& u, ObserveT* z_ptr=nullptr, ObserveT* observation_masks_ptr=nullptr, StateT* state_error_out=nullptr) const {
      StateT             x(this->state_dim), new_x(this->state_dim);
      BeliefT            new_b(this->belief_dim);
      VarianceT          sigma(this->state_dim, this->state_dim);
      StateGradT         A(this->state_dim, this->state_dim);
      StateNoiseGradT    M(this->state_dim, this->state_noise_dim);
      ObserveStateGradT  H(this->observe_dim, this->state_dim);
      ObserveNoiseGradT  N(this->observe_dim, this->observe_noise_dim);
      KalmanT            K(this->state_dim, this->observe_dim);

      StateNoiseT        zero_state_noise = StateNoiseT::Zero(this->state_noise_dim);
      ObserveNoiseT      zero_observe_noise = ObserveNoiseT::Zero(this->observe_noise_dim);

      extract_state(b, &x);
      extract_sigma(b, &sigma);

      this->f->linearize(x, u, zero_state_noise, &A, nullptr, &M);

      sigma = A*sigma*A.transpose() + M*M.transpose();

      new_x = this->f->call(x, u, zero_state_noise);

      ObserveT observation_masks;
      if (observation_masks_ptr != nullptr) {
        observation_masks = *observation_masks_ptr;
      } else {
        observation_masks = this->h->observation_masks(new_x, this->approx_factor);
      }

      this->h->linearize(new_x, zero_observe_noise, &H, &N);

      H = observation_masks.asDiagonal() * H;

      K = matrix_div((KalmanT) (sigma*H.transpose()), (ObserveMatT) (H*sigma*H.transpose() + N*N.transpose()));

      if (z_ptr != nullptr) {
        StateT state_diff = K * observation_masks.asDiagonal() * ((*z_ptr) - this->h->call(new_x, zero_observe_noise));
        new_x = new_x + state_diff;
        if (state_error_out != nullptr) {
          *state_error_out = state_diff;
        }
      }

      EigenSolver<VarianceT> es(sigma);
      sigma = sigma - K*(H*sigma);

      compose_belief(new_x, matrix_sqrt(sigma), &new_b);

      return new_b;
    }
  };
}
