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
  class UkfBeliefFunc : public BeliefFunc<_StateFuncT, _ObserveFuncT, _BeliefT> {
  public:
    define_belief_func_types();
    typedef boost::shared_ptr< UkfBeliefFunc<_StateFuncT, _ObserveFuncT, _BeliefT> > Ptr;

    UkfBeliefFunc() : BeliefFunc<_StateFuncT, _ObserveFuncT, _BeliefT>() {}

    UkfBeliefFunc(BSPProblemHelperBasePtr helper, StateFuncPtr f, ObserveFuncPtr h) :
      BeliefFunc<_StateFuncT, _ObserveFuncT, _BeliefT>(helper, f, h) {}

    virtual BeliefT operator()(const BeliefT& b, const ControlT& u, ObserveT* z_ptr=nullptr, ObserveT* observation_masks_ptr=nullptr, StateT* state_error_out=nullptr) const {
      // TODO
      int L = state_dim + state_noise_dim + observe_noise_dim;
      double ukfalpha = 0.001; //alpha;
      double ukfbeta = 2; //beta;
      double ukfkappa = 0; //kappa;

      double lambda = ukfalpha*ukfalpha*(L + ukfkappa) - L;
      double w = 1 / (2*(L + lambda));
      double mw = lambda / (L + lambda);
      double vw = mw + (1 - ukfalpha*ukfalpha + ukfbeta);

      StateT x(this->state_dim);
      VarianceT sigma(this->state_dim, this->state_dim);
      extract_state_and_sigma(b, &x, &sigma);

	    // propagate sigma points through f
      StateNoiseT state_noise = StateNoiseT::Zero(this->state_noise_dim);
      StateSigmaPointsT sigmapts(this->state_dim, this->state_sigma_points_dim);

      int idx = 0;
      sigmapts.col(idx++) = this->f->call(x, u, state_noise);

      double step = sqrt(L+lambda);
      VarianceT scaled_sigma = step * sigma;

      for (int i = 0; i < this->state_dim; ++i) {
        sigmapts.col(idx++) = this->f->call(x + scaled_sigma.col(i), u, state_noise);
        sigmapts.col(idx++) = this->f->call(x - scaled_sigma.col(i), u, state_noise);
      }

      for (int i = 0; i < this->state_noise_dim; ++i) {
        state_noise(i) = step;
        sigmapts.col(idx++) = this->f->call(x, u, state_noise);
        state_noise(i) = -step;
        sigmapts.col(idx++) = this->f->call(x, u, state_noise);
        state_noise(i) = 0;
      }

	    // calculate mean -- O(state_dim^2)
      StateT new_x = (mw + 2*observe_noise_dim*w) * sigmapts.col(0);
      for (int i = 1; i < sigmapts.cols(); ++i) {
        new_x += w * sigmapts.col(i);
      }

  	  // calculate variance -- O(state_dim^3)
      VarianceT new_sigma = (vw + 2*r_dim*w) * (sigmapts.col(0) - new_x)*(sigmapts.col(0) - new_x).transpose();
      for (int i = 1; i < sigmapts.cols(); ++i) {
        new_sigma += w * (sigmapts.col(i) - new_x)*(sigmapts.col(i) - new_x).transpose();
      }

      ObserveSigmaPointsT Z(this->observe_dim, this->observe_sigma_points_dim);
      ObserveNoiseT observe_noise = ObserveNoiseT::Zero(observe_noise_dim);

      idx = 0;
      for (int i = 0; i < sigmapts.cols(); ++i) {
        Z.col(idx++) = this->h->call(sigmapts.col(i), observe_noise); // TODO add observation masks
      }

      for (int i = 0; i < observe_noise_dim; ++i) {
        observe_noise(i) = step;
        Z.col(idx++) = this->h->call(sigmapts.col(0), observe_noise);
      }
    }
  };
}
