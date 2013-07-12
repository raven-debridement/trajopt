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
  class BeliefFunc : public ProblemState {
  public:

    define_belief_func_types;

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

    virtual BeliefT operator()(const BeliefT& b, const ControlT& u, ObserveT* z_ptr=NULL, ObserveT* observation_masks_ptr=NULL, StateT* state_error_out=NULL) const {
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

      //cout << "sigma before: " << sigma << endl;

      sigma = A*sigma*A.transpose() + M*M.transpose();
      //cout << "sigma after: " << sigma << endl;

      new_x = f->call(x, u, zero_state_noise);

      ObserveT observation_masks;
      if (observation_masks_ptr != NULL) {
        observation_masks = *observation_masks_ptr;
      } else {
        observation_masks = h->observation_masks(new_x, approx_factor);
      }
      //cout << "obs masks: " << observation_masks.transpose() << endl;

      h->linearize(new_x, zero_observe_noise, &H, &N);

      H = observation_masks.asDiagonal() * H;

      //cout << "N: " << endl << N << endl;
      //cout << "H: " << endl << H << endl;
      //cout << "sigma: " << endl << sigma << endl;
      //cout << "nominator:" << endl << (sigma*H.transpose()) << endl;
      //cout << "denominator:" << endl <<  (H*sigma*H.transpose() + N*N.transpose())<< endl;
      K = matrix_div((KalmanT) (sigma*H.transpose()), (ObserveMatT) (H*sigma*H.transpose() + N*N.transpose()));
      //cout << "K: " << endl << K << endl;

      //cout << "obs masks:\n" << observation_masks.asDiagonal() << endl;

      if (z_ptr != NULL) {
        StateT state_diff = K * observation_masks.asDiagonal() * ((*z_ptr) - h->call(new_x, zero_observe_noise));
        new_x = new_x + state_diff;
        if (state_error_out != NULL) {
          *state_error_out = state_diff;
        }
      }

      EigenSolver<VarianceT> es(sigma);
      //cout << "eigvals: " << es.eigenvalues() << endl;
      sigma = sigma - K*(H*sigma);//ensure_precision((VarianceT) (sigma - K*(H*sigma)));
      //VarianceT sqrtsigma = matrix_sqrt(sigma);
      //for (int i = 0; i < sqrtsigma.rows(); ++i) {
      //  for (int j = 0; j < sqrtsigma.cols(); ++j) {
      //    if (isnan(sqrtsigma(i, j))) {
      //      cout << "NAN!!" << endl;
      //      cout << "sigma: " << endl << sigma << endl;
      //      cout << "sqrtsigma: " << endl << sqrtsigma << endl;
      //      goto finish;
      //    }
      //  }
      //}
//finish:
      //for (int i = 0; i < new_x.rows(); ++i) {
      //  if (isnan(new_x(i))) {
      //    cout << "NANX!!" << endl;
      //  }
      //}


      //if (isnan(sigma.norm())) {
      //  cout << "NAN!!!" << endl;
      //}
      //cout << "new sigma: " << endl << sigma << endl;
      //cout << "new sqrt sigma:" << endl << matrix_sqrt(sigma) << endl;

      compose_belief(new_x, matrix_sqrt(sigma), &new_b);

      return new_b;
    }

    void linearize(const BeliefT& b
                 , const ControlT& u
                 , BeliefGradT* output_A
                 , BeliefControlGradT* output_B
                 , BeliefT* output_c
                  ) const {
      if (output_A) num_diff((boost::function<BeliefT (const BeliefT& )>) boost::bind(&BeliefFunc<StateFuncT, ObserveFuncT, BeliefT>::operator(), this, _1, u, (ObserveT*) NULL, (ObserveT*) NULL, (StateT*) NULL), b, belief_dim, this->epsilon, output_A);
      if (output_B) num_diff((boost::function<BeliefT (const ControlT& )>) boost::bind(&BeliefFunc<StateFuncT, ObserveFuncT, BeliefT>::operator(), this, b, _1, (ObserveT*) NULL, (ObserveT*) NULL, (StateT*) NULL), u, belief_dim, this->epsilon, output_B);
      if (output_c) *output_c = this->call(b, u);
    }


    //BeliefT call(const BeliefT& b, const ControlT& u) const {
    //  return operator()(b, u);
    //}

    BeliefT call(const BeliefT& b, const ControlT& u, ObserveT* z_ptr=NULL, ObserveT* observation_masks_ptr=NULL, StateT* state_error_out=NULL) const {
      return operator()(b, u, z_ptr, observation_masks_ptr, state_error_out);
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
