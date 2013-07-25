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
    define_belief_func_types();
    typedef boost::shared_ptr< BeliefFunc<_StateFuncT, _ObserveFuncT, _BeliefT> > Ptr;

    double epsilon;
    double approx_factor;
    double sigma_pts_scale;
    BSPProblemHelperBasePtr helper;
    StateFuncPtr f;
    ObserveFuncPtr h;

    BeliefFunc() : epsilon(DefaultEpsilon), approx_factor(0.5) {}

    BeliefFunc(BSPProblemHelperBasePtr helper, StateFuncPtr f, ObserveFuncPtr h) :
      ProblemState(helper), helper(helper), f(f), h(h), epsilon(DefaultEpsilon), sigma_pts_scale(0.01), approx_factor(0.5) {}

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

    virtual BeliefT operator()(const BeliefT& b, const ControlT& u, ObserveT* z_ptr=nullptr, ObserveT* observation_masks_ptr=nullptr, StateT* state_error_out=nullptr) const = 0;

    BeliefT call(const BeliefT& b, const ControlT& u, ObserveT* z_ptr=nullptr, ObserveT* observation_masks_ptr=nullptr, StateT* state_error_out=nullptr) const {
      return operator()(b, u, z_ptr, observation_masks_ptr, state_error_out);
    }

    virtual void linearize(const BeliefT& b // state
                         , const ControlT& u // control
                         , BeliefGradT* output_A // df/dx
                         , BeliefControlGradT* output_B // df/du
                         , BeliefT* output_c // df/dm
                          ) const {
      if (output_A) {
        boost::function<BeliefT (const BeliefT& )> f_b;
        f_b = boost::bind(&BeliefFunc::operator(), this, _1, u, nullptr, nullptr, nullptr);
        num_diff(f_b, b, belief_dim, this->epsilon, output_A);
      }
      if (output_B) {
        boost::function<BeliefT (const ControlT& )> f_u;
        f_u = boost::bind(&BeliefFunc::operator(), this, b, _1, nullptr, nullptr, nullptr);
        num_diff(f_u, u, belief_dim, this->epsilon, output_B);
      }
      if (output_c) {
        *output_c = this->call(b, u);
      }
    }

    virtual void extract_state(const BeliefT& belief, StateT* output_state) const {
      assert (belief.size() == belief_dim);
      assert (output_state != nullptr);
      *output_state = belief.head(state_dim);
    }

    virtual void extract_sqrt_sigma(const BeliefT& belief, VarianceT* output_sqrt_sigma) const {
      assert (belief.size() == belief_dim);
      assert (output_sqrt_sigma != nullptr);
      sqrt_sigma_vec_to_sqrt_sigma(belief.tail(sigma_dof), output_sqrt_sigma, state_dim);
    }

    void extract_state_and_sqrt_sigma(const BeliefT& belief, StateT* output_state, VarianceT* output_sqrt_sigma) const {
      extract_state(belief, output_state);
      extract_sqrt_sigma(belief, output_sqrt_sigma);
    }

    void extract_state_and_sigma(const BeliefT& belief, StateT* output_state, VarianceT* output_sigma) const {
      extract_state(belief, output_state);
      extract_sigma(belief, output_sigma);
    }

    virtual void compose_belief(const StateT& state, const VarianceT& sqrt_sigma, BeliefT* output_belief) const {
      assert (state.size() == state_dim);
      assert (output_belief != nullptr);
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
      assert (output_sigma != nullptr);
      sqrt_sigma_vec_to_sigma(belief.tail(sigma_dof), output_sigma, state_dim);
    }

    virtual StateT sigma_point(const StateT& state, const VarianceT& sqrt_sigma, int sigma_pt_ind) const {
      assert (0 <= sigma_pt_ind && sigma_pt_ind <= 2*state_dim);
      if (sigma_pt_ind == 0) {
        return state;
      } else {
        if (sigma_pt_ind % 2 == 1) {
          return state + sigma_pts_scale*sqrt_sigma.col( (sigma_pt_ind-1) / 2 );
        } else {
          return state - sigma_pts_scale*sqrt_sigma.col( (sigma_pt_ind-1) / 2 );
        }
      }
    }
    
    virtual StateT sigma_point(const BeliefT& belief, int sigma_pt_ind) const {
      StateT state;
      VarianceT sqrt_sigma;
      extract_state_and_sqrt_sigma(belief, &state, &sqrt_sigma);
      return sigma_point(state, sqrt_sigma, sigma_pt_ind);
    }

    virtual void sigma_points(const BeliefT& belief, vector<DblVec>* output_sigma_points) const {
      assert (output_sigma_points != nullptr);
      StateT state;
      VarianceT sqrt_sigma;
      extract_state_and_sqrt_sigma(belief, &state, &sqrt_sigma);
      output_sigma_points->clear();
      //output_sigma_points->push_back( toDblVec(sigma_point(state, sqrt_sigma, 0)) );
      for (int i = 0; i <= 2*state_dim; ++i) {
        output_sigma_points->push_back( toDblVec(sigma_point(state, sqrt_sigma, i)) );
      }
    }

    virtual void sigma_points_grad(const BeliefT& belief, int sigma_pt_ind, SigmaPointsGradT* output_grad) const {
      assert (output_grad != nullptr);
      output_grad->resize(state_dim, belief_dim);
      BeliefT belief_plus = belief, belief_minus = belief;
      for (int i = 0; i < belief_dim; ++i) {
        belief_plus(i) += epsilon;
        belief_minus(i) -= epsilon;
        output_grad->col(i) = (sigma_point(belief_plus, sigma_pt_ind) - sigma_point(belief_minus, sigma_pt_ind)) / (2*epsilon);
        belief_plus(i) = belief_minus(i) = belief(i);
      }
    }
  };

}
