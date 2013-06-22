#include "bsp.hpp"

namespace BSP {

  void StateFunc::linearize(const VectorXd& x // state
                          , const VectorXd& u // control
                          , const VectorXd& m // state noise
                          , MatrixXd* output_A // df/dx
                          , MatrixXd* output_B // df/du
                          , MatrixXd* output_M // df/dm
                           ) const {
    if (output_A) *output_A = num_diff(boost::bind(&StateFunc::operator(), this, _1, u, m), x, output_size(), this->epsilon);
    if (output_B) *output_B = num_diff(boost::bind(&StateFunc::operator(), this, x, _1, m), u, output_size(), this->epsilon);
    if (output_M) *output_M = num_diff(boost::bind(&StateFunc::operator(), this, x, u, _1), m, output_size(), this->epsilon);
  }

  void ObserveFunc::linearize(const VectorXd& x  // state
                            , const VectorXd& n  // observation noise
                            , MatrixXd* output_H // df/dx
                            , MatrixXd* output_N // df/dn
                             ) {
    if (output_H) *output_H = num_diff(boost::bind(&ObserveFunc::operator(), this, _1, n), x, output_size(), this->epsilon);
    if (output_N) *output_N = num_diff(boost::bind(&ObserveFunc::operator(), this, x, _1), n, output_size(), this->epsilon);
  }

  void BeliefFunc::linearize(const VectorXd& b // belief
                           , const VectorXd& u // control
                           , const StateFunc& f // state propagation function
                           , const ObserveFunc& h // observation function
                           , MatrixXd* output_A // df/db
                           , MatrixXd* output_B // df/du
                           , VectorXd* output_c // new belief at current point
                            ) {
    if (output_A) *output_A = num_diff(boost::bind(&BeliefFunc::operator(), this, _1, u, f, h), b, output_size(), this->epsilon);
    if (output_B) *output_B = num_diff(boost::bind(&BeliefFunc::operator(), this, b, _1, f, h), u, output_size(), this->epsilon);
    if (output_c) *output_c = this->call(b, u, f, h);
  }

  void BeliefFunc::extract_state(const VectorXd& belief, VectorXd* output_state) const {
    assert (belief.size() == helper->belief_dim);
    *output_state = belief.head(helper->state_dim);
  }

  void BeliefFunc::extract_sqrt_sigma(const VectorXd& belief, MatrixXd* output_sqrt_sigma) const {
    assert (belief.size() == helper->belief_dim);
    output_sqrt_sigma->resize(helper->state_dim, helper->state_dim);
    for (int index = helper->state_dim, i = 0; i < helper->state_dim; ++i) {
      for (int j = i; j < helper->state_dim; ++j) {
        // the upper diagonal entries of sqrt_sigma are stored in row order in the belief vector
        *output_sqrt_sigma(i, j) = *output_sqrt_sigma(j, i) = belief(index++);
      }
    }
  }

  void BeliefFunc::compose_belief(const VectorXd& state
                                , const MatrixXd& sqrt_sigma
                                , VectorXd* output_belief
                                 ) const {
    assert (state.size() == helper->state_dim);
    assert (output_belief->resize(helper->belief_dim));
    output_belief->head(helper->state_dim) = state;
    for (int index = helper->state_dim, i = 0; i < helper->state_dim; ++i) {
      for (int j = i; j < helper->state_dim; ++j) {
        *output_belief(index++) = 0.5 * (sqrt_sigma(i, j) + sqrt_sigma(j, i));
      }
    }
  }

  void BeliefFunc::extract_sigma(const VectorXd& belief, MatrixXd* output_sigma) {
    extract_sqrt_sigma(belief, output_sigma);
    (*output_sigma) *= output_sigma->transpose();
  }

}
