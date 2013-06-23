#pragma once

#include "common.hpp"

#define BSP_TYPEDEFS(state_dim, state_noise_dim, control_dim, observe_dim, observe_noise_dim, sigma_dof, belief_dim) \
  typedef Matrix<double, state_dim, 1> StateT; \
  typedef Matrix<double, state_noise_dim, 1> StateNoiseT; \
  typedef Matrix<double, control_dim, 1> ControlT; \
  typedef Matrix<double, observe_dim, 1> ObserveT; \
  typedef Matrix<double, observe_noise_dim, 1> ObserveNoiseT; \
  typedef Matrix<double, belief_dim, 1> BeliefT; \
  typedef Matrix<double, state_dim, state_dim> VarianceT; \
  typedef Matrix<double, state_dim, state_dim> VarianceCostT; \
  typedef Matrix<double, state_dim, state_dim> StateGradT; \
  typedef Matrix<double, control_dim, control_dim> ControlGradT; \
  typedef Matrix<double, control_dim, control_dim> ControlCostT; \
  typedef Matrix<double, state_dim, state_noise_dim> StateNoiseGradT; \
  typedef Matrix<double, observe_dim, state_dim> ObserveStateGradT; \
  typedef Matrix<double, observe_dim, observe_noise_dim> ObserveNoiseGradT; \
  typedef Matrix<double, state_dim, state_dim> KalmanT; \
  typedef Matrix<double, state_dim, state_dim> GammaT;

namespace BSP {
  template<class MatType>
  MatType matrix_div(const MatType& A, const MatType& B) {
    PartialPivLU<MatType> solver(B);
    return solver.solve(A.transpose()).transpose();
  }

  template<class MatType>
  MatType matrix_sqrt(const MatType& A) {
    SelfAdjointEigenSolver<MatType> solver(A);
    return (solver.eigenvectors().real()) // U
         * (solver.eigenvalues().real().cwiseSqrt().asDiagonal()) // sqrt(Sigma)
         * (solver.eigenvectors().real().transpose()); // U'
  }

  template<class InputType, class OutputType, class MatType>
  void num_diff(const boost::function<OutputType (const InputType& )>& f
              , const InputType& cur // current point
              , int m // output dimension
              , double epsilon
              , MatType* out
               ) {
    assert (out != NULL);
    int n = cur.size();
    InputType x = cur;
    out->resize(m, n);
    for (int i = 0; i < n; ++i) {
      x(i) = cur(i) + epsilon;
      OutputType fplus = f(x);
      x(i) = cur(i) - epsilon;
      OutputType fminus = f(x);
      out->col(i) = (fplus - fminus) / (2 * epsilon);
      x(i) = cur(i);
    }
  }

  template<class VecT, class MatT>
  void sqrt_sigma_vec_to_sqrt_sigma(const VecT& sqrt_sigma_vec, MatT* sqrt_sigma, int dim) {
    sqrt_sigma->resize(dim, dim);
    for (int index = 0, i = 0; i < dim; ++i) {
      for (int j = i; j < dim; ++j) {
        // the upper diagonal entries of sqrt_sigma are stored in row order in sqrt_sigma_vec
        (*sqrt_sigma)(i, j) = (*sqrt_sigma)(j, i) = sqrt_sigma_vec(index++);
      }
    }
  }

  template<class VecT, class MatT>
  void sqrt_sigma_vec_to_sigma(const VecT& sqrt_sigma_vec, MatT* sigma, int dim) {
    sqrt_sigma_vec_to_sqrt_sigma(sqrt_sigma_vec, sigma, dim);
    (*sigma) = (*sigma) * (sigma->transpose());
  }

  template<typename T, size_t N>
  T* end(T (&ra)[N]) {
    return ra + N;
  }

}
