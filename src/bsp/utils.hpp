#pragma once

#include <armadillo>
#include <mlpack/core.hpp>
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
  typedef Matrix<double, observe_dim, observe_dim> ObserveMatT; \
  typedef Matrix<double, state_dim, observe_dim> KalmanT;

#define ENSURE_VECTOR(T) BOOST_STATIC_ASSERT( (MatrixTraits<T>::cols == 1) );

#define BSP_CONCAT_IMPL( x, y, z ) x##y##z

#define BSP_CONCAT( x, y, z ) BSP_CONCAT_IMPL( x, y, z )

#define __ENSURE_SAME_SCALAR_TYPE(prefix, T1, T2) \
  typedef typename MatrixTraits<T1>::scalar_type BSP_CONCAT(_scalar_type_, prefix, _1); \
  typedef typename MatrixTraits<T2>::scalar_type BSP_CONCAT(_scalar_type_, prefix, _2); \
  BOOST_STATIC_ASSERT( (boost::is_same< BSP_CONCAT(_scalar_type_, prefix, _1), BSP_CONCAT(_scalar_type_, prefix, _2) >::value) );

#define ENSURE_SAME_SCALAR_TYPE(T1, T2) __ENSURE_SAME_SCALAR_TYPE(__COUNTER__, T1, T2)
  

namespace BSP {
  template<class MatT1, class MatT2>
  MatT1 matrix_div(const MatT1& A, const MatT2& B) {
    PartialPivLU<MatT2> solver(B);
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

  void AddVarArrays(OptProb& prob, int rows, const vector<int>& cols, const vector<double>& lbs, const vector<double>& ubs, const vector<string>& name_prefix, const vector<VarArray*>& newvars);

  void AddVarArrays(OptProb& prob, int rows, const vector<int>& cols, const vector<string>& name_prefix, const vector<VarArray*>& newvars);

  void AddVarArray(OptProb& prob, int rows, int cols, double lb, double ub, const string& name_prefix, VarArray& newvars);

  void AddVarArray(OptProb& prob, int rows, int cols, const string& name_prefix, VarArray& newvars);

  template<class VecT>
  arma::vec to_arma_vec(VecT x) {
    return arma::conv_to<arma::vec>::from(toDblVec(x));
  }

  template<class MatT>
  arma::mat to_arma_mat(MatT m) {
    arma::mat out(m.rows(), m.cols());
    for (int i = 0; i < m.rows(); ++i) {
      for (int j = 0; j < m.cols(); ++j) {
        out(i, j) = m(i, j);
      }
    }
    return out;
  }

  template<class VecT>
  VecT to_eigen_vec(arma::vec x) {
    return (VecT) toVectorXd(arma::conv_to<DblVec>::from(x));
  }

  template<class MeanT, class VarianceT>
  MeanT sample_gaussian(const MeanT& mean, const VarianceT& cov) {
    mlpack::distribution::GaussianDistribution distribution(to_arma_vec<MeanT>(mean), to_arma_mat<VarianceT>(cov));
    return to_eigen_vec<MeanT>(distribution.Random());
  }
}
