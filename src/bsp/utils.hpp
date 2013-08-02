#pragma once

#include "common.hpp"
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

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
    assert (out != nullptr);
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

  void AddVarArrays(OptProb& prob, int rows, const vector<int>& cols, const vector<double>& lbs, const vector<double>& ubs, const vector<string>& name_prefix, const vector<VarArray*>& newvars);

  void AddVarArrays(OptProb& prob, int rows, const vector<int>& cols, const vector<string>& name_prefix, const vector<VarArray*>& newvars);

  void AddVarArray(OptProb& prob, int rows, int cols, double lb, double ub, const string& name_prefix, VarArray& newvars);

  void AddVarArray(OptProb& prob, int rows, int cols, const string& name_prefix, VarArray& newvars);

  inline double rndnum() {
	  return ((double) random()) / RAND_MAX;
  }

  inline double normal() {
    double u_1 = 0;
    while (u_1 == 0) {
      u_1 = rndnum();
    }
    double u_2 = 0;
    while (u_2 == 0) {
      u_2 = rndnum();
    }
    return sqrt(-2*log(u_1)) * sin(2*PI*u_2);
  }

  template<class MeanT, class VarianceT>
  MeanT sample_gaussian(const MeanT& mean, const VarianceT& cov, double scaling=1.0/9.0) {
    static boost::normal_distribution<> standard_normal;
    static boost::mt19937 rng(static_cast<unsigned int>(std::time(0)));
    static boost::variate_generator<boost::mt19937&, 
                                    boost::normal_distribution<> > gen(rng, standard_normal);
    MeanT sample(mean.size());
    for (int i = 0; i < mean.size(); ++i) {
    	sample(i) = normal();
    }
    SelfAdjointEigenSolver<VarianceT> solver(cov*scaling);
    return (solver.eigenvectors().real())
         * (solver.eigenvalues().real().cwiseSqrt().asDiagonal())
         * sample + mean;
  }

  void seed_random();

  double sigmoid(double x);

  string get_current_directory();
  
}
