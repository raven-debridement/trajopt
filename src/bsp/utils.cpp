#include "bsp.hpp"

namespace BSP {
  MatrixXd matrix_div(const MatrixXd& A, const MatrixXd& B) {
    PartialPivLU<MatrixXd> solver(B);
    return solver.solve(A.transpose()).transpose();
  }

  MatrixXd matrix_sqrt(const MatrixXd& A) {
    SelfAdjointEigenSolver<MatrixXd> solver(A);
    return (solver.eigenvectors().real()) // U
         * (solver.eigenvalues().real().cwiseSqrt().asDiagonal()) // sqrt(Sigma)
         * (solver.eigenvectors().real().transpose()); // U'
  }

  MatrixXd num_diff(const boost::function<VectorXd>& f
                  , const VectorXd& cur // current point
                  , int m // output dimension
                  , double epsilon
                   ) {
    int n = cur.size();
    VectorXd x = cur;
    MatrixXd out(n, m);
    for (int i = 0; i < n; ++i) {
      x(i) = cur(i) + epsilon;
      VectorXd fplus = f(x);
      x(i) = cur(i) - epsilon;
      VectorXd fminus = f(x);
      out.col(i) = (fplus - fminus) / (2 * epsilon);
      x(i) = cur(i);
    }
    return out;
  }
}
