#include "common.hpp"
#include "costs.hpp"
#include "utils.hpp"

namespace BSP {

  VarianceCost::VarianceCost(const VarVector& sqrt_sigma_vars, const MatrixXd& Q) : Cost("variance"), sqrt_sigma_vars(sqrt_sigma_vars), Q(Q), state_dim(Q.cols()) {
    // vars ONLY contains the sqrt sigma entries in a belief vector
    // tr(Q*Sigma*Q') = sum_square(Q*SqrtSigma) 
    sigma_dof = state_dim * (state_dim + 1) / 2;
    for (int i = 0; i < state_dim; ++i) {
      for (int j = 0; j < state_dim; ++j) {
        AffExpr aff;
        for (int k = 0; k < state_dim; ++k) {
          exprInc(aff, exprMult(
            AffExpr(sqrt_sigma_vars[sigma_id(k, j)]),
            Q(i, k)
          ));
        }
        exprInc(expr, exprSquare(aff));
      }
    }
  }

  int VarianceCost::sigma_id(int i, int j) {
    return i <= j ? ((2 * state_dim - i + 1) * i / 2 + j - i)
                  : ((2 * state_dim - j + 1) * j / 2 + i - j);
  }

  double VarianceCost::value(const vector<double>& xvec) {
    VectorXd sqrt_sigma_vec = getVec(xvec, sqrt_sigma_vars);
    MatrixXd sigma;
    sqrt_sigma_vec_to_sigma(sqrt_sigma_vec, &sigma, state_dim);
    return (Q * sigma * Q.transpose()).trace();
  }

  ConvexObjectivePtr VarianceCost::convex(const vector<double>& xvec, Model* model) {
    ConvexObjectivePtr out(new ConvexObjective(model));
    out->addQuadExpr(expr);
    return out;
  }

  ControlCost::ControlCost(const VarVector& control_vars, const MatrixXd& R): Cost("control"), control_vars(control_vars), R(R), control_dim(R.cols()) {
    // lazy XD
    MatrixXd sqrt_R = matrix_sqrt(R);
    for (int i = 0; i < control_dim; ++i) {
      AffExpr aff;
      for (int j = 0; j < control_dim; ++j) {
        exprInc(aff, exprMult(control_vars[j], sqrt_R(i, j)));
      }
      exprInc(expr, exprSquare(aff));
    }
    //expr = exprSquare(aff);
  }

  double ControlCost::value(const vector<double>& xvec) {
    VectorXd u = getVec(xvec, control_vars);
    return u.transpose() * R * u;
  }

  ConvexObjectivePtr ControlCost::convex(const vector<double>& xvec, Model* model) {
    ConvexObjectivePtr out(new ConvexObjective(model));
    out->addQuadExpr(expr);
    return out;
  }
}
