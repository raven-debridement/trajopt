#pragma once

#include "common.hpp"
#include "belief_collision_cost.hpp"

namespace BSP {
  class VarianceCost : public Cost {
  public:
    VarianceCost(const VarVector& sqrt_sigma_vars, const MatrixXd& Q);
    int sigma_id(int i, int j);
    double value(const vector<double>& xvec, Model* model);
    ConvexObjectivePtr convex(const vector<double>& xvec);
  protected:
    VarVector sqrt_sigma_vars;
    QuadExpr expr;
    int state_dim;
    int sigma_dof;
    MatrixXd Q;
  };

  class ControlCost : public Cost {
  public:
    ControlCost(const VarVector& control_vars, const MatrixXd& R);
    double value(const vector<double>& xvec, Model* model);
    ConvexObjectivePtr convex(const vector<double>& xvec);
  protected:
    VarVector control_vars;
    QuadExpr expr;
    int control_dim;
    MatrixXd R;
  };
}
