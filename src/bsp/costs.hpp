#include "common.hpp"

namespace BSP {
  class VarianceCost : public Cost {
  public:
    VarianceCost(const VarVector& sqrt_sigma_vars, const MatrixXd& Q);
    int sigma_id(int i, int j);
    virtual double value(const vector<double>& xvec);
    virtual ConvexObjectivePtr convex(const vector<double>& xvec, Model* model);
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
    virtual double value(const vector<double>& xvec);
    virtual ConvexObjectivePtr convex(const vector<double>& xvec, Model* model);
  protected:
    VarVector control_vars;
    QuadExpr expr;
    int control_dim;
    MatrixXd R;
  };
}
