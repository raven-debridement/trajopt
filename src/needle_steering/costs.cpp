#include "needle_steering.hpp"

namespace Needle {

  ConstantSpeedCost::ConstantSpeedCost(const Var& var, double coeff, NeedleProblemHelperPtr helper) : Cost("Speed"), var(var), coeff(coeff), helper(helper) {
    assert (helper->speed_formulation == NeedleProblemHelper::ConstantSpeed);
    exprInc(expr, exprMult(var, coeff * helper->T));
  }

  double ConstantSpeedCost::value(const vector<double>& xvec) {
    double speed = getVec(xvec, singleton<Var>(var))[0];
    return speed * coeff * helper->T;
  }

  ConvexObjectivePtr ConstantSpeedCost::convex(const vector<double>& xvec, Model* model) {
    ConvexObjectivePtr out(new ConvexObjective(model));
    out->addAffExpr(expr);
    return out;
  }

  VariableSpeedCost::VariableSpeedCost(const VarVector& vars, double coeff, NeedleProblemHelperPtr helper) : Cost("Speed"), vars(vars), coeff(coeff), helper(helper) {
    assert (helper->speed_formulation == NeedleProblemHelper::VariableSpeed);
    for (int i = 0; i < vars.size(); ++i) {
      exprInc(expr, exprMult(vars[i], coeff));//exprSquare(exprAdd(AffExpr(vars[i]), -deviation)), coeff));
    }
  }

  double VariableSpeedCost::value(const vector<double>& xvec) {
    VectorXd speeds = getVec(xvec, vars);
    return coeff * speeds.array().sum();
  }

  ConvexObjectivePtr VariableSpeedCost::convex(const vector<double>& xvec, Model* model) {
    ConvexObjectivePtr out(new ConvexObjective(model));
    out->addAffExpr(expr);
    return out;
  }

  SpeedDeviationCost::SpeedDeviationCost(const VarVector& vars, double deviation, double coeff, NeedleProblemHelperPtr helper) : Cost("Speed"), vars(vars), deviation(deviation), coeff(coeff), helper(helper) {
    assert (helper->speed_formulation == NeedleProblemHelper::VariableSpeed);
    for (int i = 0; i < vars.size(); ++i) {
      exprInc(expr, exprMult(exprSquare(exprAdd(AffExpr(vars[i]), -deviation)), coeff));
    }
  }

  double SpeedDeviationCost::value(const vector<double>& xvec) {
    VectorXd speeds = getVec(xvec, vars);
    return coeff * (speeds.array() - deviation).square().sum();
  }

  ConvexObjectivePtr SpeedDeviationCost::convex(const vector<double>& xvec, Model* model) {
    ConvexObjectivePtr out(new ConvexObjective(model));
    out->addQuadExpr(expr);
    return out;
  }

  RotationQuadraticCost::RotationQuadraticCost(const VarVector& vars, double coeff, NeedleProblemHelperPtr helper) : Cost("Rotation"), vars(vars), coeff(coeff), helper(helper) {
    for (int i = 0; i < vars.size(); ++i) {
      exprInc(expr, exprMult(exprSquare(vars[i]), coeff));
    }
  }

  double RotationQuadraticCost::value(const vector<double>& xvec) {
    VectorXd vals = getVec(xvec, vars);
    return vals.array().square().sum() * coeff;
  }

  ConvexObjectivePtr RotationQuadraticCost::convex(const vector<double>& xvec, Model* model) {
    ConvexObjectivePtr out(new ConvexObjective(model));
    out->addQuadExpr(expr);
    return out;
  }

  RotationL1Cost::RotationL1Cost(const VarVector& vars, double coeff, NeedleProblemHelperPtr helper) : Cost("Rotation"), vars(vars), coeff(coeff), helper(helper) {}

  double RotationL1Cost::value(const vector<double>& xvec) {
    VectorXd vals = getVec(xvec, vars);
    return vals.array().abs().sum() * coeff;
  }

  ConvexObjectivePtr RotationL1Cost::convex(const vector<double>& xvec, Model* model) {
    ConvexObjectivePtr out(new ConvexObjective(model));
    for (int i = 0; i < vars.size(); ++i) {
      out->addAbs(AffExpr(vars[i]), coeff);
    }
    return out;
  }
}
