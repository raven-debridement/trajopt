#include "common.hpp"
#include "constraints.hpp"

template< class BeliefFuncT=BeliefFunc<> >
class BeliefConstraint : public EqConstraint {
public:
  typedef BeliefFuncT::Ptr BeliefFuncPtr;
  typedef BeliefFuncT::BeliefT BeliefT;
  typedef BeliefFuncT::ControlT ControlT;
  typedef BeliefFuncT::BeliefGradT BeliefGradT;
  typedef BeliefFuncT::ControlGradT ControlGradT;

  BeliefConstraint::BeliefConstraint(const VarVector& cur_belief_vars, const VarVector& cur_control_vars, const VarVector& next_belief_vars, BeliefFuncPtr f) : cur_belief_vars(cur_belief_vars), cur_control_vars(cur_control_vars), next_belief_vars(next_belief_vars), f(f), belief_dim(cur_belief_vars.size()), control_dim(cur_control_vars.size()) {
    assert (cur_belief_vars.size() == next_belief_vars.size());
  }
  virtual vector<double> BeliefConstraint::value(const vector<double>& xvec) {
    BeliefT cur_belief = getVec(xvec, cur_belief_vars);
    ControlT cur_control = getVec(xvec, cur_control_vars);
    BeliefT next_belief = getVec(xvec, next_belief_vars);

    return toDblVec(next_belief_vars - f->call(cur_belief_vars, cur_control_vars));
  }
  virtual ConvexConstraintsPtr BeliefConstraint::convex(const vector<double>& xvec, Model* model) {
    BeliefT cur_belief = getVec(xvec, cur_belief_vars);
    ControlT cur_control = getVec(xvec, cur_control_vars);

    BeliefGradT A;
    ControlGradT B;
    BeliefT c;
    f->linearize(cur_belief, cur_control, &A, &B, &c);

    ConvexConstraintsPtr out(new ConvexConstraints(model));

    for (int i = 0; i < belief_dim; ++i) {
      AffExpr aff(next_belief_vars[i]);
      for (int j = 0; j < belief_dim; ++j) {
        exprInc(aff, exprMult(-A(i, j), exprAdd(AffExpr(cur_belief_vars[j]), -cur_belief[j])));
      }
      for (int j = 0; j < control_dim; ++j) {
        exprInc(aff, exprMult(-B(i, j), exprAdd(AffExpr(cur_control_vars[j]), -cur_control[j])));
      }
      exprInc(aff, -c(i));
      out->addEqCnt(aff);
    }

    return out;
  }
protected:
  VarVector cur_belief_vars;
  VarVector cur_control_vars;
  VarVector next_belief_vars;
  BeliefFuncPtr f;
  int belief_dim;
  int control_dim;
};
