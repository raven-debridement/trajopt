#pragma once

#include "common.hpp"
#include "belief_func.hpp"
#include "constraints.hpp"

namespace BSP {
  template< class BeliefFuncT=BeliefFunc<> >
  class BeliefConstraint : public EqConstraint {
  public:
    typedef typename BeliefFuncT::Ptr BeliefFuncPtr;
    typedef typename BeliefFuncT::BeliefT BeliefT;
    typedef typename BeliefFuncT::ControlT ControlT;
    typedef typename BeliefFuncT::BeliefGradT BeliefGradT;
    typedef typename BeliefFuncT::BeliefControlGradT ControlGradT;

    typedef boost::shared_ptr< BeliefConstraint<BeliefFuncT> > Ptr;

    BeliefConstraint(const VarVector& cur_belief_vars, const VarVector& cur_control_vars, const VarVector& next_belief_vars, BeliefFuncPtr f) : cur_belief_vars(cur_belief_vars), cur_control_vars(cur_control_vars), next_belief_vars(next_belief_vars), f(f), belief_dim(cur_belief_vars.size()), control_dim(cur_control_vars.size()) {
      assert (cur_belief_vars.size() == next_belief_vars.size());
    }
    virtual vector<double> value(const vector<double>& xvec) {
      BeliefT cur_belief = getVec(xvec, cur_belief_vars);
      ControlT cur_control = getVec(xvec, cur_control_vars);
      BeliefT next_belief = getVec(xvec, next_belief_vars);

      return toDblVec(next_belief - f->call(cur_belief, cur_control));
    }
    virtual ConvexConstraintsPtr convex(const vector<double>& xvec, Model* model) {
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
          exprInc(aff, exprMult(exprSub(AffExpr(cur_belief_vars[j]), cur_belief[j]), -A(i, j)));
        }
        for (int j = 0; j < control_dim; ++j) {
          exprInc(aff, exprMult(exprSub(AffExpr(cur_control_vars[j]), cur_control[j]), -B(i, j)));
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
}
