#pragma once

#include "common.hpp"
#include "belief_collision.hpp"
#include "belief_collision_plotter_mixin.hpp"

namespace BSP {
  using namespace BSPCollision;
  template< class BeliefFuncT >
  class BeliefCollisionConstraint : public IneqConstraint, public Plotter, public BeliefCollisionPlotterMixin<BeliefFuncT> {
  public:
    typedef typename BeliefFuncT::Ptr BeliefFuncPtr;

    /* constructor for single timestep */
    BeliefCollisionConstraint(double dist_pen, double coeff, ConfigurationPtr rad, const VarVector& vars, BeliefFuncPtr belief_func, OR::KinBody::LinkPtr endeffector) : m_dist_pen(dist_pen), m_coeff(coeff), m_calc(new BeliefDiscreteCollisionEvaluator<BeliefFuncT>(rad, vars, belief_func, endeffector)) {
      name_ = "belief_collision"; 
    }

    /* constructor for continuous */
    BeliefCollisionConstraint(double dist_pen, double coeff, ConfigurationPtr rad, const VarVector& vars0, const VarVector& vars1, BeliefFuncPtr belief_func, OR::KinBody::LinkPtr endeffector) : m_dist_pen(dist_pen), m_coeff(coeff), m_calc(new BeliefContinuousCollisionEvaluator<BeliefFuncT>(rad, vars0, vars1, belief_func, endeffector)) {
      name_ = "belief_collision"; 
    }

    ConvexConstraintsPtr convex(const vector<double>& x, Model* model) {
      ConvexConstraintsPtr out(new ConvexConstraints(model));
      vector<AffExpr> exprs;
      m_calc->CalcDistExpressions(x, exprs);
      for (int i=0; i < exprs.size(); ++i) {
        AffExpr viol = exprSub(AffExpr(m_dist_pen), exprs[i]);
        out->addIneqCnt(exprMult(viol,m_coeff));
      }
      return out;
    }
    DblVec value(const vector<double>& x) {
      DblVec dists;
      m_calc->CalcDists(x, dists);
      DblVec out(dists.size());
      for (int i=0; i < dists.size(); ++i) {
        out[i] = pospart(m_dist_pen - dists[i]) * m_coeff;
      }
      return out;
    }
    void Plot(const DblVec& x, OR::EnvironmentBase& env, std::vector<OR::GraphHandlePtr>& handles) {
      BeliefCollisionPlotterMixin<BeliefFuncT>::Plot(x, env, handles, m_calc, m_dist_pen);
    }
  protected:
    boost::shared_ptr<BeliefCollisionEvaluator> m_calc;
    double m_dist_pen;
    double m_coeff;
  };
}
