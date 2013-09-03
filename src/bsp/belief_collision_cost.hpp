#pragma once

#include "common.hpp"
#include "collision/collision.hpp"
#include "belief_collision_plotter_mixin.hpp"

namespace BSP {
  using namespace BSPCollision;
  template< class BeliefFuncT >
  class BeliefCollisionCost : public Cost, public Plotter, public BeliefCollisionPlotterMixin<BeliefFuncT> {
  public:
    typedef typename BeliefFuncT::Ptr BeliefFuncPtr;

    /* constructor for single timestep */
    BeliefCollisionCost(double dist_pen, double coeff, ConfigurationPtr rad, const VarVector& vars, BeliefFuncPtr belief_func, OR::KinBody::LinkPtr endeffector) : Cost("belief_collision"), m_dist_pen(dist_pen), m_coeff(coeff), m_calc(new BeliefDiscreteCollisionEvaluator<BeliefFuncT>(rad, vars, belief_func, endeffector)) {}

    /* constructor for continuous */
    BeliefCollisionCost(double dist_pen, double coeff, ConfigurationPtr rad, const VarVector& vars0, const VarVector& vars1, BeliefFuncPtr belief_func, OR::KinBody::LinkPtr endeffector) : Cost("belief_collision"), m_dist_pen(dist_pen), m_coeff(coeff), m_calc(new BeliefContinuousCollisionEvaluator<BeliefFuncT>(rad, vars0, vars1, belief_func, endeffector)) {}

    ConvexObjectivePtr convex(const vector<double>& x) {
      ConvexObjectivePtr out(new ConvexObjective());
      vector<AffExpr> exprs;
      NamePairs bodyNames;
      m_calc->CalcDistExpressions(x, exprs, bodyNames);
      for (int i=0; i < exprs.size(); ++i) {
        AffExpr viol = exprSub(AffExpr(m_dist_pen), exprs[i]);
        out->addHinge(viol, m_coeff);
      }
      return out;
    }

    double value(const vector<double>& x, Model* model) {
      DblVec dists;
      NamePairs bodyNames;
      m_calc->CalcDists(x, dists, bodyNames);
      double out = 0;
      for (int i=0; i < dists.size(); ++i) {
        out += pospart(m_dist_pen - dists[i]) * m_coeff;
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
