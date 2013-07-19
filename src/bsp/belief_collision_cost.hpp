#pragma once

#include "common.hpp"
#include "belief_collision.hpp"
#include "belief_collision_plotter_mixin.hpp"

namespace BSP {
  using namespace BSPCollision;
  template< class BeliefFuncT >
  class BeliefCollisionCost : public Cost, public Plotter, public BeliefCollisionPlotterMixin<BeliefFuncT> {
  public:
    typedef typename BeliefFuncT::Ptr BeliefFuncPtr;

    /* constructor for single timestep */
    BeliefCollisionCost(double dist_pen, double coeff, ConfigurationPtr rad, const VarVector& vars, BeliefFuncPtr belief_func, OR::KinBody::LinkPtr endeffector) : Cost("belief_collision"), m_dist_pen(dist_pen), m_coeff(coeff), m_calc(new BeliefSingleTimestepCollisionEvaluator<BeliefFuncT>(rad, vars, belief_func, endeffector)) {}

    ConvexObjectivePtr convex(const vector<double>& x, Model* model) {
      ConvexObjectivePtr out(new ConvexObjective(model));
      vector<AffExpr> exprs;
      m_calc->CalcDistExpressions(x, exprs);
      for (int i=0; i < exprs.size(); ++i) {
        AffExpr viol = exprSub(AffExpr(m_dist_pen), exprs[i]);
        out->addHinge(viol, m_coeff);
      }
      return out;
    }

    double value(const vector<double>& x) {
      DblVec dists;
      m_calc->CalcDists(x, dists);
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
    boost::shared_ptr< BeliefSingleTimestepCollisionEvaluator<BeliefFuncT> > m_calc;
    double m_dist_pen;
    double m_coeff;

  };
}
