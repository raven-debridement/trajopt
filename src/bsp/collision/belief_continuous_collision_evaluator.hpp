#pragma once

#include "bsp/common.hpp"
#include "bsp/collision/belief_collision_evaluator.hpp"

namespace BSPCollision {

  template< class BeliefFuncT >
  class BeliefContinuousCollisionEvaluator : public BeliefDiscreteCollisionEvaluator<BeliefFuncT> {
  public:
    typedef typename BeliefFuncT::Ptr BeliefFuncPtr;
    BeliefContinuousCollisionEvaluator(ConfigurationPtr rad, const VarVector& vars0, const VarVector& vars1, BeliefFuncPtr belief_func, OR::KinBody::LinkPtr endeffector) :
        BeliefDiscreteCollisionEvaluator<BeliefFuncT>(rad, concat(vars0, vars1), belief_func, endeffector),
        m_vars0(vars0),
        m_vars1(vars1) {}

    void CalcCollisions(const DblVec& x, vector<BeliefCollision>& collisions) {
      typename BeliefFuncT::BeliefT belief0 = getVec(x, this->m_vars0);
      typename BeliefFuncT::BeliefT belief1 = getVec(x, this->m_vars1);
      vector<DblVec> dofvals0, dofvals1;
      this->belief_func->sigma_points(belief0, &dofvals0);
      this->belief_func->sigma_points(belief1, &dofvals1);
      this->m_cc->SigmaHullCastVsAll(*this->m_rad, this->m_links, dofvals0, dofvals1, collisions);
    }

    virtual void CollisionsToDistanceExpressions(const vector<BeliefCollision>& collisions, Configuration& rad,
        const Link2Int& link2ind, const VarVector& theta_vars0, const VarVector& theta_vars1, const DblVec& theta_vals0, const DblVec& theta_vals1, vector<AffExpr>& exprs, bool isTimestep1) {
      vector<AffExpr> exprs0, exprs1;
      BeliefDiscreteCollisionEvaluator<BeliefFuncT>::CollisionsToDistanceExpressions(collisions, rad, link2ind, theta_vars0, theta_vals0, exprs0, false);
      BeliefDiscreteCollisionEvaluator<BeliefFuncT>::CollisionsToDistanceExpressions(collisions, rad, link2ind, theta_vars1, theta_vals1, exprs1, true);

      exprs.resize(exprs0.size());
      for (int i=0; i < exprs0.size(); ++i) {
        exprScale(exprs0[i], (1-collisions[i].time));
        exprScale(exprs1[i], collisions[i].time);
        exprs[i] = AffExpr(0);
        exprInc(exprs[i], exprs0[i]);
        exprInc(exprs[i], exprs1[i]);
        cleanupAff(exprs[i]);
      }
    }

    virtual void CalcDistExpressions(const DblVec& x, vector<AffExpr>& exprs) {
      vector<BeliefCollision> collisions;
      this->GetCollisionsCached(x, collisions);
      DblVec theta0 = getDblVec(x, m_vars0);
      DblVec theta1 = getDblVec(x, m_vars1);
      CollisionsToDistanceExpressions(collisions, *this->m_rad, this->m_link2ind, m_vars0, m_vars1, theta0, theta1, exprs, false);
    }

    virtual void CustomPlot(const DblVec& x, std::vector<OR::GraphHandlePtr>& handles) {
      typename BeliefFuncT::BeliefT belief0 = getVec(x, this->m_vars0);
      typename BeliefFuncT::BeliefT belief1 = getVec(x, this->m_vars1);
      vector<DblVec> dofvals0, dofvals1;
      this->belief_func->sigma_points(belief0, &dofvals0);
      this->belief_func->sigma_points(belief1, &dofvals1);
      this->m_cc->PlotSigmaHullCast(*this->m_rad, this->m_links, dofvals0, dofvals1, handles);
    }
  protected:
    VarVector m_vars0, m_vars1;
  };
}
