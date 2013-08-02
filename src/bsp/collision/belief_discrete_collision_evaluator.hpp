#pragma once

#include "bsp/common.hpp"
#include "bsp/collision/belief_collision_evaluator.hpp"

namespace BSPCollision {

  template< class BeliefFuncT >
  class BeliefDiscreteCollisionEvaluator : public BeliefCollisionEvaluator {
  public:
    typedef typename BeliefFuncT::Ptr BeliefFuncPtr;
    BeliefDiscreteCollisionEvaluator(ConfigurationPtr rad, const VarVector& vars, BeliefFuncPtr belief_func, OR::KinBody::LinkPtr endeffector) :
        m_env(rad->GetEnv()),
        m_cc(BeliefCollisionChecker::GetOrCreate(*m_env)),
        m_rad(rad),
        m_vars(vars),
        m_link2ind(),
        m_links(),
        m_filterMask(-1),
        belief_func(belief_func),
        endeffector(endeffector) {
      vector<KinBody::LinkPtr> links;
      vector<int> inds;
      rad->GetAffectedLinks(m_links, true, inds);
      for (int i=0; i < m_links.size(); ++i) {
        m_link2ind[m_links[i].get()] = inds[i];
      }
    }

    virtual void GetCollisionsCached(const DblVec& x, vector<BeliefCollision>& collisions) {
      double key = vecSum(x);
      vector<BeliefCollision>* it = m_cache.get(key);
      if (it != NULL) {
        RAVELOG_DEBUG("using cached collision check\n");
        collisions = *it;
      }
      else {
        RAVELOG_DEBUG("not using cached collision check\n");
        CalcCollisions(x, collisions);
        m_cache.put(key, collisions);
      }
    }

    virtual void CalcCollisions(const DblVec& x, vector<Collision>& collisions) {
      throw std::runtime_error("not implemented");
    }

    virtual VarVector GetVars() {
      return m_vars;
    }

    void CalcCollisions(const DblVec& x, vector<BeliefCollision>& collisions) {
      typename BeliefFuncT::BeliefT belief = getVec(x, this->m_vars);
      vector<DblVec> dofvals;
      belief_func->sigma_points(belief, &dofvals);
      m_cc->SigmaHullVsAll(*m_rad, m_links, dofvals, collisions);
    }

    virtual void CollisionsToDistanceExpressions(const vector<BeliefCollision>& collisions, Configuration& rad,
        const Link2Int& link2ind, const VarVector& theta_vars, const DblVec& theta_vals, vector<AffExpr>& exprs, bool isTimestep1, NamePairs& bodyNames) {
      exprs.clear();
      exprs.reserve(collisions.size());
      typename BeliefFuncT::StateT state;
      typename BeliefFuncT::BeliefT belief = toVectorXd(theta_vals);
      belief_func->extract_state(belief, &state);
      rad.SetDOFValues(toDblVec(state));
      BOOST_FOREACH(const BeliefCollision& col, collisions) {
        Link2Int::const_iterator itA = link2ind.find(col.linkA);
        Link2Int::const_iterator itB = link2ind.find(col.linkB);
        AffExpr dist;
        bool linkAFound = itA != link2ind.end();
        bool linkBFound = itB != link2ind.end();
        for (int i=0; i<col.mi[isTimestep1].alpha.size(); i++) {
          typename BeliefFuncT::SigmaPointsGradT grad;
          belief_func->sigma_points_grad(belief, col.mi[isTimestep1].instance_ind[i], &grad);

          AffExpr dist_a(col.distance);
          if (linkAFound) {
	          MatrixXd pos_jac = rad.PositionJacobian(itA->second, col.ptA);
            VectorXd dist_grad = toVector3d(col.normalB2A).transpose()*pos_jac*grad;
            exprInc(dist_a, varDot(dist_grad, theta_vars));
            exprInc(dist_a, -dist_grad.dot(toVectorXd(theta_vals)));
          }
          if (linkBFound) {
	          MatrixXd pos_jac = rad.PositionJacobian(itB->second, (isTimestep1 && (col.cctype == CCType_Between)) ? col.ptB1 : col.ptB);
            VectorXd dist_grad = -toVector3d(col.normalB2A).transpose()*pos_jac*grad;
            exprInc(dist_a, varDot(dist_grad, theta_vars));
            exprInc(dist_a, -dist_grad.dot(toVectorXd(theta_vals)));
          }
          if (linkAFound || linkBFound) {
            exprScale(dist_a, col.mi[isTimestep1].alpha[i]);
            exprInc(dist, dist_a);
            bodyNames.push_back(pair<string, string>(col.linkA->GetParent()->GetName(), col.linkB->GetParent()->GetName()));
          }
        }
        if (dist.constant!=0 || dist.coeffs.size()!=0 || dist.vars.size()!=0) {
          exprs.push_back(dist);
        }
      }
      RAVELOG_DEBUG("%i distance expressions\n", exprs.size());
    }

    virtual void CalcDistExpressions(const DblVec& x, vector<AffExpr>& exprs, NamePairs& bodyNames) {
      vector<BeliefCollision> collisions;
      GetCollisionsCached(x, collisions);
      DblVec theta = getDblVec(x, m_vars);
      CollisionsToDistanceExpressions(collisions, *m_rad, m_link2ind, m_vars, theta, exprs, false, bodyNames);
    }
    virtual void CalcDists(const DblVec& x, DblVec& dists, NamePairs& bodyNames) {
      vector<BeliefCollision> collisions;
      GetCollisionsCached(x, collisions);
      vector<Collision> raw_collisions;
      for (int i = 0; i < collisions.size(); ++i) {
        raw_collisions.push_back(collisions[i]);
      }
      this->CollisionsToDistances(raw_collisions, m_link2ind, dists, bodyNames);
    }
    virtual void CustomPlot(const DblVec& x, std::vector<OR::GraphHandlePtr>& handles) {
      typename BeliefFuncT::BeliefT belief = getVec(x, this->m_vars);
      vector<DblVec> dofvals;
      belief_func->sigma_points(belief, &dofvals);
      m_cc->PlotSigmaHull(*m_rad, m_links, dofvals, handles);
    }
  protected:
    OR::EnvironmentBasePtr m_env;
    BeliefCollisionCheckerPtr m_cc;
    ConfigurationPtr m_rad;
    VarVector m_vars;
    Link2Int m_link2ind;
    vector<OR::KinBody::LinkPtr> m_links;
    short m_filterMask;
    BeliefFuncPtr belief_func;
    OR::KinBody::LinkPtr endeffector;
    Cache<double, vector<BeliefCollision>, 3> m_cache;
  };
}
