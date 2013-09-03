#include "trajopt/collision_terms.hpp"
#include "trajopt/collision_checker.hpp"
#include "trajopt/rave_utils.hpp"
#include "trajopt/utils.hpp"
#include "sco/expr_vec_ops.hpp"
#include "sco/expr_ops.hpp"
#include "sco/sco_common.hpp"
#include <boost/foreach.hpp>
#include "utils/eigen_conversions.hpp"
#include "sco/modeling_utils.hpp"
#include "utils/stl_to_string.hpp"
#include "utils/logging.hpp"
#include <iostream>
#include <boost/functional/hash.hpp>
using namespace OpenRAVE;
using namespace sco;
using namespace util;
using namespace std;

namespace trajopt {


void CollisionEvaluator::CollisionsToDistances(const vector<Collision>& collisions, const Link2Int& m_link2ind,
    DblVec& dists) {
  // Note: this checking (that the links are in the list we care about) is probably unnecessary
  // since we're using LinksVsAll
  dists.clear();
  dists.reserve(collisions.size());
  BOOST_FOREACH(const Collision& col, collisions) {
    Link2Int::const_iterator itA = m_link2ind.find(col.linkA);
    Link2Int::const_iterator itB = m_link2ind.find(col.linkB);
    if (itA != m_link2ind.end() || itB != m_link2ind.end()) {
      dists.push_back(col.distance);
    }
  }
}

void CollisionEvaluator::CollisionsToDistanceExpressions(const vector<Collision>& collisions, Configuration& rad,
    const Link2Int& link2ind, const VarVector& vars, const DblVec& dofvals, vector<AffExpr>& exprs, bool isTimestep1) {

  exprs.clear();
  exprs.reserve(collisions.size());
  rad.SetDOFValues(dofvals); // since we'll be calculating jacobians
  BOOST_FOREACH(const Collision& col, collisions) {
    AffExpr dist(col.distance);
    Link2Int::const_iterator itA = link2ind.find(col.linkA);
    if (itA != link2ind.end()) {
      VectorXd dist_grad = toVector3d(col.normalB2A).transpose()*rad.PositionJacobian(itA->second, col.ptA);
      exprInc(dist, varDot(dist_grad, vars));
      exprInc(dist, -dist_grad.dot(toVectorXd(dofvals)));
    }
    Link2Int::const_iterator itB = link2ind.find(col.linkB);
    if (itB != link2ind.end()) {
      VectorXd dist_grad = -toVector3d(col.normalB2A).transpose()*rad.PositionJacobian(itB->second, (isTimestep1 && (col.cctype == CCType_Between)) ? col.ptB1 : col.ptB);
      exprInc(dist, varDot(dist_grad, vars));
      exprInc(dist, -dist_grad.dot(toVectorXd(dofvals)));
    }
    if (itA != link2ind.end() || itB != link2ind.end()) {
      exprs.push_back(dist);
    }
  }
  LOG_DEBUG("%ld distance expressions\n", exprs.size());
}

void CollisionEvaluator::CollisionsToDistanceExpressions(const vector<Collision>& collisions, Configuration& rad0, Configuration& rad1, const Link2Int& link2ind,
    const VarVector& vars0, const VarVector& vars1, const DblVec& vals0, const DblVec& vals1,
    vector<AffExpr>& exprs) {
  vector<AffExpr> exprs0, exprs1;
  CollisionsToDistanceExpressions(collisions, rad0, link2ind, vars0, vals0, exprs0, false);
  CollisionsToDistanceExpressions(collisions, rad1, link2ind, vars1, vals1, exprs1, true);

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

void CollisionEvaluator::CollisionsToDistanceExpressions(const vector<Collision>& collisions, Configuration& rad00, Configuration& rad01, Configuration& rad10, Configuration& rad11, const Link2Int& link2ind0, const Link2Int& link2ind1,
    const VarVector& vars00, const VarVector& vars01, const VarVector& vars10, const VarVector& vars11,
    const DblVec& vals00, const DblVec& vals01, const DblVec& vals10, const DblVec& vals11,
    vector<AffExpr>& exprs) {

  vector<AffExpr> exprs00, exprs01, exprs10, exprs11;

  vector<OR::Vector> contact_points00, contact_points01, contact_points10, contact_points11;

  exprs00.reserve(collisions.size());
  rad00.SetDOFValues(vals00); // since we'll be calculating jacobians
  BOOST_FOREACH(const Collision& col, collisions) {
    AffExpr dist(col.distance);
    Link2Int::const_iterator itA = link2ind0.find(col.linkA);
    assert(itA != link2ind0.end());
    VectorXd dist_grad = toVector3d(col.normalB2A).transpose()*rad00.PositionJacobian(itA->second, col.pt00);
    exprInc(dist, varDot(dist_grad, vars00));
    exprInc(dist, -dist_grad.dot(toVectorXd(vals00)));
    exprs00.push_back(dist);
  }

  exprs01.reserve(collisions.size());
  rad01.SetDOFValues(vals01); // since we'll be calculating jacobians
  BOOST_FOREACH(const Collision& col, collisions) {
    AffExpr dist(col.distance);
    Link2Int::const_iterator itA = link2ind0.find(col.linkA);
    assert(itA != link2ind0.end());
    VectorXd dist_grad = toVector3d(col.normalB2A).transpose()*rad01.PositionJacobian(itA->second, col.pt01);
    exprInc(dist, varDot(dist_grad, vars01));
    exprInc(dist, -dist_grad.dot(toVectorXd(vals01)));
    exprs01.push_back(dist);
  }

  exprs10.reserve(collisions.size());
  rad10.SetDOFValues(vals10); // since we'll be calculating jacobians
  BOOST_FOREACH(const Collision& col, collisions) {
    AffExpr dist(col.distance);
    Link2Int::const_iterator itB = link2ind1.find(col.linkB);
    assert(itB != link2ind1.end());
    VectorXd dist_grad = -toVector3d(col.normalB2A).transpose()*rad10.PositionJacobian(itB->second, col.pt10);
    exprInc(dist, varDot(dist_grad, vars10));
    exprInc(dist, -dist_grad.dot(toVectorXd(vals10)));
    exprs10.push_back(dist);
  }

  exprs11.reserve(collisions.size());
  rad11.SetDOFValues(vals11); // since we'll be calculating jacobians
  BOOST_FOREACH(const Collision& col, collisions) {
    AffExpr dist(col.distance);
    Link2Int::const_iterator itB = link2ind1.find(col.linkB);
    assert(itB != link2ind1.end());
    VectorXd dist_grad = -toVector3d(col.normalB2A).transpose()*rad11.PositionJacobian(itB->second, col.pt11);
    exprInc(dist, varDot(dist_grad, vars11));
    exprInc(dist, -dist_grad.dot(toVectorXd(vals11)));
    exprs11.push_back(dist);
  }

  exprs.resize(exprs00.size());

  for (int i=0; i < exprs00.size(); ++i) {
    exprScale(exprs00[i], (1-collisions[i].timeA));
    exprScale(exprs01[i], collisions[i].timeA);
    exprScale(exprs10[i], (1-collisions[i].timeB));
    exprScale(exprs11[i], collisions[i].timeB);

    exprs[i] = AffExpr(0);
    exprInc(exprs[i], exprs00[i]);
    exprInc(exprs[i], exprs01[i]);
    exprInc(exprs[i], exprs10[i]);
    exprInc(exprs[i], exprs11[i]);
    cleanupAff(exprs[i]);
  }
  
}


inline size_t hash(const DblVec& x) {
  return boost::hash_range(x.begin(), x.end());
}

void CollisionEvaluator::GetCollisionsCached(const DblVec& x, vector<Collision>& collisions) {
  double key = hash(getDblVec(x, GetVars()));
  vector<Collision>* it = m_cache.get(key);
  if (it != NULL and false) {
    LOG_DEBUG("using cached collision check\n");
    collisions = *it;
  }
  else {
    LOG_DEBUG("not using cached collision check\n");
    CalcCollisions(x, collisions);
    m_cache.put(key, collisions);
  }
}

SingleTimestepCollisionEvaluator::SingleTimestepCollisionEvaluator(ConfigurationPtr rad, const VarVector& vars) :
  m_env(rad->GetEnv()),
  m_cc(CollisionChecker::GetOrCreate(*m_env)),
  m_rad(rad),
  m_vars(vars),
  m_link2ind(),
  m_links(),
  m_filterMask(-1) {
  vector<KinBody::LinkPtr> links;
  vector<int> inds;
  rad->GetAffectedLinks(m_links, true, inds);
  for (int i=0; i < m_links.size(); ++i) {
    m_link2ind[m_links[i].get()] = inds[i];
  }
  // TODO add argument
}


void SingleTimestepCollisionEvaluator::CalcCollisions(const DblVec& x, vector<Collision>& collisions) {
  DblVec dofvals = getDblVec(x, m_vars);
  m_rad->SetDOFValues(dofvals);
  m_cc->LinksVsAll(m_links, collisions, m_filterMask);
}

void SingleTimestepCollisionEvaluator::CalcDists(const DblVec& x, DblVec& dists) {
  vector<Collision> collisions;
  GetCollisionsCached(x, collisions);
  CollisionsToDistances(collisions, m_link2ind, dists);
}


void SingleTimestepCollisionEvaluator::CalcDistExpressions(const DblVec& x, vector<AffExpr>& exprs) {
  vector<Collision> collisions;
  GetCollisionsCached(x, collisions);
  DblVec dofvals = getDblVec(x, m_vars);
  CollisionsToDistanceExpressions(collisions, *m_rad, m_link2ind, m_vars, dofvals, exprs, false);
}

////////////////////////////////////////

CastCollisionEvaluator::CastCollisionEvaluator(ConfigurationPtr rad0, ConfigurationPtr rad1, const VarVector& vars0, const VarVector& vars1) :
  m_env(rad0->GetEnv()),
  m_cc(CollisionChecker::GetOrCreate(*m_env)),
  m_rad0(rad0),
  m_rad1(rad1),
  m_vars0(vars0),
  m_vars1(vars1),
  m_link2ind(),
  m_links() {
  vector<KinBody::LinkPtr> links;
  vector<int> inds;
  rad0->GetAffectedLinks(m_links, true, inds);
  for (int i=0; i < m_links.size(); ++i) {
    m_link2ind[m_links[i].get()] = inds[i];
  }
}

void CastCollisionEvaluator::CalcCollisions(const DblVec& x, vector<Collision>& collisions) {
  DblVec dofvals0 = getDblVec(x, m_vars0);
  DblVec dofvals1 = getDblVec(x, m_vars1);
  //m_rad0->SetDOFValues(dofvals0);
  //m_rad0->SetDOFValues(dofvals0);
  //m_rad0->SetDOFValues(dofvals0);
  //m_rad1->SetDOFValues(dofvals0);
  m_cc->CastVsAll(*m_rad0, *m_rad1, m_links, dofvals0, dofvals1, collisions);
}

void CastCollisionEvaluator::CalcDistExpressions(const DblVec& x, vector<AffExpr>& exprs) {
  vector<Collision> collisions;
  GetCollisionsCached(x, collisions);
  DblVec dofvals0 = getDblVec(x, m_vars0);
  DblVec dofvals1 = getDblVec(x, m_vars1);
  CollisionsToDistanceExpressions(collisions, *m_rad0, *m_rad1, m_link2ind, m_vars0, m_vars1, dofvals0, dofvals1, exprs);
}
void CastCollisionEvaluator::CalcDists(const DblVec& x, DblVec& dists) {
  vector<Collision> collisions;
  GetCollisionsCached(x, collisions);
  CollisionsToDistances(collisions, m_link2ind, dists);
}


//////////////////////////////////////////

CastSelfCollisionEvaluator::CastSelfCollisionEvaluator(
    ConfigurationPtr rad00, ConfigurationPtr rad01, ConfigurationPtr rad10, ConfigurationPtr rad11,
    const VarVector& vars00, const VarVector& vars01, const VarVector& vars10, const VarVector& vars11) :
  m_env(rad00->GetEnv()),
  m_cc(CollisionChecker::GetOrCreate(*m_env)),
  m_rad00(rad00),
  m_rad01(rad01),
  m_rad10(rad10),
  m_rad11(rad11),
  m_vars00(vars00),
  m_vars01(vars01),
  m_vars10(vars10),
  m_vars11(vars11),
  m_link2ind0(),
  m_links0(),
  m_link2ind1(),
  m_links1() {
  vector<KinBody::LinkPtr> links0;
  vector<int> inds0;
  rad00->GetAffectedLinks(m_links0, true, inds0);
  for (int i=0; i < m_links0.size(); ++i) {
    m_link2ind0[m_links0[i].get()] = inds0[i];
  }
  vector<KinBody::LinkPtr> links1;
  vector<int> inds1;
  rad10->GetAffectedLinks(m_links1, true, inds1);
  for (int i=0; i < m_links1.size(); ++i) {
    m_link2ind1[m_links1[i].get()] = inds1[i];
  }
}

void CastSelfCollisionEvaluator::CalcCollisions(const DblVec& x, vector<Collision>& collisions) {
  DblVec dofvals00 = getDblVec(x, m_vars00);
  DblVec dofvals01 = getDblVec(x, m_vars01);
  DblVec dofvals10 = getDblVec(x, m_vars10);
  DblVec dofvals11 = getDblVec(x, m_vars11);
  m_cc->CastVsCast(*m_rad00, *m_rad01, *m_rad10, *m_rad11, m_links0, m_links1, dofvals00, dofvals01, dofvals10, dofvals11, collisions);
}

void CastSelfCollisionEvaluator::CalcDistExpressions(const DblVec& x, vector<AffExpr>& exprs) {
  vector<Collision> collisions;
  GetCollisionsCached(x, collisions);
  DblVec dofvals00 = getDblVec(x, m_vars00);
  DblVec dofvals01 = getDblVec(x, m_vars01);
  DblVec dofvals10 = getDblVec(x, m_vars10);
  DblVec dofvals11 = getDblVec(x, m_vars11);
  CollisionsToDistanceExpressions(collisions, *m_rad00, *m_rad01, *m_rad10, *m_rad11, m_link2ind0, m_link2ind1, m_vars00, m_vars01, m_vars10, m_vars11, dofvals00, dofvals01, dofvals10, dofvals11, exprs);
}

void CastSelfCollisionsToDistances(const vector<Collision>& collisions, const Link2Int& m_link2ind0, const Link2Int& m_link2ind1,
    DblVec& dists) {
  dists.clear();
  dists.reserve(collisions.size());
  BOOST_FOREACH(const Collision& col, collisions) {
    Link2Int::const_iterator itA = m_link2ind0.find(col.linkA);
    Link2Int::const_iterator itB = m_link2ind1.find(col.linkB);
    if (itA != m_link2ind0.end() || itB != m_link2ind1.end()) {
      dists.push_back(col.distance);
    }
  }
}

void CastSelfCollisionEvaluator::CalcDists(const DblVec& x, DblVec& dists) {
  vector<Collision> collisions;
  GetCollisionsCached(x, collisions);
  CastSelfCollisionsToDistances(collisions, m_link2ind0, m_link2ind1, dists);
}

//////////////////////////////////////////

typedef OpenRAVE::RaveVector<float> RaveVectorf;

void PlotCollisions(const std::vector<Collision>& collisions, OR::EnvironmentBase& env, vector<OR::GraphHandlePtr>& handles, double safe_dist) {
  BOOST_FOREACH(const Collision& col, collisions) {
    RaveVectorf color;
    if (col.distance < 0) color = RaveVectorf(1,0,0,1);
    else if (col.distance < safe_dist) color = RaveVectorf(1,1,0,1);
    else color = RaveVectorf(0,1,0,1);
    if (col.cctype == CCType_Between) {
      handles.push_back(env.drawarrow(col.ptB, col.ptB1, .002, RaveVectorf(0,0,0,1)));
    }
    OR::Vector ptB = (col.cctype == CCType_Between)  ? ((1-col.time)* col.ptB +col.time*col.ptB1) : col.ptB;
    handles.push_back(env.drawarrow(col.ptA, ptB, .0025, color));
  }
}


CollisionCost::CollisionCost(double dist_pen, double coeff, ConfigurationPtr rad, const VarVector& vars) :
    Cost("collision"),
    m_calc(new SingleTimestepCollisionEvaluator(rad, vars)), m_dist_pen(dist_pen), m_coeff(coeff)
{}

CollisionCost::CollisionCost(double dist_pen, double coeff, ConfigurationPtr rad0, ConfigurationPtr rad1, const VarVector& vars0, const VarVector& vars1) :
    Cost("cast_collision"),
    m_calc(new CastCollisionEvaluator(rad0, rad1, vars0, vars1)), m_dist_pen(dist_pen), m_coeff(coeff)
{}
ConvexObjectivePtr CollisionCost::convex(const vector<double>& x) {
  ConvexObjectivePtr out(new ConvexObjective());
  vector<AffExpr> exprs;
  m_calc->CalcDistExpressions(x, exprs);
  for (int i=0; i < exprs.size(); ++i) {
    AffExpr viol = exprSub(AffExpr(m_dist_pen), exprs[i]);
    out->addHinge(viol, m_coeff);
  }
  return out;
}
double CollisionCost::value(const vector<double>& x, Model* model) {
  DblVec dists;
  m_calc->CalcDists(x, dists);
  double out = 0;
  for (int i=0; i < dists.size(); ++i) {
    out += pospart(m_dist_pen - dists[i]) * m_coeff;
  }
  return out;
}

void CollisionCost::Plot(const DblVec& x, OR::EnvironmentBase& env, std::vector<OR::GraphHandlePtr>& handles) {
  vector<Collision> collisions;
  m_calc->GetCollisionsCached(x, collisions);
  PlotCollisions(collisions, env, handles, m_dist_pen);
}

// ALMOST EXACTLY COPIED FROM CollisionCost

CollisionConstraint::CollisionConstraint(double dist_pen, double coeff, ConfigurationPtr rad, const VarVector& vars) :
    m_calc(new SingleTimestepCollisionEvaluator(rad, vars)), m_dist_pen(dist_pen), m_coeff(coeff)
{
  name_="collision";
}

CollisionConstraint::CollisionConstraint(double dist_pen, double coeff, ConfigurationPtr rad0, ConfigurationPtr rad1, const VarVector& vars0, const VarVector& vars1) :
    m_calc(new CastCollisionEvaluator(rad0, rad1, vars0, vars1)), m_dist_pen(dist_pen), m_coeff(coeff)
{
  name_="collision";
}

CollisionConstraint::CollisionConstraint(double dist_pen, double coeff, ConfigurationPtr rad00, ConfigurationPtr rad01, ConfigurationPtr rad10, ConfigurationPtr rad11, const VarVector& vars00, const VarVector& vars01, const VarVector& vars10, const VarVector& vars11) :
    m_calc(new CastSelfCollisionEvaluator(rad00, rad01, rad10, rad11, vars00, vars01, vars10, vars11)), m_dist_pen(dist_pen), m_coeff(coeff)
{
  name_="collision";
}

ConvexConstraintsPtr CollisionConstraint::convex(const vector<double>& x) {
  ConvexConstraintsPtr out(new ConvexConstraints());
  vector<AffExpr> exprs;
  m_calc->CalcDistExpressions(x, exprs);
  for (int i=0; i < exprs.size(); ++i) {
    AffExpr viol = exprSub(AffExpr(m_dist_pen), exprs[i]);
    out->addIneqCnt(exprMult(viol,m_coeff));
  }
  return out;
}
DblVec CollisionConstraint::value(const vector<double>& x, Model* model) {
  DblVec dists;
  m_calc->CalcDists(x, dists);
  DblVec out(dists.size());
  for (int i=0; i < dists.size(); ++i) {
    out[i] = pospart(m_dist_pen - dists[i]) * m_coeff;
  }
  return out;
}

}
