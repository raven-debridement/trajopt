#pragma once
#include "trajopt/common.hpp"
#include "trajopt/collision_checker.hpp"
#include "utils/default_map.hpp"
#include "sco/modeling.hpp"
#include "sco/sco_fwd.hpp"
#include "cache.hxx"


namespace trajopt {

typedef std::map<const OR::KinBody::Link*, int> Link2Int;
typedef std::vector< std::pair<string, string> > NamePairs;



typedef DefaultMap<string, double> Str2Dbl;

struct TRAJOPT_API CollisionEvaluator {
  virtual void CalcDistExpressions(const DblVec& x, vector<AffExpr>& exprs, NamePairs& bodyNames) = 0;
  virtual void CalcDists(const DblVec& x, DblVec& exprs, NamePairs& bodyNames) = 0;
  virtual void CalcCollisions(const DblVec& x, vector<Collision>& collisions) = 0;
  virtual void GetCollisionsCached(const DblVec& x, vector<Collision>&);
  virtual void CollisionsToDistances(const vector<Collision>& collisions, const Link2Int& m_link2ind,
    DblVec& dists, NamePairs& bodyNames);
  virtual void CollisionsToDistanceExpressions(const vector<Collision>& collisions, Configuration& rad,
    const Link2Int& link2ind, const VarVector& vars, const DblVec& dofvals, vector<AffExpr>& exprs, bool isTimestep1, NamePairs& bodyNames);
  virtual void CollisionsToDistanceExpressions(const vector<Collision>& collisions, Configuration& rad, const Link2Int& link2ind,
    const VarVector& vars0, const VarVector& vars1, const DblVec& vals0, const DblVec& vals1,
    vector<AffExpr>& exprs, NamePairs& bodyNames);
  virtual ~CollisionEvaluator() {}
  virtual VarVector GetVars()=0;

  Cache<size_t, vector<Collision>, 3> m_cache;
};
typedef boost::shared_ptr<CollisionEvaluator> CollisionEvaluatorPtr;

struct TRAJOPT_API SingleTimestepCollisionEvaluator : public CollisionEvaluator {
public:
  SingleTimestepCollisionEvaluator(ConfigurationPtr rad, const VarVector& vars);
  /**
  @brief linearize all contact distances in terms of robot dofs
  
  Do a collision check between robot and environment.
  For each contact generated, return a linearization of the signed distance function
  */
  void CalcDistExpressions(const DblVec& x, vector<AffExpr>& exprs, NamePairs& bodyNames); 
  /**
   * Same as CalcDistExpressions, but just the distances--not the expressions
   */
  void CalcDists(const DblVec& x, DblVec& exprs, NamePairs& bodyNames); 
  void CalcCollisions(const DblVec& x, vector<Collision>& collisions);
  VarVector GetVars() {return m_vars;}

  OR::EnvironmentBasePtr m_env;
  CollisionCheckerPtr m_cc;
  ConfigurationPtr m_rad;
  VarVector m_vars;
  Link2Int m_link2ind;
  vector<OR::KinBody::LinkPtr> m_links;
  short m_filterMask;
};

struct CastCollisionEvaluator : public CollisionEvaluator {
public:
  CastCollisionEvaluator(ConfigurationPtr rad, const VarVector& vars0, const VarVector& vars1);
  void CalcDistExpressions(const DblVec& x, vector<AffExpr>& exprs, NamePairs& bodyNames);
  void CalcDists(const DblVec& x, DblVec& exprs, NamePairs& bodyNames);
  void CalcCollisions(const DblVec& x, vector<Collision>& collisions);
  VarVector GetVars() {return concat(m_vars0, m_vars1);}
  

  // parameters:
  OR::EnvironmentBasePtr m_env;
  CollisionCheckerPtr m_cc;
  ConfigurationPtr m_rad;
  VarVector m_vars0;
  VarVector m_vars1;
  typedef std::map<const OR::KinBody::Link*, int> Link2Int;
  Link2Int m_link2ind;
  vector<OR::KinBody::LinkPtr> m_links;
  short m_filterMask;

};


class TRAJOPT_API CollisionCost : public Cost, public Plotter {
public:
  /* constructor for single timestep */
  CollisionCost(double dist_pen, double coeff, ConfigurationPtr rad, const VarVector& vars);
  /* constructor for cast cost */
  CollisionCost(double dist_pen, double coeff, ConfigurationPtr rad, const VarVector& vars0, const VarVector& vars1);
  virtual ConvexObjectivePtr convex(const vector<double>& x, Model* model);
  virtual double value(const vector<double>&);
  void Plot(const DblVec& x, OR::EnvironmentBase& env, std::vector<OR::GraphHandlePtr>& handles);
protected:
  CollisionEvaluatorPtr m_calc;
  double m_dist_pen;
  double m_coeff;
};
class TRAJOPT_API CollisionConstraint : public IneqConstraint {
public:
  /* constructor for single timestep */
  CollisionConstraint(double dist_pen, double coeff, ConfigurationPtr rad, const VarVector& vars);
  /* constructor for cast cost */
  CollisionConstraint(double dist_pen, double coeff, ConfigurationPtr rad, const VarVector& vars0, const VarVector& vars1);
  virtual ConvexConstraintsPtr convex(const vector<double>& x, Model* model);
  virtual DblVec value(const vector<double>&);
  void Plot(const DblVec& x, OR::EnvironmentBase& env, std::vector<OR::GraphHandlePtr>& handles);
protected:
  CollisionEvaluatorPtr m_calc;
  double m_dist_pen;
  double m_coeff;
};

class TRAJOPT_API CollisionTaggedCost : public Cost {
public:
  CollisionTaggedCost(const Str2Dbl& tag2dist_pen, const Str2Dbl& tag2coeff, ConfigurationPtr rad, const VarVector& vars);
  CollisionTaggedCost(const Str2Dbl& tag2dist_pen, const Str2Dbl& tag2coeff, ConfigurationPtr rad, const VarVector& vars0, const VarVector& vars1);
  virtual ConvexObjectivePtr convex(const vector<double>& x, Model* model);
  virtual double value(const vector<double>&);
private:
  CollisionEvaluatorPtr m_calc;
  Str2Dbl m_tag2dist_pen;
  Str2Dbl m_tag2coeff;
};

class TRAJOPT_API CollisionTaggedConstraint : public IneqConstraint {
public:
  /* constructor for single timestep */
  CollisionTaggedConstraint(const Str2Dbl& tag2dist_pen, const Str2Dbl& tag2coeff, ConfigurationPtr rad, const VarVector& vars);
  /* constructor for cast cost */
  CollisionTaggedConstraint(const Str2Dbl& tag2dist_pen, const Str2Dbl& tag2coeff, ConfigurationPtr rad, const VarVector& vars0, const VarVector& vars1);
  virtual ConvexConstraintsPtr convex(const vector<double>& x, Model* model);
  virtual DblVec value(const vector<double>&);
private:
  CollisionEvaluatorPtr m_calc;
  Str2Dbl m_tag2dist_pen;
  Str2Dbl m_tag2coeff;

};

}
