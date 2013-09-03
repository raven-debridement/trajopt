#pragma once

/*
 * Model a non-convex optimization problem by defining Cost and Constraint objects
 * which know how to generate a convex approximation
 *
 *
 */

#include <vector>
#include <boost/shared_ptr.hpp>
#include <map>
#include "sco/sco_fwd.hpp"
#include "sco/solver_interface.hpp"

namespace sco {

using std::vector;
using std::pair;
using std::map;

enum ObjectiveType {
  Objective_AffExpr,
  Objective_QuadExpr,
  Objective_Hinge, 
  Objective_Abs,
  Objective_Max
};

typedef map<Model*, vector<Var> > Model2Vars;
typedef pair<Model*, vector<Var> > ModelVarsPair;
typedef map<Model*, vector<Cnt> > Model2Cnts;
typedef pair<Model*, vector<Cnt> > ModelCntsPair;
typedef map<Model*, QuadExpr > Model2Quad;
typedef pair<Model*, QuadExpr > ModelQuadPair;

/**
Stores convex terms in a objective
For non-quadratic terms like hinge(x) and abs(x), it needs to add auxilliary variables and linear constraints to the model
Note: When this object is deleted, the constraints and variables it added to the model are removed
 */
class ConvexObjective {
public:
  
  struct ObjectiveInfo {
    ObjectiveInfo(ObjectiveType type, const QuadExpr& quadexpr, double coeff=1) : type_(type), quadexpr_(quadexpr), coeff_(coeff) {}
    ObjectiveInfo(ObjectiveType type, const AffExpr& affexpr, double coeff=1) : type_(type), affexpr_(affexpr), coeff_(coeff) {}
    ObjectiveInfo(ObjectiveType type, const AffExprVector& ev, double coeff=1) : type_(type), ev_(ev), coeff_(coeff) {}
    ObjectiveType type_;
    AffExpr affexpr_;
    QuadExpr quadexpr_;
    AffExprVector ev_;
    double coeff_;
  };


  ConvexObjective() {}
  void addAffExpr(const AffExpr&);
  void addQuadExpr(const QuadExpr&);
  void addHinge(const AffExpr&, double coeff);
  void addAbs(const AffExpr&, double coeff);
  void addHinges(const AffExprVector&);
  void addL1Norm(const AffExprVector&);
  void addL2Norm(const AffExprVector&);
  void addMax(const AffExprVector&);
  
  void addToModelAndObjective(Model* model, AffExpr& objective, bool permissive=false);
  void addToModelAndObjective(Model* model, QuadExpr& objective);
  void removeFromModel(Model* model);
  void removeFromModels();

  double value(const vector<double>& x, Model* model);
  
  ~ConvexObjective();

  Model2Vars model2vars_;
  Model2Cnts model2cnts_;
  Model2Quad model2quad_;
  vector<ObjectiveInfo> objective_infos_;
private:
  ConvexObjective(ConvexObjective&)  {}
};

/**
Stores convex inequality constraints and affine equality constraints.
Actually only affine inequality constraints are currently implemented.
*/
class ConvexConstraints {
public:
  ConvexConstraints() {}
  /** Expression that should == 0 */
  void addEqCnt(const AffExpr&);
  /** Expression that should <= 0 */
  void addIneqCnt(const AffExpr&);
  
  void addToModel(Model* model);
  void removeFromModel(Model* model);
  void removeFromModels();

  vector<double> violations(const vector<double>& x, Model* model);
  double violation(const vector<double>& x, Model* model);

  ~ConvexConstraints();
  vector<AffExpr> eqs_;
  vector<AffExpr> ineqs_;
  Model2Cnts model2cnts_;
private:
   ConvexConstraints(ConvexConstraints&) {}
};

/**
Non-convex cost function, which knows how to calculate its convex approximation (convexify() method)
*/
class Cost {
public:
  /** Evaluate at solution vector x*/
  virtual double value(const vector<double>&, Model*) = 0;
  /** Convexify at solution vector x*/
  virtual ConvexObjectivePtr convex(const vector<double>& x) = 0;

  string name() {return name_;}
  void setName(const string& name) {name_=name;}
  Cost() : name_("unnamed") {}
  Cost(const string& name) : name_(name) {}
  virtual ~Cost() {}
protected:
  string name_;
};

/**
Non-convex vector-valued constraint function, which knows how to calculate its convex approximation
*/
class Constraint {
public:

  /** inequality vs equality */
  virtual ConstraintType type() = 0;
  /** Evaluate at solution vector x*/  
  virtual vector<double> value(const vector<double>& x, Model*) = 0;
  /** Convexify at solution vector x*/  
  virtual ConvexConstraintsPtr convex(const vector<double>& x) = 0;
  /** Calculate constraint violations (positive part for inequality constraint, absolute value for inequality constraint)*/
  vector<double> violations(const vector<double>& x, Model*);
  /** Sum of violations */
  double violation(const vector<double>& x, Model*);

  string name() {return name_;}
  void setName(const string& name) {name_=name;}
  Constraint() : name_("unnamed") {}
  Constraint(const string& name) : name_(name) {}
  virtual ~Constraint() {}

protected:
  string name_;
};

class EqConstraint : public Constraint{
public:
  ConstraintType type() {return EQ;}
};

class IneqConstraint : public Constraint {
public:
  ConstraintType type() {return INEQ;}
};

/**
Non-convex optimization problem
*/
class OptProb {
public:
  OptProb();
  /** create variables with bounds [-INFINITY, INFINITY]  */
  VarVector createVariables(const vector<string>& names);
  /** create variables with bounds [lb[i], ub[i] */
  VarVector createVariables(const vector<string>& names, const vector<double>& lb, const vector<double>& ub);
  /** set the lower bounds of all the variables */
  void setLowerBounds(const vector<double>& lb);
  /** set the upper bounds of all the variables */
  void setUpperBounds(const vector<double>& ub);
  /** set lower bounds of some of the variables */
  void setLowerBounds(const vector<double>& lb, const vector<Var>& vars);
  /** set upper bounds of some of the variables */
  void setUpperBounds(const vector<double>& ub, const vector<Var>& vars);
  /** Note: in the current implementation, this function just adds the constraint to the
   * model. So if you're not careful, you might end up with an infeasible problem. */
  void addLinearConstraint(const AffExpr&, ConstraintType type);
  /** Add nonlinear cost function */
  void addCost(CostPtr);
  /** Add nonlinear constraint function */
  void addConstraint(ConstraintPtr);
  void addEqConstraint(ConstraintPtr);
  void addIneqConstraint(ConstraintPtr);
  virtual ~OptProb() {}
  /** Find closest point to solution vector x that satisfies linear inequality constraints */
  vector<double> getCentralFeasiblePoint(const vector<double>& x);
  vector<double> getClosestFeasiblePoint(const vector<double>& x);

  vector<ConstraintPtr> getConstraints() const;
  vector<CostPtr>& getCosts() {return costs_;}
  vector<ConstraintPtr>& getIneqConstraints() {return ineqcnts_;}
  vector<ConstraintPtr>& getEqConstraints() {return eqcnts_;}
  DblVec& getLowerBounds() {return lower_bounds_;}
  DblVec& getUpperBounds() {return upper_bounds_;}
  ModelPtr getModel() {return model_;}
  vector<Var>& getVars() {return vars_;}
  int getNumCosts() {return costs_.size();}
  int getNumConstraints() {return eqcnts_.size() + ineqcnts_.size();}
  int getNumVars() {return vars_.size();}

protected:
  ModelPtr model_;
  vector<Var> vars_;
  vector<double> lower_bounds_;
  vector<double> upper_bounds_;
  vector<CostPtr> costs_;
  vector<ConstraintPtr> eqcnts_;
  vector<ConstraintPtr> ineqcnts_;

  OptProb(OptProb&);
};

template <typename VecType>
inline void setVec(DblVec& x, const VarVector& vars, const VecType& vals) {
  assert(vars.size() == vals.size());
  for (int i = 0; i < vars.size(); ++i) {
    x[vars[i].var_rep->index] = vals[i];
  }
}
template <typename OutVecType>
inline OutVecType getVec1(const vector<double>& x, const VarVector& vars) {
  OutVecType out(vars.size());
  for (unsigned i=0; i < vars.size(); ++i) out[i] = x[vars[i].var_rep->index];
  return out;
}


}
