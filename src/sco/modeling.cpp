#include "modeling.hpp"
#include "utils/logging.hpp"
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <cstdio>
#include "expr_ops.hpp"
#include "sco_common.hpp"
#include "macros.h"
#include <iostream>
#include <sstream>

using namespace std;

namespace sco {

void ConvexObjective::addAffExpr(const AffExpr& affexpr) {
  objective_infos_.push_back(ObjectiveInfo(Objective_AffExpr, affexpr));
  //exprInc(quad_, affexpr);
}
void ConvexObjective::addQuadExpr(const QuadExpr& quadexpr) {
  objective_infos_.push_back(ObjectiveInfo(Objective_QuadExpr, quadexpr));
  //exprInc(quad_, quadexpr);
}
void ConvexObjective::addHinge(const AffExpr& affexpr, double coeff) {
  objective_infos_.push_back(ObjectiveInfo(Objective_Hinge, affexpr, coeff));//"hinge", 0, INFINITY
}
void ConvexObjective::addAbs(const AffExpr& affexpr, double coeff) {
  objective_infos_.push_back(ObjectiveInfo(Objective_Abs, affexpr, coeff));
  //Var neg = model_->addVar("neg", 0, INFINITY);
  //Var pos = model_->addVar("pos", 0, INFINITY);
  //vars_.push_back(neg);
  //vars_.push_back(pos);
  //AffExpr neg_plus_pos;
  //neg_plus_pos.coeffs = vector<double>(2, coeff);
  //neg_plus_pos.vars.push_back(neg);
  //neg_plus_pos.vars.push_back(pos);
  //exprInc(quad_, neg_plus_pos);
  //AffExpr affeq = affexpr;
  //affeq.vars.push_back(neg);
  //affeq.vars.push_back(pos);
  //affeq.coeffs.push_back(1);
  //affeq.coeffs.push_back(-1);
  //eqs_.push_back(affeq);
}
void ConvexObjective::addHinges(const AffExprVector& ev) {
  for (size_t i=0; i < ev.size(); ++i) addHinge(ev[i],1);
}
void ConvexObjective::addL1Norm(const AffExprVector& ev) {
  for (size_t i=0; i < ev.size(); ++i) addAbs(ev[i],1);
}
void ConvexObjective::addL2Norm(const AffExprVector& ev) {
  for (size_t i=0; i < ev.size(); ++i) addQuadExpr(exprSquare(ev[i]));
}
void ConvexObjective::addMax(const AffExprVector& ev) { // this is probably buggy
  objective_infos_.push_back(ObjectiveInfo(Objective_Max, ev));

  //Var m = model_->addVar("max", -INFINITY, INFINITY);
  //for (size_t i=0; i < ev.size(); ++i) {
  //  ineqs_.push_back(ev[i]);
  //  exprDec(ineqs_.back(), m);
  //}
}

void ConvexObjective::addToModelAndObjective(Model* model, AffExpr& objective, bool permissive) {
  assert (model != NULL);
  assert (model2vars_.count(model) == 0);
  assert (model2cnts_.count(model) == 0);
  assert (model2quad_.count(model) == 0);

  vector<AffExpr> eqs;
  vector<AffExpr> ineqs;
  vector<Var> vars;

  QuadExpr quad;

  BOOST_FOREACH(const ObjectiveInfo& info, objective_infos_) {
    switch (info.type_) {
      case Objective_AffExpr: {
        exprInc(objective, info.affexpr_);
        exprInc(quad, info.affexpr_);
        break;
      }
      case Objective_QuadExpr: {
        if (!permissive) {
          throw std::runtime_error("Error: attempt to add quadratic expression to affine objective");
        }
        break;
      }
      case Objective_Hinge: {
        Var hinge = model->addVar("hinge", 0, INFINITY);
        vars.push_back(hinge);
        ineqs.push_back(info.affexpr_);
        exprDec(ineqs.back(), hinge);
        AffExpr hinge_cost = exprMult(AffExpr(hinge), info.coeff_);
        exprInc(objective, hinge_cost);
        exprInc(quad, hinge_cost);
        break;
      }
      case Objective_Abs: {
        Var neg = model->addVar("neg", 0, INFINITY);
        Var pos = model->addVar("pos", 0, INFINITY);
        vars.push_back(neg);
        vars.push_back(pos);
        AffExpr neg_plus_pos;
        neg_plus_pos.coeffs = vector<double>(2, info.coeff_);
        neg_plus_pos.vars.push_back(neg);
        neg_plus_pos.vars.push_back(pos);
        exprInc(objective, neg_plus_pos);
        exprInc(quad, neg_plus_pos);
        AffExpr affeq = info.affexpr_;
        affeq.vars.push_back(neg);
        affeq.vars.push_back(pos);
        affeq.coeffs.push_back(1);
        affeq.coeffs.push_back(-1);
        eqs.push_back(affeq);
        break;
      }
      case Objective_Max: {
        Var m = model->addVar("max", -INFINITY, INFINITY);
        vars.push_back(m);
        for (size_t i=0; i < info.ev_.size(); ++i) {
          ineqs.push_back(info.ev_[i]);
          exprDec(ineqs.back(), m);
        }
        break;
      }
      default: {
        throw std::runtime_error("unsupported auxiliary type");
        break;
      }
    }
  }

  model->update();

  vector<Cnt> cnts;
  cnts.reserve(eqs.size() + ineqs.size());
  BOOST_FOREACH(const AffExpr& aff, eqs) {
    cnts.push_back(model->addEqCnt(aff, ""));
  }
  BOOST_FOREACH(const AffExpr& aff, ineqs) {
    cnts.push_back(model->addIneqCnt(aff, ""));
  }

  model->update();

  model2vars_.insert(ModelVarsPair(model, vars));
  model2cnts_.insert(ModelCntsPair(model, cnts));
  model2quad_.insert(ModelQuadPair(model, quad));
}

void ConvexObjective::addToModelAndObjective(Model* model, QuadExpr& objective) {
  assert (model != NULL);

  addToModelAndObjective(model, objective.affexpr, true);

  QuadExpr& quad = model2quad_.at(model);

  BOOST_FOREACH(const ObjectiveInfo& info, objective_infos_) {
    if (info.type_ == Objective_QuadExpr) {
      exprInc(objective, info.quadexpr_);
      exprInc(quad, info.quadexpr_);
    }
  }
}

//void ConvexObjective::addConstraintsToModel() {
//  cnts_.reserve(eqs_.size() + ineqs_.size());
//  BOOST_FOREACH(const AffExpr& aff, eqs_) {
//    cnts_.push_back(model_->addEqCnt(aff, ""));
//  }
//  BOOST_FOREACH(const AffExpr& aff, ineqs_) {
//    cnts_.push_back(model_->addIneqCnt(aff, ""));
//  }
//}

void ConvexObjective::removeFromModel(Model* model) {
  assert (model != NULL);
  if (model2vars_.count(model) > 0) {
    model->removeVars(model2vars_.at(model));
    model2vars_.erase(model);
  }

  if (model2cnts_.count(model) > 0) {
    model->removeCnts(model2cnts_.at(model));
    model2cnts_.erase(model);
  }

  if (model2quad_.count(model) > 0) {
    model2quad_.erase(model);
  }

}

void ConvexObjective::removeFromModels() {
  vector<Model*> models;
  BOOST_FOREACH(const Model2Vars::value_type& pair, model2vars_) {
    models.push_back(pair.first);    
  }
  BOOST_FOREACH(Model* model, models) {
    removeFromModel(model);
  }
}

ConvexObjective::~ConvexObjective() {
  removeFromModels();
}

void ConvexConstraints::addEqCnt(const AffExpr& aff) {
  eqs_.push_back(aff);
}

void ConvexConstraints::addIneqCnt(const AffExpr& aff) {
  ineqs_.push_back(aff);
}

void ConvexConstraints::addToModel(Model* model) {
  assert (model != NULL);
  assert (model2cnts_.count(model) == 0);

  vector<Cnt> cnts;
  cnts.reserve(eqs_.size() + ineqs_.size());
  BOOST_FOREACH(const AffExpr& aff, eqs_) {
    cnts.push_back(model->addEqCnt(aff, ""));
  }
  BOOST_FOREACH(const AffExpr& aff, ineqs_) {
    cnts.push_back(model->addIneqCnt(aff, ""));
  }
  model2cnts_.insert(ModelCntsPair(model, cnts));
}

void ConvexConstraints::removeFromModel(Model* model) {
  assert (model != NULL);
  if (model2cnts_.count(model) > 0) {
    model->removeCnts(model2cnts_.at(model));
    model2cnts_.erase(model);
  }
}

vector<double> ConvexConstraints::violations(const vector<double>& x, Model* model) {
  DblVec out;
  out.reserve(eqs_.size() + ineqs_.size());
  BOOST_FOREACH(const AffExpr& aff, eqs_) out.push_back(fabs(aff.value(x.data())));
  BOOST_FOREACH(const AffExpr& aff, ineqs_) out.push_back(pospart(aff.value(x.data())));
  return out;
}
double ConvexConstraints::violation(const vector<double>& x, Model* model) {
  return vecSum(violations(x, model));
}

void ConvexConstraints::removeFromModels() {
  vector<Model*> models;
  BOOST_FOREACH(const Model2Cnts::value_type& pair, model2cnts_) {
    models.push_back(pair.first);
  }
  BOOST_FOREACH(Model* model, models) {
    removeFromModel(model);
  }
}

ConvexConstraints::~ConvexConstraints() {
  removeFromModels();
}

double ConvexObjective::value(const vector<double>& x, Model* model)  {
  assert (model != NULL);
  assert (model2quad_.count(model) > 0);
  return model2quad_.at(model).value(x);
}


vector<double> Constraint::violations(const DblVec& x, Model* model) {
  DblVec val = value(x, model);
  DblVec out(val.size());

  if (type() == EQ) {
    for (size_t i=0; i < val.size(); ++i) out[i] = fabs(val[i]);
  }
  else { // type() == INEQ
    for (size_t i=0; i < val.size(); ++i) out[i] = pospart(val[i]);
  }

  return out;
}

double Constraint::violation(const DblVec& x, Model* model) {
  return vecSum(violations(x, model));
}

OptProb::OptProb() : model_(createModel()) {}

VarVector OptProb::createVariables(const vector<string>& var_names) {
  return createVariables(var_names, DblVec(var_names.size(), -INFINITY), DblVec(var_names.size(), INFINITY));
}
VarVector OptProb::createVariables(const vector<string>& var_names, const DblVec& lb, const DblVec& ub) {
  size_t n_add = var_names.size(), n_cur = vars_.size();
  assert(lb.size() == n_add);
  assert(ub.size() == n_add);
  vars_.reserve(n_cur + n_add);
  lower_bounds_.reserve(n_cur + n_add);
  upper_bounds_.reserve(n_cur + n_add);
  for (size_t i=0; i < var_names.size(); ++i) {
    vars_.push_back(model_->addVar(var_names[i], lb[i], ub[i]));
    lower_bounds_.push_back(lb[i]);
    upper_bounds_.push_back(ub[i]);
  }
  model_->update();
  return VarVector(vars_.end()-n_add, vars_.end());
}
void OptProb::setLowerBounds(const vector<double>& lb) {
  assert(lb.size() == vars_.size());
  lower_bounds_ = lb;
}
void OptProb::setUpperBounds(const vector<double>& ub) {
  assert(ub.size() == vars_.size());
  upper_bounds_ = ub;
}
void OptProb::setLowerBounds(const vector<double>& lb, const vector<Var>& vars) {
  setVec(lower_bounds_, vars, lb);
}
void OptProb::setUpperBounds(const vector<double>& ub, const vector<Var>& vars) {
  setVec(upper_bounds_, vars, ub);
}

void OptProb::addCost(CostPtr cost) {
  costs_.push_back(cost);
}
void OptProb::addConstraint(ConstraintPtr cnt) {
  if (cnt->type() == EQ) addEqConstraint(cnt);
  else addIneqConstraint(cnt);
}
void OptProb::addEqConstraint(ConstraintPtr cnt) {
  assert (cnt->type() == EQ);
  eqcnts_.push_back(cnt);
}
void OptProb::addIneqConstraint(ConstraintPtr cnt) {
  assert (cnt->type() == INEQ);
  ineqcnts_.push_back(cnt);
}
vector<ConstraintPtr> OptProb::getConstraints() const {
  vector<ConstraintPtr> out;
  out.reserve(eqcnts_.size() + ineqcnts_.size());
  out.insert(out.end(), eqcnts_.begin(), eqcnts_.end());
  out.insert(out.end(), ineqcnts_.begin(), ineqcnts_.end());
  return out;
}
void OptProb::addLinearConstraint(const AffExpr& expr, ConstraintType type) {
  if (type == EQ) model_->addEqCnt(expr, "");
  else model_->addIneqCnt(expr, "");
}

vector<double> OptProb::getCentralFeasiblePoint(const vector<double>& x) {
  assert(x.size() == lower_bounds_.size());
  DblVec center(x.size());
  for (int i=0; i < x.size(); ++i) center[i] = (lower_bounds_[i] + upper_bounds_[i])/2;
  return getClosestFeasiblePoint(center);
}
vector<double> OptProb::getClosestFeasiblePoint(const vector<double>& x) {
  LOG_DEBUG("getClosestFeasiblePoint");
  assert(vars_.size() == x.size());
  QuadExpr obj;
  for (int i=0; i < x.size(); ++i) {
    exprInc(obj, exprSquare(exprSub(AffExpr(vars_[i]),x[i])));
  }
  model_->setVarBounds(vars_, lower_bounds_, upper_bounds_);
  model_->setObjective(obj);
  CvxOptStatus status = model_->optimize();
  if(status != CVX_SOLVED) {
    model_->writeToFile("/tmp/fail2.lp");
    PRINT_AND_THROW("couldn't find a feasible point. there's probably a problem with variable bounds (e.g. joint limits). wrote to /tmp/fail2.lp");
  }
  return model_->getVarValues(vars_);
}



}
