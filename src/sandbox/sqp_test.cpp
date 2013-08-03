#include "sco/expr_ops.hpp"
#include "sco/modeling_utils.hpp"
#include "sco/modeling.hpp"
#include "osgviewer/osgviewer.hpp"
#include "trajopt/collision_checker.hpp"
#include "trajopt/collision_terms.hpp"
#include "trajopt/common.hpp"
#include "trajopt/plot_callback.hpp"
#include "trajopt/problem_description.hpp"
#include "trajopt/rave_utils.hpp"
#include "trajopt/trajectory_costs.hpp"
#include "utils/clock.hpp"
#include "utils/config.hpp"
#include "utils/eigen_conversions.hpp"
#include "utils/stl_to_string.hpp"
#include <boost/assign.hpp>
#include <boost/foreach.hpp>
#include <ctime>

using namespace std;
using namespace trajopt;
using namespace Eigen;

typedef Matrix<double, 1, 1> Vector1d;

class CostError : public ScalarOfVector {
public:
  CostError() {}
  double operator()(const VectorXd& a) const {
    return (a(1) - 1) * (a(1) - 1);
  }
};

class CntError1 : public VectorOfVector {
public:
  CntError1() {}
  VectorXd operator()(const VectorXd& a) const {
    Vector1d ans; ans << a(0)*a(0);//a(0)*a(0) + 1 - a(1);
    return ans;
  }
};

class CntError2 : public VectorOfVector {
public:
  CntError2() {}
  VectorXd operator()(const VectorXd& a) const {
    Vector1d ans; ans << a(1)*a(1)*a(1);//a(0) - 1 - a(2);
    return ans;
  }
};


void vecDiff(const DblVec& a, const DblVec& b, DblVec& c) {
  c.resize(min(a.size(), b.size()));
  for (int i = 0; i < min(a.size(), b.size()); ++i) {
    c[i] = a[i] - b[i];
  }
}

void vecScale(const DblVec& a, double scale, DblVec& c) {
  c.resize(a.size());
  for (int i = 0; i < a.size(); ++i) {
    c[i] = a[i] * scale;
  }
}


//class LineSearchSQP : public BasicTrustRegionSQP {
//public:
//	LineSearchSQP() : BasicTrustRegionSQP() {}
//	LineSearchSQP(OptProbPtr prob) : BasicTrustRegionSQP(prob) {}
//
//  bool hasViolation(const DblVec& cnt_viols) {
//    return cnt_viols.size() > 0 && vecMax(cnt_viols) < cnt_tolerance_;
//  }
//
//	OptStatus optimize() {
//
//    vector<string> cost_names = getCostNames(prob_->getCosts());
//    vector<ConstraintPtr> constraints = prob_->getConstraints();
//    vector<string> cnt_names = getCntNames(constraints);
//
//    DblVec& x_ = results_.x; // just so I don't have to rewrite code
//    if (x_.size() == 0) PRINT_AND_THROW("you forgot to initialize!");
//    if (!prob_) PRINT_AND_THROW("you forgot to set the optimization problem");    
//    
//    x_ = prob_->getClosestFeasiblePoint(x_);
//
//    // initialize LP model
//    ModelPtr lp_model = model_->cloneModel();
//    
//    assert(x_.size() == prob_->getVars().size());
//    assert(prob_->getCosts().size() > 0 || constraints.size() > 0);
//
//    OptStatus retval = INVALID;
//
//    for (int merit_increases=0; merit_increases < max_merit_coeff_increases_; ++merit_increases) { /* merit adjustment loop */
//      ++results_.n_merit_increases;
//      for (int iter=1; ; ++iter) { /* sqp loop */
//        callCallbacks(x_);
//
//        LOG_DEBUG("current iterate: %s", CSTR(x_));
//        LOG_INFO("merit increases iteration: %i; sqp iteration %i", merit_increases, iter);
//
//        // speedup: if you just evaluated the cost when doing the line search, use that
//        if (results_.cost_vals.empty() && results_.cnt_viols.empty()) { //only happens on the first iteration
//          results_.cnt_viols = evaluateConstraintViols(constraints, x_);
//          results_.cost_vals = evaluateCosts(prob_->getCosts(), x_);
//          assert(results_.n_func_evals == 0);
//          ++results_.n_func_evals;
//        }
//
//        // DblVec new_cnt_viols = evaluateConstraintViols(constraints, x_);
//        // DblVec new_cost_vals = evaluateCosts(prob_->getCosts(), x_);
//        // cout << "costs" << endl;
//        // for (int i=0; i < new_cnt_viols.size(); ++i) {
//        //   cout << cnt_names[i] << " " << new_cnt_viols[i] - results_.cnt_viols[i] << endl;
//        // }
//        // for (int i=0; i < new_cost_vals.size(); ++i) {
//        //   cout << cost_names[i] << " " << new_cost_vals[i] - results_.cost_vals[i] << endl;
//        // }
//
//
//        vector<ConvexObjectivePtr> cost_models = convexifyCosts(prob_->getCosts(),x_, model_.get());
//        vector<ConvexConstraintsPtr> cnt_models = convexifyConstraints(constraints, x_, model_.get());
//        vector<ConvexObjectivePtr> cnt_cost_models = cntsToCosts(cnt_models, merit_error_coeff_, model_.get());
//        model_->update();
//        BOOST_FOREACH(ConvexObjectivePtr& cost, cost_models)cost->addConstraintsToModel();
//        BOOST_FOREACH(ConvexObjectivePtr& cost, cnt_cost_models)cost->addConstraintsToModel();
//        model_->update();
//        QuadExpr objective;
//        BOOST_FOREACH(ConvexObjectivePtr& co, cost_models)exprInc(objective, co->quad_);
//        BOOST_FOREACH(ConvexObjectivePtr& co, cnt_cost_models){
//          exprInc(objective, co->quad_);
//        }
//  //    objective = cleanupExpr(objective);
//        model_->setObjective(objective);
//
//  //    if (logging::filter() >= IPI_LEVEL_DEBUG) {
//  //      DblVec model_cost_vals;
//  //      BOOST_FOREACH(ConvexObjectivePtr& cost, cost_models) {
//  //        model_cost_vals.push_back(cost->value(x));
//  //      }
//  //      LOG_DEBUG("model costs %s should equalcosts  %s", printer(model_cost_vals), printer(cost_vals));
//  //    }
//      // ep1 = 0.1
//      double min_cnt_improve_ratio_ = 0.1; // ep2 = 0.1
//      // tau = 0.5
//      // eta = 1e-4
//      // eta_1 = 0.25
//      // eta_2 = 0.75
//      // trust_region_min = 1e-3
//      // trust_region_max = 1e3
//      // trust_region_1 = 1
//      // merit_error_1 = 1
//      //
//
//
//      // archive the model variable values before optimization for later use
//      DblVec old_model_var_vals = model_->getVarValues(model_->getVars());
//
//
//      // 1. solve the problem constructed above _without_ a trust region size
//      CvxOptStatus status = model_->optimize(); 
//      ++results.n_qp_solves;
//      DblVec model_var_vals = model_->getVarValues(model_->getVars());
//      DblVec new_x(model_var_vals.begin(), model_var_vals.begin() + x_.size());
//      // step = new_x - x_
//      DblVec step; vecDiff(new_x, x_, step);
//
//      //  check if result is close to the previous x and if all cnt viols are satisfied for the old x
//      //    (stop in that case)
//      if (vecAbsSum(step) < OPT_EPS) {
//        if (!hasViolation(results_.cnt_viols)) {
//          if (results_.cnt_viols.size() > 0) LOG_INFO("woo-hoo! constraints are satisfied (to tolerance %.2e)", cnt_tolerance_);
//          retval = OPT_CONVERGED;
//          goto cleanup;
//        }
//      }
//
//
//      double new_merit_error_coeff = merit_error_coeff_;
//
//      // 2. if model cnt viols = 0 and old_model_merit - new_model_merit >= min_approx_improve_ratio * merit_error_coeff * old_model_viols
//      //    go to 6.
//      DblVec old_model_cost_vals = evaluateModelCosts(cost_models, old_model_var_vals);
//      DblVec old_model_cnt_viols = evaluateModelCntViols(cnt_models, old_model_var_vals);
//      DblVec new_model_cost_vals = evaluateModelCosts(cost_models, model_var_vals);
//      DblVec new_model_cnt_viols = evaluateModelCntViols(cnt_models, model_var_vals);
//
//      double old_model_merit = vecSum(old_model_cost_vals) + merit_error_coeff_ * vecSum(old_model_cnt_viols);
//      double new_model_merit = vecSum(new_model_cost_vals) + merit_error_coeff_ * vecSum(new_model_cnt_viols);
//
//      if (!hasViolation(new_model_cnt_viols) &&
//          old_model_merit - new_model_merit >= min_cnt_improve_ratio_ * merit_error_coeff * vecSum(old_model_cnt_viols)) {
//        goto linesearch; 
//      }
//
//
//      // 3. solve LP:
//      vector<ConvexObjectivePtr> lp_cnt_cost_models = cntsToCosts(cnt_models, 1, lp_model.get());
//      lp_model->update();
//      BOOST_FOREACH(ConvexObjectivePtr& cost, cnt_cost_models)cost->addConstraintsToModel();
//      lp_model->update();
//      QuadExpr objective;
//      BOOST_FOREACH(ConvexObjectivePtr& co, cnt_cost_models){
//        exprInc(objective, co->quad_);
//      }
//      lp_model->setObjective(objective);
//      //      (need to check if GUROBI successfully identify this as an LP
//      CvxOptStatus lp_status = lp_model->optimize();
//      ++results.n_lp_solves;
//      //    if old_model_viols == lp_model_viols > 0 then stop infeasible 
//      // 4. if lp_model_viols == 0:
//      //      while current_model_viols != 0
//      //        new_merit_error_coeff *= 10
//      //        resolve QP
//      //    else
//      //      while old_model_viols - current_model_viols < ep1 * (old_model_viols - lp_model_viols):
//      //        new_merit_error_coeff *= 10
//      //        resolve QP
//      // 5. while old_model_merit(new_merit_error_coeff) - current_model_merit(new_merit_error_coeff) < ep2 * new_merit_error_coeff * (old_model_viols - lp_model_viols):
//      //      new_merit_error_coeff *= 10
//      //      resolve QP
//      linesearch:
//      // 6. alpha = 1
//      //    while true_old_merit - true_new_merit(alpha*step) < eta * alpha * (old_model_merit(new_merit_coeff) - current_model_merit(new_merit_error_coeff)):
//      //      alpha = alpha * tau
//      // 7. choose trust region size:
//      //    true_improve = true_old_merit - true_new_merit
//      //    model_improve = old_model_merit - new_model_merit
//      //    if true_improve < eta_1 * model_improve:
//      //      trust_region = 0.5*vecMax(alpha*step)
//      //    else if true_improve > eta_2 * model_improve:
//      //      trust_region = 2*vecMax(alpha*step)
//      //    else:
//      //      trust_region = vecMax(alpha*step)
//      //
//      //    trust_region = max(trust_region_min, min(trust_region, trust_region_max))
//      //
//      // 8. merit_error_coeff = new_merit_error_coeff
//      //    x = x + alpha * step
//      //        
//      //        
//        while (trust_box_size_ >= min_trust_box_size_) {
//
//          setTrustBoxConstraints(x_);
//          CvxOptStatus status = model_->optimize();
//          ++results_.n_qp_solves;
//          if (status != CVX_SOLVED) {
//            LOG_ERROR("convex solver failed! set TRAJOPT_LOG_THRESH=DEBUG to see solver output. saving model to /tmp/fail2.lp");
//            model_->writeToFile("/tmp/fail2.lp");
//            retval = OPT_FAILED;
//            goto cleanup;
//          }
//          DblVec model_var_vals = model_->getVarValues(model_->getVars());
//
//          DblVec model_cost_vals = evaluateModelCosts(cost_models, model_var_vals);
//          DblVec model_cnt_viols = evaluateModelCntViols(cnt_models, model_var_vals);
//
//          // the n variables of the OptProb happen to be the first n variables in the Model
//          DblVec new_x(model_var_vals.begin(), model_var_vals.begin() + x_.size());
//
//
//
//          if (GetLogLevel() >= util::LevelDebug) {
//            DblVec cnt_costs1 = evaluateModelCosts(cnt_cost_models, model_var_vals);
//            DblVec cnt_costs2 = model_cnt_viols;
//            for (int i=0; i < cnt_costs2.size(); ++i) cnt_costs2[i] *= merit_error_coeff_;
//            LOG_DEBUG("SHOULD BE ALMOST THE SAME: %s ?= %s", CSTR(cnt_costs1), CSTR(cnt_costs2) );
//            // not exactly the same because cnt_costs1 is based on aux variables, but they might not be at EXACTLY the right value
//          }
//
//          DblVec new_cost_vals = evaluateCosts(prob_->getCosts(), new_x);
//          DblVec new_cnt_viols = evaluateConstraintViols(constraints, new_x);
//          ++results_.n_func_evals;
//
//          double old_merit = vecSum(results_.cost_vals) + merit_error_coeff_ * vecSum(results_.cnt_viols);
//          double model_merit = vecSum(model_cost_vals) + merit_error_coeff_ * vecSum(model_cnt_viols);
//          double new_merit = vecSum(new_cost_vals) + merit_error_coeff_ * vecSum(new_cnt_viols);
//          double approx_merit_improve = old_merit - model_merit;
//          double exact_merit_improve = old_merit - new_merit;
//          double merit_improve_ratio = exact_merit_improve / approx_merit_improve;
//
//          if (util::GetLogLevel() >= util::LevelInfo) {
//            LOG_INFO(" ");
//            printCostInfo(results_.cost_vals, model_cost_vals, new_cost_vals,
//                          results_.cnt_viols, model_cnt_viols, new_cnt_viols, cost_names,
//                          cnt_names, merit_error_coeff_);
//            printf("%15s | %10.3e | %10.3e | %10.3e | %10.3e\n", "TOTAL", old_merit, approx_merit_improve, exact_merit_improve, merit_improve_ratio);
//          }
//
//          if (approx_merit_improve < -1e-5) {
//            LOG_ERROR("approximate merit function got worse (%.3e). (convexification is probably wrong to zeroth order)", approx_merit_improve);
//          }
//          if (approx_merit_improve < min_approx_improve_) {
//            LOG_INFO("converged because improvement was small (%.3e < %.3e)", approx_merit_improve, min_approx_improve_);
//            retval = OPT_CONVERGED;
//            goto penaltyadjustment;
//          }
//          if (approx_merit_improve / old_merit < min_approx_improve_frac_) {
//            LOG_INFO(
//                "converged because improvement ratio was small (%.3e < %.3e)",
//                approx_merit_improve/old_merit, min_approx_improve_frac_);
//            retval = OPT_CONVERGED;
//            goto penaltyadjustment;
//          } 
//          else if (exact_merit_improve < 0 || merit_improve_ratio < improve_ratio_threshold_) {
//            adjustTrustRegion(trust_shrink_ratio_);
//            LOG_INFO("shrunk trust region. new box size: %.4f",
//                trust_box_size_);
//          } else {
//            x_ = new_x;
//            results_.cost_vals = new_cost_vals;
//            results_.cnt_viols = new_cnt_viols;
//            adjustTrustRegion(trust_expand_ratio_);
//            LOG_INFO("expanded trust region. new box size: %.4f",trust_box_size_);
//            break;
//          }
//        }
//
//        if (trust_box_size_ < min_trust_box_size_) {
//          LOG_INFO("converged because trust region is tiny");
//          retval = OPT_CONVERGED;
//          goto penaltyadjustment;
//        } else if (iter >= max_iter_) {
//          LOG_INFO("iteration limit");
//          retval = OPT_SCO_ITERATION_LIMIT;
//          goto cleanup;
//        }
//      }
//
//      penaltyadjustment:
//      if (results_.cnt_viols.empty() || vecMax(results_.cnt_viols) < cnt_tolerance_) {
//        if (results_.cnt_viols.size() > 0) LOG_INFO("woo-hoo! constraints are satisfied (to tolerance %.2e)", cnt_tolerance_);
//        goto cleanup;
//      }
//      else {
//        LOG_INFO("not all constraints are satisfied. increasing penalties");
//        merit_error_coeff_ *= merit_coeff_increase_ratio_;
//        trust_box_size_ = fmax(trust_box_size_, min_trust_box_size_ / trust_shrink_ratio_ * 1.5);
//      }
//
//
//
//    }
//    retval = OPT_PENALTY_ITERATION_LIMIT;
//    LOG_INFO("optimization couldn't satisfy all constraints");
//
//
//    cleanup:
//    assert(retval != INVALID && "should never happen");
//    results_.status = retval;
//    results_.total_cost = vecSum(results_.cost_vals);
//    LOG_INFO("\n==================\n%s==================", CSTR(results_));
//    callCallbacks(x_);
//
//    return retval;
//
//  }
//};

int main() {
  OptProbPtr prob(new OptProb());

  VarVector vars;
  {
    vector<string> var_names;
    var_names.push_back("x1");
    var_names.push_back("x2");
    //var_names.push_back("x3");
    DblVec lbs, ubs;
    lbs.push_back(-INFINITY); lbs.push_back(-INFINITY);// lbs.push_back(0);
    ubs.push_back(INFINITY); ubs.push_back(INFINITY);// ubs.push_back(INFINITY);
    vars = prob->createVariables(var_names, lbs, ubs);
  }

  {
    ScalarOfVectorPtr f(new CostError());
    prob->addCost(CostPtr(new CostFromFunc(f, vars, "cost", true)));
  }

  {
    VectorOfVectorPtr f(new CntError1());
    prob->addConstraint(ConstraintPtr(new ConstraintFromFunc(f, vars, Vector1d::Ones(), EQ, "constraint_1")));
  }

  {
    VectorOfVectorPtr f(new CntError2());
    prob->addConstraint(ConstraintPtr(new ConstraintFromFunc(f, vars, Vector1d::Ones(), EQ, "constraint_2")));
  }

  BasicTrustRegionSQP opt(prob);

  {
    DblVec initVec = toDblVec(Vector2d(1, 0));//-3, 1, 1));
    opt.initialize(initVec);
  }

  opt.optimize();
  return 0;
}
