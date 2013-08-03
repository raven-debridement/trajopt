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

#define CHECK_STATUS(status, model) {                                                                                    \
  if (status != CVX_SOLVED) {                                                                                            \
    LOG_ERROR("convex solver failed! set TRAJOPT_LOG_THRESH=DEBUG to see solver output. saving model to /tmp/fail2.lp"); \
    model->writeToFile("/tmp/fail2.lp");                                                                                 \
    retval = OPT_FAILED;                                                                                                 \
    goto cleanup;                                                                                                        \
  }                                                                                                                      \
}

#define RESOLVE_QP() {                                                                                            \
  BOOST_FOREACH(ConvexObjectivePtr& cost, cnt_cost_models) cost->removeFromModel(model_.get());                   \
  QuadExpr objective;                                                                                             \
  cnt_cost_models = cntsToCosts(cnt_models, merit_error_coeff_);                                                  \
  BOOST_FOREACH(ConvexObjectivePtr& cost, cnt_cost_models) cost->addToModelAndObjective(model_.get(), objective); \
  model_->setObjective(objective);                                                                                \
  CvxOptStatus status = model_->optimize();                                                                       \
  ++results_.n_qp_solves;                                                                                          \
  CHECK_STATUS(status, model_);                                                                                   \
  current_model_var_vals = model_->getVarValues(model_->getVars());                                               \
  current_model_cnt_viols = evaluateModelCntViols(cnt_models, current_model_var_vals, model_.get());              \
  current_model_total_cnt_viol = vecSum(current_model_cnt_viols);                                                 \
}

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
    Vector1d ans; ans << a(0)*a(0)*a(0);//a(0) - 1 - a(2);
    return ans;
  }
};


void vecDiff(const DblVec& a, const DblVec& b, DblVec& c) {
  c.resize(min(a.size(), b.size()));
  for (int i = 0; i < min(a.size(), b.size()); ++i) {
    c[i] = a[i] - b[i];
  }
}

void vecAdd(const DblVec& a, const DblVec& b, DblVec& c) {
  c.resize(min(a.size(), b.size()));
  for (int i = 0; i < min(a.size(), b.size()); ++i) {
    c[i] = a[i] + b[i];
  }
}

void vecScale(const DblVec& a, double scale, DblVec& c) {
  c.resize(a.size());
  for (int i = 0; i < a.size(); ++i) {
    c[i] = a[i] * scale;
  }
}


class LineSearchSQP : public BasicTrustRegionSQP {
public:

  double trust_shrink_ratio_;
  double trust_expand_ratio_;
  double cnt_tolerance_;
  double merit_coeff_increase_ratio_;
  double merit_error_coeff_;
  double trust_box_size_;

  double min_cnt_improve_ratio;
  double min_model_merit_improve_ratio_;
  double line_search_shrink_ratio_;
  double min_merit_improve_ratio;
  double trust_region_shrink_threshold_;
  double trust_region_expand_threshold_;
  double min_trust_box_size_;
  double max_trust_box_size_;
  double opt_eps;

	LineSearchSQP() : BasicTrustRegionSQP() {
    initParameters(); 
  }

	LineSearchSQP(OptProbPtr prob) : BasicTrustRegionSQP(prob) {
    initParameters(); 
  }

  bool hasViolation(const DblVec& cnt_viols) {
    return cnt_viols.size() > 0 && vecMax(cnt_viols) < cnt_tolerance_;
  }

  void initParameters() {
    
    trust_shrink_ratio_= 0.5;
    trust_expand_ratio_ = 2;
    cnt_tolerance_ = 1e-4;
    merit_coeff_increase_ratio_ = 10;
    merit_error_coeff_ = 1;
    trust_box_size_ = 1;

    min_cnt_improve_ratio = 0.1; // ep1
    min_model_merit_improve_ratio_ = 0.1; // ep2
    line_search_shrink_ratio_ = 0.5; // tau
    min_merit_improve_ratio = 1e-4; // eta
    trust_region_shrink_threshold_ = 0.25; // eta_1
    trust_region_expand_threshold_ = 0.75; // eta_2
    min_trust_box_size_ = 1e-3; // trust_region_min
    max_trust_box_size_ = 1e3; // trust_region_max

    opt_eps = 1e-6;

  }

	OptStatus optimize() {

    vector<string> cost_names = getCostNames(prob_->getCosts());
    vector<ConstraintPtr> constraints = prob_->getConstraints();
    vector<string> cnt_names = getCntNames(constraints);

    DblVec& x_ = results_.x; // just so I don't have to rewrite code
    if (x_.size() == 0) PRINT_AND_THROW("you forgot to initialize!");
    if (!prob_) PRINT_AND_THROW("you forgot to set the optimization problem");    
    
    x_ = prob_->getClosestFeasiblePoint(x_);

    // initialize LP model
    //ModelPtr lp_model = model_->cloneModel();
    
    assert(x_.size() == prob_->getVars().size());
    assert(prob_->getCosts().size() > 0 || constraints.size() > 0);

    OptStatus retval = INVALID;

    for (int iter=1; ; ++iter) { /* sqp loop */
      callCallbacks(x_);

      LOG_DEBUG("current iterate: %s", CSTR(x_));
      LOG_INFO("sqp iteration %i", iter);

      // speedup: if you just evaluated the cost when doing the line search, use that
      if (results_.cost_vals.empty() && results_.cnt_viols.empty()) { //only happens on the first iteration
        results_.cnt_viols = evaluateConstraintViols(constraints, x_, model_.get());
        results_.cost_vals = evaluateCosts(prob_->getCosts(), x_, model_.get());
        assert(results_.n_func_evals == 0);
        ++results_.n_func_evals;
      }

      // DblVec new_cnt_viols = evaluateConstraintViols(constraints, x_);
      // DblVec new_cost_vals = evaluateCosts(prob_->getCosts(), x_);
      // cout << "costs" << endl;
      // for (int i=0; i < new_cnt_viols.size(); ++i) {
      //   cout << cnt_names[i] << " " << new_cnt_viols[i] - results_.cnt_viols[i] << endl;
      // }
      // for (int i=0; i < new_cost_vals.size(); ++i) {
      //   cout << cost_names[i] << " " << new_cost_vals[i] - results_.cost_vals[i] << endl;
      // }


      vector<ConvexObjectivePtr> cost_models = convexifyCosts(prob_->getCosts(),x_);
      vector<ConvexConstraintsPtr> cnt_models = convexifyConstraints(constraints, x_);
      vector<ConvexObjectivePtr> cnt_cost_models = cntsToCosts(cnt_models, merit_error_coeff_);

      QuadExpr objective;

      BOOST_FOREACH(ConvexObjectivePtr& cost, cost_models) cost->addToModelAndObjective(model_.get(), objective);
      BOOST_FOREACH(ConvexObjectivePtr& cost, cnt_cost_models) cost->addToModelAndObjective(model_.get(), objective);
      
//    objective = cleanupExpr(objective);
      model_->setObjective(objective);

//    if (logging::filter() >= IPI_LEVEL_DEBUG) {
//      DblVec model_cost_vals;
//      BOOST_FOREACH(ConvexObjectivePtr& cost, cost_models) {
//        model_cost_vals.push_back(cost->value(x));
//      }
//      LOG_DEBUG("model costs %s should equalcosts  %s", printer(model_cost_vals), printer(cost_vals));
//    }
      


      // archive the model variable values before optimization for later use
      DblVec old_model_var_vals = x_;//model_->getVarValues(model_->getVars());
      for (int i = 0; i < model_->getVars().size() - x_.size(); ++i) {
        old_model_var_vals.push_back(0.); 
      }


      // 1. solve the problem constructed above _without_ a trust region size
      CvxOptStatus status = model_->optimize(); 
      ++results_.n_qp_solves;
      CHECK_STATUS(status, model_);

      DblVec current_model_var_vals = model_->getVarValues(model_->getVars());
      DblVec new_x(current_model_var_vals.begin(), current_model_var_vals.begin() + x_.size());
      // step = new_x - x_
      DblVec step; vecDiff(new_x, x_, step);

      //  check if result is close to the previous x and if all cnt viols are satisfied for the old x
      //    (stop in that case)
      if (vecAbsSum(step) < opt_eps) {
        if (!hasViolation(results_.cnt_viols)) {
          if (results_.cnt_viols.size() > 0) LOG_INFO("woo-hoo! constraints are satisfied (to tolerance %.2e)", cnt_tolerance_);
          retval = OPT_CONVERGED;
          goto cleanup;
        }
      }


      //double new_merit_error_coeff = merit_error_coeff_;

      // 2. if model cnt viols = 0 and old_model_merit - new_model_merit >= min_approx_improve_ratio * merit_error_coeff * old_model_viols
      //    go to 6.
      DblVec old_model_cost_vals = evaluateModelCosts(cost_models, old_model_var_vals, model_.get());
      DblVec old_model_cnt_viols = evaluateModelCntViols(cnt_models, old_model_var_vals, model_.get());
      DblVec current_model_cost_vals = evaluateModelCosts(cost_models, current_model_var_vals, model_.get());
      DblVec current_model_cnt_viols = evaluateModelCntViols(cnt_models, current_model_var_vals, model_.get());

      double old_model_total_cost_val = vecSum(old_model_cost_vals);
      double old_model_total_cnt_viol = vecSum(old_model_cnt_viols);
      double current_model_total_cost_val = vecSum(current_model_cost_vals);
      double current_model_total_cnt_viol = vecSum(current_model_cnt_viols);

      double old_model_merit = old_model_total_cost_val + merit_error_coeff_ * old_model_total_cnt_viol;
      double current_model_merit = current_model_total_cost_val + merit_error_coeff_ * current_model_total_cnt_viol;

      if (!hasViolation(current_model_cnt_viols) &&
          old_model_merit - current_model_merit >= min_model_merit_improve_ratio_ * merit_error_coeff_ * old_model_total_cnt_viol) {
        //goto linesearch; 
      } else {
        
        // 3. solve LP:
        vector<ConvexObjectivePtr> lp_cnt_cost_models = cntsToCosts(cnt_models, 1);
        AffExpr lp_objective;

        BOOST_FOREACH(ConvexObjectivePtr& cost, cost_models) cost->removeFromModel(model_.get());
        BOOST_FOREACH(ConvexObjectivePtr& cost, cnt_cost_models) cost->removeFromModel(model_.get());
        BOOST_FOREACH(ConvexObjectivePtr& cost, cnt_cost_models) cost->addToModelAndObjective(model_.get(), lp_objective);
        lp_model->setObjective(lp_objective);
        setTrustBoxConstraints(x_, lp_model.get());

        {
          CvxOptStatus lp_status = lp_model->optimize();
          ++results_.n_lp_solves;
          CHECK_STATUS(lp_status, lp_model);
        }

        DblVec lp_model_var_vals = lp_model->getVarValues(lp_model->getVars());
        DblVec lp_x(lp_model_var_vals.begin(), lp_model_var_vals.begin() + x_.size());

        //    if old_model_viols == lp_model_viols > 0 then stop infeasible 
        DblVec lp_model_cnt_viols = evaluateModelCntViols(cnt_models, lp_model_var_vals, lp_model.get());
        double lp_model_total_cnt_viol = vecSum(lp_model_cnt_viols);
        if (fabs(old_model_total_cnt_viol - lp_model_total_cnt_viol) < opt_eps && hasViolation(old_model_cnt_viols)) {
          retval = OPT_INFEASIBLE;
          goto cleanup;
        }

        // 4.
        if (fabs(lp_model_total_cnt_viol) < opt_eps) {
          while (fabs(current_model_total_cnt_viol) > opt_eps) {
            merit_error_coeff_ *= merit_coeff_increase_ratio_;
            RESOLVE_QP();
          }
        } else {
          while (old_model_total_cnt_viol - current_model_total_cnt_viol < min_cnt_improve_ratio * (old_model_total_cnt_viol - lp_model_total_cnt_viol)) {
            merit_error_coeff_ *= merit_coeff_increase_ratio_;
            RESOLVE_QP();
          }
        }

        current_model_cost_vals = evaluateModelCosts(cost_models, current_model_var_vals, model_.get());
        current_model_total_cost_val = vecSum(current_model_cost_vals);

        // 5.
        while ( (old_model_total_cost_val + merit_error_coeff_ * old_model_total_cnt_viol) /* old_model_merit */
               -(current_model_total_cost_val + merit_error_coeff_ * current_model_total_cnt_viol) /* current_model_merit */
              <  min_model_merit_improve_ratio_ * merit_error_coeff_ * (old_model_total_cnt_viol - lp_model_total_cnt_viol) ) {
          merit_error_coeff_ *= merit_coeff_increase_ratio_;
          RESOLVE_QP();
          current_model_cost_vals = evaluateModelCosts(cost_models, current_model_var_vals, model_.get());
          current_model_total_cost_val = vecSum(current_model_cost_vals);
        }

        old_model_merit = old_model_total_cost_val + merit_error_coeff_ * old_model_total_cnt_viol;
        current_model_merit = current_model_total_cost_val + merit_error_coeff_ * current_model_total_cnt_viol;

      }

      // 6.
      //linesearch:
      double alpha = 1;
      double old_merit = vecSum(results_.cost_vals) + merit_error_coeff_ * vecSum(results_.cnt_viols);
      DblVec current_x, current_cost_vals, current_cnt_viols, chosen_step;
      double current_merit = 0, merit_improve_bound = 0;
      do {
        vecScale(step, alpha, chosen_step);
        vecAdd(x_, chosen_step, current_x);
        //current_x = x + alpha * step;//new_x; // new_x = x + 1 * step
        current_cost_vals = evaluateCosts(prob_->getCosts(), current_x, model_.get());
        current_cnt_viols = evaluateConstraintViols(constraints, current_x, model_.get());
        current_merit = vecSum(current_cost_vals) + merit_error_coeff_ * vecSum(current_cnt_viols);
        merit_improve_bound = min_merit_improve_ratio * alpha * (old_model_merit - current_model_merit);
        alpha = alpha * line_search_shrink_ratio_;
      } while (old_merit - current_merit < merit_improve_bound);

      
      // 7. choose trust region size:

      double actual_improve = old_merit - current_merit;
      double model_improve = old_model_merit - current_model_merit;

      x_ = current_x;
      results_.cost_vals = current_cost_vals;
      results_.cnt_viols = current_cnt_viols;

      if (actual_improve < trust_region_shrink_threshold_ * model_improve) {
        trust_box_size_ = trust_shrink_ratio_ * vecMax(chosen_step);
      } else if (actual_improve > trust_region_expand_threshold_ * model_improve) {
        trust_box_size_ = trust_expand_ratio_ * vecMax(chosen_step);
      } else {
        trust_box_size_ = vecMax(chosen_step);
      }

      trust_box_size_ = max(min_trust_box_size_, min(trust_box_size_, max_trust_box_size_));
    }

    cleanup:
    assert(retval != INVALID && "should never happen");
    results_.status = retval;
    results_.total_cost = vecSum(results_.cost_vals);
    LOG_INFO("\n==================\n%s==================", CSTR(results_));
    callCallbacks(x_);

    return retval;

  }
};

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

  LineSearchSQP opt(prob);

  {
    DblVec initVec = toDblVec(Vector2d(1, 0));//-3, 1, 1));
    opt.initialize(initVec);
  }

  opt.optimize();
  return 0;
}
