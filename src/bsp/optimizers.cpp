#include "common.hpp"
#include "optimizers.hpp"

#define SEPARATE_COEFF

namespace BSP {

  BSPTrustRegionSQP::BSPTrustRegionSQP() : BasicTrustRegionSQP() {}

  BSPTrustRegionSQP::BSPTrustRegionSQP(OptProbPtr prob) : BasicTrustRegionSQP(prob) {}

  void BSPTrustRegionSQP::addMeritDoneCallback(const Callback& cb) {
    merit_done_callbacks_.push_back(cb);
  }

  void BSPTrustRegionSQP::callMeritDoneCallbacks(DblVec& x) {
    for (int i = 0; i < merit_done_callbacks_.size(); ++i) {
      merit_done_callbacks_[i](prob_.get(), x);
    }
  }

  OptStatus BSPTrustRegionSQP::optimize() {
    vector<string> cost_names = getCostNames(prob_->getCosts());
    vector<ConstraintPtr> constraints = prob_->getConstraints();
    vector<string> cnt_names = getCntNames(constraints);

    vector<ConstraintPtr> dynamics_constraints;
    vector<ConstraintPtr> collision_constraints;

    for (int i = 0; i < constraints.size(); ++i) {
      if (constraints[i]->name().find("collision") != std::string::npos) {
        collision_constraints.push_back(constraints[i]);
      } else {
        dynamics_constraints.push_back(constraints[i]);
      }
    }

    DblVec& x_ = results_.x; // just so I don't have to rewrite code
    if (x_.size() == 0) PRINT_AND_THROW("you forgot to initialize!");
    if (!prob_) PRINT_AND_THROW("you forgot to set the optimization problem");    
    
    x_ = prob_->getClosestFeasiblePoint(x_);

    assert(x_.size() == prob_->getVars().size());
    assert(prob_->getCosts().size() > 0 || constraints.size() > 0);

    double dynamics_merit_error_coeff_ = merit_error_coeff_;
    double collision_merit_error_coeff_ = merit_error_coeff_;

    //cout << "dynamics merit error coeff: " << dynamics_merit_error_coeff_ << endl;
    //cout << "collision merit error coeff: " << collision_merit_error_coeff_ << endl;


    OptStatus retval = INVALID;

    for (int merit_increases=0; merit_increases < max_merit_coeff_increases_; ) { /* merit adjustment loop */
      //++results_.n_merit_increases;

      for (int iter=1; ; ++iter) { /* sqp loop */
        callCallbacks(x_);

        LOG_DEBUG("current iterate: %s", CSTR(x_));
        LOG_INFO("merit increases iteration: %i; sqp iteration %i", merit_increases, iter);

        // speedup: if you just evaluated the cost when doing the line search, use that
        if (results_.cost_vals.empty() && results_.dynamics_cnt_viols.empty() && results_.collision_cnt_viols.empty()) { //only happens on the first iteration
          results_.dynamics_cnt_viols = evaluateConstraintViols(dynamics_constraints, x_, model_.get());
          results_.collision_cnt_viols = evaluateConstraintViols(collision_constraints, x_, model_.get());
          results_.cnt_viols = concat(results_.dynamics_cnt_viols, results_.collision_cnt_viols);
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
        vector<ConvexConstraintsPtr> dynamics_cnt_models = convexifyConstraints(dynamics_constraints, x_);
        vector<ConvexConstraintsPtr> collision_cnt_models = convexifyConstraints(collision_constraints, x_);
        vector<ConvexObjectivePtr> dynamics_cnt_cost_models = cntsToCosts(dynamics_cnt_models, dynamics_merit_error_coeff_);
        vector<ConvexObjectivePtr> collision_cnt_cost_models = cntsToCosts(collision_cnt_models, collision_merit_error_coeff_);

        QuadExpr objective;

        BOOST_FOREACH(ConvexObjectivePtr& cost, cost_models) cost->addToModelAndObjective(model_.get(), objective);
        BOOST_FOREACH(ConvexObjectivePtr& cost, dynamics_cnt_cost_models) cost->addToModelAndObjective(model_.get(), objective);
        BOOST_FOREACH(ConvexObjectivePtr& cost, collision_cnt_cost_models) cost->addToModelAndObjective(model_.get(), objective);
        //model_->update();

        //BOOST_FOREACH(ConvexObjectivePtr& co, cost_models) exprInc(objective, co->quad_);
        //BOOST_FOREACH(ConvexObjectivePtr& co, cnt_cost_models){
        //  exprInc(objective, co->quad_);
        //}
  //    objective = cleanupExpr(objective);
        model_->setObjective(objective);

  //    if (logging::filter() >= IPI_LEVEL_DEBUG) {
  //      DblVec model_cost_vals;
  //      BOOST_FOREACH(ConvexObjectivePtr& cost, cost_models) {
  //        model_cost_vals.push_back(cost->value(x));
  //      }
  //      LOG_DEBUG("model costs %s should equalcosts  %s", printer(model_cost_vals), printer(cost_vals));
  //    }
  //
        DblVec old_model_var_vals = x_;//model_->getVarValues(model_->getVars());
        for (int i = 0; i < model_->getVars().size() - x_.size(); ++i) {
          old_model_var_vals.push_back(0.); 
        }

        DblVec model_var_vals;
        DblVec model_cost_vals;
        DblVec model_dynamics_cnt_viols;
        DblVec model_collision_cnt_viols;
        DblVec new_cost_vals;
        DblVec new_dynamics_cnt_viols;
        DblVec new_collision_cnt_viols;
        DblVec new_x;

        while (trust_box_size_ >= min_trust_box_size_) {

          //cout << "setting trust box size " << trust_box_size_ << " around: ";
          //for (int i = 0; i < x_.size(); ++i) cout << x_[i] << " ";
          //cout << endl;

          //cout << "current dynamics penalty coefficient: " << dynamics_merit_error_coeff_ << endl;
          //cout << "current collision penalty coefficient: " << collision_merit_error_coeff_ << endl;

          setTrustBoxConstraints(x_, model_.get());
          CvxOptStatus status = model_->optimize();
          ++results_.n_qp_solves;
          if (status != CVX_SOLVED) {
            LOG_ERROR("convex solver failed! set TRAJOPT_LOG_THRESH=DEBUG to see solver output. saving model to /tmp/fail2.lp");
            model_->writeToFile("/tmp/fail2.lp");
            retval = OPT_FAILED;
            goto cleanup;
          }
          model_var_vals = model_->getVarValues(model_->getVars());

          model_cost_vals = evaluateModelCosts(cost_models, model_var_vals, model_.get());
          model_dynamics_cnt_viols = evaluateModelCntViols(dynamics_cnt_models, model_var_vals, model_.get());
          model_collision_cnt_viols = evaluateModelCntViols(collision_cnt_models, model_var_vals, model_.get());

          // the n variables of the OptProb happen to be the first n variables in the Model
          new_x = DblVec(model_var_vals.begin(), model_var_vals.begin() + x_.size());

          //cout << "solution: ";
          //for (int i = 0; i < new_x.size(); ++i) cout << new_x[i] << " ";
          //cout << endl;

          if (GetLogLevel() >= util::LevelDebug) {
            DblVec dynamics_cnt_costs1 = evaluateModelCosts(dynamics_cnt_cost_models, model_var_vals, model_.get());
            DblVec dynamics_cnt_costs2 = model_dynamics_cnt_viols;
            for (int i=0; i < dynamics_cnt_costs2.size(); ++i) dynamics_cnt_costs2[i] *= dynamics_merit_error_coeff_;
            DblVec collision_cnt_costs1 = evaluateModelCosts(collision_cnt_cost_models, model_var_vals, model_.get());
            DblVec collision_cnt_costs2 = model_collision_cnt_viols;
            for (int i=0; i < collision_cnt_costs2.size(); ++i) collision_cnt_costs2[i] *= collision_merit_error_coeff_;


            LOG_DEBUG("SHOULD BE ALMOST THE SAME: %s ?= %s", CSTR(concat(dynamics_cnt_costs1, collision_cnt_costs1)), CSTR(concat(dynamics_cnt_costs2, collision_cnt_costs2)) );
            // not exactly the same because cnt_costs1 is based on aux variables, but they might not be at EXACTLY the right value
          }

          new_cost_vals = evaluateCosts(prob_->getCosts(), new_x, model_.get());
          new_dynamics_cnt_viols = evaluateConstraintViols(dynamics_constraints, new_x, model_.get());
          new_collision_cnt_viols = evaluateConstraintViols(collision_constraints, new_x, model_.get());
          ++results_.n_func_evals;

          double old_merit = vecSum(results_.cost_vals) + dynamics_merit_error_coeff_ * vecSum(results_.dynamics_cnt_viols) + collision_merit_error_coeff_ * vecSum(results_.collision_cnt_viols);
          double model_merit = vecSum(model_cost_vals) + dynamics_merit_error_coeff_ * vecSum(model_dynamics_cnt_viols) + collision_merit_error_coeff_ * vecSum(model_collision_cnt_viols);
          double new_merit = vecSum(new_cost_vals) + dynamics_merit_error_coeff_ * vecSum(new_dynamics_cnt_viols) + collision_merit_error_coeff_ * vecSum(new_collision_cnt_viols);
          double approx_merit_improve = old_merit - model_merit;
          double exact_merit_improve = old_merit - new_merit;
          double merit_improve_ratio = exact_merit_improve / approx_merit_improve;

          if (util::GetLogLevel() >= util::LevelInfo) {
            LOG_INFO(" ");
            printCostInfo(results_.cost_vals, model_cost_vals, new_cost_vals,
                          results_.cnt_viols, concat(model_dynamics_cnt_viols, model_collision_cnt_viols), concat(new_dynamics_cnt_viols, new_collision_cnt_viols), cost_names,
                          cnt_names, merit_error_coeff_);
            printf("%15s | %10.3e | %10.3e | %10.3e | %10.3e\n", "TOTAL", old_merit, approx_merit_improve, exact_merit_improve, merit_improve_ratio);
          }

          if (approx_merit_improve < -1e-5) {
            LOG_ERROR("approximate merit function got worse (%.3e). (convexification is probably wrong to zeroth order)", approx_merit_improve);
          }
          if (approx_merit_improve < min_approx_improve_) {
            LOG_INFO("converged because improvement was small (%.3e < %.3e)", approx_merit_improve, min_approx_improve_);
            retval = OPT_CONVERGED;
            goto penaltyadjustment;
          }
          if (approx_merit_improve / old_merit < min_approx_improve_frac_) {
            LOG_INFO(
                "converged because improvement ratio was small (%.3e < %.3e)",
                approx_merit_improve/old_merit, min_approx_improve_frac_);
            retval = OPT_CONVERGED;
            goto penaltyadjustment;
          } 
          else if (exact_merit_improve < 0 || merit_improve_ratio < improve_ratio_threshold_) {
            adjustTrustRegion(trust_shrink_ratio_);
            LOG_INFO("shrunk trust region. new box size: %.4f",
                trust_box_size_);
          } else {
            x_ = new_x;
            //cout << "adopting new_x: ";
            //for (int i = 0; i < x_.size(); ++i) { cout << x_[i] << " "; } cout << endl;
            results_.cost_vals = new_cost_vals;
            results_.dynamics_cnt_viols = new_dynamics_cnt_viols;
            results_.collision_cnt_viols = new_collision_cnt_viols;
            results_.cnt_viols = concat(results_.dynamics_cnt_viols, results_.collision_cnt_viols);
            
            adjustTrustRegion(trust_expand_ratio_);
            LOG_INFO("expanded trust region. new box size: %.4f",trust_box_size_);
            break;
          }
        }

        if (trust_box_size_ < min_trust_box_size_) {
          LOG_INFO("converged because trust region is tiny");
          retval = OPT_CONVERGED;
          goto penaltyadjustment;
        } else if (iter >= max_iter_) {
          LOG_INFO("iteration limit");
          retval = OPT_SCO_ITERATION_LIMIT;
          goto cleanup;
        }

        continue;

        penaltyadjustment:
        if (results_.cnt_viols.empty() || vecMax(results_.cnt_viols) < cnt_tolerance_) {
          if (results_.cnt_viols.size() > 0) LOG_INFO("woo-hoo! constraints are satisfied (to tolerance %.2e)", cnt_tolerance_);
          goto cleanup;
        }
        else {
          LOG_INFO("not all constraints are satisfied. increasing penalties");
          //++merit_increases;

#ifdef SEPARATE_COEFF
          
          if (results_.dynamics_cnt_viols.empty() || vecMax(results_.dynamics_cnt_viols) < cnt_tolerance_) {
            ++merit_increases;
            ++results_.n_merit_increases;
            collision_merit_error_coeff_ *= merit_coeff_increase_ratio_;
            //dynamics_merit_error_coeff_ *= merit_coeff_increase_ratio_;
            //cout << "increasing both" << endl;
            cout << "increasing just collision" << endl;
            trust_box_size_ = fmax(trust_box_size_, min_trust_box_size_ / trust_shrink_ratio_ * 1.5);

          } else {
            ++merit_increases;
            ++results_.n_merit_increases;
            dynamics_merit_error_coeff_ *= merit_coeff_increase_ratio_;
            cout << "increasing just dynamics" << endl;
            trust_box_size_ = fmax(trust_box_size_, min_trust_box_size_ / trust_shrink_ratio_ * 1.5);
          }
#else
            ++merit_increases;
            ++results_.n_merit_increases;
            collision_merit_error_coeff_ *= merit_coeff_increase_ratio_;
            dynamics_merit_error_coeff_ *= merit_coeff_increase_ratio_;
            cout << "increasing both" << endl;

#endif

          cout << "(" << collision_merit_error_coeff_ << ", " << dynamics_merit_error_coeff_ << ")" << endl;

          callMeritDoneCallbacks(x_);
          break;
        }
      }
    }
    retval = OPT_PENALTY_ITERATION_LIMIT;
    LOG_INFO("optimization couldn't satisfy all constraints");


    cleanup:
    assert(retval != INVALID && "should never happen");
    results_.status = retval;
    results_.total_cost = vecSum(results_.cost_vals);
    LOG_INFO("\n==================\n%s==================", CSTR(results_));
    callCallbacks(x_);

    return retval;

  }

}
