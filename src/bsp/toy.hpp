#pragma once

#include "bsp.hpp"

using namespace BSP;

namespace ToyBSP {
  typedef Matrix<double, 5, 1> Vector5d;
  typedef Vector2d StateT;
  typedef Vector2d StateNoiseT;
  typedef Vector2d ControlT;
  typedef Vector2d ObserveT;
  typedef Vector2d ObserveNoiseT;
  typedef Vector5d BeliefT;
  typedef Matrix2d GammaT;

  class ToyBSPProblemHelper;
  typedef boost::shared_ptr<ToyBSPProblemHelper> ToyBSPProblemHelperPtr;

  class ToyStateFunc : public StateFunc<StateT, ControlT, StateNoiseT> {
  public:
    ToyStateFunc();
    ToyStateFunc(ToyBSPProblemHelperPtr helper);
    virtual StateT operator()(const StateT& x, const ControlT& u, const StateNoiseT& m) const ;
  protected:
    ToyBSPProblemHelperPtr helper;
  };

  class ToyObserveFunc : public ObserveFunc<StateT, ObserveT, ObserveNoiseT> {
  public:
    ToyObserveFunc();
    ToyObserveFunc(ToyBSPProblemHelperPtr helper);
    virtual ObserveT operator()(const StateT& x, const ObserveNoiseT& n) const;
  protected:
    ToyBSPProblemHelperPtr helper;
  };

  class ToyBeliefFunc : public BeliefFunc<ToyStateFunc, ToyObserveFunc> {
  public:
    ToyBeliefFunc();
    ToyBeliefFunc(ToyBSPProblemHelperPtr helper, StateFuncPtr f, ObserveFuncPtr h);
    virtual BeliefT operator()(const BeliefT& b, const ControlT& u) const;
    void set_alpha(double alpha);
    void set_tol(double tol);
  protected:
    double alpha;
    double tol;
    ToyBSPProblemHelperPtr helper;
    StateFuncPtr f;
    ObserveFuncPtr h;
    GammaT get_gamma(const StateT& x) const;
    bool sgndist(const StateT& x, StateT* dists) const;
  };

  class ToyBSPProblemHelper : public BSPProblemHelper, public boost::enable_shared_from_this<ToyBSPProblemHelper> {
  public:
    typedef typename BeliefConstraint<ToyBeliefFunc>::Ptr BeliefConstraintPtr;

    int input_dt;
    int T;
    StateT start;
    StateT goal;
    const static int n_dof = 2;

    VarArray state_vars;
    VarArray sigma_vars;
    VarArray control_vars;
    VarArray belief_vars;

    ToyStateFunc::Ptr state_func;
    ToyObserveFunc::Ptr observe_func;
    ToyBeliefFunc::Ptr belief_func;
    vector< BeliefConstraintPtr > belief_constraints;

    ToyBSPProblemHelper();
    void set_state_dim(int new_state_dim);
    void set_observe_dim(int new_observe_dim);
    void set_control_dim(int new_control_dim);
    void configure_problem(OptProb& prob);
    void create_variables(OptProb& prob);
    void add_variance_cost(OptProb& prob);
    void add_control_cost(OptProb& prob);
    void add_start_constraint(OptProb& prob);
    void add_goal_constraint(OptProb& prob);
    void add_belief_constraint(OptProb& prob);
    void configure_optimizer(OptProb& prob, BasicTrustRegionSQP& opt);
    void init_optimize_variables(OptProb& prob, BasicTrustRegionSQP& opt);
  };
}
