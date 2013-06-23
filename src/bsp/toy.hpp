#pragma once

#include "bsp.hpp"

using namespace BSP;

namespace ToyBSP {

  // This will generate a bunch of types like StateT, ControlT, etc.
  BSP_TYPEDEFS(
    2, // state_dim
    2, // state_noise_dim
    2, // control_dim
    2, // observe_dim
    2, // observe_noise_dim
    3, // sigma_dof
    5 // belief_dim
  );

  class ToyBSPProblemHelper;
  typedef boost::shared_ptr<ToyBSPProblemHelper> ToyBSPProblemHelperPtr;

  class ToyStateFunc : public StateFunc<StateT, ControlT, StateNoiseT> {
  public:
    typedef boost::shared_ptr<ToyStateFunc> Ptr;
    ToyStateFunc();
    ToyStateFunc(ToyBSPProblemHelperPtr helper);
    ToyBSPProblemHelperPtr toy_helper;
    virtual StateT operator()(const StateT& x, const ControlT& u, const StateNoiseT& m) const ;
  };

  class ToyObserveFunc : public ObserveFunc<StateT, ObserveT, ObserveNoiseT> {
  public:
    typedef boost::shared_ptr<ToyObserveFunc> Ptr;
    ToyObserveFunc();
    ToyObserveFunc(ToyBSPProblemHelperPtr helper);
    ToyBSPProblemHelperPtr toy_helper;
    virtual ObserveT operator()(const StateT& x, const ObserveNoiseT& n) const;
  };

  class ToyBeliefFunc : public BeliefFunc<ToyStateFunc, ToyObserveFunc> {
  public:
    typedef boost::shared_ptr<ToyBeliefFunc> Ptr;
    ToyBeliefFunc();
    ToyBeliefFunc(ToyBSPProblemHelperPtr helper, StateFuncPtr f, ObserveFuncPtr h);
    virtual BeliefT operator()(const BeliefT& b, const ControlT& u) const;
    void set_alpha(double alpha);
    void scale_alpha(double scale);
    void set_tol(double tol);
    ToyBSPProblemHelperPtr toy_helper;
  protected:
    double alpha;
    double tol;
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
    void configure_optimizer(OptProb& prob, BSPTrustRegionSQP& opt);
    void init_optimize_variables(OptProb& prob, BSPTrustRegionSQP& opt);
    void merit_done_callback(OptProb*, DblVec& x);
    void add_optimizer_callback(OptProb& prob, BSPTrustRegionSQP& opt);
  };
}
