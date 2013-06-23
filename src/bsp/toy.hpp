#pragma once

#include "bsp.hpp"

using namespace BSP;

namespace ToyBSP {
  typedef Matrix<double, 5, 1> Vector5d;

  class ToyBSPProblemHelper : public BSPProblemHelper {
  public:
    int input_dt;
    int T;
    Vector2d start;
    Vector2d goal;
    const static int n_dof = 2;

    ToyBSPProblemHelper();
    void set_state_dim(int new_state_dim);
    void set_observe_dim(int new_observe_dim);
    void set_control_dim(int new_control_dim);
    void configure_problem(OptProb& prob);
    void configure_optimizer(OptProb& prob, BasicTrustRegionSQP& opt);
  };

  typedef boost::shared_ptr<ToyBSPProblemHelper> ToyBSPProblemHelperPtr;


  class ToyStateFunc : public StateFunc<Vector2d, Vector2d, Vector2d> {
  public:
    ToyStateFunc(ToyBSPProblemHelperPtr helper);
    virtual Vector2d operator()(const Vector2d& x, const Vector2d& u, const Vector2d& m) const ;
  protected:
    ToyBSPProblemHelperPtr helper;
  };

  class ToyObserveFunc : public ObserveFunc<Vector2d, Vector2d, Vector2d> {
  public:
    ToyObserveFunc(ToyBSPProblemHelperPtr helper);
    virtual Vector2d operator()(const Vector2d& x, const Vector2d& n) const;
  protected:
    ToyBSPProblemHelperPtr helper;
  };

  class ToyBeliefFunc : public BeliefFunc<ToyStateFunc, ToyObserveFunc> {
  public:
    ToyBeliefFunc(ToyBSPProblemHelperPtr helper);
    ToyBeliefFunc(ToyBSPProblemHelperPtr helper, double alpha);
    virtual Vector5d operator()(const Vector5d& b, const Vector2d& u) const;
    void set_alpha(double alpha);
    void set_tol(double tol);
  protected:
    double alpha;
    double tol;
    ToyBSPProblemHelperPtr helper;
    Matrix2d get_gamma(const Vector2d& x) const;
    bool sgndist(const Vector2d& x, Vector2d* dists) const;
  };
}
