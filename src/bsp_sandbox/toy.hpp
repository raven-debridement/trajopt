#pragma once

#include "bsp/bsp.hpp"
#include "geometry_2d.hpp"
#include <QtCore>
#include <QImage>
#include <qevent.h>

using namespace BSP;
using namespace Geometry2D;

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

  typedef Matrix<double, 1, 1> Vector1d;

  // state: { x, y }
  // control: { dx, dy }
  // observation: { x, y }

  class ToyBSPProblemHelper;
  typedef boost::shared_ptr<ToyBSPProblemHelper> ToyBSPProblemHelperPtr;

  class ToyStateFunc : public StateFunc<StateT, ControlT, StateNoiseT> {
  public:
    typedef boost::shared_ptr<ToyStateFunc> Ptr;
    ToyBSPProblemHelperPtr toy_helper;

    ToyStateFunc();
    ToyStateFunc(BSPProblemHelperBasePtr helper);
    virtual StateT operator()(const StateT& x, const ControlT& u, const StateNoiseT& m) const ;
  };

  class ToyObserveFunc : public ObserveFunc<StateT, ObserveT, ObserveNoiseT> {
  public:
    typedef boost::shared_ptr<ToyObserveFunc> Ptr;
    ToyBSPProblemHelperPtr toy_helper;

    ToyObserveFunc();
    ToyObserveFunc(BSPProblemHelperBasePtr helper);
    bool sgndist(const StateT& x, Vector1d* dists) const;
    virtual ObserveT observation_masks(const StateT& x, double approx_factor) const;
    virtual ObserveT operator()(const StateT& x, const ObserveNoiseT& n) const;
  };

  class ToyBeliefFunc : public BeliefFunc<ToyStateFunc, ToyObserveFunc, BeliefT> {
  public:
    typedef boost::shared_ptr<ToyBeliefFunc> Ptr;
    ToyBSPProblemHelperPtr toy_helper;
    ToyBeliefFunc();
    ToyBeliefFunc(BSPProblemHelperBasePtr helper, StateFuncPtr f, ObserveFuncPtr h);
  };

  class ToyBSPProblemHelper : public BSPProblemHelper<ToyBeliefFunc> {
  public:
    typedef typename BeliefConstraint<ToyBeliefFunc>::Ptr BeliefConstraintPtr;
    double input_dt;
    double goaleps;

    virtual void add_goal_constraint(OptProb& prob);
    ToyBSPProblemHelper();
  };

  class ToyPlotter : public BSPQtPlotter, public ProblemState {
  protected:
    ToyBSPProblemHelperPtr toy_helper;
    ToyPlotter(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper, QWidget* parent=NULL);
    void compute_distmap(QImage* distmap, StateT* state, double approx_factor);
  };

  class ToySimulationPlotter : public ToyPlotter {
    Q_OBJECT
  public:
    ToySimulationPlotter(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper, QWidget* parent=NULL);
    virtual void update_plot_data(void* data_x, void* data_sim);
    vector<StateT> simulated_positions;
  protected:
    QImage distmap;
    vector<StateT> states;
    vector<VarianceT> sigmas;
    vector<StateT> waypoints;
    virtual void paintEvent(QPaintEvent*);
  };

  class ToyOptPlotter : public ToyPlotter {
    Q_OBJECT
  public:
    ToyOptPlotter(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper, QWidget* parent=NULL);
  public slots:
    virtual void update_plot_data(void *data);
  public:
    boost::shared_ptr<ToySimulationPlotter> simplotptr;
  protected:
    double old_approx_factor, cur_approx_factor;
    QImage distmap;
    vector<StateT> states_opt, states_actual, states_01, states_waypoints;
    vector<StateT> samples;
    vector<VarianceT> sigmas_opt, sigmas_actual, sigmas_01;
    vector<StateT> openlooppts;
    bool fileinp;
    virtual void paintEvent(QPaintEvent*);
  };

  class ToyOptimizerTask : public BSPOptimizerTask {
    Q_OBJECT
  public:
    ToyOptimizerTask(QObject* parent=NULL);
    ToyOptimizerTask(int argc, char **argv, QObject* parent=NULL);
    virtual void run();
    void stage_plot_callback(boost::shared_ptr<ToyOptPlotter> plotter, OptProb*, DblVec& x);
  };

  class ToyBSPPlanner : public BSPPlanner<ToyBSPProblemHelper> {
  public:
    ToyBSPPlanner();
    virtual void custom_simulation_update(StateT* state, VarianceT* sigma, const StateT& actual_state);
    virtual void initialize_optimizer_parameters(BSPTrustRegionSQP& opt, bool is_first_time=true);
    virtual void initialize();
  };

  typedef boost::shared_ptr<ToyBSPPlanner> ToyBSPPlannerPtr;

}

template ToyBSP::ToyOptPlotter* BSPOptimizerTask::create_plotter<ToyBSP::ToyOptPlotter>(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper);
template ToyBSP::ToySimulationPlotter* BSPOptimizerTask::create_plotter<ToyBSP::ToySimulationPlotter>(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper);
