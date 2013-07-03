#pragma once

#include "bsp/bsp.hpp"
#include <QtCore>
#include <QImage>
#include <qevent.h>

using namespace BSP;

namespace CarBSP {

  // This will generate a bunch of types like StateT, ControlT, etc.
  BSP_TYPEDEFS(
      4, // state_dim
      4, // state_noise_dim
      2, // control_dim
      4, // observe_dim
      4, // observe_noise_dim
      10, // sigma_dof
      14 // belief_dim
  );

  // state: { x, y, angle, velocity }
  // control: { theta, acceleration }

  class CarBSPProblemHelper;
  typedef boost::shared_ptr<CarBSPProblemHelper> CarBSPProblemHelperPtr;

  class CarStateFunc : public StateFunc<StateT, ControlT, StateNoiseT> {
  public:
    typedef boost::shared_ptr<CarStateFunc> Ptr;
    CarBSPProblemHelperPtr car_helper;

    CarStateFunc();
    CarStateFunc(BSPProblemHelperBasePtr helper);
    virtual StateT operator()(const StateT& x, const ControlT& u, const StateNoiseT& m) const ;
  };

  class CarObserveFunc : public ObserveFunc<StateT, ObserveT, ObserveNoiseT> {
  public:
    typedef boost::shared_ptr<CarObserveFunc> Ptr;
    CarBSPProblemHelperPtr car_helper;

    CarObserveFunc();
    CarObserveFunc(BSPProblemHelperBasePtr helper);
    bool sgndist(const Vector2d& x, Vector2d* dists) const;
    ObserveMatT compute_gamma(const StateT& x, double approx_factor) const;
    ObserveMatT compute_inverse_gamma(const StateT& x, double approx_factor) const;
    virtual ObserveT operator()(const StateT& x, const ObserveNoiseT& n) const;
    virtual ObserveT operator()(const StateT& x, const ObserveNoiseT& n, double approx_factor) const;
  };

  class CarBeliefFunc : public BeliefFunc<CarStateFunc, CarObserveFunc, BeliefT> {
  public:
    typedef boost::shared_ptr<CarBeliefFunc> Ptr;
    CarBSPProblemHelperPtr car_helper;
    CarBeliefFunc();
    CarBeliefFunc(BSPProblemHelperBasePtr helper, StateFuncPtr f, ObserveFuncPtr h);
  };

  class CarBSPProblemHelper : public BSPProblemHelper<CarBeliefFunc> {
  public:
    struct RRTNode {
      StateT x;
      ControlT u;
      int bp;
    };

    typedef typename BeliefConstraint<CarBeliefFunc>::Ptr BeliefConstraintPtr;
    double input_dt;
    double carlen;
    double goaleps;

    vector<Vector4d> rrt_edges;

    virtual void RRTplan();
    virtual void add_goal_constraint(OptProb& prob);
    CarBSPProblemHelper();
  };

  class CarPlotter : public BSPQtPlotter, public ProblemState {
  protected:
    CarBSPProblemHelperPtr car_helper;
    CarPlotter(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper, QWidget* parent=NULL);
    void compute_distmap(QImage* distmap, StateT* current_state, double approx_factor);
  };

  class CarOptPlotter : public CarPlotter {
    Q_OBJECT
  public:
    CarOptPlotter(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper, QWidget* parent=NULL);
  public slots:
    virtual void update_plot_data(void *data);
  protected:
    double old_approx_factor, cur_approx_factor;
    QImage distmap;
    vector<StateT> states_opt, states_actual;
    vector<VarianceT> sigmas_opt, sigmas_actual;
    virtual void paintEvent(QPaintEvent*);
  };

  class CarSimulationPlotter : public CarPlotter {
    Q_OBJECT
  public:
    CarSimulationPlotter(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper, QWidget* parent=NULL);
    virtual void update_plot_data(void* data_x, void* data_sim);
  protected:
    QImage distmap;
    vector<StateT> states;
    vector<VarianceT> sigmas;
    vector<StateT> simulated_positions;
    virtual void paintEvent(QPaintEvent*);
  };

  class CarOptimizerTask : public BSPOptimizerTask {
    Q_OBJECT
  public:
    CarOptimizerTask(QObject* parent=NULL);
    CarOptimizerTask(int argc, char **argv, QObject* parent=NULL);
    virtual void run();
    void stage_plot_callback(boost::shared_ptr<CarOptPlotter> plotter, OptProb*, DblVec& x);
  };

  class CarBSPPlanner : public BSPPlanner<CarBSPProblemHelper> {
  public:
    virtual void initialize();
  };

  typedef boost::shared_ptr<CarBSPPlanner> CarBSPPlannerPtr;


}

template CarBSP::CarOptPlotter* BSPOptimizerTask::create_plotter<CarBSP::CarOptPlotter>(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper);
template CarBSP::CarSimulationPlotter* BSPOptimizerTask::create_plotter<CarBSP::CarSimulationPlotter>(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper);
