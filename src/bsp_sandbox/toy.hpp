#pragma once

#include "bsp/bsp.hpp"
#include <QtCore>
#include <QImage>
#include <qevent.h>

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
    bool sgndist(const Vector2d& x, Vector2d* dists) const;
    ObserveMatT compute_gamma(const StateT& x, double approx_factor) const;
    ObserveMatT compute_inverse_gamma(const StateT& x, double approx_factor) const;
    virtual ObserveT operator()(const StateT& x, const ObserveNoiseT& n) const;
    virtual ObserveT operator()(const StateT& x, const ObserveNoiseT& n, double approx_factor) const;
  };

  class ToyBeliefFunc : public BeliefFunc<ToyStateFunc, ToyObserveFunc, BeliefT> {
  public:
    typedef boost::shared_ptr<ToyBeliefFunc> Ptr;
    ToyBSPProblemHelperPtr toy_helper;
    ToyBeliefFunc();
    ToyBeliefFunc(BSPProblemHelperBasePtr helper, StateFuncPtr f, ObserveFuncPtr h);
//<<<<<<< HEAD
//=======
//    bool sgndist(const Vector2d& x, Vector2d* dists) const;
//    virtual ObserveStateGradT sensor_constrained_observe_state_gradient(const ObserveStateGradT& H, const StateT& x) const;
//    virtual ObserveNoiseGradT sensor_constrained_observe_noise_gradient(const ObserveNoiseGradT& N, const StateT& x) const;
//    virtual VarianceT sensor_constrained_variance_reduction(const VarianceT& reduction, const StateT& x) const;
//    virtual ObserveMatT compute_gamma(const StateT& x) const;
//    virtual ObserveMatT compute_inverse_gamma(const StateT& x) const;
//>>>>>>> bsp
  };

  class ToyBSPProblemHelper : public BSPProblemHelper<ToyBeliefFunc> {
  public:
    typedef typename BeliefConstraint<ToyBeliefFunc>::Ptr BeliefConstraintPtr;
    double input_dt;
    ToyBSPProblemHelper();
    virtual void initialize();
  };

  class ToyOptPlotter : public BSPQtPlotter, public ProblemState {
    Q_OBJECT
  public:
    ToyOptPlotter(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper, QWidget* parent=NULL);
  public slots:
    virtual void update_plot_data(void* data);
  protected:
    double old_approx_factor, cur_approx_factor;
    ToyBSPProblemHelperPtr toy_helper;
    QImage distmap;
//    vector<VectorXd> states;
//    vector<MatrixXd> sigmas;
    vector<VectorXd> states_opt, states_actual;
    vector<MatrixXd> sigmas_opt, sigmas_actual;
    virtual void paintEvent(QPaintEvent*);
    virtual void keyPressEvent(QKeyEvent*);
  };

  class ToySimulationPlotter : public BSPQtPlotter, public ProblemState {
    Q_OBJECT
  public:
    ToySimulationPlotter(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper, QWidget* parent=NULL);
  //public slots:
    virtual void update_plot_data(void* data_x, void* data_sim);//DblVec& x, vector<StateT>& sim);
  protected:
    ToyBSPProblemHelperPtr toy_helper;
    QImage distmap;
//<<<<<<< HEAD
    vector<VectorXd> states;
    vector<MatrixXd> sigmas;
    vector<StateT> simulated_positions;
//=======

//>>>>>>> bsp
    virtual void paintEvent(QPaintEvent*);
  };

  class ToyOptimizerTask : public BSPOptimizerTask {
    Q_OBJECT
  public:
    ToyOptimizerTask(QObject* parent=NULL);
    ToyOptimizerTask(int argc, char **argv, QObject* parent=NULL);
    virtual void run();
    void stage_plot_callback(boost::shared_ptr<ToyOptPlotter> plotter, OptProb*, DblVec& x);
  //signals:
    //void replot_signal(DblVec& x, vector<StateT>& sim);
  };
}

template ToyBSP::ToyOptPlotter* BSPOptimizerTask::create_plotter<ToyBSP::ToyOptPlotter>(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper);
template ToyBSP::ToySimulationPlotter* BSPOptimizerTask::create_plotter<ToyBSP::ToySimulationPlotter>(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper);
