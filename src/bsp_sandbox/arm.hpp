#pragma once

#include "bsp/bsp.hpp"
#include <QtCore>
#include <QImage>
#include <qevent.h>

using namespace BSP;

namespace ArmBSP {

  // This will generate a bunch of types like StateT, ControlT, etc.
  BSP_TYPEDEFS(
    3, // state_dim
    3, // state_noise_dim
    3, // control_dim
    3, // observe_dim
    3, // observe_noise_dim
    6, // sigma_dof
    9 // belief_dim
  );

  typedef Matrix<double, 2, 3> TransT;

  class ArmBSPProblemHelper;
  typedef boost::shared_ptr<ArmBSPProblemHelper> ArmBSPProblemHelperPtr;

  class ArmStateFunc : public StateFunc<StateT, ControlT, StateNoiseT> {
  public:
    typedef boost::shared_ptr<ArmStateFunc> Ptr;
    ArmBSPProblemHelperPtr arm_helper;

    ArmStateFunc();
    ArmStateFunc(BSPProblemHelperBasePtr helper);
    virtual StateT operator()(const StateT& x, const ControlT& u, const StateNoiseT& m) const ;
  };

  class ArmObserveFunc : public ObserveFunc<StateT, ObserveT, ObserveNoiseT> {
  public:
    typedef boost::shared_ptr<ArmObserveFunc> Ptr;
    ArmBSPProblemHelperPtr arm_helper;
    ArmObserveFunc();
    ArmObserveFunc(BSPProblemHelperBasePtr helper);
    bool sgndist(const Vector2d& x, Vector2d* dists) const;
    ObserveMatT compute_gamma(const StateT& x, double approx_factor) const;
    ObserveMatT compute_inverse_gamma(const StateT& x, double approx_factor) const;
    virtual ObserveT operator()(const StateT& x, const ObserveNoiseT& n) const;
    virtual ObserveT operator()(const StateT& x, const ObserveNoiseT& n, double approx_factor) const;
  };

  class ArmBeliefFunc : public BeliefFunc<ArmStateFunc, ArmObserveFunc, BeliefT> {
  public:
    typedef boost::shared_ptr<ArmBeliefFunc> Ptr;
    ArmBSPProblemHelperPtr arm_helper;
    ArmBeliefFunc();
    ArmBeliefFunc(BSPProblemHelperBasePtr helper, StateFuncPtr f, ObserveFuncPtr h);
  };

  class ArmBSPProblemHelper : public BSPProblemHelper<ArmBeliefFunc> {
  public:
    typedef typename BeliefConstraint<ArmBeliefFunc>::Ptr BeliefConstraintPtr;
    double input_dt;
    Vector3d base_config;
    Vector2d goal_pos;
    ArmBSPProblemHelper();
    virtual void initialize();
    TransT angle_to_transform(const StateT& angle) const;
    Vector8d angle_to_endpoints(const StateT& angle) const;
    Vector2d angle_to_pos(const StateT& angle) const;
    StateT pos_to_angle(const Vector2d& pos) const;

  };

  class ArmOptPlotter : public BSPQtPlotter, public ProblemState {
    Q_OBJECT
  public:
    ArmOptPlotter(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper, QWidget* parent=NULL);
  public slots:
    virtual void update_plot_data(void* data);
  protected:
    double old_approx_factor, cur_approx_factor;
    ArmBSPProblemHelperPtr arm_helper;
    QImage distmap;
    vector<VectorXd> states;
    vector<MatrixXd> sigmas;
    virtual void paintEvent(QPaintEvent*);
    virtual void keyPressEvent(QKeyEvent*);
  };

  class ArmSimulationPlotter : public BSPQtPlotter, public ProblemState {
    Q_OBJECT
  public:
    ArmSimulationPlotter(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper, QWidget* parent=NULL);
    virtual void update_plot_data(void* data_x, void* data_sim);
  protected:
    ArmBSPProblemHelperPtr arm_helper;
    QImage distmap;
    vector<VectorXd> states;
    vector<MatrixXd> sigmas;
    vector<StateT> simulated_positions;
    virtual void paintEvent(QPaintEvent*);
  };

  class ArmOptimizerTask : public BSPOptimizerTask {
    Q_OBJECT
  public:
    ArmOptimizerTask(QObject* parent=NULL);
    ArmOptimizerTask(int argc, char **argv, QObject* parent=NULL);
    virtual void run();
    void stage_plot_callback(boost::shared_ptr<ArmOptPlotter> plotter, OptProb*, DblVec& x);
  };
}

template ArmBSP::ArmOptPlotter* BSPOptimizerTask::create_plotter<ArmBSP::ArmOptPlotter>(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper);
template ArmBSP::ArmSimulationPlotter* BSPOptimizerTask::create_plotter<ArmBSP::ArmSimulationPlotter>(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper);
