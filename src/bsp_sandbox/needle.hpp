#pragma once

#include "bsp/bsp.hpp"
#include "geometry_2d.hpp"
#include <QtCore>
#include <QImage>
#include <qevent.h>

using namespace BSP;
using namespace Geometry2D;

namespace NeedleBSP {

  // This will generate a bunch of types like StateT, ControlT, etc.
  BSP_TYPEDEFS(
      4, // state_dim
      4, // state_noise_dim
      3, // control_dim
      3, // observe_dim
      3, // observe_noise_dim
      10, // sigma_dof
      14 // belief_dim
  );

  typedef Matrix<double, 1, 1> Vector1d;
  typedef Matrix<double, 1, 1> Matrix1d;

  // state: { needle_x, needle_y, needle_angle, ultrasound_x }
  // control: { needle_theta, needle_velocity, ultrasound_velocity }
  // observation: { needle_x, needle_y, ultrasound_x }

  class NeedleBSPProblemHelper;
  typedef boost::shared_ptr<NeedleBSPProblemHelper> NeedleBSPProblemHelperPtr;

  class NeedleStateFunc : public StateFunc<StateT, ControlT, StateNoiseT> {
  public:
    typedef boost::shared_ptr<NeedleStateFunc> Ptr;
    NeedleBSPProblemHelperPtr needle_helper;

    NeedleStateFunc();
    NeedleStateFunc(BSPProblemHelperBasePtr helper);
    virtual StateT operator()(const StateT& x, const ControlT& u, const StateNoiseT& m) const ;
  };

  class NeedleObserveFunc : public ObserveFunc<StateT, ObserveT, ObserveNoiseT> {
  public:
    typedef boost::shared_ptr<NeedleObserveFunc> Ptr;
    NeedleBSPProblemHelperPtr needle_helper;

    NeedleObserveFunc();
    NeedleObserveFunc(BSPProblemHelperBasePtr helper);
    bool sgndist(double x, double ultrasound_x, Vector1d* dists) const;
    ObserveMatT compute_gamma(const StateT& x, double approx_factor) const;
    ObserveMatT compute_inverse_gamma(const StateT& x, double approx_factor) const;
    ObserveT unnoisy_observation(const StateT& x) const;
    virtual ObserveT operator()(const StateT& x, const ObserveNoiseT& n) const;
    virtual ObserveT operator()(const StateT& x, const ObserveNoiseT& n, double approx_factor) const;
  };

  class NeedleBeliefFunc : public BeliefFunc<NeedleStateFunc, NeedleObserveFunc, BeliefT> {
  public:
    typedef boost::shared_ptr<NeedleBeliefFunc> Ptr;
    NeedleBSPProblemHelperPtr needle_helper;
    NeedleBeliefFunc();
    NeedleBeliefFunc(BSPProblemHelperBasePtr helper, StateFuncPtr f, ObserveFuncPtr h);
  };

  typedef Vector3d RobotStateT;
  typedef Vector2d RobotControlT;

  class NeedleBSPProblemHelper : public BSPProblemHelper<NeedleBeliefFunc> {
  public:
    struct RRTNode {
      RobotStateT x;
      RobotControlT u;
      int bp;
    };

    typedef typename BeliefConstraint<NeedleBeliefFunc>::Ptr BeliefConstraintPtr;
    double input_dt;
    double needlelen;
    double goaleps;

    double ultrasound_span;

    double robot_state_dim;
    double robot_control_dim;

    vector<Vector4d> rrt_edges;

    virtual void RRTplan(bool);
    virtual void add_goal_constraint(OptProb& prob);
    virtual void configure_problem(OptProb& prob);
    void add_ultrasound_deviation_cost(OptProb& prob);
    NeedleBSPProblemHelper();
  };

  class NeedlePlotter : public BSPQtPlotter, public ProblemState {
  protected:
    NeedleBSPProblemHelperPtr needle_helper;
    NeedlePlotter(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper, QWidget* parent=NULL);
    void compute_distmap(QImage* distmap, StateT* state, double approx_factor);
  };

  class NeedleOptPlotter : public NeedlePlotter {
    Q_OBJECT
  public:
    NeedleOptPlotter(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper, QWidget* parent=NULL);
  public slots:
    virtual void update_plot_data(void *data);
  protected:
    double old_approx_factor, cur_approx_factor;
    QImage distmap;
    vector<StateT> states_opt, states_actual;
    vector<VarianceT> sigmas_opt, sigmas_actual;
    virtual void paintEvent(QPaintEvent*);
  };

  class NeedleSimulationPlotter : public NeedlePlotter {
    Q_OBJECT
  public:
    NeedleSimulationPlotter(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper, QWidget* parent=NULL);
    virtual void update_plot_data(void* data_x, void* data_sim);
  protected:
    QImage distmap;
    vector<StateT> states;
    vector<VarianceT> sigmas;
    vector<StateT> simulated_positions;
    virtual void paintEvent(QPaintEvent*);
  };

  class NeedleOptimizerTask : public BSPOptimizerTask {
    Q_OBJECT
  public:
    NeedleOptimizerTask(QObject* parent=NULL);
    NeedleOptimizerTask(int argc, char **argv, QObject* parent=NULL);
    virtual void run();
    void stage_plot_callback(boost::shared_ptr<NeedleOptPlotter> plotter, OptProb*, DblVec& x);
  };

  class NeedleBSPPlanner : public BSPPlanner<NeedleBSPProblemHelper> {
  public:
    virtual void initialize();
    NeedleBSPPlanner();
    virtual void custom_simulation_update(StateT* state, VarianceT* sigma, const StateT& actual_state);
  };

  typedef boost::shared_ptr<NeedleBSPPlanner> NeedleBSPPlannerPtr;

  class UltrasoundDeviationCost : public Cost {
  public:
    UltrasoundDeviationCost(const VarArray& vars, double coeff);
    virtual double value(const vector<double>& xvec);
    virtual ConvexObjectivePtr convex(const vector<double>& xvec, Model* model);
  private:
    VarArray vars;
    double coeff;
  };

}

template NeedleBSP::NeedleOptPlotter* BSPOptimizerTask::create_plotter<NeedleBSP::NeedleOptPlotter>(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper);
template NeedleBSP::NeedleSimulationPlotter* BSPOptimizerTask::create_plotter<NeedleBSP::NeedleSimulationPlotter>(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper);
