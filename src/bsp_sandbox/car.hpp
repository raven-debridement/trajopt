#pragma once

#include "bsp/bsp.hpp"
#include "geometry_2d.hpp"
#include <QtCore>
#include <QImage>
#include <qevent.h>

using namespace BSP;
using namespace Geometry2D;

namespace CarBSP {

  // This will generate a bunch of types like StateT, ControlT, etc.
  BSP_TYPEDEFS(
      7, // state_dim
      7, // state_noise_dim
      6, // control_dim
      4, // observe_dim
      4, // observe_noise_dim
      28, // sigma_dof
      35 // belief_dim
  );

  typedef Matrix<double, 7, 1> Vector7d;
  typedef Matrix<double, 7, 7> Matrix7d;

  // state: { x, y, angle, landmark1_x, landmark1_y, landmark2_x, landmark2_y }
  // control: { theta, velocity, landmark1_dx, landmark1_dy, landmark2_dx, landmark2_dy }
  // observation: { rel_landmark1_dis, rel_landmark1_angle, rel_landmark2_dis, rel_landmark2_angle }

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
    bool sgndist(const StateT& x, Vector2d* dists) const;
    ObserveT observation_masks(const StateT& x, double approx_factor=-1) const;
    virtual ObserveT operator()(const StateT& x, const ObserveNoiseT& n) const;
  };

  class CarBeliefFunc : public BeliefFunc<CarStateFunc, CarObserveFunc, BeliefT> {
  public:
    typedef boost::shared_ptr<CarBeliefFunc> Ptr;
    CarBSPProblemHelperPtr car_helper;
    CarBeliefFunc();
    CarBeliefFunc(BSPProblemHelperBasePtr helper, StateFuncPtr f, ObserveFuncPtr h);
  };

  typedef Vector3d RobotStateT;
  typedef Vector2d RobotControlT;

  class CarBSPProblemHelper : public BSPProblemHelper<CarBeliefFunc> {
  public:
    //struct RRTNode {
    //  RobotStateT x;
    //  RobotControlT u;
    //  int bp;
    //};

    typedef typename BeliefConstraint<CarBeliefFunc>::Ptr BeliefConstraintPtr;
    double input_dt;
    double carlen;
    double goaleps;

    double car_camera_depth;
    double car_camera_span_angle;

    double robot_state_dim;
    double robot_control_dim;

    vector<Vector4d> rrt_edges;

    //virtual void RRTplan(bool);
    virtual void add_goal_constraint(OptProb& prob);
    CarBSPProblemHelper();
  };

  class CarPlotter : public BSPQtPlotter, public ProblemState {
  protected:
    CarBSPProblemHelperPtr car_helper;
    CarPlotter(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper, QWidget* parent=NULL);
    void compute_distmap(QImage* distmap, StateT* state, double approx_factor);
    void draw_beam_2d(const Beam2D& beam, QPainter& painter);
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
    vector<StateT> states_opt, states_actual, states_waypoints;
    vector<StateT> samples;
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
    vector<StateT> waypoints;
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
    CarBSPPlanner();
    virtual void custom_simulation_update(StateT* state, VarianceT* sigma);
    virtual void initialize_optimizer_parameters(BSPTrustRegionSQP& opt);
  };

  typedef boost::shared_ptr<CarBSPPlanner> CarBSPPlannerPtr;


}

template CarBSP::CarOptPlotter* BSPOptimizerTask::create_plotter<CarBSP::CarOptPlotter>(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper);
template CarBSP::CarSimulationPlotter* BSPOptimizerTask::create_plotter<CarBSP::CarSimulationPlotter>(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper);
