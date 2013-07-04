#pragma once

#include "bsp/bsp.hpp"
#include "geometry_2d.hpp"
#include <QApplication>
#include <QtCore>
#include <QtCore>
#include <QImage>
#include <qevent.h>

using namespace BSP;
using namespace Geometry2D;

namespace ArmBSP {

  // This will generate a bunch of types like StateT, ControlT, etc.
  BSP_TYPEDEFS(
    5, // state_dim
    5, // state_noise_dim
    3, // control_dim
    5, // observe_dim
    5, // observe_noise_dim
    15, // sigma_dof
    20 // belief_dim
  );

  typedef Matrix<double, 8, 1> Vector8d;
  typedef Matrix<double, 5, 1> Vector5d;
  typedef Matrix<double, 5, 5> Matrix5d;
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
    virtual ObserveT real_observation(const StateT& x, const ObserveNoiseT& n) const;
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
    Vector2d real_object_pos;
    Vector3d link_lengths;
    Beam2D camera;
    vector<Beam2D> partitioned_fov;
    vector<Beam2D> cur_fov;
    int n_fov_parts;
    ArmBSPProblemHelper();
    virtual void initialize();
    virtual void add_goal_constraint(OptProb& prob);
    TransT angle_to_transform(const Vector3d& angle) const;
    Vector8d angle_to_endpoints(const Vector3d& angle) const;
    Vector2d angle_to_pos(const Vector3d& angle) const;
    vector<Segment> angle_to_segments(const Vector3d& angle) const;
    StateT pos_to_angle(const Vector2d& pos) const;

  };

  class ArmPlotter : public BSPQtPlotter, public ProblemState {
  protected:
    void draw_robot(const Vector3d& x, QPainter& painter);
    ArmBSPProblemHelperPtr arm_helper;
    ArmPlotter(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper, QWidget* parent=NULL);
    void compute_distmap(QImage* distmap, StateT* current_state, double approx_factor);
  };

  class ArmOptPlotter : public ArmPlotter {
    Q_OBJECT
  public:
    ArmOptPlotter(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper, QWidget* parent=NULL);
  public slots:
    virtual void update_plot_data(void* data);
  protected:
    double old_approx_factor, cur_approx_factor;
    QImage distmap;
    vector<StateT> states;
    vector<VarianceT> sigmas;
    virtual void paintEvent(QPaintEvent*);
  };

  class ArmSimulationPlotter : public ArmPlotter {
    Q_OBJECT
  public:
    ArmSimulationPlotter(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper, QWidget* parent=NULL);
    virtual void update_plot_data(void* data_x, void* data_sim);
  protected:
    QImage distmap;
    vector<StateT> states;
    vector<VarianceT> sigmas;
    vector<StateT> simulated_positions;
    vector<Vector2d> object_position_beliefs;
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

  class ArmGoalError : public VectorOfVector {
  public:
    ArmBSPProblemHelperPtr helper;
    ArmGoalError(ArmBSPProblemHelperPtr helper);
    VectorXd operator()(const VectorXd& a) const;
  };

  class ArmBSPPlanner : public BSPPlanner<ArmBSPProblemHelper> {
  public:
    Vector2d goal_pos;
    Vector3d link_lengths;
    Vector3d base_config;
    Vector2d real_object_pos;
    Beam2D camera;
    int n_fov_parts;
    virtual void initialize();
    virtual void custom_simulation_update(StateT* state, VarianceT* sigma);
    ArmBSPPlanner();
  };

  typedef boost::shared_ptr<ArmBSPPlanner> ArmBSPPlannerPtr;

}

template ArmBSP::ArmOptPlotter* BSPOptimizerTask::create_plotter<ArmBSP::ArmOptPlotter>(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper);
template ArmBSP::ArmSimulationPlotter* BSPOptimizerTask::create_plotter<ArmBSP::ArmSimulationPlotter>(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper);
