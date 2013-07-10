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
    6, // state_dim
    6, // state_noise_dim
    4, // control_dim
    6, // observe_dim
    6, // observe_noise_dim
    21, // sigma_dof
    27 // belief_dim
  );

  // state: { joint_1, joint_2, joint_3, camera_theta, object_x, object_y }
  // control: { djoint_1, djoint_2, djoint_3, camera_dtheta }
  // observe: { joint_1, joint_2, joint_3, camera_theta, object_x, object_y }

  typedef Matrix<double, 8, 1> Vector8d;
  typedef Matrix<double, 6, 1> Vector6d;
  typedef Matrix<double, 6, 6> Matrix6d;
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
    virtual ObserveT operator()(const StateT& x, const ObserveNoiseT& n) const;
    virtual ObserveT real_observation(const StateT& x, const ObserveNoiseT& n) const;
    virtual ObserveT observation_masks(const StateT& x, double approx_factor=-1) const;
    virtual ObserveT real_observation_masks(const StateT& x) const;
    ObserveT observe_masks_from_object_position(const StateT& x, const Vector2d& object_pos, double approx_factor) const;
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
    double link_total_lengths;
    Vector2d camera_base;
    double camera_depth;
    double camera_span_angle;
    int n_fov_parts;
    ArmBSPProblemHelper();
    virtual void initialize();
    virtual void add_goal_constraint(OptProb& prob);
    TransT angle_to_transform(const Vector3d& angle) const;
    Vector8d angle_to_endpoints(const Vector3d& angle) const;
    Vector2d angle_to_pos(const Vector3d& angle) const;
    vector<Segment> angle_to_segments(const Vector3d& angle) const;
    StateT pos_to_angle(const Vector2d& pos) const;
    void fov_from_state(const StateT& x, vector<Beam2D>* out_fov) const;
  };

  class ArmPlotter : public BSPQtPlotter, public ProblemState {
  protected:
    void draw_robot(const Vector3d& x, QPainter& painter);
    ArmBSPProblemHelperPtr arm_helper;
    ArmPlotter(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper, QWidget* parent=NULL);
    void compute_distmap(QImage* distmap, StateT* current_state, double approx_factor);
    void draw_beam_2d(const Beam2D& beam, QPainter& painter);
    void draw_beams(const vector<Beam2D>& beams, QPainter& painter);
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
    vector<Matrix2d> object_position_sigmas;
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
    Vector2d camera_base;
    double camera_depth;
    double camera_span_angle;
    double link_total_lengths;
    //Beam2D camera;
    int n_fov_parts;
    virtual void initialize();
    virtual void custom_simulation_update(StateT* state, VarianceT* sigma, const StateT& actual_state);
    virtual void initialize_optimizer_parameters(BSPTrustRegionSQP& opt, bool is_first_time=true);
    Vector2d get_feasible_pos(const Vector2d& pos) const;
    ArmBSPPlanner();
  };

  typedef boost::shared_ptr<ArmBSPPlanner> ArmBSPPlannerPtr;

}

template ArmBSP::ArmOptPlotter* BSPOptimizerTask::create_plotter<ArmBSP::ArmOptPlotter>(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper);
template ArmBSP::ArmSimulationPlotter* BSPOptimizerTask::create_plotter<ArmBSP::ArmSimulationPlotter>(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper);
