#pragma once

#include "bsp/bsp.hpp"
#include "bsp/openrave_utils.hpp"
#include "trajopt/kinematic_terms.hpp"
#include "trajopt/problem_description.hpp"
#include "geometry_2d.hpp"
#include <QtCore>
#include <QImage>
#include <qevent.h>
#include <openrave-core.h>
#include <openrave/openrave.h>
#include <array>

using namespace BSP;
using namespace Geometry2D;

namespace PointRobotBSP {

  // This will generate a bunch of types like StateT, ControlT, etc.
  BSP_TYPEDEFS(
      7, // state_dim
      7, // state_noise_dim
      7, // control_dim
      3, // observe_dim
      3, // observe_noise_dim
      28, // sigma_dof
      35 // belief_dim
  );

  typedef Matrix<double, 7, 1> Vector7d;
  typedef Matrix<double, 7, 7> Matrix7d;

  // state: { {robot_dofs} }
  // control: { {d_robot_dofs} }
  // observation: { (whatever) }

  class PointRobotBSPProblemHelper;
  typedef boost::shared_ptr<PointRobotBSPProblemHelper> PointRobotBSPProblemHelperPtr;

  class PointRobotStateFunc : public StateFunc<StateT, ControlT, StateNoiseT> {
  public:
    typedef boost::shared_ptr<PointRobotStateFunc> Ptr;
    PointRobotBSPProblemHelperPtr point_robot_helper;

    PointRobotStateFunc();
    PointRobotStateFunc(BSPProblemHelperBasePtr helper);
    StateT operator()(const StateT& x, const ControlT& u, const StateNoiseT& m) const ;
  };

  class PointRobotObserveFunc : public ObserveFunc<StateT, ObserveT, ObserveNoiseT> {
  public:
    typedef boost::shared_ptr<PointRobotObserveFunc> Ptr;
    PointRobotBSPProblemHelperPtr point_robot_helper;
    PointRobotObserveFunc();
    PointRobotObserveFunc(BSPProblemHelperBasePtr helper);
    ObserveT operator()(const StateT& x, const ObserveNoiseT& n) const;
  };

  class PointRobotBeliefFunc : public EkfBeliefFunc<PointRobotStateFunc, PointRobotObserveFunc, BeliefT> {
  public:
    typedef boost::shared_ptr<PointRobotBeliefFunc> Ptr;
    PointRobotBSPProblemHelperPtr point_robot_helper;
    PointRobotBeliefFunc();
    PointRobotBeliefFunc(BSPProblemHelperBasePtr helper, StateFuncPtr f, ObserveFuncPtr h);
  };

  class PointRobotBSPProblemHelper : public BSPProblemHelper<PointRobotBeliefFunc> {
  public:
    typedef typename BeliefConstraint<PointRobotBeliefFunc>::Ptr BeliefConstraintPtr;
    void add_goal_constraint(OptProb& prob);
    void add_collision_term(OptProb& prob);
    void configure_problem(OptProb& prob);
    PointRobotBSPProblemHelper();

    RobotBasePtr robot;
    RobotAndDOFPtr rad;
    KinBody::LinkPtr link;
    Matrix4d goal_trans;
  };

  class PointRobotOptimizerTask : public BSPOptimizerTask,
                                    public OpenRAVEPlotterMixin<PointRobotBSPProblemHelper> {
  public:
    PointRobotOptimizerTask(QObject* parent=nullptr);
    PointRobotOptimizerTask(int argc, char **argv, QObject* parent=nullptr);
    void run();
  };

  class PointRobotBSPPlanner : public BSPPlanner<PointRobotBSPProblemHelper> {
  public:
    void initialize();
    PointRobotBSPPlanner();
    void initialize_optimizer_parameters(BSPTrustRegionSQP& opt, bool is_first_time=true);
    RobotBasePtr robot;
    RobotAndDOFPtr rad;
    Matrix4d goal_trans;
    KinBody::LinkPtr link;
  };

  typedef boost::shared_ptr<PointRobotBSPPlanner> PointRobotBSPPlannerPtr;


}
