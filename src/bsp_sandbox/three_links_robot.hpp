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

namespace ThreeLinksRobotBSP {

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

  // state: { {robot_dofs} }
  // control: { {d_robot_dofs} }
  // observation: { (whatever) }

  class ThreeLinksRobotBSPProblemHelper;
  typedef boost::shared_ptr<ThreeLinksRobotBSPProblemHelper> ThreeLinksRobotBSPProblemHelperPtr;

  class ThreeLinksRobotStateFunc : public StateFunc<StateT, ControlT, StateNoiseT> {
  public:
    typedef boost::shared_ptr<ThreeLinksRobotStateFunc> Ptr;
    ThreeLinksRobotBSPProblemHelperPtr three_links_robot_helper;

    ThreeLinksRobotStateFunc();
    ThreeLinksRobotStateFunc(BSPProblemHelperBasePtr helper);
    StateT operator()(const StateT& x, const ControlT& u, const StateNoiseT& m) const ;
  };

  class ThreeLinksRobotObserveFunc : public ObserveFunc<StateT, ObserveT, ObserveNoiseT> {
  public:
    typedef boost::shared_ptr<ThreeLinksRobotObserveFunc> Ptr;
    ThreeLinksRobotBSPProblemHelperPtr three_links_robot_helper;
    ThreeLinksRobotObserveFunc();
    ThreeLinksRobotObserveFunc(BSPProblemHelperBasePtr helper);
    ObserveT operator()(const StateT& x, const ObserveNoiseT& n) const;
  };

  class ThreeLinksRobotBeliefFunc : public EkfBeliefFunc<ThreeLinksRobotStateFunc, ThreeLinksRobotObserveFunc, BeliefT> {
  public:
    typedef boost::shared_ptr<ThreeLinksRobotBeliefFunc> Ptr;
    ThreeLinksRobotBSPProblemHelperPtr three_links_robot_helper;
    ThreeLinksRobotBeliefFunc();
    ThreeLinksRobotBeliefFunc(BSPProblemHelperBasePtr helper, StateFuncPtr f, ObserveFuncPtr h);
  };

  class ThreeLinksRobotBSPProblemHelper : public BSPProblemHelper<ThreeLinksRobotBeliefFunc> {
  public:
    typedef typename BeliefConstraint<ThreeLinksRobotBeliefFunc>::Ptr BeliefConstraintPtr;
    void add_goal_constraint(OptProb& prob);
    void add_collision_term(OptProb& prob);
    void configure_problem(OptProb& prob);
    ThreeLinksRobotBSPProblemHelper();

    RobotBasePtr robot;
    RobotAndDOFPtr rad;
    KinBody::LinkPtr link;
    Matrix4d goal_trans;
  };

  class ThreeLinksRobotOptimizerTask : public BSPOptimizerTask,
                                    public OpenRAVEPlotterMixin<ThreeLinksRobotBSPProblemHelper> {
  public:
    ThreeLinksRobotOptimizerTask(QObject* parent=nullptr);
    ThreeLinksRobotOptimizerTask(int argc, char **argv, QObject* parent=nullptr);
    void run();
  };

  class ThreeLinksRobotBSPPlanner : public BSPPlanner<ThreeLinksRobotBSPProblemHelper> {
  public:
    void initialize();
    ThreeLinksRobotBSPPlanner();
    void initialize_optimizer_parameters(BSPTrustRegionSQP& opt, bool is_first_time=true);
    RobotBasePtr robot;
    RobotAndDOFPtr rad;
    Matrix4d goal_trans;
    KinBody::LinkPtr link;
  };

  typedef boost::shared_ptr<ThreeLinksRobotBSPPlanner> ThreeLinksRobotBSPPlannerPtr;


}
