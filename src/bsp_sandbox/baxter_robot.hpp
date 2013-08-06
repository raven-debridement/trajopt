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

namespace BaxterRobotBSP {

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

  class BaxterRobotBSPProblemHelper;
  typedef boost::shared_ptr<BaxterRobotBSPProblemHelper> BaxterRobotBSPProblemHelperPtr;

  class BaxterRobotStateFunc : public StateFunc<StateT, ControlT, StateNoiseT> {
  public:
    typedef boost::shared_ptr<BaxterRobotStateFunc> Ptr;
    BaxterRobotBSPProblemHelperPtr baxter_robot_helper;

    BaxterRobotStateFunc();
    BaxterRobotStateFunc(BSPProblemHelperBasePtr helper);
    StateT operator()(const StateT& x, const ControlT& u, const StateNoiseT& m) const ;
  };

  class BaxterRobotObserveFunc : public ObserveFunc<StateT, ObserveT, ObserveNoiseT> {
  public:
    typedef boost::shared_ptr<BaxterRobotObserveFunc> Ptr;
    BaxterRobotBSPProblemHelperPtr baxter_robot_helper;
    BaxterRobotObserveFunc();
    BaxterRobotObserveFunc(BSPProblemHelperBasePtr helper);
    ObserveT operator()(const StateT& x, const ObserveNoiseT& n) const;
  };

  class BaxterRobotBeliefFunc : public EkfBeliefFunc<BaxterRobotStateFunc, BaxterRobotObserveFunc, BeliefT> {
  public:
    typedef boost::shared_ptr<BaxterRobotBeliefFunc> Ptr;
    BaxterRobotBSPProblemHelperPtr baxter_robot_helper;
    BaxterRobotBeliefFunc();
    BaxterRobotBeliefFunc(BSPProblemHelperBasePtr helper, StateFuncPtr f, ObserveFuncPtr h);
  };

  class BaxterRobotBSPProblemHelper : public BSPProblemHelper<BaxterRobotBeliefFunc> {
  public:
    typedef typename BeliefConstraint<BaxterRobotBeliefFunc>::Ptr BeliefConstraintPtr;
    void add_goal_constraint(OptProb& prob);
    void add_collision_term(OptProb& prob);
    void configure_problem(OptProb& prob);
    BaxterRobotBSPProblemHelper();

    RobotBasePtr robot;
    RobotAndDOFPtr rad;
    KinBody::LinkPtr link;
    Matrix4d goal_trans;
  };

  class BaxterRobotBSPPlanner : public BSPPlanner<BaxterRobotBSPProblemHelper> {
  public:
    void initialize();
    BaxterRobotBSPPlanner();
    void initialize_optimizer_parameters(BSPTrustRegionSQP& opt, bool is_first_time=true);
    RobotBasePtr robot;
    RobotAndDOFPtr rad;
    Matrix4d goal_trans;
    KinBody::LinkPtr link;
  };

  typedef boost::shared_ptr<BaxterRobotBSPPlanner> BaxterRobotBSPPlannerPtr;


}
