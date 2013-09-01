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

namespace RavenBSP {

  // This will generate a bunch of types like StateT, ControlT, etc.
  BSP_TYPEDEFS(
      6, // state_dim
      6, // state_noise_dim
      6, // control_dim
      12, // observe_dim
      12, // observe_noise_dim
      21, // sigma_dof
      27 // belief_dim
  );

  typedef Matrix<double, 6, 1> Vector6d;
  typedef Matrix<double, 6, 6> Matrix6d;

  // state: { {robot_dofs} }
  // control: { {d_robot_dofs} }
  // observation: { (whatever) }

  class RavenBSPProblemHelper;
  typedef boost::shared_ptr<RavenBSPProblemHelper> RavenBSPProblemHelperPtr;

  class RavenStateFunc : public StateFunc<StateT, ControlT, StateNoiseT> {
  public:
    typedef boost::shared_ptr<RavenStateFunc> Ptr;
    RavenBSPProblemHelperPtr problem_helper;

    RavenStateFunc();
    RavenStateFunc(BSPProblemHelperBasePtr helper);
    StateT operator()(const StateT& x, const ControlT& u, const StateNoiseT& m) const ;
  };

  class RavenObserveFunc : public ObserveFunc<StateT, ObserveT, ObserveNoiseT> {
  public:
    typedef boost::shared_ptr<RavenObserveFunc> Ptr;
    RavenBSPProblemHelperPtr problem_helper;
    RavenObserveFunc();
    RavenObserveFunc(BSPProblemHelperBasePtr helper);
    virtual ObserveT operator()(const StateT& x, const ObserveNoiseT& n) const;
//    virtual ObserveT observation_masks(const StateT& x, double approx_factor=-1) const;
  };

  class RavenBeliefFunc : public EkfBeliefFunc<RavenStateFunc, RavenObserveFunc, BeliefT> {
  public:
    typedef boost::shared_ptr<RavenBeliefFunc> Ptr;
    RavenBSPProblemHelperPtr problem_helper;
    RavenBeliefFunc();
    RavenBeliefFunc(BSPProblemHelperBasePtr helper, StateFuncPtr f, ObserveFuncPtr h);
  };

  class RavenBSPProblemHelper : public BSPProblemHelper<RavenBeliefFunc> {
  public:
    typedef typename BeliefConstraint<RavenBeliefFunc>::Ptr BeliefConstraintPtr;
    void add_goal_constraint(OptProb& prob);
    void add_collision_term(OptProb& prob);
    void configure_problem(OptProb& prob);
    void initialize();
    RavenBSPProblemHelper();

    RobotBasePtr robot;
    RobotAndDOFPtr rad;
    KinBody::LinkPtr link;
    Matrix4d goal_pose;

    double insertion_factor;

    Matrix4d camera_pose;
  };

  class RavenBSPPlanner : public BSPPlanner<RavenBSPProblemHelper> {
  public:
    void initialize();
    RavenBSPPlanner();
    void initialize_optimizer_parameters(BSPTrustRegionSQP& opt, bool is_first_time=true);
    RobotBasePtr robot;
    RobotAndDOFPtr rad;
    Matrix4d goal_pose;
    KinBody::LinkPtr link;

    double insertion_factor;

    Matrix4d camera_pose;
  };

  typedef boost::shared_ptr<RavenBSPPlanner> RavenBSPPlannerPtr;

  class RavenBSPWrapper {
  public:
  	EnvironmentBasePtr env;
  	OSGViewerPtr viewer;

  	bool sim_plotting;
  	bool stage_plotting;

  	RavenBSPPlannerPtr planner;
  	boost::function<void(OptProb*, DblVec&)> opt_callback;
  public:
  	void setEnvironment(const EnvironmentBasePtr& env);
  	void setViewer(const OSGViewerPtr& viewer, bool sim_plotting=false, bool stage_plotting=false);

  	RavenBSP::StateT start;
  	RavenBSP::VarianceT start_sigma;
  	deque<RavenBSP::ControlT> controls;
  	Matrix4d goal_pose;
  	int T;

    double insertion_factor;

  	string manip_name;
  	string link_name;

    Matrix4d camera_pose;

  	RavenBSPWrapper();

  	void initialize();
  	bool finished();
  	void solve();
  	void simulate_execution();


  	void run();
  };
  typedef boost::shared_ptr<RavenBSPWrapper> RavenBSPWrapperPtr;

}
