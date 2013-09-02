#pragma once

#include "sco/expr_ops.hpp"
#include "sco/modeling_utils.hpp"
#include "sco/modeling.hpp"
#include "osgviewer/osgviewer.hpp"
#include "trajopt/collision_checker.hpp"
#include "trajopt/bullet_collision_checker.hpp"
#include "trajopt/collision_terms.hpp"
#include "trajopt/common.hpp"
#include "trajopt/plot_callback.hpp"
#include "trajopt/problem_description.hpp"
#include "trajopt/rave_utils.hpp"
#include "trajopt/trajectory_costs.hpp"
#include "utils/clock.hpp"
#include "utils/config.hpp"
#include "utils/eigen_conversions.hpp"
#include "utils/stl_to_string.hpp"
#include <boost/assign.hpp>
#include <boost/foreach.hpp>
#include <ctime>
#include <openrave-core.h>
#include <openrave/openrave.h>
#include "strtk.hpp"

#include <stdio.h>  /* defines FILENAME_MAX */
#ifdef WINDOWS
  #include <direct.h>
  #define GetCurrentDir _getcwd
#else
  #include <unistd.h>
  #define GetCurrentDir getcwd
#endif



//#define NEEDLE_TEST
//#define USE_CURVATURE

#define SWITCH_DEFAULT default: \
                         PRINT_AND_THROW("not implemented");

using namespace trajopt;
using namespace std;
using namespace OpenRAVE;
using namespace util;
using namespace boost::assign;
using namespace Eigen;

namespace Needle {

  typedef Matrix<double, 1, 1> Vector1d;

  typedef Matrix<double, 6, 1> Vector6d;

  typedef NeedleSQP OptimizerT;
  //typedef BasicTrustRegionSQP OptimizerT;
  //typedef LineSearchSQP OptimizerT;

  inline double bound_inf(double result, double bound);

  inline string get_current_directory() {
    char cCurrentPath[FILENAME_MAX];
    if (!GetCurrentDir(cCurrentPath, sizeof(cCurrentPath))) {
      throw std::runtime_error("cannot get current path");
    }
    return string(cCurrentPath);
  }

  void AddVarArrays(OptProb& prob, int rows, const vector<int>& cols, const vector<double>& lbs, const vector<double>& ubs, const vector<string>& name_prefix, const vector<VarArray*>& newvars);
  void AddVarArrays(OptProb& prob, int rows, const vector<int>& cols, const vector<string>& name_prefix, const vector<VarArray*>& newvars);
  void AddVarArray(OptProb& prob, int rows, int cols, double lb, double ub, const string& name_prefix, VarArray& newvars);
  void AddVarArray(OptProb& prob, int rows, int cols, const string& name_prefix, VarArray& newvars);

  Matrix3d rotMat(const Vector3d& x);
  Vector3d rotVec(const Matrix3d& X);
  Matrix3d expA(const Vector3d& w);
  Matrix3d logInvA(const Vector3d& w);
  Matrix3d expRot(const Vector3d& x);
  Vector3d logRot(const Matrix3d& X);
  Matrix4d expUp(const Vector6d& x);
  Vector6d logDown(const Matrix4d& X);

  OpenRAVE::Transform matrixToTransform(const Matrix4d& X);
  OpenRAVE::Transform vecToTransform(const Vector6d& x);

  struct NeedleProblemHelper;
  typedef boost::shared_ptr<NeedleProblemHelper> NeedleProblemHelperPtr;

  struct NeedleProblemInstance;
  typedef boost::shared_ptr<NeedleProblemInstance> NeedleProblemInstancePtr;

  class ConstantSpeedCost : public Cost {
  public:
    ConstantSpeedCost(const Var& var, double coeff, NeedleProblemHelperPtr helper, NeedleProblemInstancePtr pi);
    virtual double value(const vector<double>& xvec, Model* model);
    virtual ConvexObjectivePtr convex(const vector<double>& xvec);
  private:
    Var var;
    double coeff;
    AffExpr expr;
    NeedleProblemHelperPtr helper;
    NeedleProblemInstancePtr pi;
  };

  class VariableSpeedCost : public Cost {
  public:
    VariableSpeedCost(const VarVector& vars, double coeff, NeedleProblemHelperPtr helper);
    virtual double value(const vector<double>& xvec, Model* model);
    virtual ConvexObjectivePtr convex(const vector<double>& xvec);
  private:
    VarVector vars;
    double coeff;
    AffExpr expr;
    NeedleProblemHelperPtr helper;
  };

  class SpeedDeviationCost: public Cost {
  public:
    SpeedDeviationCost(const VarVector& vars, double deviation, double coeff, NeedleProblemHelperPtr helper);
    virtual double value(const vector<double>& xvec, Model* model);
    virtual ConvexObjectivePtr convex(const vector<double>& xvec);
  private:
    VarVector vars;
    double coeff;
    double deviation;
    QuadExpr expr;
    NeedleProblemHelperPtr helper;
  };

  class RotationQuadraticCost : public Cost {
  public:
    RotationQuadraticCost(const VarVector& vars, double coeff, NeedleProblemHelperPtr helper);
    virtual double value(const vector<double>& xvec, Model* model);
    virtual ConvexObjectivePtr convex(const vector<double>& xvec);
  private:
    VarVector vars;
    double coeff;
    QuadExpr expr;
    NeedleProblemHelperPtr helper;
  };

  class RotationL1Cost : public Cost {
  public:
    RotationL1Cost(const VarVector& vars, double coeff, NeedleProblemHelperPtr helper);
    virtual double value(const vector<double>& xvec, Model* model);
    virtual ConvexObjectivePtr convex(const vector<double>& xvec);
  private:
    VarVector vars;
    double coeff;
    AffExpr expr;
    NeedleProblemHelperPtr helper;
  };

  struct LocalConfiguration : public Configuration {
    KinBodyPtr body;
    Matrix4d pose;
    LocalConfiguration(KinBodyPtr body, const Matrix4d& pose);
    LocalConfiguration(KinBodyPtr body);
    virtual void SetDOFValues(const DblVec& dofs);
    virtual void GetDOFLimits(DblVec& lower, DblVec& upper) const;
    virtual DblVec GetDOFValues() const;
    virtual int GetDOF() const;
    virtual OpenRAVE::EnvironmentBasePtr GetEnv();
    virtual DblMatrix PositionJacobian(int link_ind, const OpenRAVE::Vector& pt) const;
    virtual DblMatrix RotationJacobian(int link_ind) const;
    virtual bool DoesAffect(const KinBody::Link& link);
    virtual std::vector<KinBody::LinkPtr> GetAffectedLinks();
    virtual void GetAffectedLinks(std::vector<KinBody::LinkPtr>& links,
        bool only_with_geom, vector<int>& link_inds);
    virtual DblVec RandomDOFValues();
    virtual vector<OpenRAVE::KinBodyPtr> GetBodies();
  };

  typedef boost::shared_ptr<LocalConfiguration> LocalConfigurationPtr;

  struct PositionError : public VectorOfVector {
    LocalConfigurationPtr cfg;
    KinBodyPtr body;
    Vector6d target_pos;
    NeedleProblemHelperPtr helper;
    PositionError(LocalConfigurationPtr cfg, const Vector6d& target_pos, NeedleProblemHelperPtr helper);
    VectorXd operator()(const VectorXd& a) const;
  };

  struct PoseError : public VectorOfVector {
    LocalConfigurationPtr cfg0;
    LocalConfigurationPtr cfg1;
    NeedleProblemHelperPtr helper;
    PoseError(LocalConfigurationPtr cfg0, LocalConfigurationPtr cfg1, NeedleProblemHelperPtr helper);
    VectorXd operator()(const VectorXd& a) const;
  };

  struct SpeedDeviationError : public VectorOfVector {
    NeedleProblemHelperPtr helper;
    double deviation;
    SpeedDeviationError(double deviation, NeedleProblemHelperPtr helper);
    VectorXd operator()(const VectorXd& a) const;
  };

  struct ControlError : public VectorOfVector {
    LocalConfigurationPtr cfg0, cfg1;
    KinBodyPtr body;
    NeedleProblemHelperPtr helper;
    ControlError(LocalConfigurationPtr cfg0, LocalConfigurationPtr cfg1, NeedleProblemHelperPtr helper);
    VectorXd operator()(const VectorXd& a) const;
    int outputSize() const;
  };

  struct TrajPlotter;

  struct NeedleProblemInstance {
    double Delta_lb;
    Vector6d start;
    Vector6d goal;
    VarArray twistvars;
    VarArray phivars;
    VarArray curvature_or_radius_vars;
    VarArray Deltavars;
    Var Deltavar;
    vector<LocalConfigurationPtr> local_configs;
    vector<ConstraintPtr> dynamics_constraints;
    vector<ConstraintPtr> collision_constraints;
    DblVec initVec;
    int T;

    VectorXd GetSolution(OptimizerT& opt);
    void SetSolution(const VectorXd& sol, OptimizerT& opt);
    VectorXd GetSolutionWithoutFirstTimestep(const VectorXd& sol);
  };


  struct NeedleProblemHelper : public boost::enable_shared_from_this<NeedleProblemHelper> {
    // Formulation flag
    enum Formulation { Form1 = 1, Form2 = 2 };
    enum CurvatureConstraint { ConstantRadius = 1, BoundedRadius = 2 };
    enum Method { Colocation = 1, Shooting = 2 };
    enum CurvatureFormulation { UseRadius = 1, UseCurvature = 2 };
    enum SpeedFormulation { ConstantSpeed = 1, VariableSpeed = 2 };
    enum RotationCost { UseRotationQuadraticCost = 1, UseRotationL1Cost = 2 };
    // Config parameters
    vector<Vector6d> starts;
    vector<Vector6d> goals;
    int n_needles;
    double coeff_rotation;
    double coeff_rotation_regularization;
    double coeff_speed;
    double coeff_orientation_error;
    double improve_ratio_threshold;
    double trust_shrink_ratio;
    double trust_expand_ratio;
    double merit_error_coeff;
    int max_merit_coeff_increases;
    bool record_trust_region_history;
    vector<int> Ts;
    int n_dof;
    int formulation;
    int curvature_constraint;
    int speed_formulation;
    int method;
    int curvature_formulation;
    int rotation_cost;
    bool use_speed_deviation_constraint;
    bool use_speed_deviation_cost;
    bool plotting;
    bool verbose;
    bool plot_final_result;
    bool explicit_controls;
    bool continuous_collision;
    bool control_constraints;
    bool goal_orientation_constraint;
    double env_transparency;
    double r_min;
    vector<string> ignored_kinbody_names;
    double collision_dist_pen;
    double collision_coeff;
    vector<KinBodyPtr> robots;

    vector<NeedleProblemInstancePtr> pis;
    vector<ConstraintPtr> self_collision_constraints;

    void ConfigureProblem(OptProb& prob);
    void InitOptimizeVariables(OptimizerT& opt);
    void OptimizerCallback(OptProb*, DblVec& x);
    void ConfigureOptimizer(OptimizerT& opt);
    vector<VectorXd> GetSolutionsWithoutFirstTimestep(const vector<VectorXd>& sol);


    void InitParameters();
    void InitParametersFromConsole(int argc, char** argv);
    void Clear();

    void CreateVariables(OptProb& prob, NeedleProblemInstancePtr pi);
    void InitLocalConfigurations(const KinBodyPtr robot, OptProb& prob, NeedleProblemInstancePtr pi);
    void InitTrajectory(OptProb& prob, NeedleProblemInstancePtr pi);
    void AddRotationCost(OptProb& prob, NeedleProblemInstancePtr pi);
    void AddSpeedCost(OptProb& prob, NeedleProblemInstancePtr pi);
    void AddSpeedConstraint(OptProb& prob, NeedleProblemInstancePtr pi);
    void AddStartConstraint(OptProb& prob, NeedleProblemInstancePtr pi);
    void AddGoalConstraint(OptProb& prob, NeedleProblemInstancePtr pi);
    void AddControlConstraint(OptProb& prob, NeedleProblemInstancePtr pi);
    void AddPoseConstraint(OptProb& prob, NeedleProblemInstancePtr pi);
    void AddCollisionConstraint(OptProb& prob, NeedleProblemInstancePtr pi);
    void AddSelfCollisionConstraint(OptProb& prob, NeedleProblemInstancePtr piA, NeedleProblemInstancePtr piB);
    //void AddCollisionCost(OptProb& prob, NeedleProblemInstancePtr pi);
    void InitializeCollisionEnvironment();

    Matrix4d TransformPose(const Matrix4d& pose, double phi, double Delta, double radius) const;
    double GetPhi(const DblVec& x, int i, NeedleProblemInstancePtr pi) const;
    double GetDelta(const DblVec& x, int i, NeedleProblemInstancePtr pi) const;
    double GetCurvatureOrRadius(const DblVec& x, int i, NeedleProblemInstancePtr pi) const;
    double GetCurvature(const DblVec& x, int i, NeedleProblemInstancePtr pi) const;
    double GetRadius(const DblVec& x, int i, NeedleProblemInstancePtr pi) const;

    vector<VectorXd> GetSolutions(OptimizerT& opt);
    void SetSolutions(const vector<VectorXd>& sol, OptimizerT& opt);

    void AddNeedlesToBullet(OptimizerT& prob);
    void AddNeedleToBullet(NeedleProblemInstancePtr pi, OptimizerT& prob);

    #ifdef NEEDLE_TEST
    void checkAlignment(DblVec& x);
    #endif
  };

  struct NeedleProblemPlanner {
    int argc;
    char **argv;
    int n_needles;
    vector<int> Ts;
    NeedleProblemHelperPtr helper;

    bool plotting;
    bool plot_final_result;
    bool verbose;
    double env_transparency;
    string data_dir;
    string env_file_path;
    string robot_file_path;

    EnvironmentBasePtr env;
    boost::shared_ptr<TrajPlotter> plotter;

    vector<Vector6d> starts;
    vector<Vector6d> goals;
    DblVec x;

    NeedleProblemPlanner(int argc, char **argv);
    ~NeedleProblemPlanner();
    Vector6d PerturbState(const Vector6d& state);
    vector<VectorXd> Solve(const vector<VectorXd>& initial = vector<VectorXd>() );
    vector<VectorXd> GetSolutionsWithoutFirstTimestep(const vector<VectorXd>& sol);
    DblVec Solve(const DblVec& x);
    vector<Vector6d> SimulateExecution(const vector<Vector6d>& current_states);//const vector<VectorXd>& sol);
    //DblVec InitializeSolutionWithoutFirstTimestepAndSolve(const DblVec& x);
    bool Finished() const;
  };

  typedef boost::shared_ptr<NeedleProblemPlanner> NeedleProblemPlannerPtr;

  struct TrajPlotter {
    OSGViewerPtr viewer;
    vector<NeedleProblemInstancePtr> pis;
    TrajPlotter(const vector<NeedleProblemInstancePtr>& pis);
    void OptimizerCallback(OptProb*, DblVec& x, NeedleProblemHelperPtr helper);
  };
}
