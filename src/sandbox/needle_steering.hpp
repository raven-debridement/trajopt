#pragma once

#include "trajopt/common.hpp"
#include "sco/modeling_utils.hpp"
#include "sco/modeling.hpp"
#include <openrave-core.h>
#include <openrave/openrave.h>
#include "osgviewer/osgviewer.hpp"

//#define NEEDLE_TEST
//#define USE_CURVATURE

#define SWITCH_DEFAULT default: \
                         PRINT_AND_THROW("not implemented");

using namespace trajopt;
using namespace std;
using namespace OpenRAVE;
using namespace util;
using namespace Eigen;

namespace Needle {

  class SpeedCost : public Cost {
  public:
    SpeedCost(const Var& var, double coeff);
    virtual double value(const vector<double>& xvec);
    virtual ConvexObjectivePtr convex(const vector<double>& xvec, Model* model);
  private:
    Var var_;
    double coeff_;
    AffExpr expr_;
  };

  class RotationCost : public Cost {
  public:
    RotationCost(const VarVector& vars, double coeff);
    virtual double value(const vector<double>& xvec);
    virtual ConvexObjectivePtr convex(const vector<double>& xvec, Model* model);
  private:
    VarVector vars_;
    double coeff_;
    QuadExpr expr_;
  };

  struct LocalConfiguration : public Configuration {
    KinBodyPtr body_;
    Matrix4d pose_;
    LocalConfiguration(KinBodyPtr body, const Matrix4d& pose);
    LocalConfiguration(KinBodyPtr body);
    virtual void SetDOFValues(const DblVec& dofs);
    virtual void GetDOFLimits(DblVec& lower, DblVec& upper) const;
    virtual DblVec GetDOFValues();
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

  struct NeedleProblemHelper;

  typedef boost::shared_ptr<LocalConfiguration> LocalConfigurationPtr;
  typedef boost::shared_ptr<NeedleProblemHelper> NeedleProblemHelperPtr;

  struct PositionError : public VectorOfVector {
    LocalConfigurationPtr cfg;
    KinBodyPtr body;
    VectorXd target_pos;
    PositionError(LocalConfigurationPtr cfg, const VectorXd& target_pos);
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

  struct NeedleProblemHelper {
    // Formulation flag
    enum Formulation { Form1 = 1, Form2 = 2, Form3 = 3 };
    enum CurvatureConstraint { ConstantRadius = 1, BoundedRadius = 2 };
    enum Method { Colocation = 1, Shooting = 2 };
    enum CurvatureFormulation { UseRadius = 1, UseCurvature = 2 };
    // Config parameters
    VectorXd start;
    VectorXd goal;
    double coeff_rotation;
    double coeff_speed;
    int T;
    int n_dof;
    int formulation;
    int curvature_constraint;
    int method;
    int curvature_formulation;
    double r_min;
    vector<string> ignored_kinbody_names;
    double collision_dist_pen;
    double collision_coeff;
    double Delta_lb;
    KinBodyPtr robot;
    // Variables
    VarArray twistvars;
    VarArray phivars;
    VarArray curvaturevars;
    VarArray radiusvars;
    Var Deltavar;
    // Local configurations
    vector<LocalConfigurationPtr> local_configs;

    void ConfigureProblem(OptProb& prob);
    void InitOptimizeVariables(BasicTrustRegionSQP& opt);
    void OptimizerCallback(OptProb*, DblVec& x);
    void ConfigureOptimizer(BasicTrustRegionSQP& opt);
    void CreateVariables(OptProb& prob);
    void InitLocalConfigurations(const KinBodyPtr robot, OptProb& prob);
    void InitTrajectory(OptProb& prob);
    void AddStartConstraint(OptProb& prob);
    void AddGoalConstraint(OptProb& prob);
    void AddControlConstraint(OptProb& prob);
    void AddCollisionConstraint(OptProb& prob);
    Matrix4d TransformPose(const Matrix4d& pose, double phi, double Delta, double radius) const;
    double GetPhi(const DblVec& x, int i) const;
    double GetDelta(const DblVec& x, int i) const;
    double GetCurvature(const DblVec& x, int i) const;
    double GetRadius(const DblVec& x, int i) const;
    #ifdef NEEDLE_TEST
    void checkAlignment(DblVec& x);
    #endif
  };

  struct TrajPlotter {
    vector<LocalConfigurationPtr> local_configs;
    VarArray vars;
    OSGViewerPtr viewer;
    TrajPlotter(const vector<LocalConfigurationPtr>& local_configs, const VarArray& vars);
    void OptimizerCallback(OptProb*, DblVec& x, NeedleProblemHelperPtr helper);
    KinBodyPtr GetGoalKinBody(NeedleProblemHelperPtr helper);
    void PlotBothTrajectories(OptProbPtr prob, const BasicTrustRegionSQP& opt, NeedleProblemHelperPtr helper);
  };
}
