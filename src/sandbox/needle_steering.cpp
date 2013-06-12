#include "o3.hpp"
#include "sco/expr_ops.hpp"
#include "sco/modeling_utils.hpp"
#include "osgviewer/osgviewer.hpp"
#include "trajopt/collision_checker.hpp"
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

using namespace trajopt;
using namespace std;
using namespace OpenRAVE;
using namespace util;
using namespace boost::assign;
using namespace Eigen;


namespace Needle {

  template<typename T, size_t N>
  T* end(T (&ra)[N]) {
    return ra + N;
  }


  void AddVarArrays(OptProb& prob, int rows, const vector<int>& cols, const vector<double>& lbs, const vector<double>& ubs, const vector<string>& name_prefix, const vector<VarArray*>& newvars) {
    int n_arr = name_prefix.size();
    assert(n_arr == newvars.size());

    vector<MatrixXi> index(n_arr);
    for (int i=0; i < n_arr; ++i) {
      newvars[i]->resize(rows, cols[i]);
      index[i].resize(rows, cols[i]);
    }

    vector<string> names;
    vector<double> all_lbs;
    vector<double> all_ubs;
    int var_idx = prob.getNumVars();
    for (int i=0; i < rows; ++i) {
      for (int k=0; k < n_arr; ++k) {
        for (int j=0; j < cols[k]; ++j) {
          index[k](i,j) = var_idx;
          names.push_back( (boost::format("%s_%i_%i")%name_prefix[k]%i%j).str() );
          all_lbs.push_back(lbs[k]);
          all_ubs.push_back(ubs[k]);
          ++var_idx;
        }
      }
    }
    prob.createVariables(names, all_lbs, all_ubs); // note that w,r, are both unbounded

    const vector<Var>& vars = prob.getVars();
    for (int k=0; k < n_arr; ++k) {
      for (int i=0; i < rows; ++i) {
        for (int j=0; j < cols[k]; ++j) {
          (*newvars[k])(i,j) = vars[index[k](i,j)];
        }
      }
    }
  }

  void AddVarArrays(OptProb& prob, int rows, const vector<int>& cols, const vector<string>& name_prefix, const vector<VarArray*>& newvars) {
    vector<double> lbs(newvars.size(), -INFINITY);
    vector<double> ubs(newvars.size(), INFINITY);
    AddVrArrays(prob, rows, cols, lbs, ubs, name_prefix, new_vars);
  }

  void AddVarArray(OptProb& prob, int rows, int cols, double lb, double ub, const string& name_prefix, VarArray& newvars) {
    vector<VarArray*> arrs(1, &newvars);
    vector<string> prefixes(1, name_prefix);
    vector<int> colss(1, cols);
    vector<double> lbs(1, lb);
    vector<double> ubs(1, ub);
    AddVarArrays(prob, rows, colss, lbs, ubs, prefixes, arrs);
  }

  void AddVarArray(OptProb& prob, int rows, int cols, const string& name_prefix, VarArray& newvars) {
    AddVarArray(prob, rows, colss, -INFINITY, INFINITY, prefixes, arrs);
  }

  struct LocalConfiguration : public Configuration {

  };

  typedef boost::shared_ptr<LocalConfiguration> LocalConfigurationPtr;

  class SpeedCost : public Cost {
  public:
    SpeedCost(const Var& var, double coeff) : Cost("Speed"), var_(var), coeff_(coeff) :
      exprInc(expr_, exprMult(var, coeff));
    }
    virtual double value(const vector<double>& xvec) {
      double speed = getVec(xvec, singleton<Var>(var_))[0];
      return speed * coeff;
    }
    virtual ConvexObjectivePtr convex(const vector<double>& xvec, Model* model) {
      ConvexObjectivePtr out(new ConvexObjective(model));
      out->addAffExpr(expr_);
      return out;
    }
  private:
    Var var_;
    double coeff_;
  };

  class RotationCost : public Cost {
  public:
    RotationCost(const VarVector& vars, double coeff) : Cost("Rotation"), vars_(vars), coeff_(coeff) :
      for (int i = 0; i < vars.size(); ++i) {
        exprInc(expr_, exprMult(exprSquare(vars[i]), coeff));
      }
    }
    virtual double value(const vector<double>& xvec) {
      VectorXd vals = getVec(xvec, vars_);
      return vals.array().square().sum() * coeff_;
    }
    virtual ConvexObjectivePtr convex(const vector<double>& xvec, Model* model) {
      ConvexObjectivePtr out(new ConvexObjective(model));
      out->addQuadExpr(expr_);
      return out;
    }
  private:
    VarVector vars_;
    double coeff_;
  };

  struct EntryError : public VectorOfVector {

  };

  struct GoalError : public VectorOfVector {

  };

  struct NeedleProblemHelper {
    // Config parameters
    VectorXd start;
    VectorXd goal;
    double coeff_rotation;
    double coeff_speed;
    int T;
    int n_dof;
    double r_min;
    vector<string> ignored_kinbody_names;
    double collision_dist_pen;
    double collision_coeff;
    // Variables
    VarArray twistvars;
    VarArray phivars;
    Var Delta;
    // Local configurations
    vector<LocalConfigurationPtr> local_configs;

    void ConfigureProblem(const KinBodyPtr robot, OptProb& prob) {
      CreateVariables(prob);
      InitLocalConfigurations(robot, prob);
      prob.addCost(CostPtr(new RotationCost(phivars.col(0), coeff_rotation)));
      prob.addCost(CostPtr(new SpeedCost(Delta, coeff_speed)));
      AddEntryConstraint(prob);
      AddGoalConstraint(prob);
      AddControlConstraint(prob);
      AddCollisionConstraint(prob);
    }

    void CreateVariables(OptProb& prob) {
      AddVarArray(prob, T+1, 6, "twist", twistvars);
      AddVarArray(prob, T, 1, -PI, PI, "phi", phivars);
      double Delta_lb = (goal.topRows(3) - start.topRows(3)).norm() / T / r_min;
      Delta = prob.createVariables(singleton<string>("Delta"), singleton<double>(Delta_lb),singleton<double>(INFINITY))[0];
      // Only the twist variables are incremental (i.e. their trust regions should be around zero)
      prob.setIncremental(twistvars);
    }

    void InitLocalConfigurations(const KinBodyPtr robot, OptProb& prob) {
      for (int i = 0; i <= T; ++i) {
        local_configs.push_back(LocalConfigurationPtr(new LocalConfiguration(robot)));
      }
    }

    void AddEntryConstraint(OptProb& prob) {
      VarVector vars = twistvars.row(0);
      VectorOfVectorPtr f(new Needle::EntryError(local_configs[0]));
      VectorXd coeffs = VectorXd::Ones(6);
      prob.addConstraint(ConstraintFromFunc(f, vars, coeffs, EQ, "entry"));
    }

    void AddGoalConstraint(OptProb& prob) {
      VarVector vars = twistvars.row(T);
      VectorOfVectorPtr f(new Needle::GoalError(local_configs[T]));
      VectorXd coeffs = VectorXd::Ones(6);
      prob.addConstraint(ConstraintFromFunc(f, vars, coeffs, EQ, "goal"));
    }

    void AddControlConstraint(OptProb& prob) {

    }

    void AddCollisionConstraint(OptProb& prob) {

    }

  };

  

  struct NeedleError : public VectorOfVector {
    ConfigurationPtr cfg0, cfg1;
    double radius;
    KinBodyPtr body;
    NeedleError(ConfigurationPtr cfg0, ConfigurationPtr cfg1, double radius) : cfg0(cfg0), cfg1(cfg1), radius(radius), body(cfg0->GetBodies()[0]) {}
    VectorXd operator()(const VectorXd& a) const {
      cfg0->SetDOFValues(toDblVec(a.topRows(6)));
      OR::Transform Tw0 = body->GetTransform();
      cfg1->SetDOFValues(toDblVec(a.middleRows(6,6)));
      OR::Transform Tw1 = body->GetTransform();
      double theta = a(12);
      OR::Transform Ttarg0(OR::geometry::quatFromAxisAngle(OR::Vector(0,0,1), theta),
          OR::Vector(radius*sin(theta), radius*(1-cos(theta)),0));
          
      OR::Transform Ttargw = Tw0 * Ttarg0;
      OR::Vector position_errA = Ttargw.trans - Tw1.trans;
      OR::Vector ori_err = Ttargw * OR::Vector(1,0,0) - Tw1 * OR::Vector(1,0,0);
      return concat(toVector3d(position_errA), toVector3d(ori_err));

    }
  };

  struct TrajPlotter {
    vector<IncrementalRBPtr> rbs;
    VarArray vars;
    OSGViewerPtr viewer;
    TrajPlotter(const vector<IncrementalRBPtr>& rbs, const VarArray& vars);
    void OptimizerCallback(OptProb*, DblVec& x);
  };

  TrajPlotter::TrajPlotter(const vector<IncrementalRBPtr>& rbs, const VarArray& vars) : rbs(rbs), vars(vars) {
    viewer = OSGViewer::GetOrCreate(rbs[0]->GetEnv());
  }
  void TrajPlotter::OptimizerCallback(OptProb*, DblVec& x) {
    vector<GraphHandlePtr> handles;
    vector<KinBodyPtr> bodies = rbs[0]->GetBodies();
    MatrixXd traj = getTraj(x,vars);
    for (int i=0; i < traj.rows(); ++i) {
      rbs[i]->SetDOFValues(toDblVec(traj.row(i)));
      BOOST_FOREACH(const KinBodyPtr& body, bodies) {
        handles.push_back(viewer->PlotKinBody(body));
        SetTransparency(handles.back(), .35);
      }
    }
    viewer->Idle();
  }  

}

int main(int argc, char** argv)
{
  bool plotting=false, verbose=false;
  double env_transparency = 0.5;

  int T = 25;
  double r_min = 2;
  int n_dof = 6;

  double improve_ratio_threshold = 0.25;
  double trust_shrink_ratio = 0.7;
  double trust_expand_ratio = 1.2;
  
  double start_vec_array[] = {-12.82092, 6.80976, 0.06844, 0, 0, 0};
  double goal_vec_array[] = {-3.21932, 6.87362, -1.21877, 0, 0, 0};

  vector<double> start_vec(start_vec_array, start_vec_array + n_dof);
  vector<double> goal_vec(goal_vec_array, goal_vec_array + n_dof);

  {
    Config config;
    config.add(new Parameter<bool>("plotting", &plotting, "plotting"));
    config.add(new Parameter<bool>("verbose", &verbose, "verbose"));
    config.add(new Parameter<double>("env_transparency", &env_transparency, "env_transparency"));
    config.add(new Parameter<int>("T", &T, "T"));
    config.add(new Parameter<double>("r_min", &r_min, "r_min"));
    config.add(new Parameter<double>("improve_ratio_threshold", &improve_ratio_threshold, "improve_ratio_threshold"));
    config.add(new Parameter<double>("trust_shrink_ratio", &trust_shrink_ratio, "trust_shrink_ratio"));
    config.add(new Parameter<double>("trust_expand_ratio", &trust_expand_ratio, "trust_expand_ratio"));
    config.add(new Parameter< vector<double> >("s", &start_vec, "s"));
    config.add(new Parameter< vector<double> >("g", &goal_vec, "g"));
    CommandParser parser(config);
    parser.read(argc, argv);
  }

  

  double coeff_rotation = 1.;
  double coeff_speed = 1.;

  RaveInitialize(false, verbose ? Level_Debug : Level_Info);
  EnvironmentBasePtr env = RaveCreateEnvironment();
  env->StopSimulation();
  OSGViewerPtr viewer = OSGViewer::GetOrCreate(env);
  assert(viewer);

  env->Load(string(DATA_DIR) + "/prostate.env.xml");//needleprob.env.xml");
  viewer->SetAllTransparency(env_transparency);
  RobotBasePtr robot = GetRobot(*env);

  VectorXd start(n_dof); for (int i = 0; i < n_dof; ++i) start[i] = start_vec[i];
  VectorXd goal(n_dof); for (int i = 0; i < n_dof; ++i) goal[i] = goal_vec[i];

  const char *ignored_kinbody_c_strs = { "KinBodyProstate", "KinBodyDermis", "KinBodyEpidermis", "KinBodyHypodermis" };
  vector<string> ignored_kinbody_names(ignored_kinbody_c_strs, end(ignored_kinbody_c_strs));

  OptProbPtr prob(new OptProb());

  Needle::NeedleProblemHelper helper;
  helper.start = start;
  helper.goal = goal;
  helper.coeff_rotation = coeff_rotation;
  helper.coeff_speed = coeff_speed;
  helper.T = T;
  helper.r_min = r_min;
  helper.n_dof = n_dof;
  helper.ignored_kinbody_names = ignored_kinbody_names;
  helper.collision_dist_pen = 0.025;
  helper.collision_coeff = 20;
  helper.ConfigureProblem(*prob);

  BasicTrustRegionSQP opt(prob);
  helper.ConfigureOptimizer(opt);
  opt.max_iter_ = 500;    
  opt.improve_ratio_threshold_ = improve_ratio_threshold;
  opt.trust_shrink_ratio_ = trust_shrink_ratio;
  opt.trust_expand_ratio_ = trust_expand_ratio;

  boost::shared_ptr<Needle::TrajPlotter> plotter;
  if (plotting) {
    plotter.reset(new Needle::TrajPlotter(helper.m_rbs, trajvars));
    opt.addCallback(boost::bind(&Needle::TrajPlotter::OptimizerCallback, boost::ref(plotter), _1, _2));
  }

  MatrixXd initTraj(n_steps, n_dof);  
  for (int idof = 0; idof < n_dof; ++idof) {
    initTraj.col(idof) = VectorXd::LinSpaced(n_steps, start[idof], goal[idof]);
  }
  DblVec initVec = trajToDblVec(initTraj);
  initVec.push_back(dtheta_lb);
  opt.initialize(initVec);
  
  opt.optimize();

  RaveDestroy();


}
