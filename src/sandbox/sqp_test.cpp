#include "sco/expr_ops.hpp"
#include "sco/modeling_utils.hpp"
#include "sco/modeling.hpp"
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



using namespace std;
using namespace trajopt;
using namespace Eigen;

typedef Matrix<double, 1, 1> Vector1d;

class CostErrorA : public ScalarOfVector {
public:
  double operator()(const VectorXd& a) const {
    return a(0);
  }
};

class CntErrorA1 : public VectorOfVector {
public:
  VectorXd operator()(const VectorXd& a) const {
    Vector1d ans; ans << a(0)*a(0) + 1 - a(1);
    return ans;
  }
};

class CntErrorA2 : public VectorOfVector {
public:
  VectorXd operator()(const VectorXd& a) const {
    Vector1d ans; ans << a(0) - 1 - a(2);
    return ans;
  }
};

class CostErrorB : public ScalarOfVector {
public:
  double operator()(const VectorXd& a) const {
    return (a(1) - 1) * (a(1) - 1);
  }
};

class CntErrorB1 : public VectorOfVector {
public:
  VectorXd operator()(const VectorXd& a) const {
    Vector1d ans; ans << a(0)*a(0);//a(0)*a(0) + 1 - a(1);
    return ans;
  }
};

class CntErrorB2 : public VectorOfVector {
public:
  VectorXd operator()(const VectorXd& a) const {
    Vector1d ans; ans << a(0)*a(0)*a(0);//a(0) - 1 - a(2);
    return ans;
  }
};

class CostErrorC : public ScalarOfVector {
public:
  double operator()(const VectorXd& a) const {
    return a(0) + a(1);
  }
};

class CntErrorC1 : public VectorOfVector {
public:
  VectorXd operator()(const VectorXd& a) const {
    Vector1d ans; ans << 1 - a(1)*a(1);
    return ans;
  }
};

class CntErrorC2 : public VectorOfVector {
public:
  VectorXd operator()(const VectorXd& a) const {
    Vector1d ans; ans << a(0)*a(1);
    return ans;
  }
};

class CostErrorD : public ScalarOfVector {
public:
  double operator()(const VectorXd& a) const {
    return 2*(a(0) + a(1));
  }
};

class CntErrorD1 : public VectorOfVector {
public:
  VectorXd operator()(const VectorXd& a) const {
    Vector1d ans; ans << -a(0)*a(1);
    return ans;
  }
};

class CostErrorE : public ScalarOfVector {
public:
  double operator()(const VectorXd& a) const {
    return a(0);
  }
};

class CntErrorE1 : public VectorOfVector {
public:
  VectorXd operator()(const VectorXd& a) const {
    Vector1d ans; ans << a(0)*a(0) + 1;
    return ans;
  }
};

int main() {
  OptProbPtr prob(new OptProb());
  LineSearchSQP opt(prob);

  VarVector vars;
  #if 0 // Settings A
    vector<string> var_names;
    var_names.push_back("x1");
    var_names.push_back("x2");
    var_names.push_back("x3");
    DblVec lbs, ubs;
    lbs.push_back(-INFINITY); lbs.push_back(0); lbs.push_back(0);
    ubs.push_back(INFINITY); ubs.push_back(INFINITY); ubs.push_back(INFINITY);
    vars = prob->createVariables(var_names, lbs, ubs);
    prob->addCost(CostPtr(new CostFromFunc(ScalarOfVectorPtr(new CostErrorA()), vars, "cost", true)));
    prob->addConstraint(ConstraintPtr(new ConstraintFromFunc(VectorOfVectorPtr(new CntErrorA1()), vars, Vector1d::Ones(), EQ, "constraint_1")));
    prob->addConstraint(ConstraintPtr(new ConstraintFromFunc(VectorOfVectorPtr(new CntErrorA2()), vars, Vector1d::Ones(), EQ, "constraint_2")));
    DblVec initVec = toDblVec(Vector3d(-3, 1, 1));
    opt.initialize(initVec);
  #endif

  #if 0 // Settings B
    vector<string> var_names;
    var_names.push_back("x1");
    var_names.push_back("x2");
    DblVec lbs, ubs;
    lbs.push_back(-INFINITY); lbs.push_back(-INFINITY);
    ubs.push_back(INFINITY); ubs.push_back(INFINITY);
    vars = prob->createVariables(var_names, lbs, ubs);
    prob->addCost(CostPtr(new CostFromFunc(ScalarOfVectorPtr(new CostErrorB()), vars, "cost", true)));
    prob->addConstraint(ConstraintPtr(new ConstraintFromFunc(VectorOfVectorPtr(new CntErrorB1()), vars, Vector1d::Ones(), EQ, "constraint_1")));
    prob->addConstraint(ConstraintPtr(new ConstraintFromFunc(VectorOfVectorPtr(new CntErrorB2()), vars, Vector1d::Ones(), EQ, "constraint_2")));
    DblVec initVec = toDblVec(Vector2d(1, 0));
    opt.initialize(initVec);
  #endif

  #if 0 // Settings C
    vector<string> var_names;
    var_names.push_back("x1");
    var_names.push_back("x2");
    DblVec lbs, ubs;
    lbs.push_back(0); lbs.push_back(0);
    ubs.push_back(INFINITY); ubs.push_back(INFINITY);
    vars = prob->createVariables(var_names, lbs, ubs);
    prob->addCost(CostPtr(new CostFromFunc(ScalarOfVectorPtr(new CostErrorC()), vars, "cost", true)));
    prob->addConstraint(ConstraintPtr(new ConstraintFromFunc(VectorOfVectorPtr(new CntErrorC1()), vars, Vector1d::Ones(), INEQ, "constraint_1")));
    prob->addConstraint(ConstraintPtr(new ConstraintFromFunc(VectorOfVectorPtr(new CntErrorC2()), vars, Vector1d::Ones(), INEQ, "constraint_2")));
    DblVec initVec = toDblVec(Vector2d(0.1, 0.9));
    opt.initialize(initVec);
  #endif

  #if 0 // Settings D
    vector<string> var_names;
    var_names.push_back("x1");
    var_names.push_back("x2");
    DblVec lbs, ubs;
    lbs.push_back(0); lbs.push_back(-1);
    ubs.push_back(INFINITY); ubs.push_back(INFINITY);
    vars = prob->createVariables(var_names, lbs, ubs);
    prob->addCost(CostPtr(new CostFromFunc(ScalarOfVectorPtr(new CostErrorD()), vars, "cost", true)));
    prob->addConstraint(ConstraintPtr(new ConstraintFromFunc(VectorOfVectorPtr(new CntErrorD1()), vars, Vector1d::Ones(), INEQ, "constraint_1")));
    DblVec initVec = toDblVec(Vector2d(0, 0));
    opt.initialize(initVec);
  #endif

  #if 1 // Settings E
    vector<string> var_names;
    var_names.push_back("x");
    DblVec lbs, ubs;
    lbs.push_back(-INFINITY);
    ubs.push_back(0);
    vars = prob->createVariables(var_names, lbs, ubs);
    prob->addCost(CostPtr(new CostFromFunc(ScalarOfVectorPtr(new CostErrorE()), vars, "cost", true)));
    prob->addConstraint(ConstraintPtr(new ConstraintFromFunc(VectorOfVectorPtr(new CntErrorE1()), vars, Vector1d::Ones(), INEQ, "constraint_1")));
    DblVec initVec; initVec.push_back(10);
    opt.initialize(initVec);
  #endif

  opt.optimize();
  return 0;
}
