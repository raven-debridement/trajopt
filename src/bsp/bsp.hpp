#pragma once

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
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <ctime>
#include <openrave-core.h>
#include <openrave/openrave.h>
#include <assert.h>

#define DEBUG
#define BSP_DEFAULT_EPSILON 0.00390625

using namespace trajopt;
using namespace std;
using namespace OpenRAVE;
using namespace util;
using namespace boost::assign;
using namespace Eigen;

namespace BSP {

  class BSPProblemHelper {
  public:
    int state_dim;
    int control_dim;
    int observe_dim;
    int state_noise_dim;
    int observe_noise_dim;
  };

  typedef boost::shared_ptr<BSPProblemHelper> BSPProblemHelperPtr;

  class StateFunc {
  public:
    StateFunc(BSPProblemHelperPtr helper);
    virtual int output_size() const = 0;
    virtual VectorXd operator()(const VectorXd& x, const VectorXd& u, const VectorXd& m) const = 0;
    virtual void linearize(const VectorXd& x, const VectorXd& u, const VectorXd& m, MatrixXd* output_A, MatrixXd* output_B, MatrixXd* output_M) const;
    double epsilon;
  protected:
    BSPProblemHelperPtr helper;
  };

  class ObserveFunc {
  public:
    ObserveFunc(BSPProblemHelperPtr helper);
    virtual int output_size() const = 0;
    virtual VectorXd operator()(const VectorXd& x, const VectorXd& n) const = 0;
    virtual void linearize(const VectorXd& x, const VectorXd& n, MatrixXd* output_H, MatrixXd* output_N) const;
    double epsilon;
  protected:
    BSPProblemHelperPtr helper;
  };

  class BeliefFunc {
  public:
    BeliefFunc(BSPProblemHelperPtr helper);
    virtual int output_size() const = 0;
    virtual VectorXd operator()(const VectorXd& b, const VectorXd& u, const StateFunc& f) const = 0;
    virtual void linearize(const VectorXd& b, const VectorXd& u, const StateFunc& f, const ObserveFunc& h, MatrixXd* output_A, MatrixXd* output_B, VectorXd* output_c) const;
    double epsilon;
  protected:
    BSPProblemHelperPtr helper;
    void extract_state(const VectorXd& belief, VectorXd* output_state) const;
    void extract_sqrt_sigma(const VectorXd& belief, MatrixXd* output_sqrt_sigma) const;
    void extract_sigma(const VectorXd& belief, MatrixXd* output_sigma) const;
    void compose_belief(const VectorXd& state, const MatrixXd& sqrt_sigma, VectorXd* output_belief) const;
  };

  
}
