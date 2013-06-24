#include "bsp.hpp"
#include "toy.hpp"
#include "utils/logging.hpp"

#define CUSTOM_PREFIX "\x1b[32m[CUSTOM] "
#define LOG_CUSTOM(msg, ...) {printf(CUSTOM_PREFIX); printf(msg, ##__VA_ARGS__); printf(LOG_SUFFIX);}

using namespace BSP;
using namespace ToyBSP;

int main(int argc, char** argv) {

  int T = 20;
  int n_dof = 2;

  bool plotting = true;
  double start_vec_array[] = {-5, 2};
  double goal_vec_array[] = {-5, 0};

  vector<double> start_vec(start_vec_array, end(start_vec_array));
  vector<double> goal_vec(goal_vec_array, end(goal_vec_array));

  {
    Config config;
    config.add(new Parameter< vector<double> >("s", &start_vec, "s"));
    config.add(new Parameter< vector<double> >("g", &goal_vec, "g"));
    config.add(new Parameter<bool>("plotting", &plotting, "plotting"));
    CommandParser parser(config);
    parser.read(argc, argv);
  }

  Vector2d start = toVectorXd(start_vec);
  Vector2d goal = toVectorXd(goal_vec);
  Matrix2d start_sigma = Matrix2d::Identity();

  OptProbPtr prob(new OptProb());

  ToyBSPProblemHelperPtr helper(new ToyBSPProblemHelper());
  helper->start = start;
  helper->goal = goal;
  helper->start_sigma = start_sigma;
  helper->T = T;
  helper->configure_problem(*prob);

  BSPTrustRegionSQP opt(prob);
  opt.max_iter_ = 500;
  opt.merit_error_coeff_ = 500;
  opt.trust_shrink_ratio_ = 0.5;
  opt.trust_expand_ratio_ = 1.25;
  opt.min_trust_box_size_ = 1e-3;
  opt.min_approx_improve_ = 1e-2;
  opt.min_approx_improve_frac_ = 1e-4;
  opt.improve_ratio_threshold_ = 0.2;
  opt.trust_box_size_ = 1;

  helper->configure_optimizer(*prob, opt);

  boost::shared_ptr<ToyPlotter> plotter;
  if (plotting) {
    plotter.reset(new ToyPlotter(helper->state_vars, helper->sqrt_sigma_vars, helper->control_vars));
    opt.addCallback(boost::bind(&ToyPlotter::plot_callback, boost::ref(plotter), _1, _2, helper));
  }

  opt.optimize();

  LOG_CUSTOM("\n==================\n%s==================", CSTR(opt.results()));

  return 0;
}
