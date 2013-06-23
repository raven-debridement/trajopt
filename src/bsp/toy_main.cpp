#include "bsp.hpp"
#include "toy.hpp"
#include "utils/logging.hpp"

#define CUSTOM_PREFIX "\x1b[32m[DEBUG] "
#define LOG_CUSTOM(msg, ...) {printf(CUSTOM_PREFIX); printf(msg, ##__VA_ARGS__); printf(LOG_SUFFIX);}

using namespace BSP;
using namespace ToyBSP;

int main(int argc, char** argv) {
  int T = 20;
  int n_dof = 2;

  double start_vec_array[] = {-5, 2};
  double goal_vec_array[] = {-5, 1};

  vector<double> start_vec(start_vec_array, end(start_vec_array));
  vector<double> goal_vec(goal_vec_array, end(goal_vec_array));

  {
    Config config;
    config.add(new Parameter< vector<double> >("s", &start_vec, "s"));
    config.add(new Parameter< vector<double> >("g", &goal_vec, "g"));
    CommandParser parser(config);
    parser.read(argc, argv);
  }

  Vector2d start = toVectorXd(start_vec);
  Vector2d goal = toVectorXd(goal_vec);

  OptProbPtr prob(new OptProb());

  ToyBSPProblemHelperPtr helper(new ToyBSPProblemHelper());
  helper->start = start;
  helper->goal = goal;
  helper->T = T;
  helper->configure_problem(*prob);

  BSPTrustRegionSQP opt(prob);
  opt.max_iter_ = 50;
  opt.merit_error_coeff_ = 50;
  opt.trust_shrink_ratio_ = 0.5;
  opt.trust_expand_ratio_ = 1.25;
  opt.min_trust_box_size_ = 1e-3;
  opt.trust_box_size_ = 1;

  helper->configure_optimizer(*prob, opt);

  opt.optimize();

  LOG_CUSTOM("\n==================\n%s==================", CSTR(opt.results()));
  return 0;  
}
