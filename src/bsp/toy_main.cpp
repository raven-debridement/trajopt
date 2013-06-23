#include "bsp.hpp"
#include "toy.hpp"

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

  BasicTrustRegionSQP opt(prob);
  opt.max_iter_ = 500;    

  helper->configure_optimizer(*prob, opt);

  //opt.optimize();

  cout << start.transpose() << endl;
  cout << goal.transpose() << endl;
  return 0;  
}
