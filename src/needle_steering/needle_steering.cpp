#include "needle_steering.hpp"

using namespace Needle;

void printVector(VectorXd x) {
  for (int i = 0; i < x.size(); ++i) {
    cout << x(i);
    if (i < x.size() - 1) cout << ", ";
  }
}

int main(int argc, char** argv)
{
  DblVec current_solution;
  bool is_first_time = true;

  double collision_coeff = 10;
  int max_merit_coeff_increases = 10; 
  double dynamics_coeff = 10;
  int merit_coeff_increases = 0;

  bool has_violations = true;
  bool has_dynamics_violations = true;

  bool grid_coefficient = true;

  NeedleProblemHelperPtr helper;
  OptProbPtr prob;
  OptimizerT opt;
  DblVec seq_collision_coeff, seq_dynamics_coeff;

  double prev_trust_box_size = .1;

  if (grid_coefficient) {
    while (has_violations && merit_coeff_increases <= max_merit_coeff_increases) {
      while (has_dynamics_violations && merit_coeff_increases <= max_merit_coeff_increases) {
        helper.reset(new NeedleProblemHelper());
        helper->Initialize(argc, argv);

        helper->collision_coeff = collision_coeff;
        helper->dynamics_coeff = dynamics_coeff;
        helper->initVec = current_solution;
        helper->max_merit_coeff_increases = 1;

        seq_collision_coeff.push_back(collision_coeff);
        seq_dynamics_coeff.push_back(dynamics_coeff);

        prob.reset(new OptProb());
        helper->ConfigureProblem(*prob);

        opt = OptimizerT(prob);
        helper->ConfigureOptimizer(opt);
        opt.trust_box_size_ = fmax(prev_trust_box_size, opt.min_trust_box_size_ / opt.trust_shrink_ratio_ * 1.5);

        OptStatus status = opt.optimize();
        current_solution = opt.x();
        has_dynamics_violations = helper->HasDynamicsViolations(current_solution, opt.getModel().get());
        prev_trust_box_size = opt.trust_box_size_;
        if (has_dynamics_violations) {
          dynamics_coeff *= 10;
          ++merit_coeff_increases;
        }
      }
      has_violations = has_dynamics_violations || helper->HasCollisionViolations(current_solution, opt.getModel().get());
      if (has_violations) {
        has_dynamics_violations = true;
        collision_coeff *= 10;
        dynamics_coeff *= 10;
        ++merit_coeff_increases;
        ++merit_coeff_increases;
      }
    }

    if (has_violations) {
      cout << "status: PENALTY_ITERATION_LIMIT" << endl;
    } else {
      cout << "status: CONVERGED" << endl;
    }

    cout << "n merit increases: " << merit_coeff_increases << endl;

    cout << "penalty coefficient sequence: ";

    for (int i = 0; i < seq_dynamics_coeff.size(); ++i) {
      cout << "(" << seq_collision_coeff[i] << ", " << seq_dynamics_coeff[i] << ")";
      if (i < seq_dynamics_coeff.size() - 1) cout << ", ";
      else cout << endl;
    }
  } else {

    helper.reset(new NeedleProblemHelper());
    helper->Initialize(argc, argv);
    helper->collision_coeff = 100;
    helper->dynamics_coeff = 10;
    helper->max_merit_coeff_increases = 10;
    prob.reset(new OptProb());
    helper->ConfigureProblem(*prob);
    opt = OptimizerT(prob);
    helper->ConfigureOptimizer(opt);
    OptStatus status = opt.optimize();
    current_solution = opt.x();
    if (!helper->HasDynamicsViolations(current_solution, opt.getModel().get()) && !helper->HasCollisionViolations(current_solution, opt.getModel().get())) {
      cout << "status: CONVERGED" << endl;
    } else {
      cout << "status: PENALTY_ITERATION_LIMIT" << endl;
    }
  }

  RaveDestroy();

  

  return 0;
  
}
