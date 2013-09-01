#include "needle_steering.hpp"

using namespace Needle;

int main(int argc, char** argv) {

  NeedleProblemPlannerPtr planner(new NeedleProblemPlanner(argc, argv));

  DblVec x;
  vector< vector<Vector6d> > states;
  states.push_back(planner->starts);
  while (!planner->Finished()) {
    x = planner->InitializeSolutionWithoutFirstTimestepAndSolve(x);
    states.push_back(planner->SimulateExecution(states.back(), x));
  }

  return 0;

}
