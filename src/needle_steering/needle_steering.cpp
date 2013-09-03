#include "needle_steering.hpp"

using namespace Needle;

int main(int argc, char** argv) {

  NeedleProblemPlannerPtr planner(new NeedleProblemPlanner(argc, argv));

  vector<VectorXd> sol;
  vector< vector<Vector6d> > states;
  states.push_back(planner->starts);
  while (!planner->Finished()) {
    cout << planner->Ts[0] << endl;
    sol = planner->GetSolutionsWithoutFirstTimestep(planner->Solve(sol));
    states.push_back(planner->SimulateExecution(states.back()));
    cout << "simulate one step" << endl;
  }

  return 0;

}
