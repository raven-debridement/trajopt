#include "needle_steering.hpp"

namespace Needle {

  TrajPlotter::TrajPlotter(const vector<NeedleProblemInstancePtr>& pis) : pis(pis) {
    viewer = OSGViewer::GetOrCreate(pis[0]->local_configs[0]->GetEnv());
  }

  void TrajPlotter::OptimizerCallback(OptProb* prob, DblVec& x, NeedleProblemHelperPtr helper) {
    vector<GraphHandlePtr> handles;
    //BOOST_FOREACH(CostPtr& cost, prob->getCosts()) {
    //  if (Plotter* plotter = dynamic_cast<Plotter*>(cost.get())) {
    //    plotter->Plot(x, *(helper->local_configs[0]->GetEnv()), handles);
    //  }
    //}
    //vector<ConstraintPtr> constraints = prob->getConstraints();
    //BOOST_FOREACH(ConstraintPtr& cnt, constraints) {
    //  if (Plotter* plotter = dynamic_cast<Plotter*>(cnt.get())) {
    //    plotter->Plot(x, *(helper->local_configs[0]->GetEnv()), handles);
    //  }
    //}
    //CollisionChecker::GetOrCreate(
    EnvironmentBasePtr env = helper->pis[0]->local_configs[0]->GetEnv();
    CollisionChecker::GetOrCreate(*env)->PlotCollisionGeometry(handles);//SetContactDistance(collision_dist_pen + 0.05);
    for (int k = 0; k < pis.size(); ++k) {
      vector<KinBodyPtr> bodies = pis[k]->local_configs[0]->GetBodies();
      MatrixXd vals = getTraj(x, pis[k]->twistvars);
      for (int i=0; i < vals.rows(); ++i) {
        pis[k]->local_configs[i]->SetDOFValues(toDblVec(vals.row(i)));
        BOOST_FOREACH(const KinBodyPtr& body, bodies) {
          handles.push_back(viewer->PlotKinBody(body));
          SetTransparency(handles.back(), .35);
        }
      }
    }
    viewer->Idle();
  }
}
