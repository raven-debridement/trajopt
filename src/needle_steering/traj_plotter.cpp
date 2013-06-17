#include "needle_steering.hpp"

namespace Needle {

  TrajPlotter::TrajPlotter(const vector<LocalConfigurationPtr>& local_configs, const VarArray& vars) : local_configs(local_configs), vars(vars) {
    viewer = OSGViewer::GetOrCreate(local_configs[0]->GetEnv());
  }

  void TrajPlotter::OptimizerCallback(OptProb*, DblVec& x, NeedleProblemHelperPtr helper) {
    vector<GraphHandlePtr> handles;
    vector<KinBodyPtr> bodies = local_configs[0]->GetBodies();
    MatrixXd vals = getTraj(x, vars);
    for (int i=0; i < vals.rows(); ++i) {
      local_configs[i]->SetDOFValues(toDblVec(vals.row(i)));
      BOOST_FOREACH(const KinBodyPtr& body, bodies) {
        handles.push_back(viewer->PlotKinBody(body));
        SetTransparency(handles.back(), .35);
      }
    }
    handles.push_back(viewer->PlotKinBody(GetGoalKinBody(helper))); 
    SetTransparency(handles.back(), 1);
    viewer->Idle();
  }

  KinBodyPtr TrajPlotter::GetGoalKinBody(NeedleProblemHelperPtr helper) {
    OpenRAVE::Transform T;
    T.trans.x = helper->goal[0];
    T.trans.y = helper->goal[1];
    T.trans.z = helper->goal[2];
    OpenRAVE::Vector rot(helper->goal[3], helper->goal[4], helper->goal[5]);
    T.rot = OpenRAVE::geometry::quatFromAxisAngle(rot);
    helper->robot->SetTransform(T);
    return helper->robot;
  }

  void TrajPlotter::PlotBothTrajectories(OptProbPtr prob, const BasicTrustRegionSQP& opt, NeedleProblemHelperPtr helper) {
    DblVec x = prob->getModel()->getVarValues(prob->getModel()->getVars());
    vector<GraphHandlePtr> handles;
    KinBodyPtr robot = helper->robot;
    // plot real trajectory
    Matrix4d current_pose = expUp(helper->start);
    robot->SetTransform(Needle::matrixToTransform(current_pose));
    handles.push_back(viewer->PlotKinBody(robot));
    SetTransparency(handles.back(), .5);
    for (int i = 0; i < helper->T; ++i) {
      double phi = helper->GetPhi(x, i);
      double Delta = helper->GetDelta(x, i);
      double curvature_or_radius = helper->GetCurvatureOrRadius(x, i);
      current_pose = helper->TransformPose(current_pose, phi, Delta, curvature_or_radius);
      robot->SetTransform(Needle::matrixToTransform(current_pose));
      handles.push_back(viewer->PlotKinBody(robot));
      SetTransparency(handles.back(), .5);
    }
    // plot ideal trajectory
    MatrixXd vals = getTraj(x, vars);
    for (int i=0; i < vals.rows(); ++i) {
      local_configs[i]->SetDOFValues(toDblVec(vals.row(i)));
      handles.push_back(viewer->PlotKinBody(robot));
      SetTransparency(handles.back(), .35);
    }
    viewer->Idle();
  }
}
