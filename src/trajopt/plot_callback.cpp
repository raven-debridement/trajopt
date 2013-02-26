#include "trajopt/plot_callback.hpp"
#include "trajopt/common.hpp"
#include "osgviewer/osgviewer.hpp"
#include "utils/eigen_conversions.hpp"
#include <boost/foreach.hpp>
#include "trajopt/problem_description.hpp"
using namespace OpenRAVE;
using namespace util;
using namespace Eigen;
using namespace std;
namespace trajopt {

void PlotTraj(OSGViewer& viewer, BeliefRobotAndDOFPtr rad, const TrajArray& traj, vector<GraphHandlePtr>& handles) {
	const int n_dof = rad->GetDOF();
	const int theta_size = n_dof + n_dof*(n_dof+1)/2;
	const int n_steps = traj.rows();

	for (int i=0; i < n_steps; ++i) {
    rad->SetDOFValues(toDblVec(traj.block(i,0,1,n_dof).transpose()));
    handles.push_back(viewer.PlotKinBody(rad->GetRobot()));
		double trans_param = (((double)i)/((double)n_steps-1.0)+0.35)/1.35;
    SetTransparency(handles.back(), trans_param);
  }

	OR::RobotBase::RobotStateSaver saver = const_cast<BeliefRobotAndDOF*>(rad.get())->Save();

	for (int i=0; i < n_steps; ++i) {
		VectorXd theta = traj.block(i,0,1,theta_size).transpose();
		VectorXd x;
		MatrixXd rt_Sigma;
		rad->decomposeBelief(theta, x, rt_Sigma);

		MatrixXd jac = rad->EndEffectorJacobian(x);

		rad->SetDOFValues(toDblVec(x));
		OR::Vector trans = rad->link->GetTransform().trans;
		Vector3d trans_eig(trans.x, trans.y, trans.z);
		MatrixXd cov = MatrixXd::Identity(3,3);
		cov.topLeftCorner(n_dof, n_dof) = jac * rt_Sigma * rt_Sigma.transpose() * jac.transpose();
		cov(2,2) = 0.00000001; // otherwise, flat ellipsoids have a coloring artifact

		handles.push_back(viewer.PlotEllipsoid(gaussianToTransform(trans_eig,cov), OR::Vector(1,0,0,1)));
    SetTransparency(handles.back(), 0.35);
	}
}

void PlotCosts(OSGViewer& viewer, vector<CostPtr>& costs, vector<ConstraintPtr>& cnts, BeliefRobotAndDOFPtr rad, const VarArray& vars, const DblVec& x) {
  vector<GraphHandlePtr> handles;
  handles.clear();
  BOOST_FOREACH(CostPtr& cost, costs) {
    if (Plotter* plotter = dynamic_cast<Plotter*>(cost.get())) {
      plotter->Plot(x, *rad->GetRobot()->GetEnv(), handles);
    }
  }
  BOOST_FOREACH(ConstraintPtr& cnt, cnts) {
    if (Plotter* plotter = dynamic_cast<Plotter*>(cnt.get())) {
      plotter->Plot(x, *rad->GetRobot()->GetEnv(), handles);
    }
  }
  TrajArray traj = getTraj(x, vars);
  PlotTraj(viewer, rad, traj, handles);
  viewer.Idle();
  rad->SetDOFValues(toDblVec(traj.row(traj.rows()-1)));
}



Optimizer::Callback PlotCallback(TrajOptProb& prob) {
  OSGViewerPtr viewer = OSGViewer::GetOrCreate(prob.GetEnv());
  vector<ConstraintPtr> cnts = prob.getConstraints();
  return boost::bind(&PlotCosts, boost::ref(*viewer),
                      boost::ref(prob.getCosts()),
                      cnts,
                      prob.GetRAD(),
                      boost::ref(prob.GetVars()),
                      _1);
}

}
