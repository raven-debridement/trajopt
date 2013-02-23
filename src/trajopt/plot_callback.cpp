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

void PlotTraj(OSGViewer& viewer, BeliefRobotAndDOFPtr rad, const TrajArray& x, vector<GraphHandlePtr>& handles) {
	for (int i=0; i < x.rows(); ++i) {
    rad->SetDOFValues(toDblVec(x.row(i)));
    handles.push_back(viewer.PlotKinBody(rad->GetRobot()));
    SetTransparency(handles.back(), .35);
  }

	OR::RobotBasePtr robot = rad->GetRobot();
	OR::KinBody::LinkPtr link = robot->GetLink("Finger");
	int n_dof = robot->GetDOF();
	int n_steps = x.rows();
	OR::RobotBase::RobotStateSaver saver = const_cast<BeliefRobotAndDOF*>(rad.get())->Save();
	MatrixXd rt_Sigma = MatrixXd::Identity(n_dof,n_dof)*0.1;
	for (int i=0; i < n_steps; ++i) {
		VectorXd x0 = x.block(i,0,1,n_dof).transpose();

		MatrixXd jac = rad->EndEffectorJacobian(x0);

		robot->SetDOFValues(toDblVec(x0), false);
		OR::Vector trans = link->GetTransform().trans;
		Vector3d trans_eig(trans.x, trans.y, trans.z);
		MatrixXd cov = jac * rt_Sigma * rt_Sigma.transpose() * jac.transpose();

		handles.push_back(viewer.PlotEllipsoid(gaussianToTransform(trans_eig,cov), OR::Vector(1,0,0,1)));

		if (i < n_steps-1) {
			VectorXd x1 = x.block(i+1,0,1,n_dof).transpose();
			VectorXd u0 = x1-x0;
			VectorXd x_unused;
			rad->ekfUpdate(u0, x0, rt_Sigma, x_unused, rt_Sigma);
		}
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
