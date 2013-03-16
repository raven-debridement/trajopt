#include "trajopt/plot_callback.hpp"
#include "trajopt/common.hpp"
#include "osgviewer/osgviewer.hpp"
#include "utils/eigen_conversions.hpp"
#include <boost/foreach.hpp>
#include "trajopt/problem_description.hpp"
#include "trajopt/collision_checker.hpp"
#include "trajopt/collision_avoidance.hpp"
#include "sco/modeling_utils.hpp"
using namespace OpenRAVE;
using namespace util;
using namespace Eigen;
using namespace std;
namespace trajopt {

void PlotTraj(OSGViewer& viewer, BeliefRobotAndDOFPtr rad, const TrajArray& traj, vector<GraphHandlePtr>& handles) {
	const int n_dof = rad->GetDOF();
	const int n_steps = traj.rows();

	vector<double> color_params;
	for (int i=0; i<n_steps; i++) {
		color_params.push_back(((double)i)/((double)n_steps-1.0));
	}

	for (int i=n_steps-1; i>=0; i--) {
    rad->SetDOFValues(toDblVec(traj.block(i,0,1,n_dof).transpose()));
    handles.push_back(viewer.PlotKinBody(rad->GetRobot()));
//		SetColor(handles.back(), osg::Vec4f(0,color_params[i],1.0-color_params[i],1));
  }
  rad->SetDOFValues(toDblVec(traj.block(0,0,1,n_dof).transpose()));

//	boost::shared_ptr<CollisionChecker> cc = CollisionChecker::GetOrCreate(*rad->GetRobot()->GetEnv());
//	cc->PlotCollisionGeometry(handles);
//	cc->PlotDebugGeometry(handles);

	if (traj.cols() > n_dof) { // if belief space?
		const int n_theta = rad->GetNTheta();
		OR::RobotBase::RobotStateSaver saver = const_cast<BeliefRobotAndDOF*>(rad.get())->Save();

		for (int i=0; i < n_steps; ++i) {
			VectorXd theta = traj.block(i,0,1,n_theta).transpose();
			VectorXd mean;
			MatrixXd cov;
			rad->GetEndEffectorNoiseAsGaussian(theta, mean, cov);

			if (cov(2,2) == 0)
				handles.push_back(viewer.PlotEllipseXYContour(gaussianAsTransform(mean,cov), OR::Vector(0,color_params[i],1.0-color_params[i],1)));
//			else // for some reason, the color is not rendered well for ellipses
//				handles.push_back(viewer.PlotEllipsoid(gaussianAsTransform(mean,cov), OR::Vector(0,color_params[i],1.0-color_params[i],1)));
		}

//		for (int i=0; i<n_steps; i++) {
//			VectorXd theta = traj.block(i,0,1,n_theta).transpose();
//			VectorXd x;
//			MatrixXd rt_Sigma;
//			rad->decomposeBelief(theta, x, rt_Sigma);
//
//			MatrixXd sigma_pts = rad->sigmaPoints(x, rt_Sigma*rt_Sigma. transpose());
//			for (int j=1; j<sigma_pts.cols(); j++) {
//				rad->SetDOFValues(toDblVec(sigma_pts.col(j)));
//				handles.push_back(viewer.PlotKinBody(rad->GetRobot()));
//				SetColor(handles.back(), osg::Vec4f(0,color_params[i],1.0-color_params[i],0.3));
//			}
//		}

		VectorXd theta = traj.block(0,0,1,n_theta).transpose();
		for (int i=0; i<n_steps; i++) {
			VectorXd mean;
			MatrixXd cov;
			rad->GetEndEffectorNoiseAsGaussian(theta, mean, cov);

			if (cov(2,2) == 0)
				handles.push_back(viewer.PlotEllipseXYContour(gaussianAsTransform(mean,cov), OR::Vector(0,color_params[i],1.0-color_params[i],1), true));
//			else
//				handles.push_back(viewer.PlotEllipsoid(gaussianAsTransform(mean,cov), OR::Vector(0,color_params[i],1.0-color_params[i],1)));

			if (i != (n_steps-1)) {
				VectorXd u = traj.block(i,n_theta,1,n_dof).transpose();
				theta = rad->BeliefDynamics(theta,u);
			}
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
                      _2);
}

}
