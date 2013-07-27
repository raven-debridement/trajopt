#include "common.hpp"
#include "bsp_planner.hpp"
#include <openrave-core.h>
#include <openrave/openrave.h>

namespace BSP {
  OpenRAVE::Transform matrixToTransform(const Matrix4d& X) {
    OpenRAVE::TransformMatrix M;
    M.trans.x = X(0, 3);
    M.trans.y = X(1, 3);
    M.trans.z = X(2, 3);
    for (int row = 0; row < 3; ++row) {
      for (int col = 0; col < 3; ++col) {
        M.m[row*4+col] = X(row, col);
      }
    }
    return OpenRAVE::Transform(M);
  }

  template< class StateT, class VarianceT >
  void belief_to_endeffector_noise(RobotAndDOFPtr rad, OR::KinBody::LinkPtr endeffector, const StateT& state, const VarianceT& sigma, Vector3d* mean, Matrix3d* cov) {
    assert (mean != nullptr);
    assert (cov != nullptr);
    Configuration::SaverPtr saver = rad->Save();
    rad->SetDOFValues(toDblVec(state));
    OR::Vector trans = endeffector->GetTransform().trans;
    MatrixXd jac = rad->PositionJacobian(endeffector->GetIndex(), trans);
    *mean = Vector3d(trans.x, trans.y, trans.z);
    *cov = jac * sigma * jac.transpose();
  }

  osg::Matrix gaussian_as_transform(const Vector3d& mean, const Matrix3d& cov) {
    EigenSolver<Matrix3d> es(cov);
    Matrix4d t = Matrix4d::Identity();
    t.block(0,0,3,3) = es.eigenvectors().real() * es.eigenvalues().real().cwiseSqrt().asDiagonal();
    t.block(0,3,3,1) = mean;
    Matrix4d t_transpose = t.transpose();

    osg::Matrix osg_t;
    osg_t.set(t_transpose.data());
    return osg_t;
  }

  template< class BSPPlannerT >
  struct OpenRAVEPlotterMixin {
    static void plot_opt_trajectory(boost::shared_ptr<BSPPlannerT> planner, RobotAndDOFPtr rad, OSGViewerPtr viewer, OptProb* prob, DblVec& x, vector<GraphHandlePtr>* handles) {
      assert (handles != nullptr);
      BOOST_FOREACH(CostPtr& cost, prob->getCosts()) {
        if (Plotter* plotter = dynamic_cast<Plotter*>(cost.get())) {
          plotter->Plot(x, *(rad->GetEnv()), *handles);
        }
      }
      vector<ConstraintPtr> constraints = prob->getConstraints();
      BOOST_FOREACH(ConstraintPtr& cnt, constraints) {
        if (Plotter* plotter = dynamic_cast<Plotter*>(cnt.get())) {
          plotter->Plot(x, *(rad->GetEnv()), *handles);
        }
      }
      for (int i = 0; i <= planner->helper->T; ++i) {
        cout << "dof: " << getVec(x, planner->helper->state_vars.row(i)).transpose() << endl;
        rad->SetDOFValues(toDblVec(getVec(x, planner->helper->state_vars.row(i))));
        handles->push_back(viewer->PlotKinBody(rad->GetRobot()));
        SetTransparency(handles->back(), .35);
      }
    }

    static void plot_sim_trajectory(boost::shared_ptr<BSPPlannerT> planner, RobotAndDOFPtr rad, OSGViewerPtr viewer, vector<GraphHandlePtr>* handles) {
      assert (handles != nullptr);
      for (int i = 0; i < planner->simulated_positions.size(); ++i) {
        rad->SetDOFValues(toDblVec(planner->simulated_positions[i]));
        handles->push_back(viewer->PlotKinBody(rad->GetRobot()));
        SetTransparency(handles->back(), .35);
      }
    }

    static void stage_plot_callback(boost::shared_ptr<BSPPlannerT> planner, RobotAndDOFPtr rad, OSGViewerPtr viewer, OptProb* prob, DblVec& x) {
      vector<GraphHandlePtr> handles;
      plot_opt_trajectory(planner, rad, viewer, prob, x, &handles);
      viewer->Idle();
    }

    static void sim_plot_callback(boost::shared_ptr<BSPPlannerT> planner, RobotAndDOFPtr rad, OSGViewerPtr viewer) {
      vector<GraphHandlePtr> handles;
      plot_sim_trajectory(planner, rad, viewer, &handles);
      viewer->Idle();
    }
  };
}
