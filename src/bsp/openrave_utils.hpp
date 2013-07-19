#include "common.hpp"
#include "bsp_problem_helper.hpp"
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

  template< class BSPProblemHelperT >
  struct OpenRAVEPlotterMixin {
    void plot(boost::shared_ptr<BSPProblemHelperT> helper, RobotAndDOFPtr rad, OSGViewerPtr viewer, OptProb* prob, DblVec& x) {
      vector<GraphHandlePtr> handles;
      handles.clear();
      BOOST_FOREACH(CostPtr& cost, prob->getCosts()) {
        if (Plotter* plotter = dynamic_cast<Plotter*>(cost.get())) {
          plotter->Plot(x, *(rad->GetEnv()), handles);
        }
      }
      vector<ConstraintPtr> constraints = prob->getConstraints();
      BOOST_FOREACH(ConstraintPtr& cnt, constraints) {
        if (Plotter* plotter = dynamic_cast<Plotter*>(cnt.get())) {
          plotter->Plot(x, *(rad->GetEnv()), handles);
        }
      }
      for (int i = 0; i <= helper->T; ++i) {
        rad->SetDOFValues(toDblVec(getVec(x, helper->state_vars.row(i))));
        handles.push_back(viewer->PlotKinBody(rad->GetRobot()));
        SetTransparency(handles.back(), .35);
      }

      viewer->Idle();
    }
  };
}
