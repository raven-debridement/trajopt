#pragma once

#include "common.hpp"
#include "collision/collision.hpp"

namespace BSP {
  using namespace BSPCollision;
  template< class BeliefFuncT >
  struct BeliefCollisionPlotterMixin {
    void PlotCollisions(const std::vector<BeliefCollision>& collisions, OR::EnvironmentBase& env, vector<OR::GraphHandlePtr>& handles, double safe_dist) {
      BOOST_FOREACH(const BeliefCollision& col, collisions) {
        RaveVectorf color;
        if (col.distance < 0) {
        	color = RaveVectorf(0,0,0,1);

        }
        else if (col.distance < safe_dist) {
        	//continue;
        	color = RaveVectorf(1,1,0,1);
        }
        //else color = RaveVectorf(0,1,0,1);
        //if (col.cctype == CCType_Between) {
        //  handles.push_back(env.drawarrow(col.ptB, col.ptB1, .0008, RaveVectorf(0,0,0,1)));
        //}
        OR::Vector ptB = (col.cctype == CCType_Between)  ? ((1-col.time)* col.ptB +col.time*col.ptB1) : col.ptB;
		handles.push_back(env.drawarrow(col.ptA, ptB, .0004, color));
      }
    }
    void Plot(const DblVec& x,
              OR::EnvironmentBase& env,
              std::vector<OR::GraphHandlePtr>& handles,
              boost::shared_ptr<BeliefCollisionEvaluator> m_calc,
              double m_dist_pen) {
      vector<BeliefCollision> collisions;
      m_calc->GetCollisionsCached(x, collisions);
      PlotCollisions(collisions, env, handles, m_dist_pen);
      m_calc->CustomPlot(x, handles);
    }
    
  };
}
