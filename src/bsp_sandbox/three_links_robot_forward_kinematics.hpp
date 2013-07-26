#pragma once

#include "three_links_robot.hpp"

namespace ThreeLinksRobotBSP {
  Vector3d forward_kinematics_test(const StateT& state) {
    double t0 = state(0),
           t1 = state(1),
           t2 = state(2);
    Vector3d eetrans;
    eetrans(0) = 0.16*cos(t0) + 0.16*cos(t0+t1) + 0.08*cos(t0+t1+t2);
    eetrans(1) = 0.16*sin(t0) + 0.16*sin(t0+t1) + 0.08*sin(t0+t1+t2);
    eetrans(2) = 0;
    return eetrans;
    
  }
  
  Vector3d forward_kinematics(const StateT& state) {
    Vector3d eetrans;
    double x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11;
		x0=cos((double)state[0]);
		x1=cos((double)state[1]);
		x2=sin((double)state[0]);
		x3=sin((double)state[1]);
		x4=sin((double)state[2]);
		x5=cos((double)state[2]);
		x6=0.160000000000000*x2;
		x7=0.0800000000000000*x2;
		x8=0.160000000000000*x0;
		x9=0.0800000000000000*x0;
		x10=((((x1)*(x9)))+((-1.00000000000000*(x3)*(x7))));
		x11=((((x1)*(x7)))+(((x3)*(x9))));
		eetrans[0]=((((x1)*(x8)))+(x8)+((-1.00000000000000*(x11)*(x4)))+((-1.00000000000000*(x3)*(x6)))+(((x10)*(x5))));
		eetrans[1]=((((x1)*(x6)))+(x6)+(((x11)*(x5)))+(((x3)*(x8)))+(((x10)*(x4))));
		eetrans[2] = 0;

    cout << "kinematics: " << eetrans.transpose() << "; " << "other: " << forward_kinematics_test(state).transpose() << endl;
    return eetrans;
  }
  
}
