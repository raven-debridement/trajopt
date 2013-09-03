#ifndef __PLANNER_3D_H__
#define __PLANNER_3D_H__

#include "planner.h"

class Planner3D : public Planner
{
public:
	Planner3D() {};
	~Planner3D() {};

	// set RRT planner defaults for 3D planner
	void setPlannerDefaults();

	// Get the best available control input (given planning time)
	// plantime: specifies the time interval in real-time that one is prepared to invest in re-planning
	bool getControlSequence(std::list<Control>& clist);
	
private:
	// Utility function to select best control input if replanning is not possible
	// Select best initial twist that would get the needle closest to the goal
	void selectBestTwist(const Matrix<4,4>& T, const std::vector<Matrix<3,1> >& uvec, Matrix<4,4>& bestT, Matrix<3,1>& bestU);

	void testControls();
};

#endif