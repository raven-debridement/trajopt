#ifndef __PLANNER_2D_H__
#define __PLANNER_2D_H__

#include "planner.h"

class Planner2D : public Planner
{
public:
	Planner2D() {};
	~Planner2D() {};

	// set RRT planner defaults for 2D planner
	void setPlannerDefaults();

	// Get the best available control input (given planning time)
	// plantime: specifies the time interval in real-time that one is prepared to invest in re-planning
	bool getControlSequence(std::list<Control>& clist);

private:
	// for testing only
	void testControls();
};

#endif