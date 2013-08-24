#ifndef __INTERFACE_H__
#define __INTERFACE_H__

#include "planner.h"

class Interface
{
public:
    Interface(PLANNER_TYPE ptype);
    ~Interface();

	// planning/control loop
    void runNeedleSteering();

private:
	
    void initPlanner2D();
	void initPlanner3D();
  
	void senseAndUpdateNeedleTipPose();
    void senseAndUpdateGoalPos();
    void senseAndUpdateObstaclePos();

    // Planner class interface
    Planner* m_pPlanner;
    
	// Needle tip pose (position and orientation encoded as a 4x4 matrix)
    Matrix<4,4> m_needlepose;
    // Goal position (in R3)
    Matrix<3,1> m_goalpos;

    // Temporary struct for describing obstacle (obstacle id and center of the obstacle where a magnetic tracker is placed)
    typedef struct Obstacle
    {
        int obs_id;
        Matrix<3,1> center;
    } Obstacle;
    std::vector<Obstacle> m_obstacles;
};

#endif