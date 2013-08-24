#ifndef __PLANNER_H__
#define __PLANNER_H__

#define _CRT_RAND_S

#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <list>
#include <stack>
#include <map>
#include <set>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <time.h>

#include "matrix.h"
#include "utils.h"

enum METRIC_TYPE {DISTANCE, CLEARANCE};
enum PLANNER_TYPE {PLANNER2D, PLANNER3D};

#define INFTY 9e9

struct Control 
{
	double m_v, m_w, m_dt;

	Control() {
		m_v = m_w = m_dt = 0.0;
	}

	Control& operator = (const Control& c) {
		m_v = c.m_v;
		m_w = c.m_w;
		m_dt = c.m_dt;
		return (*this);
	}
};

class Planner
{
public:
	Planner() {};
	~Planner() {};

	// Set needle radius
	void setNeedleRadius(double nradius);

	// Set default insertion speed
	void setDefaultInsertionSpeed(double v);
	// Set default twist speed
	void setDefaultTwistSpeed(double w);
	// Set time step
	void setTimeStep(double timestep);
	// Set duty cycle period
	void setDutyCyclePeriod(double dcp);

	// Set initial needle tip pose (4x4 matrix describing the position and orientation of the needle-tip)
	void setNeedleTipPose(const Matrix<4,4> & T);
	// Set initial goal position (the goal position is expected to change due to deformation of tissue as needle is inserted)
	void setGoalPosition(const Matrix<3,1>& g);

	// Update workspace (sensing)

	// Update obstacle position (track the center of the obstacle (sphere, box, cylinder))
	void updateObstaclePosition(int obsid, const Matrix<3,1> pt);
	
	// Update routines for needle pose and goal position
	void updateNeedleTipPose(const Matrix<4,4> & T);
	void updateGoalPosition(const Matrix<3,1>& g);

	// Reached goal test (required for termination)
	bool goalReached();

	// Clean up planner (after planning or on exit)
	void cleanup();
		
	// Initialize Callisto environment (for visualization, collision detection and planning)
	void callistoSetup();
	
	// Initialize environment bounding box
	void initBoundingBox(const Matrix<3,1>& boxmin, const Matrix<3,1>& boxmax);

	// Currently, the user is allowed to add four kinds of obstacles: sphere, cylinder, box, triangle mesh defining obstacle
	// Specify center and radius of spherical obstacle
	int addSphericalObstacle(const Matrix<3,1>& center, double radius);
	// Specify center, radius, height of cylinder (default orientation along +z axis)
	// Also specify angle of rotation around x, y and z axes (rotations composed in that order)
	int addCylindricalObstacle(const Matrix<3,1>& center, double radius, double height, double xrot, double yrot, double zrot);
	// Specify center and dimensions of box obstacle (axis-aligned)
	int addBoxObstacle(const Matrix<3,1>& center, double xsize, double ysize, double zsize);
	// input triangle mesh (anatomical obstacles)
	int addTriangleMesh(char * filename, char * label, bool colcheck);
		
	// add polygonal obstacle, only for 2D environments
	int addPolygonalObstacle(int np, float* p, char* label);

	// Execute given control inputs (assuming idealized needle kinematic model)
	Matrix<4,4> executeControls(const Matrix<4,4> & T, const std::list<Control>& clist);
	// Execute given control inputs with artificially injected Gaussian noise
	Matrix<4,4> executeControlsWithNoise(const Matrix<4,4> & T, const std::list<Control>& clist);

	// Execute closed loop (replanning)
	void closedLoopExecution(bool replan) { m_replanFlag = replan;}

public:
	// Initialize planner defaults (planner specific)
	virtual void setPlannerDefaults() = 0;

	// Get the best available control input given current needle pose
	virtual bool getControlSequence(std::list<Control>& clist) = 0;

	// Debug routine to test controls
	virtual void testControls() = 0;
	
protected:
	// Planning specific classes
	struct TreeNode 
	{
		Matrix<4,4> m_T;
		Matrix<3,1> m_u;
		bool m_flip;

		int m_bp;

		bool m_marked;
		int m_attempts;

		double m_depth;
		double m_clearance;

		std::vector<int> m_children;

		TreeNode() 
		{
			m_bp = -1;
			m_marked = false;
			m_depth = 0.0;
			m_clearance = INFTY;
			m_attempts = 0;
		}
	};

	struct PathNode 
	{
		Matrix<4,4> m_T;
		Matrix<3,1> m_u;
		bool m_flip;
	};

	// Routines to select a good path from a set of candidate paths
	// Compute clearance (distance to nearest obstacle)
	double getClearance(const Matrix<3,1>& pt);
	double evalMetric(const TreeNode& tu);

	// Dijkstra implementation (two source paths based on user-specified metric type)
	int dijkstraSearch(const std::vector<TreeNode>& tree, const std::vector<int>& paths);
	
	// Initialize RRT tree for planning
	void initTree(std::vector<TreeNode>& tree, const Matrix<4,4>& pose);

	// RRT routines
	double dist(const Matrix<4,4>& T, const Matrix<3,1>& point);
	int nearestNeighbor(std::vector<TreeNode>& tree, const Matrix<3,1>& point);
	
	// Simulate needle insertion step given control inputs and duration of simulation (duty-cycling model)
	Matrix<4,4> stepTDiscrete(const Matrix<4,4>& T, const Matrix<3,1>& u, bool flip, double tfrac);

	// Local planner for connecting nodes in RRT planner
	void localPlanner(std::vector<TreeNode>& tree, int node, const Matrix<3,1>& point, Matrix<3,1>& u, bool& flip, bool clamp);

	// collision detection routines for RRT planner
	bool checkCollision(std::vector<TreeNode>& tree, int node, const Matrix<3,1>& u, bool flip);

	// add tree node to RRT 
	int addTreeNode(std::vector<TreeNode>& tree, int node, const Matrix<3,1>& u, bool flip, bool marked);

	// one iteration of the RRT planner
	bool rrtStep(std::vector<TreeNode>& tree, std::vector<int>& paths); 
	
	// build RRT (3D)
	bool buildRRT(std::vector<TreeNode>& tree, std::vector<int>& paths, double plantime);

	// Choose best path from a set of candidate paths
	void bestPath(const std::vector<TreeNode>& tree, const std::vector<int>& paths);

	// display best path
	void displayBestPath();

	// Translate (l,r,w) inputs into low-level control inputs (duty-cycling)
	void getControls(const Matrix<3,1>& u, bool flip, std::list<Control>& clist);
	
	// Simulate needle insertion step given control inputs and duration of simulation
	Matrix<4,4> stepT(const Matrix<4,4>& T, const Control & c);
	
	// Perturb control by artificially injecting Gaussian noise (to simulate noisy execution)
	Control perturbControls(const Control& c);
	
protected:
	Matrix<4,4> m_npose;
	Matrix<3,1> m_goal;

	Matrix<3,1> m_bbmin, m_bbmax;

	double m_vmax, m_wmax, m_kmax;
	double m_nrad;

	// default insertion and rotation speeds from system
	double m_vsys, m_wsys;
	// duty cycle period 
	double m_tau;

	// flag to turn replanning on and off
	bool m_replanFlag;
	
	double m_tstep;

	// Termination condition
	bool m_reachedgoal;
	double m_rthreshold;
	
	// Callisto members
	int m_caltree;
	int m_calbox;
	int m_calobs;
	int m_calenv;
	int m_calpoint;
	
	int m_calbestpath;
	int m_caltrace;

	// planning variables
	double m_factor;
	double m_planbias;

	int m_maxchildren;
	int m_maxattempts;
	int m_maxiter;
	int m_maxNNiter;

	int m_maxEvals;
	
	double m_plantime;

	std::vector<PathNode> m_rrtpath;
	std::vector<TreeNode> m_rrttree;
	std::vector<int> m_rrtpaths;

	// For Dijkstra's search
	METRIC_TYPE m_mtype;
	std::vector<std::pair<double, int> > m_metric;

	// Planner type
	PLANNER_TYPE m_ptype;
	
	// Best path
	int m_currNode;
	std::vector<PathNode> m_solution;

	// keep track of number of planning iterations
	int m_numIter;

	// Y axis (to keep track of net twist)
	Matrix<3,1> yaxis;
	double appliedTwist, realTwist, accumTwist;
};

#endif
