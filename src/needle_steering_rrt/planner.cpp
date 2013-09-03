#include "planner.h"

// Set needle radius
void Planner::setNeedleRadius(double nradius) {
	m_nrad = nradius;
	m_kmax = 1.0/m_nrad;
}

// Set default insertion speed
void Planner::setDefaultInsertionSpeed(double v) {
	m_vsys = v;
}

// Set default twist speed
void Planner::setDefaultTwistSpeed(double w) {
	m_wsys = w;
}

// Set time step
void Planner::setTimeStep(double timestep) {
	m_tstep = timestep;
}

// Set duty cycle period
void Planner::setDutyCyclePeriod(double dcp) {
	m_tau = dcp;
}

// Set initial needle pose and goal position
void Planner::setNeedleTipPose(const Matrix<4,4> & T)
{
	m_npose = T;

	int obj_id;
	//CAL_CreateSphere(m_caltree, 0.2f, (float)m_npose(0,3), (float)m_npose(1,3), (float)m_npose(2,3), &obj_id);
	//CAL_SetObjectColor(obj_id, 0.8f, 0.1f, 0.1f);

	// Set y-axis
	yaxis[0] = m_npose(0,0); yaxis[1] = m_npose(1,0); yaxis[2] = m_npose(2,0);
	
	appliedTwist = 0;
	realTwist = 0;
	accumTwist = 0;
}

void Planner::setGoalPosition(const Matrix<3,1>& g) 
{
	// Setup goal
	m_goal = g;

	m_rthreshold = 0.1;

	int obj_id;
	//CAL_CreateSphere(m_caltree, (float)m_rthreshold, (float)m_goal[0], (float)m_goal[1], (float)m_goal[2], &obj_id);
	//CAL_SetObjectColor(obj_id, 0, 0, 1);
}


// Update obstacle position (tracker)
void Planner::updateObstaclePosition(int obsid, const Matrix<3,1> pt) 
{
	//CAL_SetObjectPosition(obsid, (float)pt[0], (float)pt[1], (float)pt[2]);
}

// Update routines for needle pose and goal position
void Planner::updateNeedleTipPose(const Matrix<4,4> & T) 
{
	// To prevent excessive windup in needle, we limit rotations to [-2pi, 2pi]
	double dp = yaxis[0]*T(0,1) + yaxis[1]*T(1,1) + yaxis[2]*T(2,1);
	if (dp > 1) {
		dp = 1;
	} else if (dp < -1) {
		dp = -1;
	}
	realTwist = acos(dp);

	// compare applied twist at each step to the real twist
	if (appliedTwist*realTwist < 0) {
		realTwist *= -1;
	}

	accumTwist += realTwist;

	double result = fmod(accumTwist, (2*M_PI));
	if (result < 0) { result += (2*M_PI); }
	accumTwist = result;

	//std::cout << "appliedTwist: " << appliedTwist << " realTwist: " << realTwist << " accumTwist: " << accumTwist << std::endl;

	yaxis[0] = T(0,1); yaxis[1] = T(1,1); yaxis[2] = T(2,1);

	//CAL_CreateSphere(m_caltrace, 0.1f, (float)T(0,3), (float)T(1,3), (float)T(2,3));

	double sc = 0.5;
	float line[6] = {(float)T(0,3), (float)T(1,3), (float)T(2,3), (float)(T(0,3) + sc*T(0,0)), (float)(T(1,3) + sc*T(1,0)), (float)(T(2,3) + sc*T(2,0))};
	int np[1] = {2};
	
	//glLineWidth(3.0);
	//CAL_CreatePolyline(m_caltrace, 1, np, line);

	// Update needle tip pose
	m_npose = T;
}

void Planner::updateGoalPosition(const Matrix<3,1>& g) {
	m_goal = g;
}

// Initialize Callisto environment (for visualization, collision detection and planning)
//void Planner::callistoSetup() 
//{
//	// Create Obstacles in Callisto
//	CAL_CreateGroup(&m_calenv, 0, true, "Environment");
//
//	CAL_CreateGroup(&m_calobs, m_calenv, true, "Obstacles");
//	CAL_SetGroupColor(m_calobs, 0.5, 0.5, 0.5);
//	
//	// Create Point in Callisto for distance calculation
//	CAL_CreateGroup(&m_calpoint, 0, true, "Point");
//	CAL_CreateSphere(m_calpoint, 0, 0, 0, 0);
//
//	// Create Box in Callisto
//	CAL_CreateGroup (&m_calbox, 0, false, "Box");
//	CAL_SetGroupColor(m_calbox, 0, 0, 0);
//
//	// Create Tree group in Callisto
//	CAL_CreateGroup (&m_caltree, 0, false, "Tree");
//	CAL_SetGroupColor(m_caltree, 0, 1, 0);
//
//	// Show the current best path
//	CAL_CreateGroup (&m_calbestpath, 0, false, "Best Path");
//	CAL_SetGroupColor(m_calbestpath, 0, 0, 1);
//
//	// Show execution trace
//	CAL_CreateGroup (&m_caltrace, 0, false, "Trace");
//	CAL_SetGroupColor(m_caltrace, 1, 0, 0);
//}

// Initialize environment bounding box
void Planner::initBoundingBox(const Matrix<3,1>& boxmin, const Matrix<3,1>& boxmax)
{
	//int np[1] = {2};

	m_bbmin = boxmin;
	m_bbmax = boxmax;
	
	//float side1[6] = {(float)m_bbmin[0], (float)m_bbmin[1], (float)m_bbmin[2], (float)m_bbmax[0], (float)m_bbmin[1], (float)m_bbmin[2]};
	//float side2[6] = {(float)m_bbmin[0], (float)m_bbmin[1], (float)m_bbmin[2], (float)m_bbmin[0], (float)m_bbmax[1], (float)m_bbmin[2]};
	//float side3[6] = {(float)m_bbmin[0], (float)m_bbmin[1], (float)m_bbmin[2], (float)m_bbmin[0], (float)m_bbmin[1], (float)m_bbmax[2]};

	//float side4[6] = {(float)m_bbmax[0], (float)m_bbmax[1], (float)m_bbmax[2], (float)m_bbmin[0], (float)m_bbmax[1], (float)m_bbmax[2]};
	//float side5[6] = {(float)m_bbmax[0], (float)m_bbmax[1], (float)m_bbmax[2], (float)m_bbmax[0], (float)m_bbmin[1], (float)m_bbmax[2]};
	//float side6[6] = {(float)m_bbmax[0], (float)m_bbmax[1], (float)m_bbmax[2], (float)m_bbmax[0], (float)m_bbmax[1], (float)m_bbmin[2]};
	//
	//float side7[6] = {(float)m_bbmin[0], (float)m_bbmax[1], (float)m_bbmax[2], (float)m_bbmin[0], (float)m_bbmin[1], (float)m_bbmax[2]};
	//float side8[6] = {(float)m_bbmin[0], (float)m_bbmax[1], (float)m_bbmax[2], (float)m_bbmin[0], (float)m_bbmax[1], (float)m_bbmin[2]};
	//float side9[6] = {(float)m_bbmax[0], (float)m_bbmin[1], (float)m_bbmax[2], (float)m_bbmin[0], (float)m_bbmin[1], (float)m_bbmax[2]};
	//
	//float side10[6] = {(float)m_bbmax[0], (float)m_bbmin[1], (float)m_bbmax[2], (float)m_bbmax[0], (float)m_bbmin[1], (float)m_bbmin[2]};
	//float side11[6] = {(float)m_bbmax[0], (float)m_bbmax[1], (float)m_bbmin[2], (float)m_bbmin[0], (float)m_bbmax[1], (float)m_bbmin[2]};
	//float side12[6] = {(float)m_bbmax[0], (float)m_bbmax[1], (float)m_bbmin[2], (float)m_bbmax[0], (float)m_bbmin[1], (float)m_bbmin[2]};

	//CAL_CreatePolyline(m_calbox, 1, np, side1);
	//CAL_CreatePolyline(m_calbox, 1, np, side2);
	//CAL_CreatePolyline(m_calbox, 1, np, side3);
	//CAL_CreatePolyline(m_calbox, 1, np, side4);
	//CAL_CreatePolyline(m_calbox, 1, np, side5);
	//CAL_CreatePolyline(m_calbox, 1, np, side6);
	//CAL_CreatePolyline(m_calbox, 1, np, side7);
	//CAL_CreatePolyline(m_calbox, 1, np, side8);
	//CAL_CreatePolyline(m_calbox, 1, np, side9);
	//CAL_CreatePolyline(m_calbox, 1, np, side10);
	//CAL_CreatePolyline(m_calbox, 1, np, side11);
	//CAL_CreatePolyline(m_calbox, 1, np, side12);

	//// Proxy obstacles
	//int m_calobsbox;
	//CAL_CreateGroup(&m_calobsbox, m_calobs, true, "", false, CAL_FALSE);
	//float xm = (float)((m_bbmin[0] + m_bbmax[0])*0.5);
	//float ym = (float)((m_bbmin[1] + m_bbmax[1])*0.5);
	//float zm = (float)((m_bbmin[2] + m_bbmax[2])*0.5);

	//float width = 0.1f;

	//float xw = (float)(m_bbmax[0] - m_bbmin[0]) + width;
	//float yw = (float)(m_bbmax[1] - m_bbmin[1]) + width;
	//float zw = (float)(m_bbmax[2] - m_bbmin[2]) + width;

	//CAL_CreateBox(m_calobsbox, width, yw, zw, (float)(m_bbmin[0]) - width, ym, zm);
	//CAL_CreateBox(m_calobsbox, width, yw, zw, (float)(m_bbmax[0]) + width, ym, zm);
	//CAL_CreateBox(m_calobsbox, xw, width, zw, xm, (float)(m_bbmin[1]) - width, zm);
	//CAL_CreateBox(m_calobsbox, xw, width, zw, xm, (float)(m_bbmax[1]) + width, zm);
	//CAL_CreateBox(m_calobsbox, xw, yw, width, xm, ym, (float)(m_bbmin[2]) - width);
	//CAL_CreateBox(m_calobsbox, xw, yw, width, xm, ym, (float)(m_bbmax[2]) + width);
}
	
// Add obstacles
// Currently, the user is allowed to add four kinds of obstacles: sphere, cylinder, box, triangle mesh

// Specify center and radius of spherical obstacle
int Planner::addSphericalObstacle(const Matrix<3,1>& center, double radius)
{
	int obs_id;
	CAL_CreateSphere(m_calobs, (float)radius, (float)center[0], (float)center[1], (float)center[2], &obs_id);
	return obs_id;
}

// Specify center, radius, height of cylinder (default orientation along +z axis)
// Also specify angle of rotation around x, y and z axes (rotations composed in that order)	
int Planner::addCylindricalObstacle(const Matrix<3,1>& center, double radius, double height, double xrot, double yrot, double zrot)
{
	int obs_id;
	CAL_CreateCylinder(m_calobs, (float)radius, (float)height, (float)center[0], (float)center[1], (float)center[2], &obs_id);
	CAL_SetObjectOrientation(obs_id, (float)xrot, (float)yrot, (float)zrot);
	return obs_id;
}

// Specify center and dimensions of box obstacle (axis-aligned)
int Planner::addBoxObstacle(const Matrix<3,1>& center, double xsize, double ysize, double zsize)
{
	int obs_id;
	CAL_CreateBox(m_calobs, (float)xsize, (float)ysize, (float)zsize, (float)center[0], (float)center[1], (float)center[2], &obs_id);
	return obs_id;
}

// Add triangle mesh (anatomical obstacles)
int Planner::addTriangleMesh(char * filename, char * label, bool colcheck)
{
	int obs_id;
	obs_id = objLoader(filename, label, m_calobs);
	return obs_id;
}

// add polygonal obstacle, only for 2D environments
int Planner::addPolygonalObstacle(int np, float* p, char* label)
{
	int obs_id;
	CAL_CreatePolygon(m_calobs, np, p, &obs_id, label);
	return obs_id;
}

// Simulate needle insertion step given control inputs and duration of simulation
Matrix<4,4> Planner::stepT(const Matrix<4,4>& T, const Control & c) 
{
	if (c.m_dt == 0) 
		return T;

	Matrix<3,1> v, w;
	Matrix<4,4> U;

	if (m_ptype == PLANNER2D) {
		w[0] = c.m_w; w[1] = 0.0; w[2] = c.m_v * m_kmax;
		v[0] = c.m_v; v[1] = 0.0; v[2] = 0.0;
	} else if (m_ptype == PLANNER3D) {
		w[0] = c.m_w; w[1] = c.m_v * m_kmax; w[2] = 0.0;
		v[0] = c.m_v; v[1] = 0.0; v[2] = 0.0;
	}

	U = zeros<4,4>();
	U.insert(0,0, cpMatrix(w));
	U.insert(0,3, v);
	
	Matrix<4,4> Tnxt = T*exp(c.m_dt*U);

	if (m_ptype == PLANNER2D) {
		Tnxt(0,2) = Tnxt(1,2) = 0;
		(Tnxt(2,2) < 0)? Tnxt(2,2) = -1 : Tnxt(2,2) = 1;
	}

	return Tnxt;
}

// Execute given control input (assuming idealized needle kinematic model)
Matrix<4,4> Planner::executeControls(const Matrix<4,4> & T, const std::list<Control>& clist)
{
	Matrix<4,4> nT = T;

	for(std::list<Control>::const_iterator cit = clist.begin(); cit != clist.end(); ++cit) {
		nT = stepT(nT, *cit);
		CAL_CreateSphere(m_caltrace, 0.07f, (float)nT(0,3), (float)nT(1,3), (float)nT(2,3));
	}

	// Reached goal? test
	Matrix<3,1> pt = nT.subMatrix<3,1>(0,3);
	double dg = tr(~(pt - m_goal) * (pt - m_goal));

	if (dg < m_rthreshold * m_rthreshold) {
		m_reachedgoal = true;
	}

	return nT;
}

Control Planner::perturbControls(const Control& c) 
{
	Control cnoise = c;
	
	if (cnoise.m_w != 0) {
		cnoise.m_w += normal()*cnoise.m_w*0.0025;
	}
	
	if (cnoise.m_v != 0) {
		cnoise.m_v += normal()*cnoise.m_v*0.0025; 
	}
	
	return cnoise;
}

// Execute given control input with artificially injected Gaussian noise
Matrix<4,4> Planner::executeControlsWithNoise(const Matrix<4,4> & T, const std::list<Control>& clist)
{
	Matrix<4,4> nT = T;
	for(std::list<Control>::const_iterator cit = clist.begin(); cit != clist.end(); ++cit) {
		nT = stepT(nT, perturbControls(*cit));
		CAL_CreateSphere(m_caltrace, 0.07f, (float)nT(0,3), (float)nT(1,3), (float)nT(2,3));
	}

	// Reached goal? test
	Matrix<3,1> pt = nT.subMatrix<3,1>(0,3);
	double dg = tr(~(pt - m_goal) * (pt - m_goal));

	if (dg < m_rthreshold * m_rthreshold) {
		m_reachedgoal = true;
	}

	return nT;
}

// Reached goal? test
bool Planner::goalReached() {
	return m_reachedgoal;
}

// Shutdown planner
void Planner::cleanup() 
{
	m_rrttree.clear();
	m_rrtpath.clear();
	m_rrtpaths.clear();
	m_solution.clear();
	m_metric.clear();
}

void Planner::initTree(std::vector<TreeNode>& tree, const Matrix<4,4>& pose) 
{
	TreeNode n;
	n.m_T = pose;

	n.m_bp = -1;

	n.m_depth = 0.0;
	n.m_clearance = getClearance(pose.subMatrix<3,1>(0,3));
	n.m_marked = false;
	n.m_attempts = 0;

	tree.push_back(n);
}

double Planner::dist(const Matrix<4,4>& T, const Matrix<3,1>& point) 
{
	Matrix<3,1> proj_point = ~T.subMatrix<3,3>(0,0) * (point - T.subMatrix<3,1>(0,3));
	double z = proj_point[2];
	if (z < 0) { 
		return INFTY;
	}
	double y = sqrt(proj_point[0]*proj_point[0] + proj_point[1]*proj_point[1]);
	if (y == 0) {
		return z;
	}
	double r = (z*z + y*y) / (2*y);
	if (r < m_nrad) {
		return INFTY;
	}
	return atan2(z, r-y) * r;
}

int Planner::nearestNeighbor(std::vector<TreeNode>& tree, const Matrix<3,1>& point) 
{
	int closest = -1;
	double mindist = INFTY;
	int col = 1;

	std::stack<int> st;
	st.push(0);
	int nc = (int)tree[0].m_children.size();

	for (int i = 0; i < nc; ++i) {
		st.push(tree[0].m_children[i]);
	}

	// DFS traversal of tree to find nearest neighbor (s.t. point lies in the reachable set of closest node)
	while (!st.empty()) 
	{
		int i = st.top();
		st.pop();

		TreeNode & tn = tree[i];
		double d = dist(tn.m_T, point);

		if ( ((int)tn.m_children.size() > m_maxchildren) || tn.m_marked ) {
			d = INFTY;
		}

		if (d < INFTY) 
		{
			double cost = m_factor*d + tn.m_depth;
			if ( cost < mindist) {
				col = 1;
				CAL_CheckLineCollision(m_calobs, (float)point(0,0), (float)point(1,0), (float)point(2,0), (float)tn.m_T(0,3), (float)tn.m_T(1,3), (float)tn.m_T(2,3), false, &col);
				if (col == 0) {
					closest = i;
					mindist = cost;
				}
			}

			int nc = (int) tn.m_children.size();
			for (int c = 0; c < nc; ++c) {
				st.push(tn.m_children[c]);
			}
		}
	}

	return closest;
}

// Local planner for connecting nodes in RRT planner
void Planner::localPlanner(std::vector<TreeNode>& tree, int node, const Matrix<3,1>& point, Matrix<3,1>& u, bool clamp) 
{
	Matrix<4,4> & T = tree[node].m_T;
	Matrix<3,1> proj_point = ~T.subMatrix<3,3>(0,0) * (point - T.subMatrix<3,1>(0,3));

	double x = proj_point[0];
	if (x < 0) { 
		//return INFTY;
		std::cerr << "Outside reachable set, dist INFTY" << std::endl;
		std::exit(-1);
	}

	double y = sqrt(proj_point[2]*proj_point[2] + proj_point[1]*proj_point[1]);
	if (y == 0) {
		//return x;
		//std::cout << "r: INFTY l: " << x << std::endl;
		return;
	}

	double r = (x*x + y*y) / (2*y);
	if (r < m_nrad) {
		//return INFTY;
		std::cerr << "Outside reachable set, dist INFTY" << std::endl;
		std::exit(-1);
	}
	
	// Output controls in terms of l,theta,r
	u[0] = atan2(x, r-y) * r;
	u[1] = atan2(proj_point[1], -proj_point[2]);
	u[2] = r;

	if (clamp && u[0] > m_vmax*m_tstep) {
		u[0] = m_vmax*m_tstep; 
	}
}

// Simulate needle insertion step given control inputs and duration of simulation (duty-cycling model)
Matrix<4,4> Planner::stepTDiscrete(const Matrix<4,4>& T, const Matrix<3,1>& u, double tfrac)
{
	double len = u[0];
	double rot = u[1];
	double rad = u[2];
	
	double phi = (len/rad)*tfrac;

	double cp = cos(phi);
	double sp = sin(phi);

	Matrix<3,1> pt, x, y, z;

  double ca = cos(rot);
  double sa = sin(rot);

  pt[0] = rad * sp;
  pt[1] = rad * sa * (1.0 - cp);
  pt[2] = -rad * ca * (1.0 - cp);

  x[0] = cp;
  x[1] = sa*sp;
  x[2] = -ca*sp;
  
  y[0] = 0;
  y[1] = ca;
  y[2] = sa;

  z[0] = sp;
  z[1] = -sa*cp;
  z[2] = ca*cp;

	Matrix<4,4> G;
	G(0,0) = x[0]; G(0,1) = y[0]; G(0,2) = z[0]; G(0,3) = pt[0];
	G(1,0) = x[1]; G(1,1) = y[1]; G(1,2) = z[1]; G(1,3) = pt[1];
	G(2,0) = x[2]; G(2,1) = y[2]; G(2,2) = z[2]; G(2,3) = pt[2];
	G(3,0) = 0.0;  G(3,1) = 0.0; G(3,2) = 0.0;   G(3,3) = 1.0;	

	return T*G;
}

// collision detection routines for RRT planner
bool Planner::checkCollision(std::vector<TreeNode>& tree, int node, const Matrix<3,1>& u, bool flip) 
{
	double len = u[0];
	double rot = u[1];
	double rad = u[2];

	double ca = cos(rot);
	double sa = sin(rot);

	double phi, cp, sp;
	double pt[3], curr[3];
	float line[6];
	int col;

	const Matrix<4,4>& T = tree[node].m_T;

	int nsegs;
	nsegs = 3 * (int)ceil(len/(m_vmax * m_tstep));

	double tfrac = 1.0/(double)nsegs;

	double (*ipts)[3] = new double[nsegs+1][3];
	ipts[0][0] = T(0,3); 
	ipts[0][1] = T(1,3);
	ipts[0][2] = T(2,3);

	bool valid = true;
	for(int si = 1; si <= nsegs && valid; ++si) 
	{
		phi	= (len/rad)*tfrac*si;

		cp = cos(phi);
		sp = sin(phi);

    pt[0] = rad * sp;
    pt[1] = rad * sa * (1.0 - cp);
    pt[2] = -rad * ca * (1.0 - cp);

		curr[0] = T(0,0) * pt[0] + T(0,1) * pt[1] + T(0,2) * pt[2] + T(0,3);
		curr[1] = T(1,0) * pt[0] + T(1,1) * pt[1] + T(1,2) * pt[2] + T(1,3);
		curr[2] = T(2,0) * pt[0] + T(2,1) * pt[1] + T(2,2) * pt[2] + T(2,3);

		line[0] = (float)ipts[si-1][0]; line[1] = (float)ipts[si-1][1]; line[2] = (float)ipts[si-1][2];
		line[3] = (float)curr[0]; line[4] = (float)curr[1]; line[5] = (float)curr[2];
		col = 1;
		CAL_CheckLineCollision(m_calobs, line[0], line[1], line[2], line[3], line[4], line[5], false, &col);
		if (col != 0) valid = false;

		ipts[si][0] = curr[0]; ipts[si][1] = curr[1]; ipts[si][2] = curr[2];
	}

	// if valid, add all pts for visualization
	/*
	if (valid) 
	{
	int np[1] = {2};
	for(int si = 1; si <= nsegs; ++si) 
	{
	line[0] = (float)ipts[si-1][0]; line[1] = (float)ipts[si-1][1]; line[2] = (float)ipts[si-1][2];
	line[3] = (float)ipts[si][0]; line[4] = (float)ipts[si][1]; line[5] = (float)ipts[si][2];
	//std::cout << line[0] << " " << line[1] << " " << line[2] << " " << line[3] << " " << line[4] << " " << line[5] << std::endl;

	//CAL_CreatePolyline(cal_tree, 1, np, line);	
	}
	}
	*/

	// cleanup
	delete[] ipts;

	return valid;
}

int Planner::addTreeNode(std::vector<TreeNode>& tree, int node, const Matrix<3,1>& u, bool flip, bool marked)
{
	TreeNode newnode;
	newnode.m_bp = node;

	newnode.m_u = u;
	newnode.m_flip = flip;

	newnode.m_T = stepTDiscrete(tree[node].m_T, newnode.m_u, flip, 1.0);

	newnode.m_depth = tree[node].m_depth + u[0];
	newnode.m_clearance = getClearance(newnode.m_T.subMatrix<3,1>(0,3));

	//std::cout << "node: " << tree.size() << " addTreeNode: " << newnode.T(0,3) << " " << newnode.T(1,3) << " " << newnode.T(2,3) << std::endl;

	newnode.m_marked = marked;
	newnode.m_attempts = 0;

	//float line[6] = {(float) tree[node].m_T(0,3), (float)tree[node].m_T(1,3), (float)tree[node].m_T(2,3), 
	//				 (float)newnode.m_T(0,3), (float)newnode.m_T(1,3), (float)newnode.m_T(2,3)};
	//int np[1] = {2};
	//CAL_CreatePolyline(m_caltree, 1, np, line);

	int newid = (int) tree.size();
	tree[node].m_children.push_back(newid);
	tree.push_back(newnode);

	return newid;
}

bool Planner::rrtStep(std::vector<TreeNode>& tree, std::vector<int>& paths) 
{
	int node;
	int tries = 0;

	// Pick Random Point
	Matrix<3,1> point;
	double rgoal = 0.5; //rthreshold;

	do {
		if (tries > 25) {
			rgoal = 1.0;
		}
		if (random() < m_planbias) {
			// Perturb goal (for reachability test)
			point[0] = m_goal[0] + rgoal*(2.0*random()-1.0);
			point[1] = m_goal[1] + rgoal*(2.0*random()-1.0);
			point[2] = m_goal[2] + rgoal*(2.0*random()-1.0);

		} else {
			point[0] = m_bbmin[0] + (m_bbmax[0] - m_bbmin[0])*random(); 
			point[1] = m_bbmin[1] + (m_bbmax[1] - m_bbmin[1])*random(); 
			point[2] = m_bbmin[2] + (m_bbmax[2] - m_bbmin[2])*random();
		}
		// Find Nearest Node
		node = nearestNeighbor(tree, point);
		//std::cout << "Point: " << point[0] << " " << point[1] << " " << point[2] << " Node: " << node << std::endl;
	} while ((node == -1) && (++tries < m_maxNNiter));

	// Handle case when no samples lies within reachable set 
	if (tries == m_maxNNiter) {
		return false;
	}

	//std::cout << "Point: " << point[0] << " " << point[1] << " " << point[2] << " Node: " << node << std::endl;

	Matrix<3,1> ctrl;
	//bool flip = false;
	localPlanner(tree, node, point, ctrl, true);

	bool valid = checkCollision(tree, node, ctrl);

	if (valid) 
	{
		int newid = addTreeNode(tree, node, ctrl, flip, false);

		// Reached goal? test
		TreeNode & newnode = tree[newid];
		Matrix<3,1> pnew = newnode.m_T.subMatrix<3,1>(0,3);
		double dg = tr(~(pnew - m_goal) * (pnew - m_goal));

		if (dg < (m_rthreshold * m_rthreshold)) 
		{
			// Debug!
			//CAL_CreateSphere(m_caltree, 0.1f, (float)pnew[0], (float)pnew[1], (float)pnew[2]);	
			//tree.back().marked = true;
			newnode.m_marked = true;
			tree[node].m_attempts++;
			tree[node].m_marked = true;
			paths.push_back(newid);
		}

		// Attempt to directly connect to the goal node
		if (dist(newnode.m_T, m_goal) < INFTY && !newnode.m_marked) 
		{

			//std::cout << "m_T\n" << newnode.m_T << std::endl;
			//std::cout << "m_goal: " << ~m_goal << std::endl;

			localPlanner(tree, newid, m_goal, ctrl, flip, false);

			//int num;
			//std::cin >> num;

			valid = checkCollision(tree, newid, ctrl, flip);

			if (valid) 
			{
				int nsegs = (int)ceil(ctrl[0]/(m_vmax * m_tstep));

				//double dfrac = ctrl[0]/(double)nsegs;
				double dfrac = ctrl[0] - (double)(nsegs-1)*(m_vmax*m_tstep);

				Matrix<3,1> ictrl;
				ictrl[0] = (m_vmax * m_tstep); ictrl[1] = ctrl[1]; ictrl[2] = ctrl[2];
				newid = addTreeNode(tree, newid, ictrl, flip, true);

				flip = false;
				ictrl[1] = 0.0;
				for(int si = 1; si < nsegs-1; ++si) {
					newid = addTreeNode(tree, newid, ictrl, flip, true);
				}

				ictrl[0] = dfrac;
				//newid = addTreeNode(newid, ictrl, true);

				tree[newid].m_attempts++;
				tree[newid].m_marked = true;
				int goalid = addTreeNode(tree, newid, ictrl, flip, true);

				paths.push_back(goalid);

				//CAL_CreateSphere(m_caltree, 0.1f, (float)tree[goalid].m_T(0,3), (float)tree[goalid].m_T(1,3), (float)tree[goalid].m_T(2,3));

				//std::cout << "goalid: " << goalid << " bp: " << tree[goalid].bp << std::endl;
				//std::cout << "pT: " << tree[tree[goalid].bp].T(0,3) << " " << tree[tree[goalid].bp].T(1,3) << " " << tree[tree[goalid].bp].T(2,3) << std::endl;
				//std::cout << "gT: " << tree[goalid].T(0,3) << " " << tree[goalid].T(1,3) << " " << tree[goalid].T(2,3) << std::endl;
			}
		}
	} // if valid

	return true;
}

bool Planner::buildRRT(std::vector<TreeNode>& tree, std::vector<int>& paths, double plantime) 
{
	// Build Tree
	clock_t startClk = clock();

	bool stepRet = true;

	// Initialize tree
	initTree(tree, m_npose);

	for (int i = 0; stepRet && (((double)(clock() - startClk) / CLOCKS_PER_SEC) < plantime); ++i) {
		stepRet = rrtStep(tree, paths);
	}

	if (stepRet && !paths.empty()) 
	{
		std::cout << "\nRRT build: " << (double) (clock() - startClk) / CLOCKS_PER_SEC << std::endl;
		int numpaths = (int)paths.size();
		std::cout << "Num paths found: " << numpaths << std::endl;
	} 
	else {
		std::cerr << "Unable to find solution, reusing previous solution" << std::endl;
	}

	//std::cout << "RRT BUILD PAUSED" << std::endl;
	//int num;
	//std::cin >> num;

	return (stepRet && !paths.empty());
}

// Control proceeds in three stages:
// Stage 1: Rotate needle tip in place
// Stage 2: Duty-cycle needle stage I
// Stage 3: Duty-cycle needle stage II

void Planner::getControls(const Matrix<3,1>& u, bool flip, std::list<Control>& clist)
{
	//std::cout << "l: " << u[0] << " theta: " << u[1] << " r: " << u[2] << std::endl;

	clist.clear();

	Control c;

	double dc = (1.0 - (m_nrad/u[2]));

	if ((dc < 0) || (dc > 1)) {
		std::cerr << "Duty cycling factor falls outside limits [0,1]" << std::endl;
		std::exit(-1);
	}

	double theta;
	
	if (m_ptype == PLANNER2D) {
		(flip)? theta = M_PI : theta = 0.0;
	} else {
		theta = u[1];
	}

	double twist = theta;
	if ((accumTwist + theta) > 2.0*M_PI) {
		twist = (theta - 2.0*M_PI);
	}
	else if ((accumTwist + theta) < -2.0*M_PI) {
		twist = (2.0*M_PI + theta);
	}
	
	// First stage: rotate in place
	if ( ((m_ptype == PLANNER2D) && flip) || m_ptype == PLANNER3D) {
		(twist > 0)? c.m_w = m_wsys : c.m_w = -m_wsys;
		//c.m_dt = theta/c.m_w; // abs needed? 
		c.m_dt = twist/c.m_w; // abs needed? 
		c.m_v = 0.0;
		clist.push_back(c);
	}

	// duty-cycle period (tau) (round off to nearest integer)
	int dcPeriods = (int)((m_tstep/m_tau) + 0.5);

	double dclen = u[0]/(double)dcPeriods;
	
	int rotCycles = 2;

	double tdel = (2.0*2.0*M_PI)/m_wsys;

	double vdc = (dclen * dc)/tdel;

	while (vdc > 0.4) {
		tdel *= 2.0;
		vdc = (dclen * dc)/tdel;
		rotCycles *= 2;
	}
	//std::cout << "vdc: " << vdc << " rotCycles: " << rotCycles << std::endl;

	double invrotcycles = 1.0/(double)rotCycles;

	double twistsgn;
	((accumTwist + twist) < 0)? twistsgn = 1 : twistsgn = -1;

	for(int i = 0; i < dcPeriods; ++i) 
	{
		c.m_v = m_vsys; 
		c.m_dt = (dclen*(1.0 - dc))/c.m_v; 
		c.m_w = 0.0;
		clist.push_back(c);

		for(int j = 0; j < rotCycles; j+=2) {
			c.m_w = twistsgn*m_wsys;
			c.m_dt = (2.0*M_PI)/m_wsys;
			c.m_v = (dclen*dc*invrotcycles)/c.m_dt;
			clist.push_back(c);

			c.m_w = -twistsgn*m_wsys;
			c.m_dt = (2.0*M_PI)/m_wsys;
			c.m_v = (dclen*dc*invrotcycles)/c.m_dt;
			clist.push_back(c);
		}
	}
}

double Planner::getClearance(const Matrix<3,1>& pt) 
{
	CAL_SetGroupPosition(m_calpoint, (float)pt[0], (float)pt[1], (float)pt[2]);

	int num_pairs;
	CAL_GetClosestPairs (m_calpoint, m_calenv, &num_pairs);
	SCALResult* results = new SCALResult[num_pairs];
	CAL_GetResults(results);
	double distance = results[0].distance;
	delete[] results;

	return distance;
}

double Planner::evalMetric(const TreeNode& tu) 
{
	double mval;
	switch(m_mtype)
	{
	case DISTANCE:
		mval = tu.m_depth;
		break;
	case CLEARANCE:
		mval = -(tu.m_clearance);
		break;
	default:
		std::cerr << "Metric type unsupported" << std::endl;
		std::exit(-1);
	}
	return mval;
}

int Planner::dijkstraSearch(const std::vector<TreeNode>& tree, const std::vector<int>& paths) 
{
	std::multimap<double, int> Q;

	int nn = (int)tree.size();
	m_metric.assign(nn, std::make_pair(INFTY, -1));
	std::vector<std::multimap<double, int>::iterator> pos_in_Q(nn, Q.end());

	const TreeNode & root = tree[0];
	int nc = (int)root.m_children.size();
	for (int j = 0; j < nc; ++j) 
	{
		int u = root.m_children[j];
		const TreeNode & tu = tree[u];
		double dnbr = evalMetric(tu);
		m_metric[u] = std::make_pair(dnbr, -1);
		pos_in_Q[u] = Q.insert(std::make_pair(dnbr, u));
	}

	int u, v;
	int end;
	bool found = false;
	int attempts = -1;

	while (!Q.empty() && !found) 
	{
		u = Q.begin()->second;

		if(std::find(paths.begin(), paths.end(), u) != paths.end()) 
		{
			++attempts;
			if (attempts == (int)(0.015*paths.size() + 0.5)) {
				end = u;
				found = true;
			}
		}

		Q.erase(Q.begin());
		pos_in_Q[u] = Q.end();

		const TreeNode & tu = tree[u];
		int nc = (int)tu.m_children.size();

		for (int j = 0; j < nc; ++j) 
		{
			v = tu.m_children[j];
			const TreeNode & tv = tree[v];
			double metric_uv;

			if (m_mtype == DISTANCE) 
			{
				metric_uv = evalMetric(tv);
				if (m_metric[v].first > (m_metric[u].first + metric_uv)) 
				{
					m_metric[v] = std::make_pair(m_metric[u].first + metric_uv, u);
					if (pos_in_Q[v] == Q.end()) {
						pos_in_Q[v] = Q.insert(std::make_pair(m_metric[v].first, v));
					} else {
						pos_in_Q[v] = Q.insert(Q.erase(pos_in_Q[v]), std::make_pair(m_metric[v].first, v));
					}
				}
			}
			else if (m_mtype == CLEARANCE) 
			{
				metric_uv = std::min(evalMetric(tu), evalMetric(tv));

				if (m_metric[v].first > std::max(m_metric[u].first, metric_uv)) 
				{
					m_metric[v] = std::make_pair(std::max(m_metric[u].first, metric_uv), u);

					if (pos_in_Q[v] == Q.end()) {
						pos_in_Q[v] = Q.insert(std::make_pair(m_metric[v].first, v));
					} else {
						pos_in_Q[v] = Q.insert(Q.erase(pos_in_Q[v]), std::make_pair(m_metric[v].first, v));
					}
				}
			}
		}
	}
	return end;
}

void Planner::bestPath(const std::vector<TreeNode>& tree, const std::vector<int>& paths) 
{
	// Search for best path
	int goalid = -1;
	if (!paths.empty()) {
		goalid = dijkstraSearch(tree, paths);
	} else {
		std::cerr << "No solution found!" << std::endl;
		std::exit(-1);
	}

	// Clear previous best path
	m_solution.clear();

	int node = goalid;
	std::vector<int> path;
	while (node != -1) 
	{
		path.push_back(node);
		node = tree[node].m_bp;
	}
	std::reverse(path.begin(), path.end());

	// Iterate over all path segments
	int ns = (int)path.size();
	m_solution.resize(ns);

	for(int si = 0; si < ns-1; ++si) 
	{
		const TreeNode & ti = tree[path[si]];
		const TreeNode & tnxt = tree[path[si+1]];

		PathNode & pn = m_solution[si];

		pn.m_T = ti.m_T;
		pn.m_u = tnxt.m_u;
		pn.m_flip = tnxt.m_flip;
	}
}

//void Planner::displayBestPath()
//{
//	// Clear out group for visualization
//	CAL_EmptyGroup(m_calbestpath);
//
//	// Iterate over all path segments
//	int ns = (int)m_solution.size();
//
//	int nsegs = 5;
//	double tfrac = 1.0/(double)nsegs;
//	int objid;
//
//	for(int i = 0; i < ns-1; ++i) 
//	{
//		PathNode & pi = m_solution[i];
//
//		Matrix<3,1> curr;
//		const Matrix<4,4>& T = pi.m_T;
//
//		double len = pi.m_u[0];
//		double rot = pi.m_u[1];
//		double rad = pi.m_u[2];
//		bool flip = pi.m_flip;
//
//		double cw = cos(rot);
//		double sw = sin(rot);
//
//		double phi, cp, sp;
//		double pt[3];
//
//		curr[0] = T(0,3); curr[1] = T(1,3); curr[2] = T(2,3);
//
//		CAL_CreateSphere(m_calbestpath, (float)0.025, (float)curr[0], (float)curr[1], (float)curr[2], &objid);
//		//std::cout << "Creating sphere at: " << curr[0] << ", " << curr[1] << ", "  << curr[2] << std::endl;
//
//		for(int si = 1; si <= nsegs; ++si)
//		{
//			phi	= (len/rad)*tfrac*si;
//
//			cp = cos(phi);
//			sp = sin(phi);
//
//			if (m_ptype == PLANNER2D) {
//				if (flip) {
//					pt[0] = rad * sp;
//					pt[1] = -rad * (1.0 - cp);
//					pt[2] = 0.0;
//				}
//				else {
//					pt[0] = rad * sp;
//					pt[1] = rad * (1.0 - cp);
//					pt[2] = 0.0;
//				}
//			}
//			else if (m_ptype == PLANNER3D) {
//				pt[0] = rad * sp;
//				pt[1] = rad * sw * (1.0 - cp);
//				pt[2] = -rad * cw * (1.0 - cp);
//			}
//			
//			curr[0] = T(0,0) * pt[0] + T(0,1) * pt[1] + T(0,2) * pt[2] + T(0,3);
//			curr[1] = T(1,0) * pt[0] + T(1,1) * pt[1] + T(1,2) * pt[2] + T(1,3);
//			curr[2] = T(2,0) * pt[0] + T(2,1) * pt[1] + T(2,2) * pt[2] + T(2,3);
//
//			CAL_CreateSphere(m_calbestpath, (float)0.025, (float)curr[0], (float)curr[1], (float)curr[2], &objid);
//			//std::cout << "Creating sphere at: " << curr[0] << ", " << curr[1] << ", "  << curr[2] << std::endl;
//		}
//	}
//}
