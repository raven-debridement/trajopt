#include "planner2D.h"

// set RRT planner defaults for 2D planner
void Planner2D::setPlannerDefaults()
{
	m_tstep = 1.0;
	
	m_factor = 2.0;

	m_plantime = 0.25;
	m_planbias = 0.2;

	m_maxiter = 1000;
	m_maxNNiter = 200;

	m_maxEvals = 10;

	m_maxchildren = 15;
	m_maxattempts = 5;

	m_reachedgoal = false;

	m_mtype = CLEARANCE;

	m_vmax = 1.0;
	m_wmax = 2.0*M_PI;

	m_ptype = PLANNER2D;
}

// Get the best available control input (given planning time)
bool Planner2D::getControlSequence(std::list<Control>& clist)
{
	m_rrttree.clear();
	m_rrtpath.clear();
	m_rrtpaths.clear();
	m_metric.clear();
	
	PathNode pn;
	bool terminate = false;
	bool retcode = false;

	Matrix<3,1> ctrl;
	bool flip;

	if (m_numIter == 0 || m_replanFlag) {
		retcode = buildRRT(m_rrttree, m_rrtpaths, m_plantime);
	}

	if (retcode)
	{
		m_solution.clear();
		m_currNode = 0;

		if (m_mtype == DISTANCE || m_mtype == CLEARANCE) {
			bestPath(m_rrttree, m_rrtpaths);
		} else {
			std::cerr << "Metric type not supported" << std::endl;
			std::exit(-1);
		}

		pn = m_solution[m_currNode];
		ctrl = pn.m_u;
		flip = pn.m_flip;

		//PathNode & pnxt = m_solution[m_currNode+1];
		//std::cout << "Expected pose:\n" << pnxt.m_T << std::endl;

		displayBestPath();
	} 
	else 
	{
		if (m_currNode < (int)m_solution.size()-2) {
			++m_currNode;
			pn = m_solution[m_currNode];

			ctrl = pn.m_u;
			flip = pn.m_flip;

		} else {
			terminate = true;
			//std::cout << "TERMINATED!" << std::endl;
		}
	}

	//std::cout << "PAUSED" << std::endl;
	//int num;
	//std::cin >> num;

	if (!terminate) {
		std::cout << "Planning matrix:\n" << pn.m_T << std::endl;
		getControls(ctrl, flip, clist);
	}

	++m_numIter;
	return terminate;
}

void Planner2D::testControls()
{
	initTree(m_rrttree, m_npose);
	
	Matrix<3,1> u;
	bool flip = false;

	Matrix<3> pt;
	
	pt[0] = 18; pt[1] = -2; pt[2] = 0;
	double d = dist(m_npose, pt);
	std::cout << "dist: " << d << std::endl;

	localPlanner(m_rrttree, 0, pt, u, flip, false); 
	std::cout << "u: " << u[0] << " " << u[1] << " " << u[2] << std::endl;
	std::cout << flip << std::endl;

	Matrix<4,4> nT = stepTDiscrete(m_npose, u, flip, 1.0);
	std::cout << "Projected transform matrix:\n" << nT << std::endl;

	std::list<Control> clist;
	
	/*
	u[0] = 3*M_PI*0.5; u[1] = 0.0; u[2] = 5;
	flip = true;

	Matrix<4,4> nT = stepTDiscrete(m_npose, u, flip, 1.0);
	std::cout << "Projected transform matrix:\n" << nT << std::endl;

	Control c;
	c.m_v = 0.0; c.m_w = M_PI; c.m_dt = 1.0;
	clist.push_back(c);
	c.m_v = 3.0*M_PI*0.5; c.m_w = 0.0; c.m_dt = 1.0;
	clist.push_back(c);

	m_npose = executeControls(m_npose, clist);
	*/

	getControls(u, flip, clist);

	//m_npose = executeControlsWithNoise(m_npose, clist);
	m_npose = executeControls(m_npose, clist);

	for(std::list<Control>::iterator cit = clist.begin(); cit != clist.end(); ++cit) {
		std::cout << "v: " << cit->m_v << " w: " << cit->m_w << " dt: " << cit->m_dt << std::endl;
	}
	
	std::cout << "New transform matrix:\n" << m_npose << std::endl;

	updateNeedleTipPose(m_npose);
}