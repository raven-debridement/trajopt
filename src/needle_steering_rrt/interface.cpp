#include "planner2D.h"
#include "planner3D.h"

#include "interface.h"

// Flag to run/pause the simulation
bool m_bRunSim = false;

// Key press handling routine
// Press 'space' to toggle between pausing and running the simulation
void keyFunc(char key, bool pressed){
  if (key == ' ' && pressed) {
    m_bRunSim = !m_bRunSim;
  }
}

Interface::Interface(PLANNER_TYPE ptype)
{
  m_bRunSim = false;

  if (ptype == PLANNER3D) {
    initPlanner3D();
  } else {
    std::cerr << "Planner type not applicable. Choose PLANNER3D" << std::endl;
    std::exit(-1);
  }
}

Interface::~Interface()
{
  m_pPlanner->cleanup();

  CAL_End();
}

void Interface::runNeedleSteering()
{
  //m_pPlanner->testControls();
  
  //int num;
  //std::cin >> num;

  // Main simulation loop
  bool terminate = false;
  while (!m_pPlanner->goalReached() && !terminate)
  {
    if (m_bRunSim)
    {
      std::cerr << "Running" << std::endl;
      senseAndUpdateNeedleTipPose();

      senseAndUpdateGoalPos();

      senseAndUpdateObstaclePos();

      std::list<Control> clist;
      terminate = m_pPlanner->getControlSequence(clist);

      m_needlepose = m_pPlanner->executeControlsWithNoise(m_needlepose, clist);
      //m_needlepose = m_pPlanner->executeControls(m_needlepose, clist);
      std::cout << m_needlepose << std::endl;

      //for(std::list<Control>::iterator cit = clist.begin(); cit != clist.end(); ++cit) {
      //  std::cout << "insert: " << cit->m_v*cit->m_dt << " rotate: " << cit->m_w*cit->m_dt << std::endl;
      //}
    }
  }
  std::cerr << "Completed";
}

// Initialize the 3D planner
void Interface::initPlanner3D()
{
  m_pPlanner = new Planner3D();
  m_pPlanner->bulletSetup();

  // Initialize callisto environment
  //CAL_Initialisation (true, true, true);
  //m_pPlanner->callistoSetup();
  //CAL_SetKeypressCallback(keyFunc);

  // Setup bounding box for environment
  Matrix<3,1> bbmin, bbmax;
  bbmin[0] = 0.0; bbmin[1] = -5.0; bbmin[2] = -3.0;
  bbmax[0] = 20.0; bbmax[1] = 5.0; bbmax[2] = 3.0;
  m_pPlanner->initBoundingBox(bbmin, bbmax);

  // Setup obstacles
  Matrix<3,1> c;

  /*
  c[0] = 0.0; c[1] = 0.0; c[2] = 3.0;
  Obstacle o1;
  o1.center = c;
  o1.obs_id = m_pPlanner->addSphericalObstacle(o1.center, 0.5);
  m_obstacles.push_back(o1);
  */

  //c[0] = 10.0; c[1] = 0.0; c[2] = 0.0;
  //Obstacle o2;
  //o2.center = c;
  //o2.obs_id = m_pPlanner->addBoxObstacle(o2.center, 2.0, 2.0, 2.0);
  //m_obstacles.push_back(o2);
  
  // Set needle radius (experimentally observed)
  m_pPlanner->setNeedleRadius(10.0); //1/cm

  m_pPlanner->setTimeStep(1.0); //sec
  
  // Set default insertion speed
  m_pPlanner->setDefaultInsertionSpeed(0.25); // cm/sec

  // Set default twist speed
  m_pPlanner->setDefaultTwistSpeed(8*M_PI); // rad/s

  // Set duty cycle period
  m_pPlanner->setDutyCyclePeriod(0.1); // s

  // Set initial needle pose
  m_needlepose = identity<4>();
  m_needlepose(0,3) = 0.0; m_needlepose(1,3) = 0.0; m_needlepose(2,3) = 0.0;
  m_pPlanner->setNeedleTipPose(m_needlepose);

  // Set initial goal location
  m_goalpos[0] = 18.0; m_goalpos[1] = -2.0; m_goalpos[2] = 0.0;
  m_pPlanner->setGoalPosition(m_goalpos);
  
  // set closed-loop execution
  m_pPlanner->closedLoopExecution(true);

  // Set planner defaults
  m_pPlanner->setPlannerDefaults();

  // Setup view params
  RaveInitialize(false, this->verbose ? Level_Debug : Level_Info);
  EnvironmentBasePtr env = RaveCreateEnvironment();
  OSGViewerPtr viewer;
  viewer = OSGViewer::GetOrCreate(env);
  assert(viewer);
  env->Load(string(DATA_DIR) + "/prostate.env.xml");
  viewer->SetAllTransparency(0.1);

  collision_checker.reset(new BulletCollisionChecker(env));

  //CAL_SetViewParams(0, (float)((bbmin[0] + bbmax[0])*0.5), (float)((bbmin[1] + bbmax[1])*0.5), 22.0f, (float)((bbmin[0] + bbmax[0])*0.5), (float)((bbmin[1] + bbmax[1])*0.5), (float)((bbmin[2] + bbmax[2])*0.5));
}

// TODO: Replace stub code with code that senses the current needle tip position and orientation and updates it for replanning
void Interface::senseAndUpdateNeedleTipPose()
{
  Matrix<4,4> needleposeNew = m_needlepose;
  m_pPlanner->updateNeedleTipPose(needleposeNew);
}

// TODO: Replace stub code with code that senses (displaced) goal position and updates it for replanning
void Interface::senseAndUpdateGoalPos()
{
  Matrix<3,1> goalposNew = m_goalpos;
  m_pPlanner->updateGoalPosition(goalposNew);
}

// TODO: Replace stub code with code that senses (displaced) obstacle centers and updates it for replanning
void Interface::senseAndUpdateObstaclePos()
{
  // Update centers for obstacles
  for(std::vector<Obstacle>::iterator it = m_obstacles.begin(); it != m_obstacles.end(); ++it)
  {
    Matrix<3,1> centerNew = it->center;
    m_pPlanner->updateObstaclePosition(it->obs_id, centerNew);
  }
}
