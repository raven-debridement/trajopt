import trajoptpy
import openravepy as rave
import numpy as np
import json
import trajoptpy.math_utils as mu
import trajoptpy.kin_utils as ku
import trajoptpy.make_kinbodies as mk

def move_arm_to_grasp(xyz_targ, quat_targ, link_name, manip_name):
    
    request = {
        "basic_info" : {
            "n_steps" : 10,
            "manip" : manip_name,
            "start_fixed" : True,
            "belief_space" : True
        },
        "costs" : [
        {
            "type" : "collision",
            "params" : {"coeffs" : [30],"dist_pen" : [0.05], "belief_space" : True}
        },
#        {
#            "type" : "joint_vel",
#            "params": {"coeffs" : [1]}
#        },
        {
            "type" : "control",
            "params": {"coeffs" : [1]}
        },
        {
            "type" : "covariance",
            "params": {
                "Q" : (np.eye(3)*1).tolist()
            }
        },
        ],
        "constraints" : [
        {
            "type" : "pose",
            "name" : "final_pose",
            "params" : {
                "pos_coeffs" : [1,1,1],
                "rot_coeffs" : [1,1,1],
                "xyz" : list(xyz_targ),
                "wxyz" : list(quat_targ),
                "link" : link_name,
            },
        },
#        {
#            "type" : "collision",
#            "params" : {"coeffs" : [1],"dist_pen" : [0], "belief_space" : True}
#        },
        {
            "type" : "control",
            "params": {
                "u_min" : -0.4,
                "u_max" : 0.4
            }
        },
        ],
        "init_info" : {
            "type" : "stationary",
            "initial_rt_sigma" : (np.eye(3)).tolist()
        }
    }
    
#    traj_init = np.load("/home/alex/Desktop/traj.npy")
#    request["init_info"]["type"] = "given_traj"
#    request["init_info"]["data"] = [x.tolist() for x in traj_init]
    
    return request



if __name__ == "__main__":
        
    ### Parameters ###
    ENV_FILE = "../data/three_links.env.xml"
    MANIP_NAME = "arm"
    LINK_NAME = "Finger"
    INTERACTIVE = True
    ##################
    
    ### Env setup ####
    env = rave.RaveGetEnvironment(1)
    if env is None:
        env = rave.Environment()
        env.StopSimulation()
        env.Load(ENV_FILE)
    robot = env.GetRobots()[0]
    robot.SetDOFValues([ 0 ,0 ,0]);
    manip = robot.GetManipulator(MANIP_NAME)
    viewer = trajoptpy.GetViewer(env)
    viewer.SetCameraTransformation([0,0,2], [0,0,0], [0,1,0])
    ##################
    T_gripper = robot.GetLink(LINK_NAME).GetTransform()
    T_grasp = T_gripper.copy()
    T_grasp[:3,3]  += np.array([-0.5,0.2,0])
    T_grasp = T_grasp.dot(rave.matrixFromAxisAngle([0,0,1],np.pi/2))
    xyz_targ = T_grasp[:3,3]
    #success = mk.create_cylinder(env, xyz_targ-np.array([.03,0,0]), .02, .5)
    quat_targ = rave.quatFromRotationMatrix(T_grasp[:3,:3])
    #success = mk.create_cylinder(env, T_gripper[:3,3]-np.array([.1,-.1,0]), .02, .5)

    mk.create_mesh_box(env, np.array([0.2,0.25,0]), np.array([0.05,0.05,0.3]), "box1")

    request = move_arm_to_grasp(xyz_targ, quat_targ, LINK_NAME, MANIP_NAME)
    
#    ##################
#    # first optimize ignoring collision costs
#    all_costs = request["costs"]
#    all_constraints = request["constraints"]
#    noncollision_costs = [cost for cost in request["costs"] if "collision" not in cost["type"]]
#    noncollision_constraints = [constraint for constraint in request["constraints"] if "collision" not in constraint["type"]]
#    
#    request["costs"] = noncollision_costs
#    request["constraints"] = noncollision_constraints
#    
#    saver = rave.Robot.RobotStateSaver(robot)
#    s = json.dumps(request)
#    print "REQUEST:",s
#    trajoptpy.SetInteractive(False);
#    prob = trajoptpy.ConstructProblem(s, env)
#    result = trajoptpy.OptimizeProblem(prob)
#    del saver # reverts the robot state
#
#    # add collision cost and constraint back again
#    request["costs"] = all_costs
#    request["constraints"] = all_constraints
#    
#    # use the resulting trajectory as initialization for the new optimization that includes collision costs
#    path_init = result.GetTraj()
#    request["init_info"]["type"] = "given_traj"
#    request["init_info"]["data"] = [x.tolist() for x in path_init]
#    ##################
    
    s = json.dumps(request)
    print "REQUEST:",s
    trajoptpy.SetInteractive(INTERACTIVE);
    prob = trajoptpy.ConstructProblem(s, env)
    result = trajoptpy.OptimizeProblem(prob)

    #raw_input('Hit ENTER to continue.')	
    
#    for i in xrange(10):
#        request["init_info"]["type"] = "given_traj"
#        request["init_info"]["data"] = [x.tolist() for x in result.GetTraj()]
#        trajoptpy.SetInteractive(INTERACTIVE);
#        prob = trajoptpy.ConstructProblem(json.dumps(request), env)
#        result = trajoptpy.OptimizeProblem(prob)
#        
#    traj = result.GetTraj()
#    dofs = traj[0,:3]
#    u = traj[0,-3:]
#    dofs1 = prob.RobotDynamics(dofs,u)

#    np.save("/home/alex/Desktop/traj.npy", result.GetTraj())
