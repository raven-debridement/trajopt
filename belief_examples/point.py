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
            "n_steps" : 20,
            "manip" : "base_point",
            "start_fixed" : True,
            "belief_space" : True
        },
        "costs" : [
        {
            "type" : "collision",
            "params" : {"coeffs" : [30],"dist_pen" : [0.2], "belief_space" : True}
        },
#        {
#            "type" : "joint_vel",
#            "params": {"coeffs" : [1]}
#        },
        {
            "type" : "control",
            "params": {"coeffs" : [0.1]}
        },
        {
            "type" : "covariance",
            "params": {
                "Q" : (np.eye(2)*1).tolist()
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
            }
        },
        {
            "type" : "control",
            "params": {
                "u_min" : -0.5,
                "u_max" : 0.5
            }
        },
        ],
        "init_info" : {
            "type" : "stationary",
            "initial_rt_sigma" : (np.eye(2)).tolist()
        }
    }
    
#    traj_init = np.load("/home/alex/Desktop/traj.npy")
#    request["init_info"]["type"] = "given_traj"
#    request["init_info"]["data"] = [x.tolist() for x in traj_init]
    
    return request



if __name__ == "__main__":
        
    ### Parameters ###
    ENV_FILE = "../data/point.env.xml"
    MANIP_NAME = "base"
    LINK_NAME = "Base"
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
    viewer.SetCameraTransformation([0,0,20], [0,0,0], [0,1,0])
    ##################
    
    T_gripper = np.eye(4)
    T_gripper[:3,3] = np.array([2,2,0])
    robot.GetLink(LINK_NAME).SetTransform(T_gripper)
    T_grasp = np.eye(4)
    xyz_targ = T_grasp[:3,3]
#    mk.create_mesh_box(env, np.array([5+0.01,0,0]), np.array([0.01,10,0.01]))
    mk.create_mesh_box(env, np.array([3.5,3,0]), np.array([1,1,2]), "box1")
    mk.create_mesh_box(env, np.array([3.5,-0.5,0]), np.array([1,1,2]), "box2")
#    mk.create_mesh_box(env, np.array([3.5,1.25,0]), np.array([0.2,0.2,1]), "box3")
#    mk.create_mesh_box(env, np.array([3.5,2,0]), np.array([0.2,0.2,1]), "box4")
    quat_targ = rave.quatFromRotationMatrix(T_grasp[:3,:3])

    request = move_arm_to_grasp(xyz_targ, quat_targ, LINK_NAME, MANIP_NAME)
    
    ##################
    """
    # first optimize ignoring collision costs
    all_costs = request["costs"]
    noncollision_costs = [cost for cost in request["costs"] if "collision" not in cost["type"]]
    
    request["costs"] = noncollision_costs
    
    saver = rave.Robot.RobotStateSaver(robot)
    s = json.dumps(request)
    print "REQUEST:",s
    trajoptpy.SetInteractive(False);
    prob = trajoptpy.ConstructProblem(s, env)
    result = trajoptpy.OptimizeProblem(prob)
    del saver # reverts the robot state

    # add collision cost back again
    request["costs"] = all_costs
    
    # use the resulting trajectory as initialization for the new optimization that includes collision costs
    path_init = result.GetTraj()
    request["init_info"]["type"] = "given_traj"
    request["init_info"]["data"] = [x.tolist() for x in path_init]
    """
    ##################

    s = json.dumps(request)
    print "REQUEST:",s

    trajoptpy.SetInteractive(INTERACTIVE);
    prob = trajoptpy.ConstructProblem(s, env)
    result = trajoptpy.OptimizeProblem(prob)
    
    print "Sum of final costs is ", sum([cost[1] for cost in result.GetCosts()])
    
#    np.save("/home/alex/Desktop/traj.npy", result.GetTraj())
