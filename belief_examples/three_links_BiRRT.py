import trajoptpy
import trajoptpy.make_kinbodies as mk
from openravepy import *
import numpy as np
import time

env = Environment()
env.Load('../data/three_links_BiRRT.xml')
env.SetViewer('qtcoin')
robot = env.GetRobots()[0]

ikmodel=databases.inversekinematics.InverseKinematicsModel(robot,iktype=IkParameterization.Type.TranslationXY2D)
if not ikmodel.load():
    ikmodel.autogenerate()

print('IK model autogenerated!')

T_gripper = robot.GetLink("Finger").GetTransform()
T_grasp = T_gripper.copy()
T_grasp[:3,3]  += np.array([-0.5,0.2,0])
T_grasp = T_grasp.dot(matrixFromAxisAngle([0,0,1],np.pi/2))
xyz_targ = T_grasp[:3,3]
quat_targ = quatFromRotationMatrix(T_grasp[:3,:3])

mk.create_mesh_box(env, np.array([0.2,0.25,0]), np.array([0.05,0.05,0.3]), "box1")

target=T_grasp[:2,3]
print(target)

h=env.plot3(np.array([target[0],target[1],ikmodel.manip.GetTransform()[2,3]]),10.0)

robot.SetDOFLimits([-10, -10, -10],[10,10,10]) # to restrict space of IK solutions
robot.SetActiveDOFs(robot.SetActiveManipulator('arm').GetArmIndices()) 

ikp=IkParameterization(target,IkParameterization.Type.TranslationXY2D) 
manip=interfaces.BaseManipulation(robot)

#traj = manip.MoveToHandPosition(ikparam=ikp,seedik=10,execute=False,outputtrajobj=True)
#print(traj)
	
iksolns=ikmodel.manip.FindIKSolutions(ikp, IkFilterOptions.CheckEnvCollisions)
print(iksolns)
print('Num IK solutions: ' + str(len(iksolns)))
for i in np.random.permutation(len(iksolns))[0:min(25,len(iksolns))]:
	#robot.SetDOFValues(iksolns[i],ikmodel.manip.GetArmIndices())
	robot.SetDOFValues([0, 0, 0])
	traj = manip.MoveActiveJoints(goal=iksolns[i],execute=False,outputtrajobj=True)
	print('Num waypoints: ' + str(traj.GetNumWaypoints()))
	for j in range(traj.GetNumWaypoints()):
    		# get the waypoint values, this holds velocites, time stamps, etc
		data=traj.GetWaypoint(j)
		# extract the robot joint values only
		dofvalues=traj.GetConfigurationSpecification().ExtractJointValues(data,robot,robot.GetActiveDOFIndices())
		raveLogInfo('waypoint %d is %s'%(j,dofvalues))
		robot.SetDOFValues(dofvalues)
		print(robot.GetActiveManipulator().GetEndEffectorTransform())
		raw_input('Hit ENTER to continue.')	

#robot.GetController().SetPath(traj)

raw_input('Hit ENTER to continue.')
env.Destroy()

