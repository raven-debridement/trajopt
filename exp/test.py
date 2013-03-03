import openravepy
import numpy, time
env = openravepy.Environment()
env.SetViewer('qtcoin')
env.Load('test.xml')
robot = env.GetRobots()[0]

#joint_start = [-1.832, -0.332, -1.011, -1.437, -1.1  , -2.106,  3.074]
#robot.SetDOFValues(joint_start, robot.GetManipulator('leftarm').GetArmIndices())

import pdb; pdb.set_trace()
