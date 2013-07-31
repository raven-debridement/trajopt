from openravepy import *
import openravepy.ikfast as ikfast

import sys
print sys.path

env = Environment()
kinbody = env.ReadRobotXMLFile("./baxter/baxter-reduced.dae")
env.Add(kinbody)
solver = ikfast.IKFastSolver(kinbody=kinbody)
chaintree = solver.generateIkSolver(baselink=0, eelink=8, freeindices=[2,3,4,5,6,7,8], solvefn=ikfast.IKFastSolver.solveFullIK_6D)
code = solver.writeIkSolver(chaintree)
open('ik.cpp','w').write(code)
