from re import I
import numpy as np
from ModernRobotics import FKinSpace, Adjoint, JacobianSpace, TransInv, se3ToVec, MatrixLog6
from scipy.linalg import logm

# Implements the Newton-Raphson method to find all possible configurations of the robot
def IKinSpace(Slist, M, T, thetalist0, eomg, ev):
    thetalist = np.array(thetalist0)
    i = 0
    maxiterations = 20
    Tsb = FKinSpace(M,Slist,thetalist)
    temp = se3ToVec(MatrixLog(np.matmul(TransInv(Tsb),T)))
    Vs = np.matmul(Adjoint(Tsb),temp)
    err = (np.linalg.norm(Vs[0:3]) > eomg) or (np.linalg.norm(Vs[3:6]) > ev)
    while err and i < maxiterations:
        thetalist = thetalist + np.matmul(np.linalg.pinv(JacobianSpace(Slist,thetalist)), Vs)
        i = i + 1
        Tsb = FKinSpace(M, Slist, thetalist)
        #temp1 = se3ToVec(MatrixLog(np.matmul(InverseT(Tsb),T)))
        Vs = np.matmul(Adjoint(Tsb),se3ToVec(MatrixLog(np.matmul(TransInv(Tsb),T))))
        err = (np.linalg.norm(Vs[0:3]) > eomg) or (np.linalg.norm(Vs[3:6]) > ev)
    return (thetalist % (np.pi*2), not err)

# Returns the Log of a matrix
def MatrixLog(T):
    return logm(T)


Slist = np.array([[0, 0,  1,  0, 0,    0],
                [0, 1,  0,  -89.159, 0,    0],
                [0, 1,  0,  -89.159, 0,    425],
                [0, 1,  0,  -89.159, 0,    817.25],
                [0, 0,  -1,  -135.85, 817.25,    0],
                [0, 1, 0, 5.491, 0, 817.25]]).T
M = np.array([[-1, 0,  0, 817.25],
            [ 0, 0,  1, 218.15],    
            [ 0, 1, 0, -5.491],
            [ 0, 0,  0, 1]])
T = np.array([[0, 1,  0,     -0.5],
            [0, 0,  -1,      0.1],
            [-1, 0, 0, 0.1],
            [0, 0,  0,      1]])
thetalist0 = np.array([0.1, 0.1, 0.1, 0.1, 0.1, 0.1])
eomg = 0.001
ev = 0.01

#temp = se3ToVec(MatrixLog(np.matmul(InverseT(T),T)))
#Vs = np.matmul(Adjoint(T),temp)
#SxT = np.dot(np.array(Slist)[:, 0],thetalist[0])
print(IKinSpace(Slist,M,T,thetalist0,eomg,ev))






    


