import math
import numpy as np
from scipy.sparse.linalg import expm

# Defines the forward kinematics from the Space Frame
def FKinSpace(M, Slist, thetalist):
        # Convert into arrays to do multiplication with matrices
        S1 = np.array(Slist[0])
        S2 = np.array(Slist[1])
        S3 = np.array(Slist[2])
        S4 = np.array(Slist[3])
        S5 = np.array(Slist[4])
        S6 = np.array(Slist[5])
        # Raise the matrix with base exponential
        T1 = expm(S1*thetalist[0])
        T2 = expm(S2*thetalist[1])
        T3 = expm(S3*thetalist[2])
        T4 = expm(S4*thetalist[3])
        T5 = expm(S5*thetalist[4])
        T6 = expm(S6*thetalist[5])
        #Return the final transformation
        T11 = np.dot(T1,T2)
        T22 = np.dot(T11,T3)
        T33 = np.dot(T22,T4)
        T44 = np.dot(T33,T5)
        T55 = np.dot(T44,T6)
        return np.dot(T55, M)
        
        

# Defines the forward kinematics from the Body Frame
def FKinBody(M, Blist, thetalist):
        # Convert into arrays to do multiplication with matrices
        B1 = np.array(Blist[0])
        B2 = np.array(Blist[1])
        B3 = np.array(Blist[2])
        B4 = np.array(Blist[3])
        B5 = np.array(Blist[4])
        B6 = np.array(Blist[5])
        # Raise the matrix with base exponential
        T1 = expm(B1*thetalist[0])
        T2 = expm(B2*thetalist[1])
        T3 = expm(B3*thetalist[2])
        T4 = expm(B4*thetalist[3])
        T5 = expm(B5*thetalist[4])
        T6 = expm(B6*thetalist[5])
        # Calculate and return the final transformation
        T11 = np.dot(T1,T2)
        T22 = np.dot(T11,T3)
        T33 = np.dot(T22,T4)
        T44 = np.dot(T33,T5)
        T55 = np.dot(T44,T6)
        return np.dot(T55, M)

# Make the matrices
M = [[1, 0, 0, 0],
    [0, 1, 0, 1066.9],
    [0, 0, 1, 118.5],
    [0, 0, 0, 1]]

Slist6 = [[0, 0, 1, 0],
        [0, 0, 0, 0],
        [-1, 0, 0, 0],
        [0, 0, 0, 0]]
Slist5 = [[0, -1, 0, 275.5],
        [1, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0]]
Slist4 = [[0, -1, 0, 685.5],
        [1, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0]]
Slist3 = [[0, 0, 1, 9.8],
        [0, 0, 0, 0],
        [-1, 0, 0, 0],
        [0, 0, 0, 0]]
Slist2 = [[0, -0.866, 0.5, 778.07],
        [0.866, 0, 0, 0],
        [-0.5, 0, 0, 0],
        [0, 0, 0, 0]]
Slist1 = [[0, 0, 1, -118.54],
        [0, 0, 0, 0],
        [-1, 0, 0, 0],
        [0, 0, 0, 0]]   

BList1 = [[0, 0, 1, 0],
        [0, 0, 0, 0],
        [-1, 0, 0, 0],
        [0, 0, 0, 0]]
BList2 = [[0, -0.866, 0.5, 0],
        [0.866, 0, 0, 0],
        [-0.5, 0, 0, 0],
        [0, 0, 0, 0]]
BList3 = [[0, 0, 1, 108.74],
        [0, 0, 0, 0],
        [-1, 0, 0, 0],
        [0, 0, 0, 0]]
BList4 = [[0, -1, 0, -381.5],
        [1, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0]]
BList5 = [[0, -1, 0, -791.5],
        [1, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0]]
BList6 = [[0, 0, 1, 98.94],
        [0, 0, 0, 0],
        [-1, 0, 0, 0],
        [0, 0, 0, 0]]

newB = [BList1, BList2, BList3, BList4, BList5, BList6]
newS = [Slist1, Slist2, Slist3, Slist4, Slist5, Slist6]
thetalist = [math.pi/8, math.pi/4, -math.pi/4, math.pi, math.pi/2, -math.pi/4]

if __name__ == "__main__":
    print(FKinSpace(M,newS,thetalist))    
    print(FKinBody(M,newB,thetalist))
