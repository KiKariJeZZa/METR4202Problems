import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np
from numpy import cos
from numpy import sin
import sympy as sym
import math

def Mellipse(L1, L2, thetalist):
    # Initiate the values and parameters needed
    theta1 = thetalist[0]
    theta = thetalist[0] + thetalist[1]
    Jv = np.array([[-L1*sin(theta1) - L2*sin(theta), -L2*sin(theta)],
                [L1*cos(theta1) + L2*cos(theta), L2*cos(theta)]])
    Jv_t = np.transpose(Jv)
    A = np.matmul(Jv,Jv_t)
    eigenvalues = np.linalg.eig(A)

    # Calculate the manipulability features
    eigen_max = np.amax(eigenvalues[0])
    eigen_min = np.amin(eigenvalues[0])
    mu_1 = np.sqrt(eigen_max) / np.sqrt(eigen_min)
    mu_2 = eigen_max / eigen_min
    mu_3 = np.sqrt(np.linalg.det(A))
    mu_test = np.sqrt(eigen_max*eigen_min)
    print('Mu_1:',mu_1)
    print('Mu_2:',mu_2)
    print('Mu_3:',mu_3)
    print('Mu_test', mu_test)

    # Get the features of the ellipse
    width = np.sqrt(eigenvalues[0][0])
    height = np.sqrt(eigenvalues[0][1])
    angle_rad = math.atan(eigenvalues[1][0][1]/eigenvalues[1][0][0])
    angle_deg = angle_rad*180/np.pi

    # Plot the arms of the robot
    link1_posx = [0, L1*cos(theta1)]
    link1_posy = [0, L1*sin(theta1)]
    link2_posx = [L1*cos(theta1), L1*cos(theta1)+L2*cos(theta)]
    link2_posy = [L1*sin(theta1), L1*sin(theta1)+L2*sin(theta)]
    base = plt.Circle((0,0), 0.05, color='black')
    ellipse = Ellipse(xy=(L1*cos(theta1)+L2*cos(theta), L1*sin(theta1)+L2*sin(theta)), 
                        width=width, height=height, angle = -angle_deg, 
                        edgecolor='r', fc='None', lw=2)
    fig, ax = plt.subplots()
    line1, = ax.plot(link1_posx, link1_posy, linewidth=5)
    line2, = ax.plot(link2_posx, link2_posy, linewidth=5)
    ax.legend([line1,line2],['Link 1', 'Link 2'])
    plt.xlim(-2, 4)
    plt.ylim(-2, 4)
    ax.add_patch(base)
    ax.add_patch(ellipse)
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.title('2R Planar Robot configuration')
    plt.show()

    
# Change arm lengths and angles
L1 = 1
L2 = 1
theta1 = 60
theta2 = 60
thetalist = (theta1*np.pi/180, theta2*np.pi/180)

print(Mellipse(L1, L2, thetalist))
    
    
   