import numpy as np
from math import cos,sin,sqrt,acos
from numpy import matrix,mat
#rotate counterclockwise at z+ axis, and then y+
def Rotation_angle(phi:float,theta:float)->matrix:
    A=mat([[cos(phi),-sin(phi),0,0],[sin(phi),cos(phi),0,0],[0,0,1,0],[0,0,0,1]])#yaw angle
    B=mat([[cos(theta),0,sin(theta),0],[0,1,0,0],[-sin(theta),0,cos(theta),0],[0,0,0,1]])#pitch angle
    C=A*B
    return C
def Translation(x:float,y:float,z:float)->matrix:
    A=np.mat([[1,0,0,x],[0,1,0,y],[0,0,1,z],[0,0,0,1]])
    return A
#input a line vector, not column for convenient. but calculate at column vector firstly
#this formulate have a name.
def Rotation(u:matrix,v:matrix)->matrix:
    p:matrix=u.T
    q:matrix=v.T
    #unitised
    normp=sqrt(p.T*p)
    normq=sqrt(q.T*q)
    p=p/normp
    q=q/normq
    #p cross q 
    p_ref=mat([[0,-p[2,0],p[1,0]],[p[2,0],0,-p[0,0]],[-p[1,0],p[0,0],0]])
    n:matrix=p_ref*q
    normn=sqrt(n.T*n)
    n=n/normn
    n_ref=mat([[0,-n[2,0],n[1,0]],[n[2,0],0,-n[0,0]],[-n[1,0],n[0,0],0]])
    costheta=(p.T*q)[0,0]
    I=mat([[1,0,0],[0,1,0],[0,0,1]])
    #the formulate itself
    R:matrix=costheta*I+(1-costheta)*n*n.T+sqrt(1-costheta**2)*n_ref
    return R

if __name__ !="__main__":
    pass