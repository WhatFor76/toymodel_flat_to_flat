'''
Copyright @WhatFor, contact with 332449502@qq.com
Used in step 3.
This is the main of step 3 for the thickening of a plane model.
to control the max-min rate of the height of shafts, modify the 
value of q, and the target function c cdot x is a random vector.

'''

import pickle
from Unit import *
import math
from scipy.linalg import null_space
from plotclasses import PlotModel
import matplotlib.pyplot as plt
import random
from scipy.optimize import linprog
class GlobalParameter:
    def __init__(self):
        ax=None
        bx=None
        with open('outState1.pickle','rb') as f:
           ax=pickle.load(f)
        with open('outState2.pickle','rb') as f:
           bx=pickle.load(f)
        self.Surfaces1:list[Surface]=ax[0]
        self.Shafts1:list[Edge]=ax[1]
        self.Vertice1:list[Vertex]=ax[2]
        self.Surfaces2:list[Surface]=bx[0]
        self.Shafts2:list[Edge]=bx[1]
        self.Vertice2:list[Vertex]=bx[2]
        self.V=len(self.Vertice1)
        self.E=len(self.Shafts1)   
        self.S=len(self.Surfaces1)
eps=0.000001
GV=GlobalParameter()

def get_Constraint()->np.matrix:
    #constraints of state 1, at 2D plate
    Con=np.zeros(shape=(6*GV.E,GV.S+GV.E+6*GV.S))
    con_addr=0
    def write_to_Con1(a:int,b:int,c:int,Con0:np.matrix)->None:
        for i in range(3):
            Con[con_addr+i,a]=Con0[i,0]
            Con[con_addr+i,b]=Con0[i,1]
            Con[con_addr+i,GV.S+c]=Con0[i,2]
            Con[con_addr+i,GV.S+GV.E+a*3+i]=1
            Con[con_addr+i,GV.S+GV.E+b*3+i]=-1
    #though state 1's constraints only a rank 0 matrix in fact, we can simply ignore it. but for the future construction, we write it to Con.
    #function: x1-Sn1(s1-h)=x2-Sn2(s2-h), V=x+V0, so [Sn1*D1-Sn2*D2,T][s,h,x_state1,x_state2].T=0
    #this can provide a condition that [s,h] will be in the column space of [x] for get the real position can get a legal solution.
    for e in GV.Shafts1:
        if(e.sur2==-1):
            continue
        a=GV.Surfaces1[e.sur1].uppersurface()
        b=GV.Surfaces1[e.sur2].uppersurface()
        x1=a/math.sqrt(a*a.T)
        x2=b/math.sqrt(b*b.T)
        diffm1=np.mat([1,0,-1])
        diffm2=np.mat([0,1,-1])
        B=x1.T*diffm1-x2.T*diffm2
        s1=e.sur1
        s2=e.sur2
        h=e.id
        write_to_Con1(s1,s2,h,B)#write this con to the proper location#
        con_addr+=3
    #combine with state 2.
    def write_to_Con2(a:int,b:int,c:int,Con0:np.matrix)->None:
        for i in range(3):
            Con[con_addr+i,a]=Con0[i,0]
            Con[con_addr+i,b]=Con0[i,1]
            Con[con_addr+i,GV.S+c]=Con0[i,2]
            Con[con_addr+i,GV.S+GV.E+3*GV.S+a*3+i]=1
            Con[con_addr+i,GV.S+GV.E+3*GV.S+b*3+i]=-1
    for e in GV.Shafts2:
        if(e.sur2==-1):
            continue
        a=GV.Surfaces2[e.sur1].uppersurface()
        b=GV.Surfaces2[e.sur2].uppersurface()
        x1=a/math.sqrt(a*a.T)
        x2=b/math.sqrt(b*b.T)
        diffm1=np.mat([1,0,-1])
        diffm2=np.mat([0,1,-1])
        B=x1.T*diffm1-x2.T*diffm2
        s1=e.sur1
        s2=e.sur2
        h=e.id
        write_to_Con2(s1,s2,h,B)
        con_addr+=3
    
    Con=np.mat(Con)
    return Con


#a is the index of surface and b is the index of shaft, shaft is negative to surface: [:,a]-T[:,b]>=0, Xs' parameter is 0 for X can be any.
def get_Diff()->np.matrix:
    T:list[float]=np.zeros(shape=(2*GV.E,GV.S+GV.E+6*GV.S))
    T_addr=0
    q_state:list[int]=np.zeros(shape=2*GV.E)#to mark if this element of q is 0 or positive.
    def write_to_T(a:int,b:int,type:int)->None:
        T[T_addr,a]=1
        T[T_addr,b]=-1
        if(type==1):
            q_state[T_addr]=0
        if(type==2):
            q_state[T_addr]=1
        return 
    for e in GV.Shafts1:
        if(e.sur2==-1):
            #force the side shafts to zero thickness, the same as type1 shafts
            write_to_T(e.sur1,GV.S+e.id,1)
            T_addr+=1
        else:
            write_to_T(e.sur1,GV.S+e.id,e.type)
            T_addr+=1
            write_to_T(e.sur2,GV.S+e.id,e.type)
            T_addr+=1
    #set surface0 's height the baseline 0.
    T[T_addr,0]=1
    T_addr+=1
    T=np.mat(T)
    return T,q_state
#only state 2, for state 1, GV. imported _1 items.
def get_real_position(S:list[Surface],x:np.matrix)->np.matrix:
    T:list[float]=np.zeros(shape=(GV.E*3,GV.S*3))
    T_addr=0
    b:list[float]=np.zeros(shape=(GV.E*3,1))
    #TV=b
    def write_to_T(a:int,b:int):
        for i in range(3):
            T[T_addr+i,a*3+i]=1
            T[T_addr+i,b*3+i]=-1
    def write_to_b(e:Edge):
        s1=x[e.sur1,0]
        s2=x[e.sur2,0]
        h=x[GV.S+e.id,0]
        sn1=S[e.sur1].uppersurface()
        sn2=S[e.sur2].uppersurface()
        _b:np.matrix=(sn1*(s1-h)-sn2*(s2-h)).T
        for i in range(3):
            b[T_addr+i,0]=_b[i,0]
    for e in GV.Shafts1:
        if(e.sur2==-1):
            continue
        write_to_T(e.sur1,e.sur2)
        write_to_b(e)
        T_addr+=3
    # eliminate 3 DOF from forcing surface1 have a offset 0
    for i in range(3):
        T[T_addr+i,i]=1
        b[T_addr+i,0]=0
    T_addr+=3
    T:np.matrix=np.mat(T)
    T_M=np.linalg.inv(T.T*T)
    Vertex_position=T_M*T.T*b
    return Vertex_position

def get_dir(A:np.matrix):
        
    q=np.zeros(shape=(A.shape[0]*2))
    c=np.zeros(shape=(A.shape[1]))
    random.seed(10086)
    for i in range(A.shape[1]):
        c[i]=-random.random()
    '''the control rate of max and min length of the shafts'''
    for i in range(A.shape[0]):
            q[i]=0.3
            q[i+A.shape[0]]=-0.09
       
    A1=A.A
    A2=(-A).A
    A=np.append(A1,A2,axis=0)
    d=linprog(c=c,A_ub=A,b_ub=q,bounds=(-10,10))
    print(d.success)
    d=np.mat(d.x).T
    return d

def flaten_surfaces(x:np.matrix,h:float=-1)->np.matrix:
    if(h==-1):
        return x
    for s in GV.Surfaces1:
        if(abs(s.e[0].L-s.e[1].L)<eps):
            for e in s.e:
                if(h<=x[GV.S+e.id,0]):
                    print("error!")
                    return x 
            x[s.id,0]=h
    print("flaten_surface Succeed!")
    return x
def main():
    
    Con=get_Constraint()
    Diff,q_state=get_Diff()
    Null=np.mat(null_space(Con,rcond=1e-6))
    ACR=Diff*Null #A have two parts, ACRx=q this x=ker(R)=Nc, CNc=q>0
    #null space of A will be the num of the side edge plus one, which is the undefined baseline height.
    #so we believe that the code have no error. and we just ignore the null, for baseline is 0 and all side edge have no thickness.
    
    
# the method to get q is temporary, will be replaced by a function or even a simulator further.
    
    C=np.mat(np.zeros(shape=ACR.shape))
    indexC=0
    R=np.mat(np.zeros(shape=ACR.shape))
    indexR=0
    for i in range(q_state.shape[0]):
        if(q_state[i]==1):
            #q[i] is any positive num
            C[indexC,:]=ACR[i,:]
            indexC+=1
        if(q_state[i]==0):
            R[indexR,:]=ACR[i,:]
            indexR+=1
    N=np.mat(null_space(R,rcond=1e-6))
    C=C*N   
    C=C[0:indexC]
    d=get_dir(C)
    x=Null*N*d
    
    x=flaten_surfaces(x[0:GV.S+GV.E,0],h=0.08)# this is not so necessary, but to get a flat of all equilateral triangle plates in two states, it's necessary. 
    vertex_offset1=get_real_position(GV.Surfaces1,x=x[0:GV.S+GV.E,0])
    vertex_offset2=get_real_position(GV.Surfaces2,x=x[0:GV.S+GV.E,0])



    ploter=PlotModel()
    plt.ion()
    ploter.plot2(S=GV.Surfaces2,v0=vertex_offset2)  
    with open("outHingeShift.pickle","wb") as f:
        x0=np.array(x[0:GV.S+GV.E])
        x0=x0.flatten()
        pickle.dump(x0,f)
    with open("outSurfacesOffset1.pickle","wb") as f:
        x0=np.array(vertex_offset1[0:3*GV.S])
        x0=x0.flatten()
        pickle.dump(x0,f)
    with open("outSurfacesOffset2.pickle","wb") as f:
        x0=np.array(vertex_offset2[0:3*GV.S])
        x0=x0.flatten()
        pickle.dump(x0,f)
    plt.ioff()
    ploter.plotflat(S=GV.Surfaces1,x=x[0:GV.S+GV.E,0])
    plt.show()
if __name__=="__main__":
    main()



