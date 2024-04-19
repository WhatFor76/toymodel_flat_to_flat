
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
        with open('outSurfacesOffset1.pickle','rb') as f:
            self.SurfacesOff1:list[float]=pickle.load(f)
        with open('outSurfacesOffset2.pickle','rb') as f:
            self.SurfacesOff2:list[float]=pickle.load(f)
        with open('outHingeShift.pickle','rb') as f:
            self.shift:list[float]=pickle.load(f)
        self.V=len(self.Vertice1)
        self.E=len(self.Shafts1)   
        self.S=len(self.Surfaces1)
eps=0.0001
GV=GlobalParameter()
#to build the surfaces and shafts using hinges in 3d space, 
#additionally, type1 surfaces are specialized dealt with, share the vertice in surfaces with its nearby surface with no shafts.
#and it leads to two same constraint, but has no influence to the final result, so we allow these useless hinges' existance.
def buildTruss(Surfaces:list[Surface],Shafts:list[Edge])->tuple((np.matrix,list)):
    A=np.zeros((3*GV.S+12*GV.E+9,3*(3*GV.S+2*GV.E)))    
    indextag=np.zeros(A.shape[1],dtype=int)
    for i in range(indextag.shape[0]):
        indextag[i]=-1
    for s in Surfaces:
        col_id=[[0,0,0],[0,0,0],[0,0,0]]
        for i in range(3):
            for j in range(3):
                col_id[i][j]=3*s.v_shifted[i].id+j
                indextag[3*(s.id*3+i)+j]=3*(s.id*3+i)+j
                if(s.v_shifted[i].type==3):
                    A[0,3*(s.id*3+i)+j]=-999#(s.id*3+i) is the ref to the common independent v_shift's id, not used in the special ones, so it is marked to zero
                    col_id[i][j]=9*GV.S+3*s.v_shifted[i].id+j    
                    indextag[9*GV.S+3*s.v_shifted[i].id+j]=3*(s.id*3+i)+j 
        for j in range(3):
            #surfaces    
            A[s.id*3+0,col_id[0][j]]=s.v_shifted[0].x[j]-s.v_shifted[1].x[j]
            A[s.id*3+0,col_id[1][j]]=s.v_shifted[1].x[j]-s.v_shifted[0].x[j]
            A[s.id*3+1,col_id[1][j]]=s.v_shifted[1].x[j]-s.v_shifted[2].x[j]
            A[s.id*3+1,col_id[2][j]]=s.v_shifted[2].x[j]-s.v_shifted[1].x[j]            
            A[s.id*3+2,col_id[2][j]]=s.v_shifted[2].x[j]-s.v_shifted[0].x[j]
            A[s.id*3+2,col_id[0][j]]=s.v_shifted[0].x[j]-s.v_shifted[2].x[j]
    #hinge of shafts
    for e in Shafts:
        #for sur1!=-1 is predefined, so it is no necessary to check again.
        if(e.sur2==-1):
            for j in range(3):
                A[0,9*GV.S+3*e.u_shifted.id+j]=-999
                A[0,9*GV.S+3*e.v_shifted.id+j]=-999
        if(e.type==2 and e.sur2!=-1):
            s=Surfaces[e.sur1]
            for i in range(3):
                sur_col_id=[0,0,0]
                for j in range(3):
                    sur_col_id[j]=3*s.v_shifted[i].id+j
                    if(s.v_shifted[i].type==3):
                        sur_col_id[j]=9*GV.S+3*s.v_shifted[i].id+j                    
                for j in range(3):
                    A[GV.S*3+e.u_shifted.id*3+i,            sur_col_id[j]]=s.v_shifted[i].x[j]-e.u_shifted.x[j]
                    A[GV.S*3+e.u_shifted.id*3+i,9*GV.S+3*e.u_shifted.id+j]=e.u_shifted.x[j]-s.v_shifted[i].x[j]
                    A[GV.S*3+e.v_shifted.id*3+i,            sur_col_id[j]]=s.v_shifted[i].x[j]-e.v_shifted.x[j]
                    A[GV.S*3+e.v_shifted.id*3+i,9*GV.S+3*e.v_shifted.id+j]=e.v_shifted.x[j]-s.v_shifted[i].x[j]
            s=Surfaces[e.sur2]
            for i in range(3):
                sur_col_id=[0,0,0]
                for j in range(3):
                    sur_col_id[j]=3*s.v_shifted[i].id+j
                    if(s.v_shifted[i].type==3):
                        sur_col_id[j]=9*GV.S+3*s.v_shifted[i].id+j    
                for j in range(3):
                    A[GV.S*3+6*GV.E+e.u_shifted.id*3+i,            sur_col_id[j]]=s.v_shifted[i].x[j]-e.u_shifted.x[j]
                    A[GV.S*3+6*GV.E+e.u_shifted.id*3+i,9*GV.S+3*e.u_shifted.id+j]=e.u_shifted.x[j]-s.v_shifted[i].x[j]
                    A[GV.S*3+6*GV.E+e.v_shifted.id*3+i,            sur_col_id[j]]=s.v_shifted[i].x[j]-e.v_shifted.x[j]
                    A[GV.S*3+6*GV.E+e.v_shifted.id*3+i,9*GV.S+3*e.v_shifted.id+j]=e.v_shifted.x[j]-s.v_shifted[i].x[j]
    #then fix a equilant triange surface to ground. these surface have no type3 vertex, so col_id neednot to be specially judged 
    #the fix constraint is dx=0, the referenced matrix is a unit matrix.
    for s in Surfaces:
        if(abs(s.e[0].L-s.e[1].L)<eps):
            for i in range(3):
                for j in range(3):
                    A[3*GV.S+12*GV.E+3*i+j,3*s.v_shifted[i].id+j]=1
            break


    indextag=indextag[[A[0,i]!=-999 for i in range(A.shape[1])]]
    A=A[[not np.all(abs(A[i,:])<eps) for i in range(A.shape[0])],:]
    A=A[:,[A[0,i]!=-999 for i in range(A.shape[1])]]

    return (A,indextag)
def main():
    v0=GV.SurfacesOff1
    for s in GV.Surfaces1:
        s.set_v_shifted(v0[s.id*3:(s.id+1)*3])
    for e in GV.Shafts1:
        e.set_e_shifted(GV.Surfaces1,GV.shift[GV.Surfaces1[e.sur1].id]-GV.shift[GV.S+e.id],GV.shift[GV.Surfaces1[e.sur2].id]-GV.shift[GV.S+e.id])
    A,idref=buildTruss(GV.Surfaces1,GV.Shafts1)
    N=np.mat(null_space(A,rcond=eps))
    nullrank=N.shape[1]
    print("State1DOF:",nullrank)
    print(A.shape)
    #for state 1 needs a plot of modium, so idref is necessary to mark which line is not deleted.
    #and also for a consideration of visual, only does we use the surface to plot.
    #move_vector use the order of a vertex in the surface, so type9 vertex's col in matrix A ought to be specially judged. 
    move_vector=np.zeros(shape=(N.shape[1],GV.S*9))
    for k in range(N.shape[1]):
        for i in range(N.shape[0]):
            if(idref[i]!=-1):
                move_vector[k,idref[i]]=-N[i,k]
    ploter=PlotModel()
    print(move_vector.T)
    plt.ion()
    #0: the first main dof direction, and 1:, the second... et al.
    mv=move_vector[3,:].flatten()
    print(mv)
    ploter.plotmovement(GV.Surfaces1,mv)  
    plt.ioff()
    plt.show()
    v0=GV.SurfacesOff2
    for s in GV.Surfaces2:
        s.set_v_shifted(v0[s.id*3:(s.id+1)*3])

    for e in GV.Shafts2:
        e.set_e_shifted(GV.Surfaces2,GV.shift[GV.Surfaces2[e.sur1].id]-GV.shift[GV.S+e.id])
    A,idref=buildTruss(GV.Surfaces2,GV.Shafts2)
    N=np.mat(null_space(A,rcond=eps))
    nullrank=N.shape[1]
    print("State2DOF:",nullrank)
    print(A.shape)
if __name__=="__main__":
    main()