import pickle
from DHMatrix import Rotation,Translation
from Unit import *
import math
from scipy.linalg import null_space
from scipy.linalg import svd,pinv
from plotclasses import PlotModel
import matplotlib.pyplot as plt
import random
from scipy.optimize import linprog
class GlobalParameter:
    def __init__(self):
        with open('outSurfaceState1.pickle','rb') as f:
            self.Surfaces1:list[Surface]=pickle.load(f)
        with open('outShaftState1.pickle','rb') as f:
            self.Shafts1:list[Edge]=pickle.load(f)
        with open('outVertexState1.pickle','rb') as f:
            self.Vertice1:list[Vertex]=pickle.load(f)
        with open('outSurfaceState2.pickle','rb') as f:
            self.Surfaces2:list[Surface]=pickle.load(f)
        with open('outShaftState2.pickle','rb') as f:
            self.Shafts2:list[Edge]=pickle.load(f)
        with open('outVertexState2.pickle','rb') as f:
            self.Vertice2:list[Vertex]=pickle.load(f)
        self.V=len(self.Vertice1)
        self.E=len(self.Shafts1)   
        self.M=len(self.Surfaces1)

GV=GlobalParameter()

def get_Constraint()->np.matrix:
    #constraints of state 1, at 2D plate
    Con=np.zeros(shape=(6*GV.V,GV.M+GV.E))
    con_addr=0
    def write_to_Con(a:int,b:int,c:int,Con0:np.matrix)->None:
        for i in range(3):
            Con[con_addr+i,a]+=Con0[i,0]
            Con[con_addr+i,b]+=Con0[i,1]
            Con[con_addr+i,GV.M+c]+=Con0[i,2]
    #though state 1's constraints only a rank 0 matrix in fact, we can simply ignore it. but for the future construction, we write it to Con.
    for v in GV.Vertice1:
        if(v.type==-1):
            continue
        T:np.matrix=np.mat([[1,0,0],[0,1,0],[0,0,1]])

        #(H0=H6)+T5H5+T5T4H4+T5T4T3H3+T5T4T3T2H2+T5T4T3T2T1H1
        for i in range(v.adj_num-1,-1,-1):
            a=v.adj_surface[i].uppersurface()
            b=v.adj_surface[(i+1)%v.adj_num].uppersurface()
            x1=a/math.sqrt(a*a.T)
            x2=b/math.sqrt(b*b.T)
            diffm1=np.mat([1,0,-1])
            diffm2=np.mat([0,1,-1])
            B=-x1.T*diffm1+x2.T*diffm2
            Con0=T*B
            s1=v.adj_surface[i].id
            s2=v.adj_surface[(i+1)%v.adj_num].id
            h=v.adj_edge[(i+1)%v.adj_num].id
            write_to_Con(s1,s2,h,Con0)#write this con to the proper location#
            #get rotate matrix for state e0 to state e1 in a surface
            e0=v.adj_edge[i].get_vector()
            e1=v.adj_edge[(i+1)%v.adj_num].get_vector()
            #though Rotation has unitised, for rubust and readability, e0,e1 we also make a unitisation.
            e0=e0/math.sqrt(e0*e0.T)
            e1=e1/math.sqrt(e1*e1.T)
            R=Rotation(e0,e1)
            T=T*R
        con_addr+=3
    #combine with state 2.
    for v in GV.Vertice2:
        if(v.type==-1):
            continue
        T:np.matrix=np.mat([[1,0,0],[0,1,0],[0,0,1]])
        #(H0=H6)+T5H5+T5T4H4+T5T4T3H3+T5T4T3T2H2+T5T4T3T2T1H1
        for i in range(v.adj_num-1,-1,-1):
            a=v.adj_surface[i].uppersurface()
            b=v.adj_surface[(i+1)%v.adj_num].uppersurface()
            x1=a/math.sqrt(a*a.T)
            x2=b/math.sqrt(b*b.T)
            diffm1=np.mat([1,0,-1])
            diffm2=np.mat([0,1,-1])
            B=-x1.T*diffm1+x2.T*diffm2
            Con0=T*B
            s1=v.adj_surface[i].id
            s2=v.adj_surface[(i+1)%v.adj_num].id
            h=v.adj_edge[(i+1)%v.adj_num].id
            write_to_Con(s1,s2,h,Con0)#write this con to the proper location#
            #get rotate matrix for state e0 to state e1 in a surface
            e0=v.adj_edge[i].get_vector(v)
            e1=v.adj_edge[(i+1)%v.adj_num].get_vector(v)
            #though Rotation has unitised, for rubust and readability, e0,e1 we also make a unitisation.
            e0=e0/math.sqrt(e0*e0.T)
            e1=e1/math.sqrt(e1*e1.T)
            R=Rotation(e0,e1)
            T=T*R
        con_addr+=3
    Con=np.mat(Con)
    return Con



#then SVD to get
#here have 3 types of shaft


#a is the index of surface and b is the index of shaft, shaft is negative to surface: [:,a]-T[:,b]>=0
def get_Diff()->np.matrix:
    T:list[float]=np.zeros(shape=(2*GV.E,GV.M+GV.E))
    T_addr=0
    q_state:list[int]=np.zeros(shape=2*GV.E)#to mark if this element of q is 0 or positive.
    def write_to_T(a:int,b:int,type:int)->None:
        T[T_addr,a]=1
        T[T_addr,b]=-1
        if(type==2):
            q_state[T_addr]=1
        return 
    for e in GV.Shafts1:
        if(e.sur2==-1):
            #force the side shafts to zero thickness, the same as type1 shafts
            write_to_T(e.sur1,GV.M+e.id,1)
            T_addr+=1
        else:
            write_to_T(e.sur1,GV.M+e.id,e.type)
            T_addr+=1
            write_to_T(e.sur2,GV.M+e.id,e.type)
            T_addr+=1
    #set surface0 's height the baseline 0.
    T[T_addr,0]=1
    T_addr+=1
    T=np.mat(T)
    return T,q_state
#only state 2, for state 1, GV. imported _1 items.
def get_real_position(S:list[Surface],x:np.matrix)->np.matrix:
    T:list[float]=np.zeros(shape=(GV.E*3,GV.M*3))
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
        h=x[GV.M+e.id,0]
        sn1=S[e.sur1].uppersurface()
        sn2=S[e.sur2].uppersurface()
        _b:np.matrix=(sn1*(s1-h)-sn2*(s2-h))
        _b=_b.T
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
    '''
    U,sigma,VT=svd(T)    
    b:np.matrix=np.mat(b)
    U:np.matrix=np.mat(U)
    VT:np.matrix=np.mat(VT)
    #null rank of T is 3, for the v+v0, v0 is any.
    Sigma_inv:np.matrix=np.mat(np.zeros(shape=(VT.shape[0],U.shape[1])))
    for i in range(sigma.shape[0]):
        if(abs(sigma[i])>1e-8):
            Sigma_inv[i,i]=1/sigma[i]
    T_inv=VT.T*Sigma_inv*U.T
    Vertex_position=T_inv*b
    '''
    T_M=np.linalg.inv(T.T*T)
    Vertex_position=T_M*T.T*b
    print(T.shape)
    #linprog instead.
    for i in range(GV.M):
        print(i,Vertex_position[i*3:i*3+3,0])
    return Vertex_position
def main():
    
    Con=get_Constraint()
    Diff,q_state=get_Diff()
    Null=np.mat(null_space(Con,rcond=1e-6))
    A=Diff*Null
    #null space of A will be the num of the side edge plus one, which is the undefined baseline height.
    #so we believe that the code have no error. and we just ignore the null, for baseline is 0 and all side edge have no thickness.
    
    (U,sigma,VT)=svd(A)
    U:np.matrix=np.mat(U)
    VT:np.matrix=np.mat(VT)
    #to get the pesudo-inverse matrix (Sigma-1)^T
    Sigma_inv=np.mat(np.zeros(shape=(VT.shape[0],U.shape[1])))
    for i in range(sigma.shape[0]):
        if(abs(sigma[i])>1e-8):
            Sigma_inv[i,i]=1/sigma[i]
    A_inv=VT.T*Sigma_inv*U.T
    
    
# the method to get q is temporary, will be replaced by a function or even a simulator further.
    
    q=np.zeros(shape=(q_state.shape[0]*2))
    c=np.zeros(shape=(A.shape[1]))
    for i in range(A.shape[1]):
        c[i]=-random.random()
    for i in range(q_state.shape[0]):
        if(q_state[i]==1):
            #q[i] is any positive num
            q[i]=0.1
       
    A1=A.A
    A2=(-A).A
    A=np.append(A1,A2,axis=0)
    d=linprog(c=c,A_ub=A,b_ub=q,bounds=(-1,1),)
    d=np.mat(d.x).T
    
    '''
    q=np.zeros(shape=(q_state.shape[0]))
    for i in range(q_state.shape[0]):
        if(q_state[i]==1):
            #q[i] is any positive num
            q[i]=0.1
    q:np.matrix=np.mat(q).T
    d=A_inv*q
    print(A.shape)
    '''
    x=Null*d
    for e in GV.Shafts1:
        if(e.type==1 and e.sur2!=-1):
            print(e.id,x[e.sur1],x[e.sur2],x[GV.M+e.id])
            #print(e.id,GV.Surfaces2[e.sur1].uppersurface(),GV.Surfaces2[e.sur2].uppersurface())
    #to return to the real space, which is the height of surfaces and shafts.
    vertex_position=get_real_position(GV.Surfaces1,x)
    ploter=PlotModel()
    plt.ion()
    ploter.plot2(S=GV.Surfaces1,v0=vertex_position)
    plt.ioff()
    plt.show()
if __name__=="__main__":
    main()



