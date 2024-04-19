import numpy as np
from numpy import matrix
from scipy.linalg import null_space,svd
import random
from Unit import *
import math
from math import acos,sqrt
eps=0.0001
dxrate=1e-6
sample_rate=0.5
YoungModium=500
Kappa=1000000
StressE=100000
HyperKappa=10000
bound=0.005
torrent_len_error=0.001
'''
v in e and v in s tag can be gorbal to accerate.
'''
def fequ(a,b):
    return abs(a-b)<eps
class OTPSimu:
    #this share the same memory with V,E which used to initize the class, just a light copy, not a deep one
    def __init__(self,numE:int,numS:int,S:list[Surface],shafts:list[Edge],gamma:Edge=None) -> None:
        self.numEdge=numE
        self.numSurface=numS
        self.Surfaces=S
        self.shafts=shafts
        self.gamma=gamma
        #loss to plot in figure to evaluate the quality of a simulate, to help to determine when to quit and output.
        self.loss:list[float]=[0]*10000000
        self.ittime=0
        self.v_in_s_tag=np.zeros(3*(3*self.numSurface+2*self.numEdge),dtype=int)
        self.v_in_e_tag=np.zeros(3*(3*self.numSurface+2*self.numEdge),dtype=int)
        self.bar_tag=np.zeros(3*self.numSurface+(8+1)*self.numEdge,dtype=int)
        for i in range(self.v_in_s_tag.shape[0]):
            self.v_in_s_tag[i]=-1
        for i in range(self.v_in_e_tag.shape[0]):
            self.v_in_e_tag[i]=-1
        for i in range(self.bar_tag.shape[0]):
            self.bar_tag[i]=i
        for s in self.Surfaces:
            for i in range(3):
                for j in range(3):
                    self.v_in_s_tag[3*(s.id*3+i)+j]=s.id*3+i
                    if(s.v[i].type==3):
                        self.v_in_s_tag[3*(s.id*3+i)+j]=-999
                        self.v_in_e_tag[3*(s.id*3+i)+j]=-999

        for e in self.shafts:
            #for sur1!=-1 is predefined, so it is no necessary to check again.
            #e.u.id =2e.id+0 and e.v.id =2e.id+1, nothing special.
            if(e.sur2==-1):
                for j in range(3):
                    self.v_in_e_tag[9*self.numSurface+3*e.u.id+j]=-999
                    self.v_in_e_tag[9*self.numSurface+3*e.v.id+j]=-999
                    self.v_in_s_tag[9*self.numSurface+3*e.u.id+j]=-999
                    self.v_in_s_tag[9*self.numSurface+3*e.v.id+j]=-999
            #type1 and side edge have been updated in surface.
            else:
                for j in range(3):
                    self.v_in_e_tag[9*self.numSurface+3*e.u.id+j]=e.id
                    self.v_in_e_tag[9*self.numSurface+3*e.v.id+j]=e.id
        self.v_in_s_tag=self.v_in_s_tag[[self.v_in_s_tag[i]!=-999  for i in range(self.v_in_s_tag.shape[0])]]
        self.v_in_e_tag=self.v_in_e_tag[[self.v_in_e_tag[i]!=-999  for i in range(self.v_in_e_tag.shape[0])]]
        
        for e in self.shafts:
            if( e.sur2==-1 or e.type==1):
                for j in range(4):
                    self.bar_tag[self.numSurface*3+e.u.id*4+j]=-999
                    self.bar_tag[self.numSurface*3+e.v.id*4+j]=-999
                self.bar_tag[self.numSurface*3+8*self.numEdge+e.id]=-999
        self.bar_tag=self.bar_tag[[self.bar_tag[i]!=-999  for i in range(self.bar_tag.shape[0])]]
    def SGD_getdir(self,t:int=0,Gamma_Tar:float=0,lenerr:float=0)->matrix:
            #calculate A
            #u.dx,u.dy,u.dz at line 3*u.id+0..3*u.id+2
            #f0 is the loss of current loss!
            f0=self.loss[self.ittime-1]
            #to save the loss current, it is a legal, and have happended state.

            A=self._getConstraint()
            #A*delta=0
            Null=np.mat(null_space(A,rcond=1e-6))
            nullrank=Null.shape[1]
            #found direction in ker(A) with P, andN a ker
            #finally need a target function and a punishment fucntion
            sampledir=random.sample(range(0,nullrank),round(nullrank*sample_rate))
            omega=np.mat(np.zeros((nullrank,1)))
            for i in sampledir:
                dx=dxrate*Null[:,i]
                self.update(dx)
                f1=self._relax(t=t,Gamma_Tar=Gamma_Tar,lenerr=lenerr)
                omega[i,0]=(f1-f0)/dxrate
                self.update(-dx)
            #return direction as dir
            dx=-Null*omega
            return dx

    #modify need    
    def update(self,dx:matrix)->None:
        Len=self.v_in_s_tag.shape[0]
        for i in range(0,Len,3):
            id=self.v_in_s_tag[i]
            if(id!=-1):
                for j in range(3):
                    self.Surfaces[id//3].v[id%3].x[j]+=dx[i+j,0]
        #type 3 vertex is ignore in surface and update in edge.
        Len=self.v_in_e_tag.shape[0]
        for i in range(0,Len,6):
            id=self.v_in_e_tag[i]
            if(id!=-1):
                for j in range(3):
                    self.shafts[id].u.x[j]+=dx[i+j,0]
                    self.shafts[id].v.x[j]+=dx[i+3+j,0]
    def stepping(self,step_length:float,Gamma_Tar:float)->float:
        #P(x)=a^2
        # x=x+ldx, where dx have just a mod |ldx|=length
        current_error=self._relax(Gamma_Tar=Gamma_Tar)
        self.loss[self.ittime]=current_error
        self.ittime+=1
        dx=self.SGD_getdir(Gamma_Tar=Gamma_Tar)
        normdx:float=np.linalg.norm(dx)
        if(normdx>eps):
        # if dx is not a 0 num
        # update x=x+dx
        # maybe need to fix the error
        # update the primer equation f(x)=a^2, which is quitely identical to the first part of get_dir.
        # WE USE a ANNEALING METHOD to CONTROL if THIS can BE ACCEPTED!!!
            dx0=dx/normdx*step_length

            self.update(dx0)
            #we accept this update
            return current_error
        else:
            print("wronganswer\n")
            return -1
  
    #x+dx dot y+dy=0 means xy+xdy+ydx+dxdy=0, the ignore of dxdy leads to a error in tensselation. 100*Young's*m**2, m is the error from pi/2.
    
    #theta=arccos(b**2-a**2-c**2/ 2ac) and which we want is abs(bata-pi/2-theta), where bata is the angle of surface 
    # and rectangle, c is the width of rectangle and a,b is the length of shifted hinges.  
    # and E=YoungModium beta**2
    #the index tag is not the same as stable comfirm for the edge needs to be calculated
    #the the type 3 vertex will be updated when edge upd.
    def _getConstraint(self)->matrix:
        A=np.zeros((3*self.numSurface+(8+1)*self.numEdge,3*(3*self.numSurface+2*self.numEdge)))    
        for s in self.Surfaces:
            col_id=[[0,0,0],[0,0,0],[0,0,0]]
            for i in range(3):
                for j in range(3):
                    col_id[i][j]=3*s.v[i].id+j
                    if(s.v[i].type==3):
                        A[0,3*(s.id*3+i)+j]=-999#(s.id*3+i) is the ref to the common independent v_shift's id, not used in the special ones, so it is marked to zero
                        col_id[i][j]=9*self.numSurface+3*s.v[i].id+j    
            for j in range(3):
                #surfaces    
                A[s.id*3+0,col_id[0][j]]=s.v[0].x[j]-s.v[1].x[j]
                A[s.id*3+0,col_id[1][j]]=s.v[1].x[j]-s.v[0].x[j]
                A[s.id*3+1,col_id[1][j]]=s.v[1].x[j]-s.v[2].x[j]
                A[s.id*3+1,col_id[2][j]]=s.v[2].x[j]-s.v[1].x[j]            
                A[s.id*3+2,col_id[2][j]]=s.v[2].x[j]-s.v[0].x[j]
                A[s.id*3+2,col_id[0][j]]=s.v[0].x[j]-s.v[2].x[j]
        #hinge of shafts
        for e in self.shafts:
            #for sur1!=-1 is predefined, so it is no necessary to check again.
            #e.u.id =2e.id+0 and e.v.id =2e.id+1, nothing special.
            
            if(e.sur2==-1):
                for j in range(3):
                    A[0,9*self.numSurface+3*e.u.id+j]=-999
                    A[0,9*self.numSurface+3*e.v.id+j]=-999
            if(e.type==2 and e.sur2!=-1):               
                s=self.Surfaces[e.sur1]
                sur_col_idu=[0,0,0]
                sur_col_idv=[0,0,0]
                for j in range(3):
                    sur_col_idu[j]=3*s.v[e.uid1].id+j
                    if(s.v[e.uid1].type==3):
                        sur_col_idu[j]=9*self.numSurface+3*s.v[e.uid1].id+j    
                    sur_col_idv[j]=3*s.v[e.vid1].id+j
                    if(s.v[e.vid1].type==3):
                        sur_col_idv[j]=9*self.numSurface+3*s.v[e.vid1].id+j      
                for j in range(3):               
                    A[self.numSurface*3+e.u.id*4+0,              sur_col_idu[j]]=s.v[e.uid1].x[j]-e.u.x[j]
                    A[self.numSurface*3+e.u.id*4+0,9*self.numSurface+3*e.u.id+j]=e.u.x[j]-s.v[e.uid1].x[j]
                    A[self.numSurface*3+e.u.id*4+1,              sur_col_idv[j]]=s.v[e.vid1].x[j]-e.u.x[j]
                    A[self.numSurface*3+e.u.id*4+1,9*self.numSurface+3*e.u.id+j]=e.u.x[j]-s.v[e.vid1].x[j]
                    A[self.numSurface*3+e.v.id*4+0,              sur_col_idu[j]]=s.v[e.uid1].x[j]-e.v.x[j]
                    A[self.numSurface*3+e.v.id*4+0,9*self.numSurface+3*e.v.id+j]=e.v.x[j]-s.v[e.uid1].x[j]
                    A[self.numSurface*3+e.v.id*4+1,              sur_col_idv[j]]=s.v[e.vid1].x[j]-e.v.x[j]
                    A[self.numSurface*3+e.v.id*4+1,9*self.numSurface+3*e.v.id+j]=e.v.x[j]-s.v[e.vid1].x[j]
                    
                s=self.Surfaces[e.sur2]
                sur_col_idu=[0,0,0]
                sur_col_idv=[0,0,0]
                for j in range(3):
                    sur_col_idu[j]=3*s.v[e.uid2].id+j
                    if(s.v[e.uid2].type==3):
                        sur_col_idu[j]=9*self.numSurface+3*s.v[e.uid2].id+j    
                    sur_col_idv[j]=3*s.v[e.vid2].id+j
                    if(s.v[e.vid2].type==3):
                        sur_col_idv[j]=9*self.numSurface+3*s.v[e.vid2].id+j      
                for j in range(3):               
                    A[self.numSurface*3+e.u.id*4+2,              sur_col_idu[j]]=s.v[e.uid2].x[j]-e.u.x[j]
                    A[self.numSurface*3+e.u.id*4+2,9*self.numSurface+3*e.u.id+j]=e.u.x[j]-s.v[e.uid2].x[j]
                    A[self.numSurface*3+e.u.id*4+3,              sur_col_idv[j]]=s.v[e.vid2].x[j]-e.u.x[j]
                    A[self.numSurface*3+e.u.id*4+3,9*self.numSurface+3*e.u.id+j]=e.u.x[j]-s.v[e.vid2].x[j]
                    A[self.numSurface*3+e.v.id*4+2,              sur_col_idu[j]]=s.v[e.uid2].x[j]-e.v.x[j]
                    A[self.numSurface*3+e.v.id*4+2,9*self.numSurface+3*e.v.id+j]=e.v.x[j]-s.v[e.uid2].x[j]
                    A[self.numSurface*3+e.v.id*4+3,              sur_col_idv[j]]=s.v[e.vid2].x[j]-e.v.x[j]
                    A[self.numSurface*3+e.v.id*4+3,9*self.numSurface+3*e.v.id+j]=e.v.x[j]-s.v[e.vid2].x[j]
                for j in range(3):
                    A[self.numSurface*3+8*self.numEdge+e.id,9*self.numSurface+3*e.u.id+j]=e.u.x[j]-e.v.x[j]
                    A[self.numSurface*3+8*self.numEdge+e.id,9*self.numSurface+3*e.v.id+j]=e.v.x[j]-e.u.x[j]

        A=A[:,[A[0,i]!=-999 for i in range(A.shape[1])]]
        A=A[self.bar_tag,:]
        return A






        
    def _relax(self,t:int=0,Gamma_Tar:float=0,lenerr:float=0)->float:
        if(t==0):
            energy=0
            energy+=self.get_gamma(Gamma_Tar)[0]         
            energy+=self.get_energy()
            return energy
        else:
            energy=0
            energy+=StressE*lenerr
            return energy

    #this is to get the vertical vector of the shaft in the triangle. line vector. relations differ from the plate simu.
    def _get_edge_and_counterpart(self,s:Surface,uid:int,vid:int)->matrix:
        #the edge the same as the shaft
        tid=[x for x in [0,1,2] if x not in [uid,vid]]
        tid=tid[0]
        v1:matrix=Edge(u=s.v[uid],v=s.v[tid]).get_vector()
        v2:matrix=Edge(u=s.v[uid],v=s.v[vid]).get_vector()
        v=v1-v1*v2.T/(v2*v2.T)*v2
        v=v/sqrt(v*v.T)
        return v
    
    #if return anglesig >=0 a valley, and <=0 a mount.
    def _get_dihedral(self,e:Edge)->tuple[float,float]:
        a:Surface=self.Surfaces[e.sur1]
        b:Surface=self.Surfaces[e.sur2]
        v1=self._get_edge_and_counterpart(a,e.uid1,e.vid1)
        v2=self._get_edge_and_counterpart(b,e.uid2,e.vid2)
        p=a.uppersurface()
        normp=np.linalg.norm(p)
        q=b.uppersurface()
        normq=np.linalg.norm(q)
        v3=p/normp+q/normq
        costheta=v1*v2.T
        #to fix the float error such as 1.0+1e-16
        if(costheta>1):
            costheta=1
        if(costheta<-1):
            costheta=-1
        value=acos(costheta)
        anglesig=v1*v3.T
        return (value,anglesig)
    def get_gamma(self,Gamma_Tar)->tuple[float,bool]:
        energy=0
        flag=True
        for g in self.gamma:
            g=self._get_dihedral(g)[0]
            if(not(Gamma_Tar-bound<g and g<Gamma_Tar+bound)):
                energy+=HyperKappa*(math.exp(abs(g-Gamma_Tar))-1)#the force from mota should be firstly satisfied.
                flag=False 
            else:
                energy+=Kappa*(g-Gamma_Tar)**2
        return energy,flag
    #get the surfaces made by the two edges in the surfaces related to the shafts, returns two pairs of vectors (type:Edge) AD,AB and CB,CD of ABCD
    #and ab in sur1, cd in sur2 of E
    #attention, ad is always from sur1 to sur2. This is used to calculate the bias
    def _getFoldedShaft(self,e:Edge)->tuple[Edge,Edge]:
        s1=self.Surfaces[e.sur1]
        s2=self.Surfaces[e.sur2]
        return Edge(u=s1.v[e.uid1],v=e.u),Edge(u=s2.v[e.vid2],v=e.v)
    '''
    @return two bias angle of two shafts.
    @return error: no need to calculate
    '''
    def _getBiasAngle(self,e:Edge)->tuple[float,float]:
        if(e.type==1 or e.sur2==-1):
            print("neednot to calc")
            return 0,0
        e1,e2=self._getFoldedShaft(e)
        e1=e1.get_vector()/e1.L
        e2=e2.get_vector()/e2.L
         #unitized e
        s1=self.Surfaces[e.sur1]
        ns1=self._get_edge_and_counterpart(s1,e.uid1,e.vid1)
        beta1=acos(ns1*e1.T)
        s2=self.Surfaces[e.sur2]
        ns2=self._get_edge_and_counterpart(s2,e.uid2,e.vid2)
        beta2=acos(ns2*e2.T)
        angle1=abs(beta1-0.5*math.pi)
        angle2=abs(beta2-0.5*math.pi)
        if(angle1>0.17 or angle2>0.17):
            print("pause",e.id,e.H1,self._getFoldedShaft(e)[0].L,self._getFoldedShaft(e)[1].L)
        return angle1,angle2
    def get_energy(self)->float:

        energy=0
        for e in self.shafts:
            if(e.sur2==-1):
                continue
            value,anglesig=self._get_dihedral(e)
            if((anglesig>0 and e.type==2)or(anglesig<0 and e.type==1)):
                energy+=10*(math.exp(10*min(value,math.pi-value))-1)
                continue
            if(e.type==1):
                continue
            theta1,theta2=self._getBiasAngle(e)
            energy+=YoungModium/e.H1*theta1**2+YoungModium/e.H2*theta2**2
            #E=0.5*EI(\theta)^2/L=k\theta^2/L
        return energy


    def length_repair(self)->tuple[matrix,matrix]:
        A=self._getConstraint()
        #Adx=0
        U,sigma,VT=svd(2*A.T)
        U=np.mat(U)
        VT=np.mat(VT)
        r=0
        dc=np.mat(np.zeros(shape=(U.shape[0],1)))
        e=self.get_length_error()
        while(r<self.numEdge and sigma[r]>1e-12):
            dc-=(VT[r,:]*e)[0,0]*U[:,r]/sigma[r]
            r+=1
        return dc,e
    def get_length_error(self)->matrix:
        err=np.zeros(shape=(3*self.numSurface+(8+1)*self.numEdge))    
        for s in self.Surfaces:
            for j in range(3):
                #surfaces
                ec=s.e[j].get_vector()
                ec=sqrt(ec*ec.T)
                err[s.id*3+j]=ec-s.e[j].L
        #hinge of shafts
        for e in self.shafts:
            #for sur1!=-1 is predefined, so it is no necessary to check again.
            #e.u.id =2e.id+0 and e.v.id =2e.id+1, nothing special.     
            if(e.sur2==-1 or e.type==1):
                continue
            s1=self.Surfaces[e.sur1]
            s2=self.Surfaces[e.sur2]         
            L0=Edge(u=e.u,v=s1.v[e.uid1]).L
            L1=Edge(u=e.u,v=s1.v[e.vid1]).L
            L2=Edge(u=e.v,v=s1.v[e.uid1]).L
            L3=Edge(u=e.v,v=s1.v[e.vid1]).L
            err[self.numSurface*3+e.u.id*4+0]=L0-e.H1
            err[self.numSurface*3+e.u.id*4+1]=L1-sqrt(e.H1**2+e.L**2)
            err[self.numSurface*3+e.v.id*4+0]=L2-sqrt(e.H1**2+e.L**2)
            err[self.numSurface*3+e.v.id*4+1]=L3-e.H1
            L4=Edge(u=e.u,v=s2.v[e.uid2]).L
            L5=Edge(u=e.u,v=s2.v[e.vid2]).L
            L6=Edge(u=e.v,v=s2.v[e.uid2]).L
            L7=Edge(u=e.v,v=s2.v[e.vid2]).L
            err[self.numSurface*3+e.u.id*4+2]=L4-e.H2
            err[self.numSurface*3+e.u.id*4+3]=L5-sqrt(e.H2**2+e.L**2)
            err[self.numSurface*3+e.v.id*4+2]=L6-sqrt(e.H2**2+e.L**2)
            err[self.numSurface*3+e.v.id*4+3]=L7-e.H2
            ec=e.get_vector()
            ec=sqrt(ec*ec.T)
            err[self.numSurface*3+8*self.numEdge+e.id]=ec-e.L
        err=err[self.bar_tag]          
        return np.mat(err).T