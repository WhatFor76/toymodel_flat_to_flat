import numpy as np
from numpy import matrix
from scipy.linalg import null_space
import random
from Unit import *
import math
from math import acos,sqrt
eps=0.0001
dxrate=0.001
sample_rate=0.7
YoungModium=100
bound=0.002
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
        self.loss:list[float]=[0]*1000000
        self.ittime=0
    
    def SGD_getdir(self,Gamma_Tar)->matrix:
        #calculate A
        #u.dx,u.dy,u.dz at line 3*u.id+0..3*u.id+2
        #f0 is the loss of current loss!
        f0=self.loss[self.ittime-1]
        #to save the loss current, it is a legal, and have happended state.

        A,N=self._getConstarint()
        #A*delta=0
        Null=np.mat(null_space(A,rcond=1e-6))
        nullrank=Null.shape[1]
        #found direction in ker(A) with P, andN a ker
        #finally need a target function and a punishment fucntion
        sampledir=random.sample(range(0,nullrank),round(nullrank*sample_rate))
        omega=np.mat(np.zeros((nullrank,1)))
        for i in sampledir:
            dq=N*dxrate*Null[:,i]
            dx=self._getdx(dq)
            self.update(dx)
            f1=self._relax(Gamma_Tar)
            omega[i,0]=(f1-f0)/dxrate
            self.update(-dx)
        #return direction as dir
        delta=-N*Null*omega
        dx=self._getdx(delta)
        return dx
    def _getdx(self,delta:matrix)->matrix:
        dx=np.mat(np.zeros(shape=(3*3*self.numSurface,1)))
        #for type 3, we calculate the dA for twice.
        for s in self.Surfaces:
            for i,v in enumerate(s.v):
                PA=self._getP_A(v)
                delta_a=delta[6*s.id:6*(s.id+1),0]
                dA=PA*delta_a
                for j in range(3):
                    dx[3*(s.id+i)+j,0]=dA[j,0]
        return dx
    def update(self,dx:matrix)->None:
        for s in self.Surfaces:
            for i,v in enumerate(s.v):
                for j in range(3):
                    if(v.type==3):
                        #  v_shifted's(v's) id is too complex, so we use s.id+i instand, with a div fact to deal with repeca for type3
                        #  type3 vertice will be met twice, so we add half for one time. type3 not in the side naturally.
                        v.x[j]+=dx[3*(s.id+i)+j,0]/2
                    else:
                        v.x[j]+=dx[3*(s.id+i)+j,0]

    def stepping(self,step_length:float,Gamma_Tar:float)->int:
        #P(x)=a^2
        # x=x+ldx, where dx have just a mod |ldx|=length
        current_error=self._relax(Gamma_Tar)
        self.loss[self.ittime]=current_error
        self.ittime+=1
        dx=self.SGD_getdir(Gamma_Tar)
        normdx:float=np.linalg.norm(dx)
        if(normdx>eps):
        # if dx is not a 0 num
        # update x=x+dx
        # maybe need to fix the error
        # update the primer equation f(x)=a^2, which is quitely identical to the first part of get_dir.
        # WE USE a ANNEALING METHOD to CONTROL if THIS can BE ACCEPTED!!!
            dx0=dx/normdx*step_length
            self.update(dx0)
            print(":totE:GE:E:TE:",self._relax(Gamma_Tar),self.get_gamma(Gamma_Tar),self.get_energy(),self.get_tensselation_error())
        else:
            print("wronganswer\n")
            return -1
    


    #x+dx dot y+dy=0 means xy+xdy+ydx+dxdy=0, the ignore of dxdy leads to a error in tensselation. 100*Young's*m**2, m is the error from pi/2.
    
    #theta=arccos(b**2-a**2-c**2/ 2ac) and which we want is abs(bata-pi/2-theta), where bata is the angle of surface 
    # and rectangle, c is the width of rectangle and a,b is the length of shifted hinges.  
    # and E=YoungModium beta**2
    def get_tensselation_error(self)->float:
        ad:matrix=None
        ab:matrix=None
        cb:matrix=None
        cd:matrix=None
        energy=0
        for e in self.shafts:
            if(e.sur2==-1 or e.type==1):
                continue
            ad,ab,cb,cd=self._getFoldedShaft(e)
            if(fequ(ad.L,0)):
                continue
            ad=ad.get_vector()
            ab=ab.get_vector()
            cb=cb.get_vector()
            cd=cd.get_vector()
            ad=ad/sqrt(ad*ad.T)
            ab=ab/sqrt(ab*ab.T)
            cb=cb/sqrt(cb*cb.T)
            cd=cd/sqrt(cd*cd.T)
            err1=abs(acos(ad*ab.T)-math.pi/2)
            err2=abs(acos(ad*cd.T)-math.pi/2)
            err3=abs(acos(cb*cd.T)-math.pi/2)
            err4=abs(acos(cb*ab.T)-math.pi/2)
            energy+=10*YoungModium*err1+10*YoungModium*err2+10*YoungModium*err3+10*YoungModium*err4
        return energy

    def _getConstarint(self)->None:

        def _setpa(P:matrix,Pa:matrix,a:int,b:int):
            for i in range(3):
                for j in range(6):
                    P[a+i,b+j]=Pa[i,j]
        P0=np.mat(np.zeros(shape=(3*2*self.numEdge,6*self.numSurface)))
        for e in self.shafts:
            if(e.type==1 and e.sur2!=-1):
                PA=self._getP_A(e.u)
                _setpa(P0,PA,3*(e.id*2+0),6*e.sur1)
                _setpa(P0,-PA,3*(e.id*2+0),6*e.sur2)
                PA=self._getP_A(e.v)
                _setpa(P0,PA,3*(e.id*2+1),6*e.sur1)
                _setpa(P0,-PA,3*(e.id*2+1),6*e.sur2)
        N=np.mat(null_space(P0,rcond=1e-6))

        ad:Edge=None
        ab:Edge=None
        cb:Edge=None
        cd:Edge=None
        U=np.mat(np.zeros(shape=(3*4*self.numEdge,6*self.numSurface)))
        for e in self.shafts:
            if(e.type==2 and e.sur2!=-1):
                ad,ab,cb,cd=self._getFoldedShaft(e)
                PA=self._getP_A(ad.v)
                PB=self._getP_A(ad.u)
                _setpa(U,PA,3*(e.id*4+0),6*e.sur2)
                _setpa(U,-PB,3*(e.id*4+0),6*e.sur1)
                PA=self._getP_A(ab.v)
                PB=self._getP_A(ab.u)
                _setpa(U,PA-PB,3*(e.id*4+1),6*e.sur1)
                PA=self._getP_A(cb.v)
                PB=self._getP_A(cb.u)
                _setpa(U,PA,3*(e.id*4+2),6*e.sur1)
                _setpa(U,-PB,3*(e.id*4+2),6*e.sur2)                
                PA=self._getP_A(cd.v)
                PB=self._getP_A(cd.u)
                _setpa(U,PA-PB,3*(e.id*4+3),6*e.sur2)
        

        M=np.mat(np.zeros(shape=(4*self.numEdge,3*4*self.numEdge)))
        for e in self.shafts:
            if(e.type==2 and e.sur2!=-1):
                ad,ab,cb,cd=self._getFoldedShaft(e)
                #a cross xdy+ydx, so ad's plot at 4i+1
                Ei=ad.get_vector()
                for i in range(3):
                    M[e.id*4+0,3*(4*e.id+1)+i]=Ei[0,i]
                Ei=ab.get_vector()
                for i in range(3):
                    M[e.id*4+0,3*(4*e.id+0)+i]=Ei[0,i]

                Ei=cb.get_vector()
                for i in range(3):
                    M[e.id*4+1,3*(4*e.id+3)+i]=Ei[0,i]
                Ei=cd.get_vector()
                for i in range(3):
                    M[e.id*4+1,3*(4*e.id+2)+i]=Ei[0,i]
                
                Ei=ab.get_vector()
                for i in range(3):
                    M[e.id*4+2,3*(4*e.id+2)+i]=Ei[0,i]
                Ei=cb.get_vector()
                for i in range(3):
                    M[e.id*4+2,3*(4*e.id+1)+i]=Ei[0,i]

                Ei=ad.get_vector()
                for i in range(3):
                    M[e.id*4+3,3*(4*e.id+3)+i]=Ei[0,i]
                Ei=cd.get_vector()
                for i in range(3):
                    M[e.id*4+3,3*(4*e.id+0)+i]=Ei[0,i]
        
        null=M*U*N
        return null,N






        
    def _relax(self,Gamma_Tar:float)->float:
        
        energy=0
       
        energy+=self.get_gamma(Gamma_Tar)         
        #energy+=self.get_energy()
        energy+=self.get_tensselation_error()
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
    def get_gamma(self,Gamma_Tar)->tuple[float,float]:
        energy=0
        gamma=self._get_dihedral(self.gamma)[0]
        if(not(Gamma_Tar-bound<gamma and gamma<Gamma_Tar+bound)):
            energy=10*YoungModium*abs(gamma-Gamma_Tar)#the force from mota should be firstly satisfied.
        return energy
    #get the surfaces made by the two edges in the surfaces related to the shafts, returns two pairs of vectors (type:Edge) AD,AB and CB,CD of ABCD
    #and ab in sur1, cd in sur2 of E
    #attention, ad is always from sur1 to sur2. This is used to calculate the bias
    def _getFoldedShaft(self,e:Edge)->tuple[Edge,Edge,Edge,Edge]:
        s1=self.Surfaces[e.sur1]
        s2=self.Surfaces[e.sur2]
        u1=s1.v[e.uid1]
        v1=s1.v[e.vid1]
        u2=s2.v[e.uid2]
        v2=s2.v[e.vid2]
        return Edge(u=u1,v=u2),Edge(u=u1,v=v1),Edge(u=v2,v=v1),Edge(u=v2,v=u2)
    '''
    @return two bias angle of two shafts.
    @return error: no need to calculate
    '''
    def _getBiasAngle(self,e:Edge)->tuple[float,float]:
        if(e.type==1 or e.sur2==-1):
            print("neednot to calc")
            return 0,0
        wide_e,length_e,ax,bx=self._getFoldedShaft(e)
        c=wide_e.L
        a=e.H1
        b=e.H2
        if(fequ(b-a-c,0)or fequ(a-b-c,0)):#in a line
            return 0,0
        costheta1=(a*a+c*c-b*b)/(2*a*c)
        if(costheta1>1):
            costheta1=1
        if(costheta1<-1):
            costheta1=-1
        costheta2=(b*b+c*c-a*a)/(2*b*c)
        if(costheta2>1):
            costheta2=1
        if(costheta2<-1):
            costheta2=-1
        theta1=acos(costheta1)
        theta2=acos(costheta2)
         #unitized ne
        ne=wide_e.get_vector()/wide_e.L

        s1=self.Surfaces[e.sur1]
        ns=self._get_edge_and_counterpart(s1,e.uid1,e.vid1)
        beta1=acos(ns*ne.T)
        angle1=abs(1.5*math.pi-beta1-theta1)
        if(angle1>0.175):#which means a error calculate for angle1 is a small angle less than 10 deg.
            angle1=abs(beta1-theta1-0.5*math.pi)
        
        ne=-ne
        s2=self.Surfaces[e.sur2]
        ns=self._get_edge_and_counterpart(s2,e.uid2,e.vid2)
        beta2=acos(ns*ne.T)
        angle2=abs(1.5*math.pi-beta2-theta2)
        if(angle2>0.175):#which means a error calculate for angle1 is a small angle less than 10 deg.
            angle2=abs(beta2-theta2-0.5*math.pi)
        if(angle1>0.17):
            print("pause")
        return angle1,angle2
    def get_energy(self)->float:

        energy=0
        for e in self.shafts:
            if(e.sur2==-1):
                continue
            value,anglesig=self._get_dihedral(e)
            if((anglesig<0 and e.type==2)or(anglesig>0 and e.type==1)):
                energy+=10*(math.exp(10*min(value,math.pi-value))-1)
                continue
            if(e.type==1):
                continue
            theta1,theta2=self._getBiasAngle(e)
            energy+=YoungModium*theta1**2+YoungModium*theta2**2
        
        return energy

    def _mount_valley_error(self)->bool:
        energy=0
        for shaft in self.shafts:
            if(shaft.sur2==-1):
                continue
            value,anglesig=self._get_dihedral(shaft)
            if(anglesig<0 and shaft.type==2):
                energy+=10*(math.exp(10*min(value,math.pi-value))-1)
            if(anglesig>0 and shaft.type==1):
                energy+=10*(math.exp(10*min(value,math.pi-value))-1)
        return energy
    def _getP_A(self,A:Vertex)->matrix:
        T=np.zeros(shape=(4,3,6))
        T[0,1,2]=1
        T[0,2,1]=-1
        T[1,0,2]=-1
        T[1,2,0]=1
        T[2,0,1]=1
        T[2,1,0]=-1
        T[3,0,3]=1
        T[3,1,4]=1
        T[3,2,5]=1
        X=[A.x[0],A.x[1],A.x[2],1]
        P=np.zeros(shape=(3,6))
        for i in range(4):
            for j in range(3):
                for k in range(6):
                    P[j,k]+=T[i,j,k]*X[i]
        return P