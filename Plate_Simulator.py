'''
Copyright @WhatFor, contact with 332449502@qq.com
Used in step 2.
for outGraph.out is only a plate x,y position, this model can only initize a model that 2D counterparter is a plate, not a curved surface.
all of the detail can be found in each methed.
THIS BUFFERED SIMULATOR ONLY USED TO BUILD THE flat to flat model, a not so good simulator, for the len error is not controlled,
the only difference to the new one.
'''

import numpy as np
from scipy.linalg import null_space
import random
from Unit import *
import math
E=30#Young's modium
sample_rate=0.8
dxrate=1e-8
eps=1e-6
class SGD_simulate:
    #this share the same memory with V,E which used to initize the class, just a light copy, not a deep one
    def __init__(self,N:int,M:int,P:int,V:list[Vertex],S:list[Surface],shafts:list[Edge]) -> None:
        self.numEdge=M
        self.numVextex=N
        self.numSurface=P
        self.Surfaces=S
        self.Vertice=V
        self.shafts=shafts
        self.A=np.zeros((self.numEdge,3*self.numVextex))#3 for x,y,z each, that in 3D space
        #loss to plot in figure to evaluate the quality of a simulate, to help to determine when to quit and output.
        self.loss:list[float]=[0]*10000
        self.ittime=0
    #
    def SGD_getdir(self)->np.matrix:
        #calculate A
        #u.dx,u.dy,u.dz at line 3*u.id+0..3*u.id+2
        #f0 is the loss of current loss!
        f0=self.loss[self.ittime-1]
        #to save the loss current, it is a legal, and have happended state.


        for e in self.shafts:
            for i in range(3):
                self.A[e.id,3*e.u.id+i]=e.u.x[i]-e.v.x[i]
                self.A[e.id,3*e.v.id+i]=e.v.x[i]-e.u.x[i]

        #Adx=0
        B=np.mat(null_space(self.A,rcond=1e-6))
        nullrank=B.shape[1]
        #found direction in ker(A) with P, and a ker
        #finally need a target function and a punishment fucntion
        sampledir=random.sample(range(0,nullrank),round(nullrank*sample_rate))
        omega=np.mat(np.zeros((nullrank,1)))
        for i in sampledir:
            dq=dxrate*B[:,i]
            self.update(dq)
            f1,fx=self._punishmentfunc()
            omega[i,0]=(f1-f0)/dxrate
            self.update(-dq)
        #return direction as dir
        dir=-B*omega
        return dir
    def update(self,dx:np.matrix)->None:
        for v in self.Vertice:
            for i in range(3):
                v.x[i]+=dx[3*v.id+i,0]
    '''
    stepping will return a state code, that -1, wronganswer(div0),0 normal,1 out of torrent
    '''
    def stepping(self,torrent_temperature:float,step_length:float)->int:
        #P(x)=a^2
        # x=x+ldx, where dx have just a mod |ldx|=length
        current_error,loss_error=self._punishmentfunc()
        print(loss_error)
        self.loss[self.ittime]=current_error
        self.ittime+=1
        dx=self.SGD_getdir()
        normdx:float=np.linalg.norm(dx)
        if(normdx>eps):
        # if dx is not a 0 num
        # update x=x+dx
        # maybe need to fix the error
        # update the primer equation f(x)=a^2, which is quitely identical to the first part of get_dir.
        # WE USE a ANNEALING METHOD to CONTROL if THIS can BE ACCEPTED!!!
            dx0=dx/normdx*step_length
            self.update(dx0)
            next_error,ex=self._punishmentfunc()
            if(next_error-current_error<torrent_temperature):
                #we accept this update
                return 0
            else:
                dx0=-dx/normdx*step_length
                self.update(dx0)
                return 1
        else:
            print("wronganswer\n")
            return -1
    
    def _mount_valley_error(self)->bool:
        energy=0
        for shaft in self.shafts:
            if(shaft.sur2==-1):
                continue
            value,anglesig=self._get_dihedral(shaft)
            if(anglesig>0 and shaft.type==2):
                energy+=10*(math.exp(10*min(value,math.pi-value))-1)
            if(anglesig<0 and shaft.type==1):
                energy+=10*(math.exp(10*min(value,math.pi-value))-1)
        return energy
    #return the cos value of the upper surface, if >0 type2 and <0 type1  
    def _get_edge_and_counterpart(self,a:Surface,shaft:Edge):
        #the edge the same as the shaft
        t:np.matrix=None
        for e in a.e:
            if(e!=shaft):
                if(e.u==shaft.u):
                    t=e.get_vector(e.u)
                if(e.v==shaft.u):
                    t=e.get_vector(e.v)
        v1=t
        v2=shaft.get_vector()
        #as the method we build a shaft's surface edge, (u->v)-(x->y)=0
        v=v1-v1*np.transpose(v2)/(v2*np.transpose(v2))*v2
        return v
        
    #if return >=0 a valley, and <=0 a mount.
    def _get_dihedral(self,shaft:Edge):
        a:Surface=self.Surfaces[shaft.sur1]
        b:Surface=self.Surfaces[shaft.sur2]
        v1=self._get_edge_and_counterpart(a,shaft)
        v2=self._get_edge_and_counterpart(b,shaft)
        p=a.uppersurface()
        normp=np.linalg.norm(p)
        q=b.uppersurface()
        normq=np.linalg.norm(q)
        v3=p/normp+q/normq
        costheta=v1*np.transpose(v2)/np.linalg.norm(v1)/np.linalg.norm(v2)
        #to fix the float error such as 1.0+1e-16
        if(costheta>1):
            costheta=1
        if(costheta<-1):
            costheta=-1
        value=math.acos(costheta)
        anglesig=v1*np.transpose(v3)
        return (value,anglesig)
   
    #in this data scale, len_error is too low to fix, but not even a bad improve, rather.
    def get_length_error(self)->float:
        error=0
        for shaft in self.shafts:
            error+=E*(shaft.L-np.linalg.norm(shaft.get_vector()))**2
        return error

    def _punishmentfunc(self)->float:
        #if deviation angle is not precise enough, 2-order offset may be considered as the shaft is compressed or stretched. 

        #energy=self._mount_valley_error()#a nearly infinity value
        energy=0
        energy0=0
        #the following two items are illegal, so a large punishment weight
        energy+=self._mount_valley_error()
        energy+=self.get_length_error()
        '''
        for s in self.Surfaces:
            if(abs(s.e[0].L-s.e[1].L)<eps):
                a=s.getGravity()
                energy0+=(math.sqrt(a[0,0]**2+a[0,1]**2+(a[0,2]-1.85)**2)-1.85)**2
        energy+=energy0
        '''
        #if legal, the punishment function is linear.
        for shaft in self.shafts:
            if(shaft.type==1 and shaft.sur2!=-1):
                value,sig=self._get_dihedral(shaft)
                energy+=(value)**2
        return energy,energy0
    # need to be write as a file after step 3 is tested.
    def get_all_dihedral(self):
        for shaft in self.shafts:
            if(shaft.sur2!=-1):
                value,sig=self._get_dihedral(shaft)
                if(sig>0):
                    sig=1
                else:
                    sig=-1
                print(shaft.id,sig,value)
        

        
if __name__=="__main__":
    None