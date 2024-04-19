import numpy as np
from Unit import *
import OTPSimu
import pickle
from plotclasses import PlotModel
import matplotlib.pyplot as plt
from math import cos,exp
class GlobalParameter:
    def __init__(self):
        ax=None
        bx=None
        with open('outState1.pickle','rb') as f:
           ax=pickle.load(f)
        with open('outState2.pickle','rb') as f:
           bx=pickle.load(f)
        '''
        self.Surfaces:list[Surface]=ax[0]
        self.Shafts:list[Edge]=ax[1]
        '''
        
        self.Surfaces:list[Surface]=bx[0]
        self.Shafts:list[Edge]=bx[1]
        with open('outSurfacesOffset2.pickle','rb') as f:
            self.SurfacesOff2:list[float]=pickle.load(f)
        
        with open('outSurfacesOffset1.pickle','rb') as f:
            self.SurfacesOff1:list[float]=pickle.load(f)
        with open('outHingeShift.pickle','rb') as f:
            self.shift:list[float]=pickle.load(f)
        self.E=len(self.Shafts)   
        self.S=len(self.Surfaces)
eps=0.0001
GV=GlobalParameter()
#self.S , V, E will all be changalbe in this simulation
def main():
    v0=GV.SurfacesOff2
    for s in GV.Surfaces:
        s.set_v_shifted(v0[s.id*3:(s.id+1)*3])
    for e in GV.Shafts:
        e.set_e_shifted(GV.Surfaces,GV.shift[GV.Surfaces[e.sur1].id]-GV.shift[GV.S+e.id],GV.shift[GV.Surfaces[e.sur2].id]-GV.shift[GV.S+e.id])
        #get H1,H2 of sur1, sur2 for all edges
        if(e.type==2):
            e.set_H(GV.shift[GV.Surfaces[e.sur1].id]-GV.shift[GV.S+e.id],GV.shift[GV.Surfaces[e.sur2].id]-GV.shift[GV.S+e.id])
    #reset GV.Surfaces by v_shifted as v, so that we can get a real diheral, with only the id and position changed, no relationship changing.

    for e in GV.Shafts:
        e.u=e.u_shifted
        e.v=e.v_shifted
    for s in GV.Surfaces:
        s.v[0]=s.v_shifted[0]
        s.v[1]=s.v_shifted[1]
        s.v[2]=s.v_shifted[2]
        s.resetuppersurface_edges()
        


    solver=OTPSimu.OTPSimu(GV.E,GV.S,GV.Surfaces,GV.Shafts,[GV.Shafts[0],GV.Shafts[4],GV.Shafts[8]])
    ploter=PlotModel()
    f=open("outloss","w")
    for i in range(0,314,2):
        dc,e=solver.length_repair()
        solver.update(dc)
        print((e.T*e)[0,0])
        error0=0
        p=0
        step=0.001
        q=0
        insurance=100
        while step>1e-6:
            if(insurance==0):
                step/=10
                insurance=100
            error1=solver.stepping(step,0.01*i)
            diserr1,flag=solver.get_gamma(0.01*i)
            print(i,error1,diserr1,file=f,flush=True)
            if(flag==True):
                if(error1>error0):
                    p=p+1
                if(p>=10):
                    step/=2
                if(abs(error1-error0)<0.0001):
                    q=q+1
                else:
                    q=0
                if(q>=3):
                    break
            else:
                p=0
                insurance-=1
            error0=error1
        #if(i%10==0):
            ploter.plot_with_hinge(solver.Surfaces,solver.shafts)
        
    

    plt.ioff()
    plt.figure()
    plt.plot([i for i in range(solver.ittime)],solver.loss[0:solver.ittime])
    plt.show()
            

if __name__=="__main__":
    main()