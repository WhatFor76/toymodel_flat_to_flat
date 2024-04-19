'''
Copyright @WhatFor, contact with 332449502@qq.com
Used in step 2.
just a class to plot the model while simulating.
'''
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from math import sqrt
import labellines
import random
import numpy as np
from Unit import *
class PlotModel:
    def __init__(self):
        fig=plt.figure()
        self.ax=Axes3D(fig)
        self.t=0
        self.ax.view_init(elev=80,azim=-53)

        plt.ion()
    def plot(self,S:list[Surface]):
        plt.cla()
        for s in S:
            verts=[list(zip(((s.v[0].x[0]+1.5)*0.3,(s.v[1].x[0]+1.5)*0.3,(s.v[2].x[0]+1.5)*0.3),((s.v[0].x[1]+1.5)*0.3,(s.v[1].x[1]+1.5)*0.3,(s.v[2].x[1]+1.5)*0.3),(s.v[0].x[2]*0.3,s.v[1].x[2]*0.3,s.v[2].x[2]*0.3)))]
            self.ax.add_collection3d(Poly3DCollection(verts)) 
        plt.pause(0.01)
    def plot_with_hinge(self,S:list[Surface],E:list[Edge]):
        plt.cla()
        for s in S:
            verts=[list(zip(((s.v[0].x[0]+1.5)*0.3,(s.v[1].x[0]+1.5)*0.3,(s.v[2].x[0]+1.5)*0.3),((s.v[0].x[1]+1.5)*0.3,(s.v[1].x[1]+1.5)*0.3,(s.v[2].x[1]+1.5)*0.3),(s.v[0].x[2]*0.3,s.v[1].x[2]*0.3,s.v[2].x[2]*0.3)))]
            self.ax.add_collection3d(Poly3DCollection(verts)) 
        for e in E:
            if(e.sur2==-1 or e.type==1):
                continue
            n1=S[e.sur1].uppersurface()
            n2=S[e.sur2].uppersurface()
            ve=e.get_vector()
            ve=ve/sqrt(ve*ve.T)
            mid=np.mat([(e.u.x[0]+e.v.x[0])*0.5,(e.u.x[1]+e.v.x[1])*0.5,(e.u.x[2]+e.v.x[2])*0.5])
            x1=mid+ve*0.05
            x2=mid-ve*0.05
            x1s1=x1+n1*e.H1
            x2s1=x2+n1*e.H1
            x1s2=x1+n2*e.H2
            x2s2=x2+n2*e.H2
            verts1=[list(zip(((x1[0,0]+1.5)*0.3,(x2[0,0]+1.5)*0.3,(x2s1[0,0]+1.5)*0.3,(x1s1[0,0]+1.5)*0.3),((x1[0,1]+1.5)*0.3,(x2[0,1]+1.5)*0.3,(x2s1[0,1]+1.5)*0.3,(x1s1[0,1]+1.5)*0.3),((x1[0,2])*0.3,(x2[0,2])*0.3,(x2s1[0,2])*0.3,(x1s1[0,2])*0.3)))]
            verts2=[list(zip(((x1[0,0]+1.5)*0.3,(x2[0,0]+1.5)*0.3,(x2s2[0,0]+1.5)*0.3,(x1s2[0,0]+1.5)*0.3),((x1[0,1]+1.5)*0.3,(x2[0,1]+1.5)*0.3,(x2s2[0,1]+1.5)*0.3,(x1s2[0,1]+1.5)*0.3),((x1[0,2])*0.3,(x2[0,2])*0.3,(x2s2[0,2])*0.3,(x1s2[0,2])*0.3)))]
            self.ax.add_collection3d(Poly3DCollection(verts1,fc='grey')) 
            self.ax.add_collection3d(Poly3DCollection(verts2,fc='grey')) 
        plt.pause(0.01)
    def plotshifted(self,S:list[Surface]):
        plt.cla()
        for s in S:
            verts=[list(zip(((s.v_shifted[0].x[0]+1.5)*0.3,(s.v_shifted[1].x[0]+1.5)*0.3,(s.v_shifted[2].x[0]+1.5)*0.3),((s.v_shifted[0].x[1]+1.5)*0.3,(s.v_shifted[1].x[1]+1.5)*0.3,(s.v_shifted[2].x[1]+1.5)*0.3),(s.v_shifted[0].x[2]*0.3,s.v_shifted[1].x[2]*0.3,s.v_shifted[2].x[2]*0.3)))]
            self.ax.add_collection3d(Poly3DCollection(verts)) 
        plt.pause(0.01)
    def _randcolor(self):
        colorArray=['1','2','3','4','5','6','7','8','9','A','B','C','D','E','F']
        color="#"
        for i in range(6):
            color+=colorArray[random.randint(0,14)]
        return color
    #this is the plot of the model shifted after linear programming, where the v_shifted is not initinalized.
    def plot2(self,S:list[Surface],v0:list[float]):
        plt.cla()
        for i in range(len(S)):
            s=S[i]
            verts=[list(zip(((s.v[0].x[0]+v0[i*3+0,0]+1.5)*0.3,(s.v[1].x[0]+v0[i*3+0,0]+1.5)*0.3,(s.v[2].x[0]+v0[i*3+0,0]+1.5)*0.3),((s.v[0].x[1]+v0[i*3+1,0]+1.5)*0.3,(s.v[1].x[1]+v0[i*3+1,0]+1.5)*0.3,(s.v[2].x[1]+v0[i*3+1,0]+1.5)*0.3),((s.v[0].x[2]+v0[i*3+2,0])*0.3,(s.v[1].x[2]+v0[i*3+2,0])*0.3,(s.v[2].x[2]+v0[i*3+2,0])*0.3)))]
            self.ax.add_collection3d(Poly3DCollection(verts,fc=self._randcolor())) 
        plt.pause(0.01)
    #to plot a 2D figure with height of shafts included.
    def plotflat(self, S:list[Surface],x:list[float]):
        plt.figure(figsize=(8, 8))
        #plt.scatter(X, Y, s=100, color="red")
        for s in S:
            g=s.getGravity()
            for e in s.e:        
                plt.plot([e.u.x[0],e.v.x[0]],[e.u.x[1],e.v.x[1]], linewidth=1)
                x0=(e.u.x[0]+e.v.x[0]+g[0,0])/3
                y0=(e.u.x[1]+e.v.x[1]+g[0,1])/3 
                if(e.sur2!=-1 and e.type==2):
                    plt.text(x0,y0,str(int((x[s.id,0]-x[e.id+len(S),0])*200+0.5)))
    def plotmovement(self,S:list[Surface],vec:list[float]):
        plt.cla()
        for s in S:
            verts=[list(zip((s.v_shifted[0].x[0],s.v_shifted[1].x[0],s.v_shifted[2].x[0]),(s.v_shifted[0].x[1],s.v_shifted[1].x[1],s.v_shifted[2].x[1]),(s.v_shifted[0].x[2],s.v_shifted[1].x[2],s.v_shifted[2].x[2])))]
            self.ax.add_collection3d(Poly3DCollection(verts,fc=self._randcolor()))
        for i,s in enumerate(S):
            for j in range(3):
                L=np.linalg.norm(vec[3*(3*i+j):3*(3*i+j+1)])
                self.ax.quiver(s.v_shifted[j].x[0],s.v_shifted[j].x[1],s.v_shifted[j].x[2],vec[3*(3*i+j)+0],vec[3*(3*i+j)+1],vec[3*(3*i+j)+2],length=L)
        plt.pause(0.01)

#for i, label in enumerate(annotations):
#    plt.annotate(label, (X[i], Y[i]))
        
if __name__=="__main__":
    None