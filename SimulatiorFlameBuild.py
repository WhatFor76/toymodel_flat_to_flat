'''
Copyright @WhatFor, contact with 332449502@qq.com
Used in step 2.
this is the main function of step 2, simulate the plate model to the targe state and get its dihedral angle between each two surfaces. 
read from outGraph.out and write to dihedral.in
dihedral.in have only E-th lines, each line 3 members: index of edge, type of the dihedral(the same as the type of edge but use signal to describe
it, and type -1 edge's signal and angle value is meaningless), angle value.
in a dual graph, it can be useful to describe constraints.

'''
import numpy as np
import math
from Unit import *
import Plate_Simulator
import matplotlib.pyplot as plt
import random
from plotclasses import PlotModel
import pickle
f=open("outGraph.out","r")
item=f.readline()
(N,P,M)=item.split(maxsplit=2)
N=int(N)
P=int(P)
M=int(M)
random.seed(3)
V0:list[Vertex]=[None]*N
Shaft:list[Edge]=[None]*M
S:list[Surface]=[None]*P
for i in range(N):
    item=f.readline()
    (id,t,x,y)=item.split(maxsplit=3)
    V0[i]=Vertex(int(id),int(t),float(x),float(y),0)
#not used
#surface2 must be -1 if this edge is in the surrounding 
v_surface=[[None,None,None,0] for i in range(P)]
for i in range(M):
    item=f.readline()
    (id,t,u,v,surface1,surface2)=item.split(maxsplit=6)
    u=int(u)
    v=int(v)
    sur1=int(surface1)
    sur2=int(surface2)
    Shaft[i]=Edge(int(id),int(t),V0[u],V0[v],sur1,sur2)
    if(V0[u] not in v_surface[sur1]):
        v_surface[sur1][v_surface[sur1][3]]=V0[u]
        v_surface[sur1][3]+=1
    if(V0[v] not in v_surface[sur1]):
        v_surface[sur1][v_surface[sur1][3]]=V0[v]
        v_surface[sur1][3]+=1
    if(sur2!=-1 and V0[u] not in v_surface[sur2]):
        v_surface[sur2][v_surface[sur2][3]]=V0[u]
        v_surface[sur2][3]+=1
    if(sur2!=-1 and V0[v] not in v_surface[sur2]):
        v_surface[sur2][v_surface[sur2][3]]=V0[v]
        v_surface[sur2][3]+=1
#cannot set u,v when reading, for surface need the prime V as it's position(x,y,z)data.
for i in range(P):
    S[i]=Surface(i,v_surface[i][0],v_surface[i][1],v_surface[i][2])
for i in range(M):
    S[Shaft[i].sur1].e[S[Shaft[i].sur1].numEdge]=Shaft[i]
    S[Shaft[i].sur1].numEdge+=1
    if(Shaft[i].sur2!=-1):
        S[Shaft[i].sur2].e[S[Shaft[i].sur2].numEdge]=Shaft[i]
        S[Shaft[i].sur2].numEdge+=1
for e in Shaft:
    e.u.adj_edge[e.u.adj_num]=e
    e.v.adj_edge[e.v.adj_num]=e
    e.u.adj_num+=1
    e.v.adj_num+=1
for v in V0:
    v.sort_edges()
    v.sort_surfaces(sur=S)
for i in range(P):
    S[i].setuppersurface()
with open("outState1.pickle","wb") as f:
    pickle.dump((S,Shaft,V0),f)

solver=Plate_Simulator.SGD_simulate(N,M,P,V0,S,Shaft)
ploter=PlotModel()




#here, we use a likely simulated annealing to avoid a big shift in loss, which leads to a error when near the target point.
#but this is not a classic Temperature choice, so we may update it.
def sigmond(x:float)->float:
    return 1/(1+math.exp(-x))
init_temperature:float=2000
iterate_time=1000
step=0.01
a=-10
b=30
q=0
errortime_torrent=5
coff_step=2
coff_torr=3
p=0
T=init_temperature
for i in range(iterate_time):    
    state=solver.stepping(T,step)
    if(state==1):
        q+=1
        if(q>errortime_torrent):
            step/=coff_step
            q=0
            errortime_torrent*=coff_torr
    else:
        p+=1
        T=(1-sigmond(p/iterate_time*b+a))*init_temperature
    if(i%10==0):
        ploter.plot(solver.Surfaces)

#to print the simulated imformation for the next step. all necessary imformation is saved.
with open("outState2.pickle","wb") as f:
    pickle.dump((S,Shaft,V0),f)


plt.ioff()
plt.figure()
solver.get_all_dihedral()
plt.plot([i for i in range(solver.ittime)],solver.loss[0:solver.ittime])
print(solver.get_length_error())
plt.show()

None 


        