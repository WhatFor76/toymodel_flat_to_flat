'''
Copyright @WhatFor, contact with 332449502@qq.com
Used in step 2.
this code is to support the simulate flamebuild.py, three basic types are defined:
Vertex, Edge and Surface, Surface is the most basic unit of the simulate, detailed description are attached in front of each class.
and all edges and surfaces can be presented by vectors for a deeper calculate
'''

import numpy as np
import math
#vertex is a struct with members: index, type, x=[x,y,z]
#index is the index in ideal topology
class Vertex:
    
    def __init__(self,index:int=0,Type:int=0,x:float=0,y:float=0,z:float=0):
        self.id=index
        self.type=Type
        self.x=[x,y,z]
        #to get the edge connect itself and other vertice, 6 edges max. used in step 3.
        self.adj_edge:list[Edge]=[None]*6
        self.adj_num=0
        self.adj_surface:list[Surface]=[None]*6
    #2D_position currently is read from outGraph,'2D'(!), sorted by pole angle, and in the future may be reset specifically.
    #used to sort the edges clockwise to write the Equal Constraints in step3, and side vertex this will be not dealt
    def sort_edges(self):
        if(self.adj_num<6):
            return 
        angle=[0,0,0,0,0,0]
        for i in range(self.adj_num):
            y=self.adj_edge[i].v.x[1]-self.adj_edge[i].u.x[1]
            x=self.adj_edge[i].v.x[0]-self.adj_edge[i].u.x[0]
            if(self.adj_edge[i].u.id==self.id):
                angle[i]=math.atan2(y,x)
            else:
                angle[i]=math.atan2(-y,-x)
            if(angle[i]<0):
                angle[i]+=2*math.pi
        for i in range(self.adj_num):
            for j in range(i+1,self.adj_num):
                if(angle[i]>angle[j]):
                    t=angle[i]
                    angle[i]=angle[j]
                    angle[j]=t
                    t=self.adj_edge[i]
                    self.adj_edge[i]=self.adj_edge[j]
                    self.adj_edge[j]=t
        return
    #this need a parameter sur:list[Surface]  for that Surface is not defined yet when we define class Vertex
    def sort_surfaces(self,sur):
        if(self.adj_num<6):
            return
        for i in range(self.adj_num-1):
            if(self.adj_edge[i].sur1==self.adj_edge[i+1].sur1 or self.adj_edge[i].sur1==self.adj_edge[i+1].sur2):
                self.adj_surface[i]=sur[self.adj_edge[i].sur1]
            #a collaps confliction will not happen when edges are correctly defined.
            if(self.adj_edge[i].sur2==self.adj_edge[i+1].sur1 or self.adj_edge[i].sur2==self.adj_edge[i+1].sur2):
                self.adj_surface[i]=sur[self.adj_edge[i].sur2]
        #to meet a loop finally
        if(self.adj_edge[self.adj_num-1].sur1==self.adj_edge[0].sur1 or self.adj_edge[self.adj_num-1].sur1==self.adj_edge[0].sur2):
            self.adj_surface[self.adj_num-1]=sur[self.adj_edge[self.adj_num-1].sur1]
        if(self.adj_edge[self.adj_num-1].sur2==self.adj_edge[0].sur1 or self.adj_edge[self.adj_num-1].sur2==self.adj_edge[0].sur2):
            self.adj_surface[self.adj_num-1]=sur[self.adj_edge[self.adj_num-1].sur2]
#edge is a struct, with members: index, type, (u,v): pointer of vertex type, (sur1,sur2):index of surface
#L is the riger length, will not change once a edge object is created.
#if sur1, sur2 are not existed or meaningless, be tagged as -1 and sur2 have a priority to be -1.
#and vector can be temprorily calculated by u.x-v.x as a 1*n matrix.
class Edge:
    #type: #0 notshaft, #1 valleyshaft(thickness=0) #2mountshaft
    def __init__(self,index:int=0,Type:int=0,u:Vertex=None,v:Vertex=None,sur1:int=-1,sur2:int=-1):
        (self.id,self.type,self.u,self.v)=(index,Type,u,v)
        #this three parameters are not defined before all read, we reset it
        self.sur1=sur1
        self.sur2=sur2
        self.L=0
        for i in range(3):
            self.L+=(self.u.x[i]-self.v.x[i])**2 
        self.L=math.sqrt(self.L)
        self.u_shifted=None
        self.v_shifted=None
        #H1,H2 the height of shafts in sur1 and sur2
        self.uid1:int=None
        self.uid2:int=None
        self.vid1:int=None
        self.vid2:int=None
        self.H1=0
        self.H2=0
    #u is the start vertex, for that v can be None if direction can be ignored. Default: u->v
    def get_vector(self,u:Vertex=None)->np.matrix:
        if(self.u==u or u==None):
            t=np.mat([self.v.x[i]-self.u.x[i] for i in range(3)])
        else:
            t=np.mat([self.u.x[i]-self.v.x[i] for i in range(3)])
        return t
    #this need a parameter sur:list[Surface] Shifted  for that Surface is not defined yet when we define class Vertex
    #if we want a more precise shifted poistion of the shaft, surface1 and surface2 should be both used
    #but we now no need to do so, only surface1 is enough, and sur2 inexisted will not lead to the disappearance of this shaft.
    #you should set v and then set e s so that the mount-vally relationships can be correctly set.
    def set_e_shifted(self,sur,offset:float,off2)->None:

        s:Surface=sur[self.sur1]
        n=s.uppersurface()*offset
        u:Vertex=None
        v:Vertex=None
        #both type 1 and 2 edges shares the same position rule. and for type1 edges, the vertex is shared, so it is not the real, only has a correct position. 
        for i,vx in enumerate(s.v):
            if(self.u==vx):
                u=s.v_shifted[i]
        self.u_shifted=Vertex(2*self.id+0,x=u.x[0]-n[0,0],y=u.x[1]-n[0,1],z=u.x[2]-n[0,2])

        for i,vx in enumerate(s.v):
            if(self.v==vx):
                v=s.v_shifted[i]
        self.v_shifted=Vertex(2*self.id+1,x=v.x[0]-n[0,0],y=v.x[1]-n[0,1],z=v.x[2]-n[0,2])
        #type1 vertex needs a shared address, so surfaces' v_shifted list use shafts, and not to build hinge when comfirm its stability. 
        #these vertice are marked as type 3 vertice, and the same id as edges. type3 Vertice not in the side.
        if(self.sur2==-1):
            return
        if(self.type==1):
            self.u_shifted.type=3
            self.v_shifted.type=3
            s:Surface=sur[self.sur1]
            for i,v in enumerate(s.v):
                if(v==self.u):
                    s.v_shifted[i]=self.u_shifted
                if(v==self.v):
                    s.v_shifted[i]=self.v_shifted
            s:Surface=sur[self.sur2]
            for i,v in enumerate(s.v):
                if(v==self.u):
                    s.v_shifted[i]=self.u_shifted
                if(v==self.v):
                    s.v_shifted[i]=self.v_shifted
        
        for i,v in enumerate(sur[self.sur1].v):
            if(self.u==v):
                self.uid1=i
            if(self.v==v):
                self.vid1=i
        for i,v in enumerate(sur[self.sur2].v):
            if(self.u==v):
                self.uid2=i
            if(self.v==v):
                self.vid2=i
    def set_H(self,off1:float,off2:float)->None:
        self.H1=off1
        self.H2=off2
# surface is the basic class.
# when creating, u,v,w is set.
# uppersurface should be set after 3 edge was set by reading the file "outGraph.out", it is two vector pointer n0 and n1=e0 and e1, 
# make sure that  (n:n0 cross n1) dot (v:a vector appointed)>0, means they are in the same half space,
# and calculate temprorily as a 1*n matrix.

class Surface:
    #when height existed, this model have a initial in 2D plate
    def __init__(self,index:int=0,u:Vertex=None,v:Vertex=None,w:Vertex=None):
        self.id=index
        self.v:list[Vertex]=[u,v,w]
        self.e:list[Edge]=[None,None,None]
        self.v_shifted:list[Vertex]=[None,None,None]
        self.numEdge=0
    #current, this vertex is only the ones that face to "z+"direction
    #if not the default, n dot v>0
    #it is a line vertex not a column
    def setuppersurface(self,v:np.matrix=np.mat([[0,0,1]]))->None:
        n1=Edge(u=self.v[0],v=self.v[1])
        n2=Edge(u=self.v[0],v=self.v[2])
        a=np.mat(np.cross(n1.get_vector(),n2.get_vector()))
        t=a*v.T
        if(t>0):
            self.n=[n1,n2]
            self.nord=1
        else:
            self.n=[n2,n1]
            self.nord=2

    def resetuppersurface_edges(self)->None:
        n1=Edge(u=self.v[0],v=self.v[1])
        n2=Edge(u=self.v[0],v=self.v[2])
        if(self.nord==1):
            self.n=[n1,n2]
        elif(self.nord==2):
            self.n=[n2,n1] 
        self.e[0]=Edge(u=self.v[0],v=self.v[1])        
        self.e[1]=Edge(u=self.v[1],v=self.v[2])  
        self.e[2]=Edge(u=self.v[2],v=self.v[0])  
    #return it's unit vector(line vector)
    #return it's unit vector(line vector)
    def uppersurface(self)->np.matrix:
        a=np.cross(self.n[0].get_vector(),self.n[1].get_vector())
        a=np.mat(a)/np.linalg.norm(a)
        return a    
    #it is a line vertex not a column, two side and the bisector itself
    def getAngleBisector(self,v:Vertex)->tuple[np.matrix,np.matrix,np.matrix]:
        n=0
        bis=[None,None]
        for e in self.e:
            if(e.u==v):
                e0=e.get_vector(v)
                bis[n]=e0/math.sqrt(e0*e0.T)
                n+=1
            elif(e.v==v):
                e0=e.get_vector(v)
                bis[n]=e0/math.sqrt(e0*e0.T)
                n+=1
        b=(bis[0]+bis[1])
        b=b/math.sqrt(b*b.T)
        return bis[0],bis[1],b
    #also a line vector
    def getGravity(self)->np.matrix:
        G=np.mat([(self.v[0].x[i]+self.v[1].x[i]+self.v[2].x[i])/3 for i in range(3)])
        return G
    #you should set v and then set e s so that the mount-vally relationships can be correctly set.
    def set_v_shifted(self,offset:list[float])->None:
        #a type1 edges here will only have a correct position, the id is not true.
        for i in range(3):
            self.v_shifted[i]=Vertex(3*self.id+i,x=self.v[i].x[0],y=self.v[i].x[1],z=self.v[i].x[2])
            self.v_shifted[i].x[0]+=offset[0]
            self.v_shifted[i].x[1]+=offset[1]
            self.v_shifted[i].x[2]+=offset[2]
            #only the vertex not in the flag edge needs a new independent vertex. the other two is set when set the edge.
            #but still in the order of the original vertex list of the surface. 




if __name__=="__main__":
    None