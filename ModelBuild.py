'''
Copyright @WhatFor, contact with 332449502@qq.com
Used in step 4.
this is the core of the model build, a class of
the modify of the prime model to get a hole, the tolerance of dimension and 
joint or any other that is necessary to the 3D printer.
the out put of the model is seperately a single unit for one reference.
the joint need to be dealt in Solidworks for its complex structure. 
'''

import vtk
from Unit import *
import pickle
from math import sqrt
import copy
class GlobalParameter:
    def __init__(self):
        with open('outState1.pickle','rb') as f:
            a=pickle.load(f)
        self.Surfaces:list[Surface]=a[0]
        self.Shafts:list[Edge]=a[1]
        self.Vertice:list[Vertex]=a[2]
        #attention that surfacei z at zoff[i], and shafti z at zoff[GV.S+i]
        with open('outHingeShift.pickle','rb') as f:
            self.z_off:list[float]=pickle.load(f)
        self.V=len(self.Vertice)
        self.E=len(self.Shafts)   
        self.S=len(self.Surfaces)
        #parameters of the 3d model.
        self.realrate=100#1 unit is 100mm in the real world 
        self.t=5.5/self.realrate
        self.mh=6/self.realrate#This should be the minimum height of a type2 edge, the same as thicknize settings.
        self.epsil=0.5/self.realrate#0.5mm, half of the width of the joint
        self.w=12/self.realrate
        self.tj=0.5/self.realrate
        self.fw=5.5/self.realrate#flame width of the surface, if used, 
        #AND,! if use planebuild, fw=t
        self.hC=(5/self.realrate,0.75/self.realrate)#the center of a hole from the midline of shafts and side to the edges
        self.hr=0.3/self.realrate#the radius of the hole
class ExtrudeBuilder:
    def __init__(self,GV:GlobalParameter) -> None:
        self.GV=GV     
    def _getExclusion(self,a:list[np.matrix],n:list[float])->vtk.vtkAlgorithmOutput:
        points=vtk.vtkPoints()
        numpoints=len(a)
        for i in a:
            x=np.array(i)
            x=x.flatten()
            points.InsertNextPoint(x)

        polys=vtk.vtkCellArray()
        polys.InsertNextCell(numpoints)
        for i in range(numpoints):
            polys.InsertCellPoint(i)
        polygon=vtk.vtkPolyData()
        polygon.SetPoints(points)
        polygon.SetPolys(polys)
        extrude=vtk.vtkLinearExtrusionFilter()
        extrude.SetInputData(polygon)
        extrude.SetVector(n)
        extrude.Update()

        #normals to get the normal vector of eachpoints of this tri-prism
        normals=vtk.vtkPolyDataNormals()
        normals.SetInputConnection(extrude.GetOutputPort())
        #if the angle larger than 80 deg, use laplace method to smooth it.
        normals.SetFeatureAngle(80)
        mapper=vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(normals.GetOutputPort())
        prism=vtk.vtkActor()
        prism.SetMapper(mapper)
        prism.GetProperty().SetColor(0.8,0.5,0.7)

        #and then build the joint of the shaft
        return prism
    def refix(self,v0:list[np.matrix],sur:Surface,width:float)->list[np.matrix]:
        v=copy.deepcopy(v0)
        for i in range(3):
            a,b,bis=sur.getAngleBisector(sur.v[i])
            sintheta=sqrt(1-((a*bis.T)[0,0])**2)
            len=width/sintheta
            bis=bis*len
            v[i]=v[i]+bis
        return v
    #the Flame is a triangle flame of the surface with width GV.fw, and can be updated
    def getFlame(self,a:list[np.matrix],b:list[np.matrix],sur:Surface)->tuple[vtk.vtkAlgorithmOutput,vtk.vtkAlgorithmOutput]:
        n0=-sur.uppersurface()*self.GV.t
        n=[n0[0,0],n0[0,1],n0[0,2]]
        v1=[a[0],a[1],b[1],b[0]]
        v2=[a[0],a[2],a[1],b[1],b[2],b[0]]
        prism1=self._getExclusion(v1,n)
        prism2=self._getExclusion(v2,n)

        return prism1,prism2
    def getTrianglePrism(self,v:list[np.matrix],sur:Surface)->vtk.vtkActor:
        n0=-sur.uppersurface()*self.GV.t
        n=[n0[0,0],n0[0,1],n0[0,2]]
        prism=self._getExclusion(v,n)
        return prism
    #to get the cross-section of the edge, with it's real position, four points in counterclockwise in z+.
    def getShaftCrossSection(self,e:Edge,G:np.matrix,height:float,w:float=-1)->list[np.matrix]:
        if(w==-1):
            w=self.GV.w
        p=[None]*4
        u=np.mat(e.u.x)
        v=np.mat(e.v.x)
        a=e.get_vector(e.u)
        a=a*w/2/math.sqrt(a*a.T)
        mid=np.mat([(e.u.x[i]+e.v.x[i])/2 for i in range(3)])
        mid[0,2]=height
        uG=G-u
        uv=v-u
        proj=(uG*uv.T)[0,0]#dot product, but infact deal as projection latter
        inner_vertical:np.matrix=uG-uv*proj/((uv*uv.T)[0,0])
        inner_vertical=inner_vertical/math.sqrt(inner_vertical*inner_vertical.T)
        mid=mid+inner_vertical*self.GV.epsil
        #p[0]->p[1]have the same direction as e
        p[0]=mid-a
        p[1]=mid+a

        mid=mid+inner_vertical*self.GV.t
        p[2]=mid+a
        p[3]=mid-a
        return p
    def getJointCrossSection(self,e:Edge,sur:Surface)->list[np.matrix]:
        p=[None]*4
        G=sur.getGravity()
        u=np.mat(e.u.x)
        v=np.mat(e.v.x)
        a=e.get_vector(e.u)
        a=a*self.GV.w/2/math.sqrt(a*a.T)
        mid=np.mat([(e.u.x[i]+e.v.x[i])/2 for i in range(3)])
        uG=G-u
        uv=v-u
        proj=(uG*uv.T)[0,0]#dot product, but infact deal as projection latter
        inner_vertical:np.matrix=uG-uv*proj/((uv*uv.T)[0,0])
        inner_vertical=inner_vertical/math.sqrt(inner_vertical*inner_vertical.T)
        if(e.type==1):
            mid[0,2]=self.GV.z_off[sur.id]
        else:
            mid[0,2]=self.GV.z_off[self.GV.S+e.id]
        p[0]=mid+a
        p[1]=mid-a
        mid=mid+inner_vertical*self.GV.epsil
        p[2]=mid-a
        p[3]=mid+a
        return p
#for simplifying, we ignore epsilon of the joint when determine the location of the shaft.


    def getShaftPrism(self,e:Edge,sur:Surface)->vtk.vtkActor:
        G=sur.getGravity()
        #to build the shaft
        a=self.getShaftCrossSection(e,G,self.GV.z_off[sur.id]-self.GV.t)
        exc_dir=-sur.uppersurface()
        n0=exc_dir*(self.GV.z_off[sur.id]-self.GV.z_off[self.GV.S+e.id]-self.GV.t)
        n=[n0[0,0],n0[0,1],n0[0,2]]
        prism=self._getExclusion(a,n)
        prism.GetProperty().SetColor(0,0.5,0.7)#may delete later
        return prism
    ###! lofter is to drop each surfaces to the same zero plant for a printing 
    def getShaftPrism_lofter(self,e:Edge,sur:Surface)->vtk.vtkActor:
        G=sur.getGravity()
        #to build the shaft
        a=self.getShaftCrossSection(e,G,-self.GV.t)
        exc_dir=-sur.uppersurface()
        n0=exc_dir*(self.GV.z_off[sur.id]-self.GV.z_off[self.GV.S+e.id]-self.GV.t)
        n=[n0[0,0],n0[0,1],n0[0,2]]
        prism=self._getExclusion(a,n)
        prism.GetProperty().SetColor(0,0.5,0.7)#may delete later
        return prism
    def getJointPrism(self,e:Edge,sur:Surface)->vtk.vtkActor:
        a=self.getJointCrossSection(e,sur)
        if(e.type==1):
            exc_dir=-sur.uppersurface()
        else:
            exc_dir=sur.uppersurface()
        n0=exc_dir*(self.GV.tj)
        n=[n0[0,0],n0[0,1],n0[0,2]]
        prism=self._getExclusion(a,n)
        prism.GetProperty().SetColor(0,0.5,0.7)#may delete later
        return prism
class PlaneBuilder:
    def __init__(self,GV:GlobalParameter):
        self.GV=GV
        self._calc=ExtrudeBuilder(GV)
    def refix(self,e:Edge,sur:Surface,width:float)->tuple[np.matrix,np.matrix]:
        a,b,bis=sur.getAngleBisector(e.u)
        sintheta=sqrt(1-((a*bis.T)[0,0])**2)
        len=width/sintheta
        bis=bis*len
        u=np.mat(e.u.x)
        u=u+bis
        a,b,bis=sur.getAngleBisector(e.v)
        sintheta=sqrt(1-((a*bis.T)[0,0])**2)
        len=width/sintheta
        bis=bis*len
        v=np.mat(e.v.x)
        v=v+bis
        return (u,v)
    def getUnit(self,sur:Surface,lofter:bool=False)->vtk.vtkCleanPolyData:
        #infact, this poly data is better to input from a file, but for simplification, write directly here.
        polygonset:list[vtk.vtkPolyData]=[None,None,None]
        i=0
        for e in sur.e:

            points=vtk.vtkPoints()
             #vertex of u,v, refixed, outer and inner #vertex of u
            self.getSurfacePoly(points,e,sur,lofter)
            self.getShaftPoly(points,e,sur,lofter)
            self.getHolePoly(points,e,sur,lofter)
            
            polys=vtk.vtkCellArray()
            #upperSurface
            polys.InsertNextCell(4)
            a=[0,1,2,3]
            self._insertcellpoints(polys,a)
            #lowerSurface
            polys.InsertNextCell(4)
            a=[4,8,11,7]
            self._insertcellpoints(polys,a)
            polys.InsertNextCell(4)
            a=[9,5,6,10]
            self._insertcellpoints(polys,a)
            #SideofShaft
            polys.InsertNextCell(4)
            a=[8,12,15,11]
            self._insertcellpoints(polys,a)
            polys.InsertNextCell(4)
            a=[9,13,14,10]
            self._insertcellpoints(polys,a)
            #lowerShaft
            polys.InsertNextCell(4)
            a=[12,13,14,15]
            self._insertcellpoints(polys,a)
            polys.InsertNextCell(4)
            #outer shaft side
            a=[8,24,28,12]
            self._insertcellpoints(polys,a)
            polys.InsertNextCell(4)
            a=[9,25,29,13]
            self._insertcellpoints(polys,a)
            polys.InsertNextCell(6)
            a=[12,28,20,21,29,13]
            self._insertcellpoints(polys,a)
            polys.InsertNextCell(4)
            a=[16,17,21,20]
            self._insertcellpoints(polys,a)
            polys.InsertNextCell(10)
            a=[0,1,5,9,25,17,16,24,8,4]
            self._insertcellpoints(polys,a)          
            #inner shaft side
            polys.InsertNextCell(4)
            a=[11,27,31,15]
            self._insertcellpoints(polys,a)
            polys.InsertNextCell(4)
            a=[10,26,30,14]
            self._insertcellpoints(polys,a)
            polys.InsertNextCell(6)
            a=[15,14,30,22,23,31]
            self._insertcellpoints(polys,a)
            polys.InsertNextCell(4)
            a=[19,18,22,23]
            self._insertcellpoints(polys,a)
            polys.InsertNextCell(10)
            a=[2,3,7,11,27,19,18,26,10,6]
            self._insertcellpoints(polys,a) 
            #holeSide1
            polys.InsertNextCell(4)
            a=[16,19,23,20]
            self._insertcellpoints(polys,a)
            polys.InsertNextCell(4)
            a=[16,24,27,19]
            self._insertcellpoints(polys,a)
            polys.InsertNextCell(4)
            a=[24,27,31,28]
            self._insertcellpoints(polys,a)
            polys.InsertNextCell(4)
            a=[28,20,23,31]
            self._insertcellpoints(polys,a)
            #holeSide2
            polys.InsertNextCell(4)
            a=[17,18,22,21]
            self._insertcellpoints(polys,a)
            polys.InsertNextCell(4)
            a=[17,25,26,18]
            self._insertcellpoints(polys,a)
            polys.InsertNextCell(4)
            a=[25,26,30,29]
            self._insertcellpoints(polys,a)
            polys.InsertNextCell(4)
            a=[21,29,30,22]
            self._insertcellpoints(polys,a)          
            polygonset[i]=vtk.vtkPolyData()
            polygonset[i].SetPoints(points)
            polygonset[i].SetPolys(polys)
            i+=1
        #to get the flame polys with out the hole
        #id of the flame work is outer at 0 inner at 1 of upper Surface and 2, 3 of lower. mod4
        appendfilter=vtk.vtkAppendPolyData()
        appendfilter.AddInputData(polygonset[0])
        appendfilter.AddInputData(polygonset[1])
        appendfilter.AddInputData(polygonset[2])
        cleanfileter=vtk.vtkCleanPolyData()
        cleanfileter.SetInputConnection(appendfilter.GetOutputPort())
        triangle=vtk.vtkTriangleFilter()
        triangle.SetInputConnection(cleanfileter.GetOutputPort())
        return triangle

    def getSurfacePoly(self,p:vtk.vtkPoints,e:Edge,sur:Surface,lofter:bool=False)->None:
        if(lofter==True):
            transform=self.GV.z_off[sur.id]
        else:
            transform=0
        uo,vo=self.refix(e,sur,self.GV.epsil)    #vertex of u, refixed, outer and inner
        ui,vi=self.refix(e,sur,self.GV.epsil+self.GV.fw)   #vertex of u, refixed, inner
        uo[0,2]=self.GV.z_off[sur.id]-transform
        ui[0,2]=self.GV.z_off[sur.id]-transform
        vi[0,2]=self.GV.z_off[sur.id]-transform
        vo[0,2]=self.GV.z_off[sur.id]-transform
        print(uo)
        x=np.array(uo).flatten()
        p.InsertNextPoint(x)
        x=np.array(vo).flatten()
        p.InsertNextPoint(x)
        x=np.array(vi).flatten()
        p.InsertNextPoint(x)
        x=np.array(ui).flatten()
        p.InsertNextPoint(x)
        n=-sur.uppersurface()*self.GV.t
        (uo,vo,ui,vi)=(uo+n,vo+n,ui+n,vi+n)
        print(uo)
        x=np.array(uo).flatten()
        p.InsertNextPoint(x)
        x=np.array(vo).flatten()
        p.InsertNextPoint(x)
        x=np.array(vi).flatten()
        p.InsertNextPoint(x)
        x=np.array(ui).flatten()
        p.InsertNextPoint(x)
    def getShaftPoly(self,p:vtk.vtkPoints,e:Edge,sur:Surface,lofter:bool=False)->None:
        if(lofter==True):
            transform=self.GV.z_off[sur.id]
        else:
            transform=0
        G=sur.getGravity()
        height=self.GV.z_off[sur.id]-self.GV.t-transform
        a=self._calc.getShaftCrossSection(e,G,height)

        #for a simplized consideration, type 1 edges need to have a shaft. 
        if(e.type==1 or e.sur2==-1):
            length=self.GV.mh-self.GV.t
        elif(e.type==2):
            length=self.GV.z_off[sur.id]-self.GV.z_off[self.GV.S+e.id]-self.GV.t
        exc_dir=-sur.uppersurface()
        n=exc_dir*length
        b=[a[i]+n for i in range(4)]
        if(b[1][0,2]>0):
            print("err")
        for i in range(4):
            x=np.array(a[i]).flatten()
            p.InsertNextPoint(x)
        for i in range(4):
            x=np.array(b[i]).flatten()
            p.InsertNextPoint(x)
        
    #infact, it only cares points that reflect to the midline for a extrude, so have a strange order, from near to to far to the mid line in counterclockwise
    #if lofter = True, then the whole modle at the same plane as the buttom
    def getHolePoly(self,p:vtk.vtkPoints,e:Edge,sur:Surface,lofter:bool=False)->None:
        if(lofter==True):
            transform=self.GV.z_off[sur.id]
        else:
            transform=0
        G=sur.getGravity()
        if(e.type==1 or e.sur2==-1):
            height=self.GV.z_off[sur.id]-self.GV.hC[1]+self.GV.hr-transform
        elif(e.type==2):
            height=self.GV.z_off[self.GV.S+e.id]+self.GV.hC[1]+self.GV.hr-transform
        length=self.GV.hr*2
        exc_dir=-sur.uppersurface()
        n=exc_dir*length
        a=self._calc.getShaftCrossSection(e,G,height,2*(self.GV.hC[0]-self.GV.hr))
        print(a)
        b=[a[i]+n for i in range(4)]
        print(b)
        for i in range(4):
            x=np.array(a[i]).flatten()
            p.InsertNextPoint(x)
        for i in range(4):
            x=np.array(b[i]).flatten()
            p.InsertNextPoint(x)
        a=self._calc.getShaftCrossSection(e,G,height,2*(self.GV.hC[0]+self.GV.hr))
        b=[a[i]+n for i in range(4)]
        for i in range(4):
            x=np.array(a[i]).flatten()
            p.InsertNextPoint(x)
        for i in range(4):
            x=np.array(b[i]).flatten()
            p.InsertNextPoint(x)
    def _insertcellpoints(self,polys:vtk.vtkCellArray,a:list[float]):
        for i in a:
            polys.InsertCellPoint(i)
 