'''
Copyright @WhatFor, contact with 332449502@qq.com
Used in step 4. a out of date version, 
too much renderer in a renderwindow may cause some error with may lead to a 
fault in slicing on 3D printer. 


'''
import vtk
from Unit import *
from math import sqrt
from ModelBuild import *

GV=GlobalParameter()

def main():
    builder=ExtrudeBuilder(GV)
    colors = vtk.vtkNamedColors()
    background = map(lambda x: x / 255.0, [26, 51, 102, 255])
    colors.SetColor("BkgColor", *background)
    renderer=vtk.vtkRenderer()
    renderer.SetBackground(colors.GetColor3d("BkgColor"))
    ## ATTENTION, all of the model must be not crossed
    for sur in GV.Surfaces:
        v=[None,None,None]
        for i in range(3):
            v[i]:list[np.matrix]=np.mat(sur.v[i].x)
            v[i][0,2]=GV.z_off[sur.id]
        v_outer=builder.refix(v,sur,GV.epsil)    
        v_inner=builder.refix(v,sur,GV.epsil+GV.fw)
        #temperally ignore the shafts.
        prism1,prism2=builder.getFlame(v_outer,v_inner,sur)
        renderer.AddActor(prism1)
        renderer.AddActor(prism2)
        for e in sur.e:
            if(e.sur2==-1):
                continue
            if(e.type==2):
                prism=builder.getShaftPrism(e,sur)
                renderer.AddActor(prism)
            joint=builder.getJointPrism(e,sur)   
            renderer.AddActor(joint)
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(renderer)
    renWin.SetSize(300, 300)
    renWin.SetWindowName('ReschModel')
    ObjOuter=vtk.vtkOBJExporter()
    ObjOuter.SetFilePrefix("model")
    ObjOuter.SetInput(renWin)
    ObjOuter.Write()

    irenderer = vtk.vtkRenderWindowInteractor()
    irenderer.SetRenderWindow(renWin)

    irenderer.Initialize()

    renderer.ResetCamera()
    renderer.GetActiveCamera().Zoom(1.5)
    renWin.Render()
     
    irenderer.Start()
if __name__=="__main__":
    main()
