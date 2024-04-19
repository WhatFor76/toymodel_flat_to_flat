'''
Copyright @WhatFor, contact with 332449502@qq.com

'''
import vtk
from Unit import *
from math import sqrt
from ModelBuild import *

GV=GlobalParameter()

def main():
    builder=PlaneBuilder(GV)# only used to calculate the geometry that used, not to build the flamework
    colors = vtk.vtkNamedColors()
    background = map(lambda x: x / 255.0, [26, 51, 102, 255])
    colors.SetColor("BkgColor", *background)
    renderer=vtk.vtkRenderer()
    renderer.SetBackground(colors.GetColor3d("BkgColor"))
    ## ATTENTION, all of the model must be not crossed
    for sur in GV.Surfaces:
        unitpoly=builder.getUnit(sur,lofter=False)
        normals=vtk.vtkPolyDataNormals()
        normals.SetInputConnection(unitpoly.GetOutputPort())
        #if the angle larger than 80 deg, use laplace method to smooth it.
        normals.SetFeatureAngle(80)
        mapper=vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(normals.GetOutputPort())
        unit=vtk.vtkActor()
        unit.SetMapper(mapper)
        unit.GetProperty().SetColor(0.8,0.5,0.7)
        renderer.AddActor(unit)
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
