import vtk
from vtk.util import numpy_support
import config

DX = config.DX
DY = config.DY
NX = config.NX
NY = config.NY

def write_vti_flow(filename, u, v, p):
    imageData = vtk.vtkImageData()
    imageData.SetDimensions(NX, NY, 1)  
    imageData.SetSpacing(DX, DY, 1.0)
    imageData.SetOrigin(0.0, 0.0, 0.0)
    u_vtk = numpy_support.numpy_to_vtk(u.T.ravel(), deep=True, array_type=vtk.VTK_DOUBLE)
    v_vtk = numpy_support.numpy_to_vtk(v.T.ravel(), deep=True, array_type=vtk.VTK_DOUBLE)
    p_vtk = numpy_support.numpy_to_vtk(p.T.ravel(), deep=True, array_type=vtk.VTK_DOUBLE)
    imageData.GetPointData().AddArray(u_vtk)
    imageData.GetPointData().AddArray(v_vtk)
    imageData.GetPointData().AddArray(p_vtk)
    u_vtk.SetName("u")
    v_vtk.SetName("v")
    p_vtk.SetName("p")
    writer = vtk.vtkXMLImageDataWriter()
    writer.SetFileName(filename)
    writer.SetInputData(imageData)
    writer.Write()

def write_vti_particle(filename, particles):
    with open(filename, 'w') as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write("Data.vtk\n")
        f.write("ASCII\n")
        f.write("DATASET UNSTRUCTURED_GRID\n")
        f.write(f"POINTS {len(particles)} float\n")
        for particle in particles:
            f.write(f"{particle.x} {particle.y} 0.0\n")
        f.write(f"CELL_TYPES {len(particles)}\n")
        for _ in range(len(particles)):
            f.write("1\n")
        f.write(f"POINT_DATA {len(particles)}\n")
        f.write("SCALARS radius float\n")
        f.write("LOOKUP_TABLE default\n")
        for _ in range(len(particles)):
            f.write("0.1\n")