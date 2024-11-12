import adios2
import numpy as np
import vtk, sys, math
from vtk.util.numpy_support import numpy_to_vtk

def readArray(f, nm, nameIt=False) :
    adiosVar = f.read(nm)
    if len(adiosVar) > 1 :
        adiosVar = np.ravel(adiosVar)
    #print(nm, adiosVar.shape)
    arr = vtk.util.numpy_support.numpy_to_vtk(adiosVar, deep=True, array_type=vtk.VTK_FLOAT)
    if nameIt : arr.SetName(nm)

    return arr

def makeGrid(coords, dims, cellVars, ptVars) :
    grid = vtk.vtkRectilinearGrid()
    grid.SetDimensions(dims[0], dims[1], dims[2])
    grid.SetXCoordinates(coords[0])
    grid.SetYCoordinates(coords[1])
    grid.SetZCoordinates(coords[2])

    for var in cellVars :
        grid.GetCellData().AddArray(var)
    for var in ptVars :
        grid.GetPointData().AddArray(var)
    return grid

def dumpGrid(grid, step) :
    fname = 'grid.%02d.vtk' % step
    writer = vtk.vtkDataSetWriter()
    writer.SetFileVersion(42)
    writer.SetFileName(fname)
    writer.SetInputData(grid)
    writer.Write()

xcoords, ycoords, zcoords = (None,None,None)
f = adios2.FileReader('./output.bp')

x = f.read('coordsX')
y = f.read('coordsY')
z = f.read('coordsZ')
nx,ny,nz = (x.shape[0], y.shape[0], z.shape[0])
xcoords = readArray(f, 'coordsX')
ycoords = readArray(f, 'coordsY')
zcoords = readArray(f, 'coordsZ')
f.close()

with adios2.Stream('./output.bp', 'r') as f :
    numSteps = f.num_steps()
    for step in range(numSteps) :
        f.begin_step()
        print('step= ', step)
        varDensity = readArray(f, 'density', True)
        varEnergy = readArray(f, 'energy', True)
        varPressure = readArray(f, 'pressure', True)
        varVelX = readArray(f, 'velocityX', True)
        varVelY = readArray(f, 'velocityY', True)
        varVelZ = readArray(f, 'velocityZ', True)

        grid = makeGrid((xcoords, ycoords, zcoords), (nx,ny,nz), (varDensity, varEnergy, varPressure), (varVelX, varVelY, varVelZ))
        dumpGrid(grid, step)
        f.end_step()
