import adios2
import numpy as np
import vtk, sys, math
from vtk.util.numpy_support import numpy_to_vtk

def readArray(f, nm, nameIt=False, blockId=-1) :
    if blockId < 0 :
        adiosVar = f.read(nm)
    else :
        adiosVar = f.read(nm, blockId=blockId)

    if adiosVar.dtype == np.float64 :
        adiosVar = adiosVar.astype(np.float32)

    if len(adiosVar) > 1 :
        adiosVar = np.ravel(adiosVar)
    #print(nm, adiosVar.shape, adiosVar.dtype, adiosVar.min(), adiosVar.max())
    adiosVar[np.abs(adiosVar) < 1e-8] = 0
    #print(nm, adiosVar.shape, adiosVar.dtype, adiosVar.min(), adiosVar.max())
    #print(nm, adiosVar.shape)
    arr = vtk.util.numpy_support.numpy_to_vtk(adiosVar, deep=True, array_type=vtk.VTK_FLOAT)
    if nameIt : arr.SetName(nm)

    return arr

def makeGrid(coords, dims, cellVars, ptVars, ghostZoneVar) :
    grid = vtk.vtkRectilinearGrid()
    grid.SetDimensions(dims[0], dims[1], dims[2])
    grid.SetXCoordinates(coords[0])
    grid.SetYCoordinates(coords[1])
    grid.SetZCoordinates(coords[2])

    grid.GetCellData().AddArray(ghostZoneVar)

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
gz = f.read('ghost_zones')
nx,ny,nz = (x.shape[0], y.shape[0], z.shape[0])
xcoords = readArray(f, 'coordsX')
ycoords = readArray(f, 'coordsY')
zcoords = readArray(f, 'coordsZ')
ghostZones = readArray(f, 'ghost_zones', True)
f.close()

f = adios2.Stream('./output.bp', 'r')
numSteps = f.num_steps()
numBlocks = len(f.all_blocks_info('density')[0])
print('numBlocks= ', numBlocks)

for step in range(numSteps) :
    f.begin_step()
    print('step= ', step)
    varDensity = readArray(f, 'density', True)
    varEnergy = readArray(f, 'energy', True)
    varPressure = readArray(f, 'pressure', True)
    varVelX = readArray(f, 'velocityX', True)
    varVelY = readArray(f, 'velocityY', True)
    varVelZ = readArray(f, 'velocityZ', True)

    grid = makeGrid((xcoords, ycoords, zcoords), (nx,ny,nz), (varDensity, varEnergy, varPressure), (varVelX, varVelY, varVelZ), ghostZones)
    dumpGrid(grid, step)
    f.end_step()
