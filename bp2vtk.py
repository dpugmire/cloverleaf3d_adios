import adios2
import numpy as np
import vtk, sys, math
from vtk.util.numpy_support import numpy_to_vtk

def readArray(f, nm, nameIt=False, blockId=-1) :
    #print('read: ', nm, blockId)

    adiosVar = f.inquire_variable(nm)
    if blockId >= 0 :
        adiosVar.set_block_selection(blockId)
    adiosVar = f.read(nm)

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

    ghostZoneVar.SetName('vtkGhostCells')
    grid.GetCellData().AddArray(ghostZoneVar)

    for var in cellVars :
        grid.GetCellData().AddArray(var)
    for var in ptVars :
        grid.GetPointData().AddArray(var)
    return grid

def getGridFileName(step, blockId) :
    fname = 'vtk/grid_%03d.%03d.vtk' % (step, blockId)
    return fname

def dumpGrid(grid, step, blockId) :
    fname = getGridFileName(step, blockId)
    writer = vtk.vtkDataSetWriter()
    writer.SetFileVersion(42)
    writer.SetFileName(fname)
    writer.SetInputData(grid)
    writer.Write()

def createVisitFile(numSteps, numBlocks) :
    visitFileName = 'grid.visit'
    f = open(visitFileName, 'w')
    f.write('!NBLOCKS %d\n' % numBlocks)
    for step in range(numSteps) :
        for b in range(numBlocks) :
            gridFileName = getGridFileName(step, b)
            f.write('%s\n' % gridFileName)


xcoords, ycoords, zcoords = (None,None,None)

f = adios2.Stream('./output.bp', 'r')
numSteps = f.num_steps()
numBlocks = len(f.all_blocks_info('density')[0])
print('numBlocks= ', numBlocks)

xcoords, ycoords, zcoords, ghostZones = ([], [], [], [])

for step in range(numSteps) :
    f.begin_step()
    print('step= ', step)

    #first step, read in coords and ghost zones.
    if (step == 0) :
        createVisitFile(numSteps, numBlocks)

        for bi in range(numBlocks) :
            xcoords.append(readArray(f, 'coordsX', False, bi))
            ycoords.append(readArray(f, 'coordsY', False, bi))
            zcoords.append(readArray(f, 'coordsZ', False, bi))
            ghostZones.append(readArray(f, 'ghost_zones', True, bi))


    for bi in range(numBlocks) :
        nx = xcoords[bi].GetNumberOfValues()
        ny = ycoords[bi].GetNumberOfValues()
        nz = zcoords[bi].GetNumberOfValues()

        varDensity = readArray(f, 'density', True, bi)
        varEnergy = readArray(f, 'energy', True, bi)
        varPressure = readArray(f, 'pressure', True, bi)
        varVelX = readArray(f, 'velocityX', True, bi)
        varVelY = readArray(f, 'velocityY', True, bi)
        varVelZ = readArray(f, 'velocityZ', True, bi)

        grid = makeGrid((xcoords[bi], ycoords[bi], zcoords[bi]), (nx,ny,nz), (varDensity, varEnergy, varPressure), (varVelX, varVelY, varVelZ), ghostZones[bi])
        dumpGrid(grid, step, bi)
    f.end_step()
