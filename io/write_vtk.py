from evtk.vtk import VtkFile, VtkRectilinearGrid
import numpy as np
def writeVTK(filename, field,xdim, ydim):
    """writes linear array to vtk file"""
    nx, ny, nz = xdim, ydim, 0
    lx, ly, lz = xdim/10.0, ydim/10.0, 0
    dx, dy, dz = lx/nx, ly/ny, 0
    ncells = nx * ny
    npoints = (nx + 1) * (ny + 1) * (nz + 1)
    x = np.arange(0, lx + 0.1*dx, dx, dtype='float64')
    y = np.arange(0, ly + 0.1*dy, dy, dtype='float64')
    z = np.arange(0,0, dtype='float64')
    start, end = (0,0,0), (nx, ny, nz)

    w = VtkFile(filename, VtkRectilinearGrid)
    w.openGrid(start = start, end = end)
    w.openPiece( start = start, end = end)

# Point data
    #temp = np.random.rand(npoints)
    field=field.real.flatten()
    w.openData("Point", scalars = "Density")
    w.addData("Density", field)
    w.closeData("Point")

# Coordinates of cell vertices
    w.openElement("Coordinates")
    w.addData("x_coordinates", x);
    w.addData("y_coordinates", y);
    w.addData("z_coordinates", z);
    w.closeElement("Coordinates");

    w.closePiece()
    w.closeGrid()

    w.appendData(data = field)
    w.appendData(x).appendData(y).appendData(z)
    w.save()
    return w.getFileName()
