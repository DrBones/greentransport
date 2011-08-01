from evtk.vtk import VtkFile, VtkRectilinearGrid
import numpy as np
from numpy import flipud
# =================================
#       Helper functions      
# =================================
def __addDataToFile(vtkFile, cellData, pointData):
    # Point data
    if pointData <> None:
        keys = pointData.keys()
        vtkFile.openData("Point", scalars = keys[0])
        for key in keys:
            data = flipud(pointData[key]).real.flatten()
            vtkFile.addData(key, data)
        vtkFile.closeData("Point")

    # Cell data
    if cellData <> None:
        keys = cellData.keys()
        vtkFile.openData("Cell", scalars = keys[0])
        for key in keys:
            data = flipud(cellData[key]).real.flatten()
            vtkFile.addData(key, data)
        vtkFile.closeData("Cell")

def __appendDataToFile(vtkFile, cellData, pointData):
    # Append data to binary section
    if pointData <> None:
        keys = pointData.keys()
        for key in keys:
            data = flipud(pointData[key]).real.flatten()
            vtkFile.appendData(data)

    if cellData <> None:
        keys = cellData.keys()
        for key in keys:
            data = flipud(cellData[key]).real.flatten()
            vtkFile.appendData(data)


def writeVTK(filename,xdim, ydim, pointData=None, cellData=None):
    """
        Writes data values as a rectilinear or rectangular grid.

        PARAMETERS:
            path: name of the file without extension where data should be saved.
            cellData: dictionary containing arrays with cell centered data.
                      Keys should be the names of the data arrays.
                      Arrays must have the same dimensions in all directions and must contain
                      only scalar data.
            nodeData: dictionary containing arrays with node centered data.
                      Keys should be the names of the data arrays.
                      Arrays must have same dimension in each direction and
                      they should be equal to the dimensions of the cell data plus one and
                      must contain only scalar data.

        RETURNS:
            Full path to saved file.

    """
    nx, ny, nz = xdim, ydim, 0
    lx, ly, lz = xdim/10.0, ydim/10.0, 0
    dx, dy, dz = lx/nx, ly/ny, 0
    ncells = nx * ny
    npoints = (nx + 1) * (ny + 1) * (nz + 1)
    x = np.arange(0, lx + 0.1*dx, dx, dtype='float64')
    y = np.arange(0, ly + 0.1*dy, dy, dtype='float64')
    z = np.arange(0,0, dtype='float64')
    start, end = (0,0,0), (nx, ny, nz)
# Set up object
    w = VtkFile(filename, VtkRectilinearGrid)
#Open XML tags
    w.openGrid(start = start, end = end)
    w.openPiece( start = start, end = end)

# Coordinates of cell vertices
    w.openElement("Coordinates")
    w.addData("x_coordinates", x);
    w.addData("y_coordinates", y);
    w.addData("z_coordinates", z);
    w.closeElement("Coordinates");

# Add data from the dictionary
    #temp = np.random.rand(npoints)
    __addDataToFile(w, cellData, pointData)

#Close XML tags
    w.closePiece()
    w.closeGrid()
#Append Coordinate Data to file in binary form
    w.appendData(x).appendData(y).appendData(z)
#Append Cell and Point Data to file in binary form
    __appendDataToFile(w, cellData, pointData)
    w.save()
    return w.getFileName()
    #w.appendData(data = field)

def legacyVTK(self,value, suffix='', mode='normal'):
    """ Input values are the serialized values of the grid, either
    custom or naive (rectangular) serialization. Also accepts
    rectangular arrays """
    if value.ndim >1:
        if mode=='normal':
            xdim = value.shape[1]
            ydim = value.shape[0]
            def linetowrite(value, row, column):
                return file.write(str(value[row,column]) + "\n")
        elif mode=='spin':
            xdim = value.shape[1]
            ydim = value.shape[0]
            def linetowrite(value, row, column):
                return file.write(str(value[row,column]) + "\n")

    else:
        ydim = self.wafer.shape[1]
        xdim = self.wafer.shape[0]
        if len(value) == xdim*ydim:
            if mode=='normal':
                value = value.reshape(xdim, ydim)
                def linetowrite(value,x, y):
                    return file.write(str(value[x,y]) + "\n")
            elif mode=='spin':
                value = value.reshape(xdim, ydim)
                def linetowrite(value,x, y):
                    return file.write(str(value[x,y]) + "\n")

        else:
            def linetowrite(value,x, y):
                return file.write(str(value[self.nodes[x_pixel,y_pixel][0]]) + "\n")
    header = []
    header.append("# vtk DataFile Version 2.0\nVTK Data of Device\nASCII\nDATASET STRUCTURED_POINTS\n")
    header.append("DIMENSIONS {0} {1} 1".format(ydim, xdim))
    header.append("\nSPACING 1 1 1\nORIGIN 0 0 0\nPOINT_DATA {0}".format(xdim * ydim))
    header.append("\nSCALARS EDensity double 1\nLOOKUP_TABLE default\n")
    with open(Model.atlas + suffix + '.vtk', 'w') as file:
        for line in header:
            file.write(line)
        for x_pixel in range(xdim):
            for y_pixel in range(ydim):
                try:
                    print x_pixel, y_pixel
                    linetowrite(value,x_pixel, y_pixel)
                except KeyError:
                    print 'Error: Index (',x_pixel,',',y_pixel,') not in array'
                    file.write('0\n')

