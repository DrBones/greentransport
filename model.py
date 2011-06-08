class Model:

    def __init__(self, world):
        from scipy import linspace
        Model.canvas = world.canvas
        Model.atlas = world.atlas
        Model.hbar = world.hbar
        Model.m0 = world.m0

        self.nodes = world.nodes
        self.potential_drop = [-0.01, 0.01] #Potential Drop over legth of device
        self.a = 2e-10
        self.mass = 0.25*Model.m0                                                      #effective mass in eV
        self.t = (Model.hbar**2)/(2*self.mass*(self.a**2))
        self.fkT = 0.025                                                          #Temperature * k_boltzmann in eV, 0.0025ev~30K
        self.lambdaf = 10
        self.BField = 0
        self.zplus = 1j*1e-12
        self.Egrid = linspace(0.4,1.0,100)+self.zplus
        #Efermi = 2*t*(scipy.cos(scipy.pi/lambdaf))
        self.Efermi = 0.1
        self.dE = self.Egrid[1].real-self.Egrid[0].real

        self.__generate_potential_grid()
        self.grid2serialized(self.potential_grid)


    def writetovtk(self,value):
        xdim = Model.canvas[1]
        ydim = Model.canvas[0]
        header = []
        header.append("# vtk DataFile Version 2.0\nVTK Data of Device\nASCII\nDATASET STRUCTURED_POINTS\n")
        header.append("DIMENSIONS {0} {1} 1".format(xdim, ydim))
        header.append("\nSPACING 1 1 1\nORIGIN 0 0 0\nPOINT_DATA {0}".format(xdim * ydim))
        header.append("\nSCALARS EDensity double 1\nLOOKUP_TABLE default\n")
        with open(Model.atlas + '.vtk', 'w') as file:
            for line in header:
                file.write(line)
            for x_pixel in range(ydim):
                for y_pixel in range(xdim):
                    try:
                        file.write(str(value[self.nodes[x_pixel,y_pixel][0]]) + "\n")
                    except KeyError:
                        file.write('0\n')


    def __generate_potential_grid(self):
        from scipy import vstack, hstack, r_
        xdim = Model.canvas[1]
        ydim = Model.canvas[0]
        potential = []
        for i in range(xdim):
            potential.append(vstack(r_[self.potential_drop[0]:self.potential_drop[1]:ydim*1j]))
        potential = hstack(tuple(potential))
        self.potential_grid = potential

    def grid2serialized(self, grid):
        from scipy import concatenate, array
        potential_serialized = array([])
        for key, value in self.nodes.items():
           potential_serialized = concatenate((potential_serialized, [grid[key[0],key[1]]]))
        self.potential_serialized = potential_serialized
