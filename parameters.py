class Parameters(object):
    """ module to import which holds all the parameters for the 
    simulation, hopefully only references but i will see """
    def __init__(self):
        from numpy import pi,cos
        self.q = 1.6e-19 #Coulomb
        self.hbar = 6.58211928e-16 #eV * s
        self.c = 299792458 #m/s
        self.m0 = 0.510e6/(self.c**2) #eV*s**2/m**2
        self.eps0 = 8.854e-12 # Vacuum permittivity C/V*m
        self.epsr = 12.85
        self.kb = 8.6173324e-5 # ev /K
        self.a = 2e-9 # in meter
        self.alpha = 20e-12 # eV m
        #effective mass in eV real in GaAs 0.063
        self.mass = 0.026*self.m0
        # self.mass = 0.036*self.m0
        self.t0 = (self.hbar**2)/(2*self.mass*(self.a**2))
        self.tso = (self.alpha*1.0)/(2 * self.a)
        #tso = 0.01*t0
        #tso = 0
        #Temperature * k_boltzmann in eV, 0.0025ev~30K
        self.Temp = 2 #in Kelvin
        self.kT = self.kb * self.Temp
        # self.lambdaf = 35 # i believe in nanometer, 35 more realistic?
        self.lambdaf = 10 # i believe in nanometer, 35 more realistic?
        self.BField = 0 # in Tesla, from 0-~10
        self.zplus = 1j*1e-12
        self.band_bottom = 0
        #band_bottom = -4*t0
        self.Efermi = self.band_bottom + 2*self.t0*(1-cos(2*pi/self.lambdaf))
        #Efermi = -3.8 * t0 # close to the bottom of the band at -4.0 t0, what bottom and band in what material ?
        #self.Efermi = 0.15*self.t0
        #Potential Drop over legth of device
        self.potential_drop = [0,0]
        #electro-chemical potential in eV
        self.mu_l = self.Efermi - (self.potential_drop[1] - self.potential_drop[0])/2
        self.mu_r = self.Efermi + (self.potential_drop[1] - self.potential_drop[0])/2
        #potential_drop = [0.004*t0/2, -0.004* t0/2]# in eV
        #Egrid = linspace(Efermi-0.4*t0,Efermi +0.4*t0,100)+zplus # in eV ?
        self.mu = self.Efermi
        #dE = Egrid[1].real-Egrid[0].real
        self.raw_coords = None
        self.tuple_canvas_coordinates = None
        self.canvas = None
        self.contacts = None

    def initialize_grid(self):
        from scipy import ogrid
        size_x = self.canvas.shape[0]
        size_y = self.canvas.shape[1]
        x,y = ogrid[0:size_x:size_x*1j,0:size_y:size_y*1j]
        return x,y

    def empty_potential_grid(self):
        from numpy import zeros
        self.potential_grid = zeros(self.canvas.shape)

    def stepgrid(self, step1,step2):
        from scipy import r_
        self.potential_grid = (r_[[self.potential_drop[0]]*self.canvas.shape[1]*step1,
            [0]*self.canvas.shape[1]*(self.canvas.shape[0]-step1-step2),
            [self.potential_drop[1]]*self.canvas.shape[1]*step2].reshape(self.canvas.shape))

    def circular_qpc(self,shift=0,radius=1,scale=0):
        import aux
        x,y = self.initialize_grid()
        self.potential_grid = aux.sphericalPot(x,y,shift,radius,scale)

    def pointcharge_qpc(self,charge=0,scale=0):
        import aux
        x,y = self.initialize_grid()
        self.potential_grid = aux.pointchargePot(x,y,charge,scale)

    def linearsmooth_qpc(self,slope=1,scale=0,xi=1):
        import aux
        x,y = self.initialize_grid()
        self.potential_grid = aux.linearsmoothPot(x,y,slope,scale,xi)

    def triangular_qpc(self,shift=0,width=1,radius=1,scale=0):
        import aux
        x,y = self.initialize_grid()
        self.potential_grid = aux.triangularPot(x,y,shift,width,radius,scale)

    def rectangular_qpc(self,shift=0,width=1,scale=0):
        import aux
        x,y = self.initialize_grid()
        self.potential_grid = aux.rectangularPot(x,y,shift,width,scale)
    def __generate_potential_grid(self):
        from scipy import linspace, tile
        rowdim = self.canvas.shape[0]
        coldim = self.canvas.shape[1]
        pot_slice = linspace(p.potential_drop[0],p.potential_drop[1],rowdim)
        potential = tile(pot_slice,(coldim,1)).T
        self.potential_grid = potential

    def grid2serialized(self, grid):
        from scipy import concatenate, array
        potential_serialized = array([])
        for key, value in self.nodes.items():
           potential_serialized = concatenate((potential_serialized, [grid[key[0],key[1]]]))
        self.potential_serialized = potential_serialized

p = Parameters()
