class Model:

    def __init__(self, world):
        from scipy import linspace
        Model.canvas = world.canvas
        Model.atlas = world.atlas
        Model.hbar = world.hbar
        Model.m0 = world.m0

        self.nodes = world.nodes
        self.contacts = world.contacts
        self.active_coords = world.active_coords
        self.block_sizes = world.block_sizes
        #Potential Drop over legth of device
        self.potential_drop = [-0.01, 0.01]
        self.a = 2e-10
        #effective mass in eV
        self.mass = 0.067*Model.m0
        self.t = (Model.hbar**2)/(2*self.mass*(self.a**2))
        #Temperature * k_boltzmann in eV, 0.0025ev~30K
        self.kT = 0.025
        self.lambdaf = 10
        self.BField = 0
        self.zplus = 1j*1e-12
        self.Egrid = linspace(0.01,0.2,500)+self.zplus
        #Efermi = 2*t*(scipy.cos(scipy.pi/lambdaf))
        self.Efermi = 0.1
        self.mu = self.Efermi
        self.dE = self.Egrid[1].real-self.Egrid[0].real
        #electro-chemical potential in eV
        self.mu_l = self.Efermi + (self.potential_drop[1] - self.potential_drop[0])/2
        self.mu_r = self.Efermi - (self.potential_drop[1] - self.potential_drop[0])/2

        self.__generate_potential_grid()
        self.grid2serialized(self.potential_grid)
        self.__build_H()


    def writetovtk(self,value,suffix=''):
        xdim = Model.canvas[1]
        ydim = Model.canvas[0]
        header = []
        header.append("# vtk DataFile Version 2.0\nVTK Data of Device\nASCII\nDATASET STRUCTURED_POINTS\n")
        header.append("DIMENSIONS {0} {1} 1".format(xdim, ydim))
        header.append("\nSPACING 1 1 1\nORIGIN 0 0 0\nPOINT_DATA {0}".format(xdim * ydim))
        header.append("\nSCALARS EDensity double 1\nLOOKUP_TABLE default\n")
        with open(Model.atlas + suffix + '.vtk', 'w') as file:
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

    def __build_H(self):
        from scipy.sparse import lil_matrix
        from scipy import complex128, exp
        H = lil_matrix((len(self.nodes),len(self.nodes)),dtype=complex128)
        for item in self.nodes:
            if self.nodes[item][1] != None:
                H[self.nodes[item][0],self.nodes[item][1]] = -exp(1j*self.BField*self.nodes[item][0])*self.t
            else: pass
            if self.nodes[item][2] == None: continue
            H[self.nodes[item][0],self.nodes[item][2]] = -self.t
        H = (H.tocsr() + H.tocsr().conjugate().T).tolil()# might be faster
        H.setdiag(self.potential_serialized+4*self.t+self.zplus)
        self.H = H

    def fermifunction(self, E, mu):
        from scipy import exp
        """ Simple Fermifunction """
        fermifnc = 1/(exp((E.real-mu)/self.kT)+1)
        return fermifnc

    def __contact_greensfunction(self, ind_contact, contact_node, E):
        from scipy import arccos, exp, sin, pi
        """calculates the value of the transverse mode's square at the point
        of the contact_node, essentially sin(pi)sin(pi) multiplied with
        appropriate amplitude and constats """
        Length = (ind_contact.shape[1]-1.0)
        #Amplitude = scipy.sqrt(1/Length)
        Amplitude = 1
        ka = arccos(1-E.real/(2*self.t))
        Phase = exp(1j*ka)
        xi = Amplitude * sin(pi * contact_node/Length)
        #xi = 1
        greensfunction = (xi**2) * Phase
        return greensfunction

    def build_sigma(self, ind_contact, E):
        from scipy.sparse import lil_matrix
        from scipy import complex128
        """Build the self-energy matrix SIGMA by determining the nodes
        adjacent to a Contact and inserting the greensfunction times t**2
        (the greensfunction comes with 1/t)"""
        sigma = lil_matrix((len(self.nodes), len(self.nodes)), dtype=complex128)
        for contact_node in range(ind_contact.shape[1]):
            index = self.nodes[tuple(ind_contact.T[contact_node])][0]
            sigma[index, index] = - self.t * self.__contact_greensfunction(ind_contact, contact_node, E)
        return sigma

    def build_A(self, E):
        from scipy.sparse import eye
        from scipy import complex128
        number_of_nodes = len(self.nodes)
        sigma_l = self.build_sigma(self.contacts[0],E - self.potential_drop[0])
        sigma_r =self.build_sigma(self.contacts[1], E - self.potential_drop[1])
        sigma_in_l = -2* sigma_l.imag[0:self.block_sizes[0],0:self.block_sizes[0]] * self.fermifunction(E, mu=self.mu_l)
        sigma_in_r = -2* sigma_r.imag[-self.block_sizes[-1]:,-self.block_sizes[-1]:] * self.fermifunction(E, mu=self.mu_r)
        A = E*eye(number_of_nodes,number_of_nodes,dtype=complex128, format='lil') - self.H - sigma_l  - sigma_r
        return A, sigma_in_l, sigma_in_r

    def RRGM(self,Ablock):
        """ Performs recursive algorithm (Svizhenko et. al) to calculate
        the retarded green's function, uses views on A, i.e. the block
        matrices of A, by Ablock """
        from scipy import array, hstack
        number_of_blocks = len(self.block_sizes)
        grl = [Ablock[0,0].todense().I]
        prev_greensfnc = grl[0]
        for i in range(1,number_of_blocks):
            prev_greensfnc = (Ablock[i,i]-Ablock[i, i-1] * prev_greensfnc * Ablock[i-1,i]).I
            grl.append(prev_greensfnc)
        Gr = {}
        Gr[number_of_blocks-1,number_of_blocks-1] = grl[-1]
        rev_iterator = reversed(range(1,number_of_blocks))
        for i in rev_iterator:
            Gr[i, i-1] = -Gr[i,i] * Ablock[i,i-1] * grl[i-1]
            Gr[i-1, i] = -grl[i-1] * Ablock[i-1,i] * Gr[i,i]
            Gr[i-1, i-1] = grl[i-1]-grl[i-1] * Ablock[i-1,i] * Gr[i, i-1]
        green_diag = array([])
        for i in range(number_of_blocks):
            diag = array(Gr[i,i].diagonal()).reshape(-1)
            green_diag = hstack((green_diag, diag))
        return green_diag, grl, Gr

    def LRGM(self, Ablock, sigma_in_l, sigma_in_r):
        """ Performs recursive algorithm (Svizhenko et.al) to calculate
        lesser green's function by the use of the diagonal and offdiagonal
        matrix-elements of the retarded green's function """
        from scipy import array, hstack
        ignored, grl, Gr = self.RRGM(Ablock)
        number_of_blocks = len(self.block_sizes)
        gll = [grl[0] * sigma_in_l * grl[0].getH()]
        #len(block)-1 equals N-2 because of pythons way of building the range exclusive the endpoint
        for i in range(1,number_of_blocks-1):
            prev_lesser_greenfnc = grl[i] * (Ablock[i,i-1] * gll[i-1] * Ablock[i,i-1].getH()) * grl[i].getH()
            gll.append(prev_lesser_greenfnc)
        #N-1 step extra, belongs to for loop
        i +=1
        gll.append(grl[i] * (sigma_in_r + Ablock[i,i-1] * gll[i-1] * Ablock[i-1,i].conj()) * grl[i].getH())
        Gl = {}
        Gl[number_of_blocks-1,number_of_blocks-1] = gll[-1]
        for i in reversed(range(1,number_of_blocks)):
            Gl[i,i-1] = -Gr[i,i] * Ablock[i,i-1] * gll[i-1] - Gl[i,i] * Ablock[i-1,i].getH() * grl[i-1].getH()
            Gl[i-1,i-1] = (gll[i-1] + grl[i-1] * (Ablock[i-1,i] * Gl[i,i] * Ablock[i-1,i].getH()) * grl[i-1].getH() - (gll[i-1] * Ablock[i,i-1].getH() * Gr[i-1,i].getH() + Gr[i-1,i] * Ablock[i,i-1] * gll[i-1]))
        less_green_diag = array([])
        for i in range(number_of_blocks):
            less_diag = array(Gl[i,i].diagonal()).reshape(-1)
            less_green_diag = hstack((less_green_diag, less_diag))
        return less_green_diag

    def energyintegrate(self,integrand,sigma_in_l=None,sigma_in_r=None):
        from scipy import pi, array
        from scipy.sparse import lil_matrix
        from sparseblockslice import SparseBlocks
        #if integrand == self.RRGM:
        #    value, foo = -integrand.__call__(Ablock).imag
        #elif integrand == self.LRGM and sigma_in_l is not None and sigma_in_r is not None:
        #    value = integrand.__call__(Ablock,sigma_in_l, sigma_in_r,).real
        #else:
            #print 'Please insert supported functions RRGM or LRGM(not none sigmas)'
        integral = array([0]*len(self.nodes))
        #A = lil_matrix((len(self.nodes), len(self.nodes)))
        #print A
        max_density = []
        i=0
        for energy in self.Egrid:
            print energy
            A, sigma_in_l, sigma_in_r = self.build_A(energy)
            Ablock = SparseBlocks(A, self.block_sizes)
            integral = integrand.__call__(Ablock, sigma_in_l, sigma_in_r)*self.dE/(pi*self.a)
            print integral
            self.writetovtk(integral.real, str(i))
            i+=1
            max_density.append(integral.real.max())
        return integral, max_density

