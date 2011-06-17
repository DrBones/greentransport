class Model:

    def __init__(self, world):
        from scipy import linspace
        Model.canvas = world.canvas
        Model.atlas = world.atlas
        Model.hbar = world.hbar
        Model.m0 = world.m0

        self.nodes = world.nodes
        self.wafer = world.wafer
        self.contacts = world.contacts
        self.active_coords = world.active_coords
        self.block_sizes = world.block_sizes
        self.eps0 = world.eps0
        self.q = world.q
        #Potential Drop over legth of device
        self.potential_drop = [-0.01, 0.01]
        self.a = 2e-10
        #effective mass in eV real in GaAs 0.063
        self.mass = 0.25*Model.m0
        self.t = (Model.hbar**2)/(2*self.mass*(self.a**2))
        #Temperature * k_boltzmann in eV, 0.0025ev~30K
        self.epsr = 12.9
        self.kT = 0.025
        self.lambdaf = 10
        self.BField = 0
        self.zplus = 1j*1e-12
        self.Egrid = linspace(0.01,0.1,10)
        #Efermi = 2*t*(1-scipy.cos(2*scipy.pi/lambdaf))
        self.Efermi = 0.1
        self.mu = self.Efermi
        self.dE = self.Egrid[1].real-self.Egrid[0].real
        #electro-chemical potential in eV
        self.mu_l = self.Efermi + (self.potential_drop[1] - self.potential_drop[0])/2
        self.mu_r = self.Efermi - (self.potential_drop[1] - self.potential_drop[0])/2

        self.__generate_potential_grid()
        self.grid2serialized(self.potential_grid)
        #self.__build_H()

    def writetovtk(self,value,suffix=''):
        """ Input values are the serialized values of the grid, either
        custom or naive (rectangular) serialization """
        xdim = self.wafer.shape[1]
        ydim = self.wafer.shape[0]
        header = []
        header.append("# vtk DataFile Version 2.0\nVTK Data of Device\nASCII\nDATASET STRUCTURED_POINTS\n")
        header.append("DIMENSIONS {0} {1} 1".format(xdim, ydim))
        header.append("\nSPACING 1 1 1\nORIGIN 0 0 0\nPOINT_DATA {0}".format(xdim * ydim))
        header.append("\nSCALARS EDensity double 1\nLOOKUP_TABLE default\n")
        def linetowrite(value,x, y):
            if len(value) == xdim*ydim:
                value = value.reshape(ydim, xdim)
                return file.write(str(value[x,y]) + "\n")
            else:
                return file.write(str(value[self.nodes[x_pixel,y_pixel][0]]) + "\n")
        with open(Model.atlas + suffix + '.vtk', 'w') as file:
            for line in header:
                file.write(line)
            for x_pixel in range(ydim):
                for y_pixel in range(xdim):
                    try:
                        linetowrite(value,x_pixel, y_pixel)
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
        H.setdiag(self.potential_serialized+4*self.t)
        self.H = H

    def simpleH(self):
        from scipy.sparse import lil_matrix, triu
        from scipy import complex128, exp
        H = lil_matrix((self.wafer.shape[0]*self.wafer.shape[1],self.wafer.shape[0]*self.wafer.shape[1]),dtype=complex128)
        for row in range(self.wafer.shape[0]):
            for column in range(self.wafer.shape[1]):
                i = column+self.wafer.shape[1]*row
                if self.wafer[row,column] == 0:
                    H[i,i] = 100000
                elif self.wafer[row,column] >0:
                    H[i,i] = 4*self.t+self.potential_grid[row,column]
                    if row+1 == self.wafer.shape[0] or column+1 == self.wafer.shape[1]: continue
                    if self.wafer[row,column+1] > 0 :
                        H[i,i+1] = -self.t
                    if self.wafer[row+1,column] >0 and row+1 < self.wafer.shape[0]:
                        H[i,i+self.wafer.shape[1]] = -exp(1j*self.BField*i)*self.t
        Hupper = triu(H, 1)
        H = (Hupper.tocsr() + H.tocsr().conjugate().T).tolil()# might be faster
        self.H = H


    def simplesigma(self, ind_contact, E):
        from scipy.sparse import lil_matrix
        from scipy import complex128
        sigma = lil_matrix((self.canvas[0]*self.canvas[1],self.canvas[0]*self.canvas[1]), dtype=complex128)
        contact_node = 0
        for xypair in ind_contact.T:
            index = xypair[1] + self.wafer.shape[1]*xypair[0]
            sigma[index, index] = - self.t * self.__contact_greensfunction(ind_contact, contact_node, E)
            contact_node +=1
        return sigma

    def fermifunction(self, E, mu):
        """ Simple Fermifunction """
        from scipy import exp
        fermifnc = 1/(exp((E.real-mu)/self.kT)+1)
        return fermifnc

    def __contact_greensfunction(self, ind_contact, contact_node, E):
        """calculates the value of the transverse mode's square at the point
        of the contact_node, essentially sin(pi)sin(pi) multiplied with
        appropriate amplitude and constats """
        from scipy import arccos, exp, sin, pi
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
        """Build the self-energy matrix SIGMA by determining the nodes
        adjacent to a Contact and inserting the greensfunction times t**2
        (the greensfunction comes with 1/t)"""
        from scipy.sparse import lil_matrix
        from scipy import complex128
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
        A = (E+self.zplus)*eye(number_of_nodes,number_of_nodes,dtype=complex128, format='lil') - self.H - sigma_l  - sigma_r
        return A, sigma_in_l, sigma_in_r

    def simpleA(self, E):
        from scipy.sparse import eye
        from scipy import complex128
        number_of_nodes = self.wafer.shape[0]*self.wafer.shape[1]
        sigma_l = self.simplesigma(self.contacts[0],E - self.potential_drop[0])
        sigma_r =self.simplesigma(self.contacts[1], E - self.potential_drop[1])
        sigma_in_l = -2* sigma_l.imag[0:50, 0:50] * self.fermifunction(E, mu=self.mu_l)
        sigma_in_r = -2* sigma_r.imag[4950:5000, 4950:5000] * self.fermifunction(E, mu=self.mu_r)
        A = (E+self.zplus)*eye(number_of_nodes,number_of_nodes,dtype=complex128, format='lil') - self.H - sigma_l  - sigma_r
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
            #A, sigma_in_l, sigma_in_r = self.simpleA(energy)
            Ablock = SparseBlocks(A, self.block_sizes)
            #Ablock = SparseBlocks(A,[self.wafer.shape[1]]*self.wafer.shape[0] )
            integral = integrand.__call__(Ablock, sigma_in_l, sigma_in_r)*self.dE/(pi*self.a)
            print integral
            #self.simplewritetovtk(integral.real, str(i))
            self.writetovtk(integral.real, str(i))
            i+=1
            max_density.append(integral.real.max())
        return integral, max_density

    def simpleenergyintegrate(self,integrand,sigma_in_l=None,sigma_in_r=None):
        from scipy import pi, array
        from scipy.sparse import lil_matrix
        from sparseblockslice import SparseBlocks
        #if integrand == self.RRGM:
        #    value, foo = -integrand.__call__(Ablock).imag
        #elif integrand == self.LRGM and sigma_in_l is not None and sigma_in_r is not None:
        #    value = integrand.__call__(Ablock,sigma_in_l, sigma_in_r,).real
        #else:
            #print 'Please insert supported functions RRGM or LRGM(not none sigmas)'
        integral = array([0]*self.wafer.shape[0]*self.wafer.shape[1])
        #A = lil_matrix((len(self.nodes), len(self.nodes)))
        #print A
        max_density = []
        i=0
        for energy in self.Egrid:
            #A, sigma_in_l, sigma_in_r = self.build_A(energy)
            A, sigma_in_l, sigma_in_r = self.simpleA(energy)
            #Ablock = SparseBlocks(A, self.block_sizes)
            Ablock = SparseBlocks(A,[self.wafer.shape[1]]*self.wafer.shape[0] )
            integral = integral + integrand.__call__(Ablock, sigma_in_l, sigma_in_r)*self.dE/(pi*self.a)
            #self.writetovtk(integral.real, str(i))
            #self.writetovtk(integral.real, str(i))
            i+=1
            max_density.append(integral.real.max())
        return integral, max_density

    def build_convolutionkernel(self):
        from scipy import zeros, hstack, vstack
        from scipy.linalg import norm
        plusv_dim = self.wafer.shape[0]
        plush_dim = self.wafer.shape[1]
        kernel = zeros((plusv_dim, plush_dim))
        for i in range(plusv_dim):
            for j in range(plush_dim):
                if i==0 and j == 0: continue
                kernel[i,j] = 1/norm((i,j))
        self.kernel = kernel
        kernel = hstack((kernel[:,:0:-1], kernel))
        self.kernel = vstack((kernel[:0:-1,:], kernel))

    def hartree_from_density(self, density):
        from scipy.signal import fftconvolve
        from scipy import pi
        target = density.reshape(self.wafer.shape[0],self.wafer.shape[1])
        factor = (self.q**2)/(4* pi* self.eps0 * self.epsr)
        hartree = factor * fftconvolve(target, self.kernel, mode='valid')
        return hartree

    def selfconsistent(self):
        from scipy import ones
        initial_dens = ones((self.wafer.shape[0] * self.wafer.shape[1]))
        initial_phi = self.hartree_from_density(initial_dens)
        density = 

    def summatrix(self, mat):
        from summon import matrix
        from scipy.io import mmwrite
        from os import system
        mmwrite('tempmatrix', mat.real)
        system("sed '1,2d' tempmatrix.mtx > newfile.txt")
        m = matrix.Matrix()
        matrix.open_matrix('newfile.txt', m, format='imat')
        viewer = matrix.MatrixViewer(m, title="Sparsity")
        viewer.show()

