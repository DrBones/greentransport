class Model:

    def __init__(self, world):
        from scipy import linspace, cos, pi
        Model.canvas = world.canvas
        Model.atlas = world.atlas
        Model.hbar = world.hbar
        Model.m0 = world.m0
        Model.kb = world.kb

        self.nodes = world.nodes
        self.wafer = world.wafer
        self.contacts = world.contacts
        self.active_coords = world.active_coords
        self.block_sizes = world.block_sizes
        self.eps0 = world.eps0
        self.q = world.q
        #Potential Drop over legth of device
        self.a = 3e-9 # in meter
        self.alpha = 20e-12 # eV m
        #effective mass in eV real in GaAs 0.063
        self.mass = 0.063*Model.m0
        self.t0 = (Model.hbar**2)/(2*self.mass*(self.a**2))
        self.tso = self.alpha/(2 * self.a)
        #Temperature * k_boltzmann in eV, 0.0025ev~30K
        #self.potential_drop = [0.1*self.t0/2, -0.1* self.t0/2]# in eV
        self.potential_drop = [0,0]
        self.epsr = 12.85
        self.Temp = 2 #in Kelvin
        self.kT = Model.kb * self.Temp
        self.lambdaf = 10 # i believe in nanometer, 35 more realistic?
        self.BField = 0 # in Tesla, from 0-~10
        self.Balpha = self.BField * self.a**2 /(2 * pi *self.hbar) # without the leading q because of hbar in eV
        self.zplus = 1j*1e-12
        #self.band_bottom = 0
        self.band_bottom = -4*self.t0
        self.Efermi = self.band_bottom + 2*self.t0*(1-cos(2*pi/self.lambdaf))
        #self.Efermi = 0.1
        #self.Efermi = -3.8 * self.t0 # close to the bottom of the band at -4.0 t0, what bottom and band in what material ?
        self.Egrid = linspace(self.Efermi-0.4*self.t0,self.Efermi +0.4*self.t0,100)+self.zplus # in eV ?
        self.mu = self.Efermi
        self.dE = self.Egrid[1].real-self.Egrid[0].real
        #electro-chemical potential in eV
        self.mu_l = self.Efermi - (self.potential_drop[1] - self.potential_drop[0])/2
        self.mu_r = self.Efermi + (self.potential_drop[1] - self.potential_drop[0])/2

        self.__generate_potential_grid()
        self.grid2serialized(self.potential_grid)
        #self.__build_H()

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
                H[self.nodes[item][0],self.nodes[item][1]] = -exp(1j*self.BField*self.nodes[item][0])*self.t0
            else: pass
            if self.nodes[item][2] == None: continue
            H[self.nodes[item][0],self.nodes[item][2]] = -self.t0
        H = (H.tocsr() + H.tocsr().conjugate().T).tolil()# might be faster
        H.setdiag(self.potential_serialized+4*self.t0)
        self.H = H

    def simpleH(self):
        from scipy.sparse import lil_matrix, triu
        from scipy import complex128, exp, pi
        Ndim = self.wafer.shape[0]*self.wafer.shape[1]
        ystride = self.wafer.shape[1]
        H = lil_matrix((Ndim,Ndim),dtype=complex128)
        for row in range(self.wafer.shape[0]):
            for column in range(self.wafer.shape[1]):
                i = column+ystride*row
                if self.wafer[row,column] == 0:
                    H[i,i] = 1e4
                elif self.wafer[row,column] >0:
                    H[i,i] = 4*self.t0+self.potential_grid[row,column]
                    if column+1 < self.wafer.shape[1]:
                        if self.wafer[row,column+1] > 0 :
                            H[i,i+1] = -self.t0
                    if row+1 == self.wafer.shape[0]: continue
                    if self.wafer[row+1,column] >0:
                        H[i,i+ystride] = -exp(2 * pi*1j*self.Balpha*column%self.wafer.shape[1])*self.t0
        Hupper = triu(H, 1)
        H = (Hupper.tocsr() + H.tocsr().conjugate().T).tolil()
        self.H = H

    def spinH(self):
        from scipy.sparse import lil_matrix, triu
        from scipy import complex128, exp, pi
        Ndim = 2*self.wafer.shape[0]*self.wafer.shape[1]
        ystride = self.wafer.shape[1]
        H = lil_matrix((Ndim,Ndim),dtype=complex128)
        multi = 0
        for row in range(self.wafer.shape[0]):
            for column in range(self.wafer.shape[1]):
                i = column+ystride*row+multi*ystride
                if self.wafer[row,column] == 0:
                    H[i,i] = H[i+ystride,i+ystride] =  1e4
                elif self.wafer[row,column] >0:
                    H[i,i] = H[i+ystride,i+ystride] = 4*self.t0+self.potential_grid[row,column]
                    if column+1 < self.wafer.shape[1]:
                        if self.wafer[row,column+1] > 0 :
                            H[i,i+1] = H[i+ystride,i+1+ystride] = -self.t0
                            H[i,i+1+ystride] = self.tso
                            H[i+1,i+ystride] = -self.tso
                    if row+1 == self.wafer.shape[0]: continue
                    if self.wafer[row+1,column] >0:
                        H[i,i+2*ystride] = H[i+ystride,i+ystride+2*ystride] = -exp(2 *
                                           pi*1j*self.Balpha*column%self.wafer.shape[1])*self.t0
                        H[i,i+3*ystride] = H[i+ystride,i+2*ystride] = -1j*self.tso
            multi +=1
        Hupper = triu(H, 1)
        H = (Hupper.tocsr() + H.tocsr().conjugate().T).tolil()
        self.H = H


    def fermifunction(self, E, mu):
        """ Simple Fermifunction """
        from scipy import exp
        fermifnc = 1/(exp((E.real-mu)/self.kT)+1)
        return fermifnc

    def __contact_greensfunction(self, ind_contact, node_i, node_j, E):
        """calculates the value of the transverse mode's square at the point
        of the contact_node, essentially sin(pi)sin(pi) multiplied with
        appropriate amplitude and constats """
        from scipy import arccos, exp, sin, pi, sqrt
        length = (ind_contact.shape[1]-1.0)
        amplitude = 1/sqrt(length)
        mode_energy = self.hbar**2 * pi**2/(2*self.mass * (length*self.a)**2)
        #Amplitude = 1
        ka = arccos(1-(E+mode_energy-self.band_bottom)/(2*self.t0))
        Phase = exp(1j*ka)
        xi_i = amplitude * sin(pi * node_i/length)
        xi_j = amplitude * sin(pi * node_j/length)
        greensfunction =xi_j * xi_i *  Phase
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
            sigma[index, index] = - self.t0 * self.__contact_greensfunction(ind_contact, contact_node, E)
        return sigma

    def build_A(self, E):
        from scipy.sparse import eye
        from scipy import complex128
        number_of_nodes = len(self.nodes)
        sigma_l = self.build_sigma(self.contacts[0],E - self.potential_drop[0])
        sigma_r =self.build_sigma(self.contacts[1], E - self.potential_drop[1])
        sigma_in_l = -2* sigma_l.imag[0:self.block_sizes[0],0:self.block_sizes[0]] * self.fermifunction(E, mu=self.mu_l)
        sigma_in_r = -2* sigma_r.imag[-self.block_sizes[-1]:,-self.block_sizes[-1]:] * self.fermifunction(E, mu=self.mu_r)
        A = (E)*eye(number_of_nodes,number_of_nodes,dtype=complex128, format='lil') - self.H - sigma_l  - sigma_r
        return A, sigma_in_l, sigma_in_r

    def spinsigma(self, ind_contact, E,mode='normal'):
        from scipy.sparse import lil_matrix
        from scipy import complex128
        if mode == 'normal':
            multi = 1
        elif mode == 'spin':
            multi = 2
        ystride = self.wafer.shape[1]
        Ndim = self.canvas[0]*self.canvas[1]
        sigma = lil_matrix((multi*Ndim,multi*Ndim), dtype=complex128)
        contact_node = 0
        for xypair in ind_contact.T:
            index = xypair[1] + ystride*xypair[0]
            sigma[index, index] = - self.t0 * self.__contact_greensfunction(ind_contact, contact_node, contact_node, E)
            if mode=='normal':contact_node+=1; continue
            sigma[index+ystride, index+ystride] = - self.t0 * self.__contact_greensfunction(ind_contact, contact_node, contact_node, E)
            contact_node +=1
        return sigma

    def simpleA(self, E):
        from scipy.sparse import eye
        from scipy import complex128
        number_of_nodes = self.wafer.shape[0]*self.wafer.shape[1]
        sigma_l = self.simplesigma(self.contacts[0],E - self.potential_drop[0])
        sigma_r =self.simplesigma(self.contacts[1], E - self.potential_drop[1])
        sigma_in_l = -2* sigma_l.imag[0:50, 0:50] * self.fermifunction(E, mu=self.mu_l)
        sigma_in_r = -2* sigma_r.imag[4950:5000, 4950:5000] * self.fermifunction(E, mu=self.mu_r)
        A = (E)*eye(number_of_nodes,number_of_nodes,dtype=complex128, format='lil') - self.H - sigma_l  - sigma_r
        return A, sigma_in_l, sigma_in_r

    def spinA(self,E, mode='normal'):
        from scipy.sparse import eye
        from scipy import complex128
        if mode == 'normal':
            multi = 1
        elif mode == 'spin':
            multi = 2
        number_of_nodes = self.wafer.shape[0]*self.wafer.shape[1]
        sigma_l = self.spinsigma(self.contacts[0],E - self.potential_drop[0],mode)
        sigma_r =self.spinsigma(self.contacts[1], E - self.potential_drop[1],mode)
        print sigma_l
        start = (number_of_nodes-self.wafer.shape[1])
        end = number_of_nodes
        sigma_in_l = -2* sigma_l.imag[0:multi*self.wafer.shape[1], 0:multi*self.wafer.shape[1]] * self.fermifunction(E, mu=self.mu_l)
        sigma_in_r = -2* sigma_r.imag[multi*start:multi*end,start:multi*end] * self.fermifunction(E, mu=self.mu_r)
        A = (E+self.band_bottom)*eye(multi*number_of_nodes,multi*number_of_nodes,dtype=complex128, format='lil') - self.H - sigma_l  - sigma_r
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
        #integral = array([0]*len(self.nodes))
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
            fncvalue = -integrand.__call__(Ablock)[0]*self.fermifunction(energy, self.mu)*self.dE/(pi*self.a)
            self.writetovtk(fncvalue.imag, str(i))
            summi = integral + fncvalue
            i+=1
            max_density.append(integral.imag.min())
            self.writetovtk(summi.imag, 'summi')

        return integral, max_density

    def adaptiveenergy(self):
        pass


    def spinenergyintegrate(self,integrand,sigma_in_l=None,sigma_in_r=None, mode='normal'):
        from scipy import pi, array
        #from scipy.sparse import lil_matrix
        from sparseblockslice import SparseBlocks
        #from scipy import vstack
        #if integrand == self.RRGM:
        #    value, foo = -integrand.__call__(Ablock).imag
        #elif integrand == self.LRGM and sigma_in_l is not None and sigma_in_r is not None:
        #    value = integrand.__call__(Ablock,sigma_in_l, sigma_in_r,).real
        #else:
            #print 'Please insert supported functions RRGM or LRGM(not none sigmas)'
        #hills = array([0]*self.wafer.shape[0]*self.wafer.shape[1])
        #A = lil_matrix((len(self.nodes), len(self.nodes)))
        #print A
        if mode == 'normal':
            multiplier = 1
        elif mode == 'spin':
            multiplier = 2
        integral = array([0]*multiplier*self.wafer.shape[0]*self.wafer.shape[1])
        max_density = []
        i=0
        print "Current Energy:     ", "Left Occupation:     ", "Right Occupation:     ", "Maximum Density:"
        for energy in self.Egrid:
            #A, sigma_in_l, sigma_in_r = self.build_A(energy)
            A, sigma_in_l, sigma_in_r = self.spinA(energy,mode)
            #Ablock = SparseBlocks(A, self.block_sizes)
            Ablock = SparseBlocks(A,[self.wafer.shape[1]*multiplier]*self.wafer.shape[0] )
            fncvalue = integrand.__call__(Ablock, sigma_in_l, sigma_in_r)*self.dE/(pi*self.a**2)
            print i, energy,"            ", self.fermifunction(energy, mu=self.mu_l),"               ", self.fermifunction(energy, mu=self.mu_r),"          ", fncvalue.real.max()
            #self.writetovtk(fncvalue.real, str(i))
            integral = integral + fncvalue
            #hills = vstack((hills,integral))
            i+=1
            max_density.append(fncvalue.real.max())
        #self.writetovtk(integral.real, 'integrated')
        return integral, max_density

    def dorgm(self,name,energy, mode='normal'):
        from sparseblockslice import SparseBlocks
        from io import writeVTK
        multiplier = 1
        #energy = self.Efermi
        A, sigma_in_l, sigma_in_r = self.spinA(energy,mode)
        Ablock = SparseBlocks(A,[self.wafer.shape[1]*multiplier]*self.wafer.shape[0] )
        densi, temp1, temp2 = self.RRGM(Ablock)
        dens = -densi.imag/(self.a**2)*self.fermifunction(energy, self.mu)
        writeVTK(name, 49, 99, pointData={"Density":dens})
        return dens

    def build_convolutionkernel(self):
        from scipy import zeros, hstack, vstack
        from scipy.linalg import norm
        plusv_dim = self.wafer.shape[0]
        plush_dim = self.wafer.shape[1]
        kernel = zeros((plusv_dim, plush_dim))
        for i in range(plusv_dim):
            for j in range(plush_dim):
                if i==0 and j == 0: continue
                kernel[i,j] = 1/(self.a*norm((i,j)))
        self.kernel = kernel
        kernel = hstack((kernel[:,:0:-1], kernel))
        self.kernel = vstack((kernel[:0:-1,:], kernel))

    def hartree_from_density(self, density):
        """ Somehow, a self.a**2 is missing because of the transition
        from intgral to sum. The result has the right magnitude without though"""
        from scipy.signal import fftconvolve
        from scipy import pi
        target = density.reshape(self.wafer.shape[0],self.wafer.shape[1])
        factor = (self.q**2)/(4* pi* self.eps0 * self.epsr)
        hartree = factor * fftconvolve(target, self.kernel, mode='valid')
        return hartree

    def selfconsistent(self,E):
        from scipy import ones, pi
        from sparseblockslice import SparseBlocks
        #energy = self.Egrid[50]
        initial_dens = 1e14*ones((self.wafer.shape[0] * self.wafer.shape[1]))
        initial_phi = self.hartree_from_density(initial_dens)
        A, sigma_in_l, sigma_in_r = self.simpleA(E)
        Ablock = SparseBlocks(A,[self.wafer.shape[1]]*self.wafer.shape[0])
        density = self.LRGM(Ablock, sigma_in_l, sigma_in_r)*self.dE/(pi*self.a)
        return density

    def cells_from_points(self, array):
        from scipy.ndimage.interpolation import geometric_transform
        def shift_func(output_coords):
            return (output_coords[0] - 0.5, output_coords[1] - 0.5)
        cells = geometric_transform(array, shift_func, output_shape=(array.shape[0]-1,array.shape[1]-1))
        return cells

    def setup(self):
        from sparseblockslice import SparseBlocks
        from scipy import pi
        energy = self.Efermi
        A, sigma_in_l, sigma_in_r = self.simpleA(energy)
        Ablock = SparseBlocks(A,[self.wafer.shape[1]]*self.wafer.shape[0] )
        fncvalue = self.LRGM(Ablock, sigma_in_l, sigma_in_r)*self.dE/(pi*self.a**2)
        return fncvalue
