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
        self.multi = 1
        self.a = 3e-9 # in meter
        self.alpha = 20e-12 # eV m
        #effective mass in eV real in GaAs 0.063
        self.mass = 0.063*Model.m0
        self.t0 = (Model.hbar**2)/(2*self.mass*(self.a**2))
        #self.tso = self.alpha/(2 * self.a)
        self.tso = 0.2*self.t0
        #Temperature * k_boltzmann in eV, 0.0025ev~30K
        self.epsr = 12.85
        self.Temp = 2 #in Kelvin
        self.kT = Model.kb * self.Temp
        self.lambdaf = 10 # i believe in nanometer, 35 more realistic?
        self.BField = 0 # in Tesla, from 0-~10
        self.Balpha = self.BField * self.a**2 /(2 * pi *self.hbar) # without the leading q because of hbar in eV
        self.zplus = 1j*1e-12
        self.band_bottom = 0
        #self.band_bottom = -4*self.t0
        #self.Efermi = self.band_bottom + 2*self.t0*(1-cos(2*pi/self.lambdaf))
        self.Efermi = 0.2*self.t0
        #self.potential_drop = [0,0]
        self.potential_drop = [0.004*self.t0/2, -0.004* self.t0/2]# in eV
        #self.Efermi = -3.8 * self.t0 # close to the bottom of the band at -4.0 t0, what bottom and band in what material ?
        #self.Egrid = linspace(self.Efermi-0.4*self.t0,self.Efermi +0.4*self.t0,100)+self.zplus # in eV ?
        self.mu = self.Efermi
        #self.dE = self.Egrid[1].real-self.Egrid[0].real
        #electro-chemical potential in eV
        self.mu_l = self.Efermi - (self.potential_drop[1] - self.potential_drop[0])/2
        self.mu_r = self.Efermi + (self.potential_drop[1] - self.potential_drop[0])/2
        from scipy import r_
        self.potential_grid = r_[[self.potential_drop[0]]*30*20,[0]*30*30,[self.potential_drop[1]]*30*20].reshape(70,30)
        #self.__generate_potential_grid()
        #self.grid2serialized(self.potential_grid)
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
                            if self.wafer[row,column+1] == 239:
                                H[i,i+1+ystride] = self.tso
                                H[i+1,i+ystride] = -self.tso
                    if row+1 == self.wafer.shape[0]: continue
                    if self.wafer[row+1,column] >0:
                        H[i,i+2*ystride] = H[i+ystride,i+ystride+2*ystride] = -exp(2 *
                                           pi*1j*self.Balpha*column%self.wafer.shape[1])*self.t0
                        if self.wafer[row,column] == 239 and self.wafer[row+1,column] == 239:
                            H[i,i+3*ystride] = H[i+ystride,i+2*ystride] = -1j*self.tso
            multi +=1
        Hupper = triu(H, 1)
        H = (Hupper.tocsr() + H.tocsr().conjugate().T).tolil()
        self.H = H

    def fermifunction(self, E_tot, mu):
        """ Simple Fermifunction """
        from scipy import exp
        fermifnc = 1/(exp((E_tot.real-mu)/self.kT)+1)
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

    def eigensigma(self):
        from scipy.linalg import eig
        from scipy.sparse import lil_matrix
        from scipy import argsort,where
        #from scipy.sparse.linalg import eigen
        transverseH = lil_matrix((self.wafer.shape[1],self.wafer.shape[1]))
        transverseH.setdiag([2*self.t0]*self.wafer.shape[1])
        transverseH.setdiag([-self.t0]*self.wafer.shape[1],1)
        transverseH.setdiag([-self.t0]*self.wafer.shape[1],-1)
        #from pudb import set_trace; set_trace()
        v,d = eig(transverseH.todense())
        ndx = argsort(v)
        d=d[:,ndx]
        v=v[ndx]
        self.v = v
        self.d = d
        self.maxmode = where(self.v < self.Efermi-self.band_bottom)[0].max()
        if v.max() > self.Efermi-self.band_bottom:
            print 'Some mode energies larger than fermi energy, only up to mode {0} will fit'.format(self.maxmode)
            print 'Argument num_modes="all" takes only modes low enough'
            print ''

    def sigma(self, ind_contact, E, num_modes='all'):
        """
        Takes abolute energies from band_bottom to around Efermi and further 
        until the fermifunction puts and end to this
        """
        from scipy import sqrt,exp,asarray,dot,diag
        from scipy.sparse import lil_matrix
        from scipy import complex128
        if num_modes == 'analytical' :
            return self.analytic_sigma(ind_contact,E)
        elif num_modes == 'all':
            num_modes =self.maxmode
        Ndim = self.canvas[0]*self.canvas[1]
        #dd = diag(-self.t0*exp(1j*sqrt((E-self.band_bottom-self.v[:num_modes])/self.t0)))
        dd = diag(-self.t0*exp(1j*sqrt((E-self.v[:num_modes])/self.t0)))
        print 'Energy in Sigma used: ',E-self.v[:num_modes]
        sigma = lil_matrix((self.multi*Ndim,self.multi*Ndim), dtype=complex128)
        if ind_contact[0].min() == 0:
            sigma[0:self.wafer.shape[1], 0:self.wafer.shape[1]] = asarray(dot(dot(self.d[:,:num_modes],dd),self.d[:,:num_modes].T))
            if self.multi==2:
                sigma[self.wafer.shape[1]:self.wafer.shape[1]*self.multi, self.wafer.shape[1]:self.wafer.shape[1]*self.multi] = asarray(dot(dot(self.d[:,:num_modes],dd),self.d[:,:num_modes].T))
        elif ind_contact[0].max() == self.wafer.shape[0]-1:
            #import pudb; pudb.set_trace()
            sigma[-self.block_sizes[-1]:, -self.block_sizes[-1]:] = asarray(dot(dot(self.d[:,:num_modes],dd),self.d[:,:num_modes].T))
            if self.multi == 2:
                sigma[-self.block_sizes[-1]*self.multi:-self.block_sizes[-1], -self.block_sizes[-1]*self.multi:-self.block_sizes[-1]] = asarray(dot(dot(self.d[:,:num_modes],dd),self.d[:,:num_modes].T))
        return sigma

    def analytic_sigma(self, ind_contact, E):
        from scipy.sparse import lil_matrix
        from scipy import complex128
        ystride = self.wafer.shape[1]
        Ndim = self.canvas[0]*self.canvas[1]
        sigma = lil_matrix((self.multi*Ndim,self.multi*Ndim), dtype=complex128)
        contact_node = 0
        for xypair in ind_contact.T:
            index = xypair[1] + ystride*xypair[0]
            sigma[index, index] = - self.t0 * self.__contact_greensfunction(ind_contact, contact_node, contact_node, E)
            if self.multi==1:contact_node+=1; continue
            sigma[index+ystride, index+ystride] = - self.t0 * self.__contact_greensfunction(ind_contact, contact_node, contact_node, E)
            contact_node +=1
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

    def spinA(self,E_rel):
        from scipy.sparse import eye
        from scipy import complex128
        number_of_nodes = self.block_sizes[1]*len(self.block_sizes)
        #E_tot=self.Efermi+E_rel
        E_tot=E_rel
        sigma_l = self.sigma(self.contacts[0],E_tot- self.potential_drop[0],num_modes=1)
        sigma_r =self.sigma(self.contacts[1], E_tot - self.potential_drop[1],num_modes=1)
        sigma_in_l = -2* sigma_l.imag[0:self.multi*self.block_sizes[1], 0:self.multi*self.block_sizes[1]] * self.fermifunction(E_tot, mu=self.mu_l)
        sigma_in_r = -2* sigma_r.imag[-self.block_sizes[-1]*self.multi:,-self.block_sizes[-1]*self.multi:] * self.fermifunction(E_tot, mu=self.mu_r)
        I =eye(self.multi*number_of_nodes,self.multi*number_of_nodes,dtype=complex128, format='lil')
        A = (E_tot+self.zplus)*I - self.H - sigma_l  - sigma_r
        return A, sigma_in_l, sigma_in_r

    #def energyintegrate(self,integrand,sigma_in_l=None,sigma_in_r=None):
    #    from scipy import pi, array
    #    from scipy.sparse import lil_matrix
    #    from sparseblockslice import SparseBlocks
    #    #if integrand == self.RRGM:
    #    #    value, foo = -integrand.__call__(Ablock).imag
    #    #elif integrand == self.LRGM and sigma_in_l is not None and sigma_in_r is not None:
    #    #    value = integrand.__call__(Ablock,sigma_in_l, sigma_in_r,).real
    #    #else:
    #        #print 'Please insert supported functions RRGM or LRGM(not none sigmas)'
    #    #integral = array([0]*len(self.nodes))
    #    integral = array([0]*self.wafer.shape[0]*self.wafer.shape[1])
    #    #A = lil_matrix((len(self.nodes), len(self.nodes)))
    #    #print A
    #    max_density = []
    #    i=0
    #    for energy in self.Egrid:
    #        #A, sigma_in_l, sigma_in_r = self.build_A(energy)
    #        A, sigma_in_l, sigma_in_r = self.simpleA(energy)
    #        #Ablock = SparseBlocks(A, self.block_sizes)
    #        Ablock = SparseBlocks(A,[self.wafer.shape[1]]*self.wafer.shape[0] )
    #        fncvalue = -integrand.__call__(Ablock)[0]*self.fermifunction(energy, self.mu)*self.dE/(pi*self.a)
    #        self.writetovtk(fncvalue.imag, str(i))
    #        summi = integral + fncvalue
    #        i+=1
    #        max_density.append(integral.imag.min())
    #        self.writetovtk(summi.imag, 'summi')

    #    return integral, max_density

    #def spinenergyintegrate(self,integrand,sigma_in_l=None,sigma_in_r=None):
    #    from scipy import pi, array
    #    #from scipy.sparse import lil_matrix
    #    from sparseblockslice import SparseBlocks
    #    #from scipy import vstack
    #    #if integrand == self.RRGM:
    #    #    value, foo = -integrand.__call__(Ablock).imag
    #    #elif integrand == self.LRGM and sigma_in_l is not None and sigma_in_r is not None:
    #    #    value = integrand.__call__(Ablock,sigma_in_l, sigma_in_r,).real
    #    #else:
    #        #print 'Please insert supported functions RRGM or LRGM(not none sigmas)'
    #    #hills = array([0]*self.wafer.shape[0]*self.wafer.shape[1])
    #    #A = lil_matrix((len(self.nodes), len(self.nodes)))
    #    #print A
    #    integral = array([0]*self.multi*self.wafer.shape[0]*self.wafer.shape[1])
    #    max_density = []
    #    i=0
    #    print "Current Energy:     ", "Left Occupation:     ", "Right Occupation:     ", "Maximum Density:"
    #    for energy in self.Egrid:
    #        #A, sigma_in_l, sigma_in_r = self.build_A(energy)
    #        A, sigma_in_l, sigma_in_r = self.spinA(energy)
    #        #Ablock = SparseBlocks(A, self.block_sizes)
    #        Ablock = SparseBlocks(A,[self.wafer.shape[1]*self.multi]*self.wafer.shape[0] )
    #        fncvalue = integrand.__call__(Ablock, sigma_in_l, sigma_in_r)*self.dE/(pi*self.a**2)
    #        print i, energy,"            ", self.fermifunction(energy, mu=self.mu_l),"               ", self.fermifunction(energy, mu=self.mu_r),"          ", fncvalue.real.max()
    #        #self.writetovtk(fncvalue.real, str(i))
    #        integral = integral + fncvalue
    #        #hills = vstack((hills,integral))
    #        i+=1
    #        max_density.append(fncvalue.real.max())
    #    #self.writetovtk(integral.real, 'integrated')
    #    return integral, max_density

    def adaptiveenergy(self):
        pass

    def dorrgm(self,E_rel):
        from aux import SparseBlocks
        from greensolver import rrgm
        #from io import writeVTK
        #energy = self.Efermi
        A, sigma_in_l, sigma_in_r = self.spinA(E_rel)
        Ablock = SparseBlocks(A,self.block_sizes*self.multi )
        dens, temp1, temp2 = rrgm(Ablock)
        #dens = -densi.imag/(self.a**2)*self.fermifunction(energy, self.mu)
        #writeVTK(name, 49, 99, pointData={"Density":dens})
        return dens

    def dolrgm(self,energy):
        from aux import SparseBlocks
        from greensolver import lrgm
        #from io import writeVTK
        A, sigma_in_l, sigma_in_r = self.spinA(energy)
        Ablock = SparseBlocks(A,self.block_sizes*self.multi)
        dens = lrgm(Ablock, sigma_in_l, sigma_in_r)
        #dens = densi.real/(self.a**2)
        #name = str(energy)
        #writeVTK(name, 49, 99, pointData={"Density":dens})
        return dens

    def spindens(self,energy):
        from scipy import split,pi
        if self.multi == 1:
            self.setmode('spin')
        dens = self.dolrgm(energy)
        Gup, Gdown = split(dens.reshape(self.wafer.shape[0],self.wafer.shape[1]*2),2,axis=1)
        Sz = self.hbar/(4*pi*1j*self.a**2)*(Gup-Gdown)
        print 'max Spin Split: ', Sz.imag.max()-Sz.imag.min()
        return Sz.imag

    def edens(self,energy):
        from scipy import split,pi
        lrgm_out = self.dolrgm(energy)
        if self.multi ==1:
            edensity = lrgm_out.reshape(self.wafer.shape)*2/(2*pi*self.a**2) #times 2 for spin
        if self.multi ==2:
            Gup, Gdown = split(lrgm_out.reshape(self.wafer.shape[0],self.wafer.shape[1]*2),2,axis=1)
            edensity = 1/(2*pi*self.a**2)*(Gup+Gdown)
        return edensity

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

    def setmode(self,mode='normal'):
        if mode == 'normal':
            self.simpleH()
            self.multi = 1
        elif mode == 'spin':
            self.spinH()
            self.multi = 2
