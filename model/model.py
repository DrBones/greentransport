class Model:
    from aux import edens,spindens,transmission
    import networkx as nx

    def __init__(self, world):
        from scipy import linspace, cos, pi
        from parameters import p
        self.p = p
        Model.canvas = p.canvas
        Model.atlas = world.atlas

        #self.nodes = world.nodes
        self.contacts = p.contacts
        # self.raw_coords = p.raw_coords
        self.tuple_canvas_coordinates = p.tuple_canvas_coordinates
        self.eps0 = p.eps0
        self.p.multi = 1
        self.p.empty_potential_grid()
        #self.stepgrid(2,2)
        # self.circular_qpc()
        #self.__generate_potential_grid()
        #self.grid2serialized(self.potential_grid)
        #self.__build_H()

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
        from scipy.sparse import lil_matrix, triu,bmat,eye,hstack
        from scipy import complex128, exp, pi
        Balpha = self.BField * self.a**2 /(2 * pi *self.p.hbar) # without the leading q because of hbar in eV
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
                        H[i,i+ystride] = -exp(2 * pi*1j*Balpha*column%self.wafer.shape[1])*self.t0
        Hupper = triu(H, 1)
        H = (Hupper.tocsr() + H.tocsr().conjugate().T).tolil()
        self.H = H
        pad = lil_matrix((self.wafer.shape[1],self.wafer.shape[1]))
        hopping = -self.t0*eye(self.wafer.shape[1],self.wafer.shape[1])
        hopping_upcont = hstack((hopping,lil_matrix((hopping.shape[0],Ndim-hopping.shape[1]))))
        hopping_lowcont = hstack((lil_matrix((hopping.shape[0],Ndim-hopping.shape[1])),hopping))
        self.Hpad = bmat([[pad,hopping_upcont,None],[hopping_upcont.T,H,hopping_lowcont.T],[None,hopping_lowcont,pad]])

    def spinH(self):
        from scipy.sparse import lil_matrix, triu,bmat
        from scipy import complex128, exp, pi
# TODO just changed offdiagonal tso sign, this seems more correct
        Balpha = self.BField * self.a**2 /(2 * pi *self.p.hbar) # without the leading q because of hbar in eV
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
                        if row == 0 and self.contacts[0].SO == True:
                            H[i,i+1+ystride] = self.tso
                            H[i+1,i+ystride] = -self.tso
                        if row == self.wafer.shape[0]-1 and self.contacts[1].SO == True:
                            H[i,i+1+ystride] = self.tso
                            H[i+1,i+ystride] = -self.tso
                    if row == self.wafer.shape[0]-2 and self.contacts[1].SO == True:
                        H[i,i+3*ystride] = H[i+ystride,i+2*ystride] = -1j*self.tso
                    if row+1 == self.wafer.shape[0]: continue
                    if self.wafer[row+1,column] >0:
                        H[i,i+2*ystride] = H[i+ystride,i+ystride+2*ystride] = -exp(2 *
                                           pi*1j*Balpha*column%self.wafer.shape[1])*self.t0
                        if self.wafer[row,column] == 239 and self.wafer[row+1,column] == 239:
                            H[i,i+3*ystride] = H[i+ystride,i+2*ystride] = -1j*self.tso
                    if row == 0 and self.wafer[row+1,column] == 239:
                        H[i,i+3*ystride] = H[i+ystride,i+2*ystride] = -1j*self.tso
            multi +=1
        Hupper = triu(H, 1)
        H = (Hupper.tocsr() + H.tocsr().conjugate().T).tolil()
        self.H = H
        pad = lil_matrix((self.wafer.shape[1],self.wafer.shape[1]))
        self.Hpad = bmat([[pad,None,None],[None,H,None],[None,None,pad]])

    def generate_graph(self):
        # TODO make lazy initialized
        from aux import digraph_from_tuple_coords
        print 'Generating graph of device'
        self.graph = digraph_from_tuple_coords(self.p.tuple_canvas_coordinates,self.canvas.shape[1])
        for contact in self.contacts:
            print 'Generating graph of interface: ', contact.index
            contact.recreate_graph()

    def generate_balanced_levelstructure(self):
        from aux import BreadthFirstLevels, bisect
        #import pudb; pudb.set_trace()
        nodes_left  = self.contacts[0].node_names
        nodes_right = self.contacts[1].node_names
        BFS_levelstructure = BreadthFirstLevels(self.graph,root=nodes_left,end=nodes_right)
        N = len(list(BFS_levelstructure))
        nodes_to_bisect = set(self.graph.nodes())-nodes_left-nodes_right
        levelstructure = [nodes_left] + bisect(self.graph,N-2,nodes_left,nodes_to_bisect,nodes_right) + [nodes_right]
        self.levelstructure = levelstructure

    def hamiltonian_from_graph(self):
        from numpy import complex128
        from networkx import to_scipy_sparse_matrix
        if ('levelstructure' in dir(self)):
            print 'using levelstructure information'
            self.nodelist =[item for level in self.levelstructure for item in sorted(level)]
            self.block_sizes = [ len(level) for level in self.levelstructure]
            self.H = to_scipy_sparse_matrix(self.graph,nodelist=self.nodelist,dtype=complex128)
        else:
            self.H = to_scipy_sparse_matrix(self.graph,dtype=complex128)

    def update_hamil_diag(self):
        from numpy import repeat,array
        def serial_pot(self):
            """ Generator expression yielding the base potential
            4*t0 plus the onsite potential taken from self.potenital_grid.
            No offsetting is done! Adapt Potential of canvas to fit device """
            if self.p.multi == 2:
                coordinate_array = repeat(self.tuple_canvas_coordinates,2,axis=0)
            else:
                coordinate_array = array(self.tuple_canvas_coordinates)

            #takes levelstructure permutations into account
            coordinates = coordinate_array[self.nodelist]
            for xy in coordinates:
                yield 4*self.p.t0+self.p.potential_grid[tuple(xy)]
        self.H.setdiag(list(serial_pot(self)))

    def fermifunction(self, E, mu):
        """ Simple Fermifunction """
        from scipy import exp
        fermifnc = 1/(exp((E.real-mu)/self.p.kT)+1)
        print fermifnc
        return fermifnc

    def build_A(self, E):
        from scipy.sparse import eye
        from scipy import complex128
        number_of_nodes = len(self.nodes)
        sigma_l = self.build_sigma(self.contacts[0],E - self.potential_drop[0])
        sigma_r =self.build_sigma(self.contacts[1], E - self.potential_drop[1])
        sigma_in_l = -2* sigma_l.imag[0:self.block_sizes[0],0:self.block_sizes[0]] * self.fermifunction(E, mu=self.p.mu_l)
        sigma_in_r = -2* sigma_r.imag[-self.block_sizes[-1]:,-self.block_sizes[-1]:] * self.fermifunction(E, mu=self.p.mu_r)
        A = (E)*eye(number_of_nodes,number_of_nodes,dtype=complex128, format='lil') - self.H - sigma_l  - sigma_r
        return A, sigma_in_l, sigma_in_r

    def spinA(self,E_rel):
        from scipy.sparse import eye
        from scipy import complex128
        number_of_nodes = len(self.p.tuple_canvas_coordinates)
        #global sigma_l,sigma_r,sigma_in_l,sigma_in_r
        #E_tot=self.Efermi+E_rel
        self.p.E = E_rel
        if (not ('graph' in dir(self))):
            print "Using eigenvalue decomp of transfermatrix of blocks to build selfenergy"
            if (not ('lastenergy' in dir(self)) or self.lastenergy != E_rel):
                print 'Generating new selfenergy matrices'
                sigma_l = self.transfersigma(self.contacts[0],E_rel- self.p.potential_drop[0])
                self.sigma_l = sigma_l
                #sigma_l = self.transfersigma(self.contacts[0], E_tot)
                sigma_r =self.transfersigma(self.contacts[1], E_rel- self.p.potential_drop[1])
                self.sigma_r = sigma_r
                #sigma_r =self.transfersigma(self.contacts[1], E_tot)
                self.gamma_l = self.gamma(sigma_l)
                self.gamma_r = self.gamma(sigma_r)
                sigma_in_l = -2* sigma_l.imag[0:self.block_sizes[0], 0:self.block_sizes[0]] * self.fermifunction(E_rel, mu=self.p.mu_l)
                self.sigma_in_l = sigma_in_l
                sigma_in_r = -2* sigma_r.imag[-self.block_sizes[-1]:,-self.block_sizes[-1]:] * self.fermifunction(E_rel, mu=self.p.mu_r)
                self.sigma_in_r = sigma_in_r
            else:
                sigma_l = self.sigma_l
                sigma_r = self.sigma_r
                sigma_in_l = self.sigma_in_l
                sigma_in_r = self.sigma_in_r
        else:
            print "Using schur decomp of lead_graphs to build selfenergy"
            # if (not ('lastenergy' in dir(self)) or self.lastenergy != E_rel):
            print 'Generating new selfenergy matrices'
            sigma_l = self.contacts[0].generate_sigma()
            self.sigma_l = sigma_l
            sigma_r = self.contacts[1].generate_sigma()
            self.sigma_r = sigma_r
            # self.gamma_l = self.gamma(sigma_l)
            # self.gamma_r = self.gamma(sigma_r)
            sigma_in_l = -2* sigma_l.imag[0:self.block_sizes[0], 0:self.block_sizes[0]] * self.fermifunction(E_rel, mu=self.p.mu_l)
            self.sigma_in_l = sigma_in_l
            sigma_in_r = -2* sigma_r.imag[-self.block_sizes[-1]:,-self.block_sizes[-1]:] * self.fermifunction(E_rel, mu=self.p.mu_r)
            self.sigma_in_r = sigma_in_r
            # else:
            #     sigma_l = self.sigma_l
            #     sigma_r = self.sigma_r
            #     sigma_in_l = self.sigma_in_l
            #     sigma_in_r = self.sigma_in_r

        I =eye(number_of_nodes*self.p.multi,number_of_nodes*self.p.multi,dtype=complex128, format='lil')
        self.lastenergy = E_rel
        print 'sigmainl',sigma_in_l.shape
        print 'sigmainr',sigma_in_r.shape
        print 'I', I.shape
        print 'H', self.H.shape
        A = (E_rel+self.p.zplus)*I - self.H - sigma_l  - sigma_r
        return A, sigma_in_l, sigma_in_r

    def adaptiveenergy(self):
        pass

    def dorrgm(self,energy):
        from aux import SparseBlocks
        from greensolver import rrgm
        #from io_spinr import writeVTK
        #energy = self.Efermi
        A, sigma_in_l, sigma_in_r = self.spinA(energy)
        Ablock = SparseBlocks(A,self.block_sizes )
        diag, grl, Gr= rrgm(Ablock)
        self.grl=grl
        #dens = -densi.imag/(self.a**2)*self.fermifunction(energy, self.mu)
        #writeVTK(name, 49, 99, pointData={"Density":dens})
        return diag, grl, Gr

    def conductance(self,energy):
        """ WRRRROOOOOOONNNNNGGGGG """
        from aux import transmission
        from scipy import pi
        rrgm_out = self.dorrgm(energy)
        G_pq = self.e**2/(self.h*2*pi)*transmission(rrgm_out)
        return G_pq

    def pdfit(self, energy):
        import matplotlib
        matplotlib.use('pdf')
        import matplotlib.pylab as mp
        lrgm_out = self.dolrgm(energy)
        dens = self.edens(lrgm_out)
        dens_spin = self.spindens(lrgm_out)
        mp.imshow(dens[0])
        mp.colorbar()
        mp.savefig('dummy',transparent=True)
        mp.close('all')
        mp.imshow(dens_spin,cmap='RdBu_r')
        mp.colorbar()
        mp.savefig('dummy'+'spindens',transparent=True)
        


    def dolrgm(self,energy):
        #import pudb; pudb.set_trace()
        from aux import SparseBlocks
        from greensolver import lrgm
        #from io_spinr import writeVTK
        A, sigma_in_l, sigma_in_r = self.spinA(energy)
        Ablock = SparseBlocks(A,self.block_sizes)
        lrgm_value= lrgm(Ablock, sigma_in_l, sigma_in_r)
        self.grl = lrgm_value[1]
        #dens = densi.real/(self.a**2)
        #name = str(energy)
        #writeVTK(name, 49, 99, pointData={"Density":dens})
        return lrgm_value[0]

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
        factor = (self.p.q**2)/(4* pi* self.eps0 * self.epsr)
        hartree = factor * fftconvolve(target, self.kernel, mode='valid')
        return hartree

    def expand_to_spinspace(self):
        from aux import spingraph_from_graph
        self.graph = spingraph_from_graph(self.graph)
        for contact in self.contacts:
            contact.graph = spingraph_from_graph(contact.graph)
            contact.expand_to_spin_space()

    def set_current(self,mode='unpolarized'):
        for contact in self.contacts:
            contact.current = mode

    def setmode(self,mode='normal'):
        if mode == 'normal':
            if ('graph' in dir(self)): del self.graph
            self.p.multi = 1
            self.order = 'even'
            self.block_sizes=[self.wafer.shape[1]*self.p.multi]*self.wafer.shape[0]
            self.simpleH()
        elif mode == 'spin':
            if ('graph' in dir(self)): del self.graph
            self.p.multi = 2
            self.order = 'even'
            self.block_sizes=[self.wafer.shape[1]*self.p.multi]*self.wafer.shape[0]
            self.spinH()
        elif mode == 'graph':
            self.p.multi = 1
            self.order = 'even'
            self.generate_graph()
            print 'Generate balanced levelstructure'
            self.generate_balanced_levelstructure()
            print 'Generate Hamiltonian from Graph'
            self.hamiltonian_from_graph()
            print 'Update Hamiltonian Diagonal'
            self.update_hamil_diag()
        elif mode == 'spin_graph':
            self.p.multi = 2
            self.order = 'odd'
            self.generate_graph()
            print 'Expanding to Spinspace'
            self.expand_to_spinspace()
            print 'Generate balanced levelstructure'
            self.generate_balanced_levelstructure()
            print 'Generate Hamiltonian from Graph'
            self.hamiltonian_from_graph()
            print 'Update Hamiltonian Diagonal'
            self.update_hamil_diag()
