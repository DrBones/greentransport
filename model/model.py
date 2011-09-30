class Model:
    from aux import edens,spindens,transmission
    import networkx as nx

    def __init__(self, world):
        from scipy import linspace, cos, pi
        Model.canvas = world.canvas
        Model.atlas = world.atlas
        Model.hbar = world.hbar
        Model.m0 = world.m0
        Model.kb = world.kb

        #self.nodes = world.nodes
        self.wafer = world.wafer
        self.leads = world.leads
        self.contacts = world.contacts
        self.raw_coords = world.raw_coords
        self.active_coords = world.active_coords
        self.block_sizes = world.block_sizes
        self.eps0 = world.eps0
        self.q = world.q
        #Potential Drop over legth of device
        self.multi = 1
        self.a = 3e-9 # in meter
        self.alpha = 20e-12 # eV m
        #effective mass in eV real in GaAs 0.063
        self.mass = 0.036*Model.m0
        self.t0 = (Model.hbar**2)/(2*self.mass*(self.a**2))
        self.tso = self.alpha/(2 * self.a)
        #self.tso = 0.01*self.t0
        #self.tso = 0
        #Temperature * k_boltzmann in eV, 0.0025ev~30K
        self.epsr = 12.85
        self.Temp = 2 #in Kelvin
        self.kT = Model.kb * self.Temp
        self.lambdaf = 10 # i believe in nanometer, 35 more realistic?
        self.BField = 0 # in Tesla, from 0-~10
        self.zplus = 1j*1e-12
        self.band_bottom = 0
        #self.band_bottom = -4*self.t0
        #self.Efermi = self.band_bottom + 2*self.t0*(1-cos(2*pi/self.lambdaf))
        self.Efermi = 0.15*self.t0
        self.potential_drop = [0,0]
        #self.potential_drop = [0.004*self.t0/2, -0.004* self.t0/2]# in eV
        #self.Efermi = -3.8 * self.t0 # close to the bottom of the band at -4.0 t0, what bottom and band in what material ?
        #self.Egrid = linspace(self.Efermi-0.4*self.t0,self.Efermi +0.4*self.t0,100)+self.zplus # in eV ?
        self.mu = self.Efermi
        #self.dE = self.Egrid[1].real-self.Egrid[0].real
        #electro-chemical potential in eV
        self.mu_l = self.Efermi - (self.potential_drop[1] - self.potential_drop[0])/2
        self.mu_r = self.Efermi + (self.potential_drop[1] - self.potential_drop[0])/2
        #self.stepgrid(2,2)
        self.circular_qpc()
        #self.__generate_potential_grid()
        #self.grid2serialized(self.potential_grid)
        #self.__build_H()

    def stepgrid(self, step1,step2):
        from scipy import r_
        self.potential_grid = (r_[[self.potential_drop[0]]*self.wafer.shape[1]*step1,
            [0]*self.wafer.shape[1]*(self.wafer.shape[0]-step1-step2),
            [self.potential_drop[1]]*self.wafer.shape[1]*step2].reshape(self.wafer.shape))

    def circular_qpc(self,shift=0,radius=1,scale=0):
        from scipy import ogrid
        import aux
        size_x = self.wafer.shape[0]
        size_y = self.wafer.shape[1]
        x,y = ogrid[0:size_x:size_x*1j,0:size_y:size_y*1j]
        self.potential_grid = aux.sphericalPot(x,y,shift,radius,scale)

    def __generate_potential_grid(self):
        from scipy import linspace, tile
        rowdim = self.wafer.shape[0]
        coldim = self.wafer.shape[1]
        pot_slice = linspace(self.potential_drop[0],self.potential_drop[1],rowdim)
        potential = tile(pot_slice,(coldim,1)).T
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
        from scipy.sparse import lil_matrix, triu,bmat,eye,hstack
        from scipy import complex128, exp, pi
        Balpha = self.BField * self.a**2 /(2 * pi *self.hbar) # without the leading q because of hbar in eV
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
        Balpha = self.BField * self.a**2 /(2 * pi *self.hbar) # without the leading q because of hbar in eV
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

    def generate_graphs(self):
        # TODO make lazy initialized
        from aux import digraph_from_coords
        print 'Generating graph of device'
        self.graph,self.tuple_of_coords = digraph_from_coords(self,self.raw_coords)
        self.lead_graphs = []
        for lead in self.leads:
            print 'Generating graph of lead: ',lead.index
            lead_list = list(digraph_from_coords(self,lead))
            lead_list.append(lead.index)
            self.lead_graphs.append(lead_list)
        self.add_nodenames_to_contacts()

    def add_nodenames_to_contacts(self):
        for contact in self.contacts:
            contact.names = set()
            contact_tuple= tuple(zip(contact[0],contact[1]))
            for initial_tuple in contact_tuple:
                initial_name = self.tuple_of_coords.index(initial_tuple)
                contact.names.add(initial_name)
                #for partner_tuple in contact_tuple:
                #    partner_name = self.tuple_of_coords.index(partner_tuple)
                #    self.graph.add_edge(initial_name,partner_name)

    def expand_contacts_to_spin_space(self):
        for contact in self.contacts:
            list_of_nodenames = list(contact.names)
            contact.names=set()
            for contact_node in list_of_nodenames:
                contact.names.add(contact_node*2)
                contact.names.add(2*contact_node+1)

    def generate_balanced_levelstructure(self):
        from aux import BreadthFirstLevels, bisect
        #import pudb; pudb.set_trace()
        BFS_levelstructure = BreadthFirstLevels(self.graph,root=self.contacts[0].names,end=self.contacts[1].names)
        N = len(list(BFS_levelstructure))
        nodes_left  = self.contacts[0].names
        nodes_right = self.contacts[1].names
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
            if self.multi == 2:
                coordinate_array = repeat(self.tuple_of_coords,2,axis=0)
            else:
                coordinate_array = array(self.tuple_of_coords)

            #takes levelstructure permutations into account
            coordinates = coordinate_array[self.nodelist]
            for xy in coordinates:
                yield 4*self.t0+self.potential_grid[tuple(xy)]
        self.H.setdiag(list(serial_pot(self)))

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
        from scipy.sparse import lil_matrix,bmat,eye
        from scipy import argsort,where
        #from scipy.sparse.linalg import eigen
        transverseH = lil_matrix((self.wafer.shape[1],self.wafer.shape[1]))
        transverseH.setdiag([2*self.t0]*self.wafer.shape[1])
        transverseH.setdiag([-self.t0]*self.wafer.shape[1],1)
        transverseH.setdiag([-self.t0]*self.wafer.shape[1],-1)
#following is wrong
        #SO=eye(self.wafer.shape[1],self.wafer.shape[1],1)*self.tso-eye(self.wafer.shape[1],self.wafer.shape[1],-1)*self.tso
        #transverseHspin = bmat([[transverseH, SO],[SO,transverseH]])
        #self.HH = transverseHspin
        #from pudb import set_trace; set_trace()
        v,d = eig(transverseH.todense())
        ndx = argsort(v)
        d=d[:,ndx]
        v=v[ndx]
        self.v = v
        self.d = d
        try:
            self.maxmode = where(self.v < self.Efermi-self.band_bottom)[0].max()+1
        except ValueError:
            print "ValueError probably no modes will fit at that energy"
        if v.max() > self.Efermi-self.band_bottom:
            print 'Some mode energies larger than fermi energy, only up to mode {0} will fit'.format(self.maxmode)
            print 'Argument num_modes="all" takes only modes low enough'
            print ''

    def gamma(self, sigma):
        gamma = 1j*(sigma-sigma.conj())
        return gamma

    def generate_transverse_hamil(self,ind_contact):
        from numpy import zeros, eye
        contact_length = ind_contact.shape[1]

    def transfersigma(self,ind_contact,E):
        from scipy.linalg import eig,inv
        from scipy.sparse import lil_matrix
        from aux import SparseBlocks
        from scipy import argsort,dot,eye,hstack,vstack,zeros,complex128,asarray,split
        E= E+self.zplus
        Ndim = len(self.active_coords)
        block=ind_contact.shape[1]
        Hblock = SparseBlocks(self.H,self.block_sizes)
        self.Hblock = Hblock
        I=eye(block*self.multi)
        Zeros = zeros((block*self.multi,block*self.multi))
        #if ind_contact.SO is False:
        #    H00 = 4*self.t0*eye(block) -self.t0*eye(block,k=1) -self.t0*eye(block,k=-1)
        #    Hhop = -self.t0*I #no Bfield as of now
        #    inv_Hhop = -1/self.t0*I #no Bfield as of now
        if ind_contact.index == 0:
            H00 = asarray(Hblock[0,0].todense())
            #H10 = asarray(Hblock[1,0].todense())
            H01 = asarray(Hblock[0,1].todense())
            inv_H01 = inv(H01)
        if ind_contact.index == 1:
            H00 = asarray(Hblock[-1,-1].todense())
            #H10 = asarray(Hblock[-1,-2].todense())
            H01 = asarray(Hblock[-2,-1].todense())
            """indices switch because hopping matrices go the other
            direction x --> -x, better results this way altough not much difference"""
            inv_H01 = inv(H01)
        TransMatrix =vstack((hstack((dot(inv_H01,E*I-H00),dot(-inv_H01,H01.conj().T))),hstack((I,Zeros))))
        v,S = eig(TransMatrix)
        ndx = argsort(abs(v))
        S=S[:,ndx]
        v=v[ndx]
        self.S=S
        self.v =v
        Sleft,Sright = split(S,2,axis=1)
        S4,S3 = split(Sright,2,axis=0)
        S2,S1 = split(Sleft,2,axis=0)
        #S2 =S[:block*self.multi,:block*self.multi]
        #S1= S[block*self.multi:,:block*self.multi]
        self.S2 = S2
        self.S1 = S1
        self.S4 = S4
        self.S3 = S3
        print 'S1 shape: ',S1.shape
        print 'S2 shape: ',S2.shape
        if ind_contact.index == 0:
            dotted = dot(S2,inv(S1))
        if ind_contact.index == 1:
            dotted = dot(S3,inv(S4))
        invBracket =inv(E*I-H00-dot(H01,dotted))
        SigmaRet=self.t0**2*invBracket
        if ind_contact.index == 0:
            self.SigmaRet1 = SigmaRet
            #temp = zeros((60,60),dtype=complex128)
            #temp[30:60,30:60] =SigmaRet[:30,:30]
            #SigmaRet=temp
        else :
            self.SigmaRet2 = SigmaRet
        print 'SigaRet shape: ',SigmaRet.shape
        sigma = lil_matrix((self.multi*Ndim,self.multi*Ndim), dtype=complex128)
        print 'sigma shape: ',sigma.shape
        if ind_contact.index == 0:
            sigma[0:SigmaRet.shape[0], 0:SigmaRet.shape[1]] = SigmaRet
            self.sigma1=sigma
        elif ind_contact.index == 1:
            sigma[-SigmaRet.shape[0]:, -SigmaRet.shape[1]:] = SigmaRet
            self.sigma2=sigma
        import pudb; pudb.set_trace()
        return sigma

    def sigma_from_lead_graph(self,lead_graph,E):
        from scipy.linalg import inv,schur
        from scipy.sparse import lil_matrix
        from aux import SparseBlocks, eigenvector_from_eigenvalue, all_elements_are_unique
        from scipy import argsort,dot,eye,hstack,vstack,zeros,complex128,split,asarray,diag,array
        # should be able to turn zplus off, but then i need better eigenvalue comparison
        E= E+self.zplus
        graph, coords, index = lead_graph
        Ndim = len(self.active_coords)
        block=graph.order()/2
        I=eye(block*self.multi)
        Zeros = zeros((block*self.multi,block*self.multi))
        Hlead = self.nx.to_numpy_matrix(graph,dtype=complex128)
        # add 4*self.t0 because the matrix from graph lacks diagonal (no self-loops)
        H00 = asarray(Hlead[:block,:block])+4*self.t0 * I
        H01 = asarray(Hlead[:block,block:])
        inv_H01 = inv(H01)
        CompanionMatrix_array =vstack((hstack((Zeros,I)),hstack((dot(-inv_H01,H01.conj().T),dot(inv_H01,E*I-H00)))))
        #CompanionMatrix_array =vstack((hstack((dot(inv_H01,E*I-H00),dot(-inv_H01,H01.conj().T))),hstack((I,Zeros))))
        # the following 'complex' might be superfluous and only for real input matrices.
        T,Z = schur(CompanionMatrix_array,output='complex')
        eigenvalues = diag(T)
        propagating_eigenvalues = []
        propagating_eigenvectors = []
        for eigenvalue in eigenvalues:
            if abs(abs(eigenvalue)-1) < 0.01:
                propagating_eigenvalues.append(eigenvalue)
                eigenvector = eigenvector_from_eigenvalue(CompanionMatrix_array, eigenvalue)
                propagating_eigenvectors.append(eigenvector)
        prop_eig_array = array(propagating_eigenvectors).T
        if not all_elements_are_unique:
            print "--------------WARNING!!!!!---------------"
            print "One or more eigenvalues are identical, please rotate eigenvectors, I don't know how to do that"
        # sort eigenvalues and Z according to acending abs(eigenvalue), TODO: do better sorting, now it
        # depends on luck. sort using the calulated eigenvectors above
        sorting_indices = abs(eigenvalues).argsort()
        #T = T[:,sorting_indices][sorting_indices,:]
        #Z = Z[:,sorting_indices][sorting_indices,:]
        Zleft,Zright = split(Z,2,axis=1)
        #S4,S3 = split(Sright,2,axis=0)
        Z11,Z21 = split(Zleft,2,axis=0)
        SigmaRet = dot(H01,dot(Z21,inv(Z11)))
        import pudb; pudb.set_trace()
        if index == 0:
            self.SigmaRet1 = SigmaRet
        else :
            self.SigmaRet2 = SigmaRet
        print 'SigaRet shape: ',SigmaRet.shape
        sigma = lil_matrix((self.multi*Ndim,self.multi*Ndim), dtype=complex128)
        print 'sigma shape: ',sigma.shape
        if index == 0:
            sigma[0:SigmaRet.shape[0], 0:SigmaRet.shape[1]] = SigmaRet
            self.sigma1=sigma
        elif index == 1:
            sigma[-SigmaRet.shape[0]:, -SigmaRet.shape[1]:] = SigmaRet
            self.sigma2=sigma
        return sigma

    def sigma(self, ind_contact, E, num_modes='all'):
        """
        Takes abolute energies from band_bottom to around Efermi and further
        until the fermifunction puts and end to this
        """
        from scipy import sqrt,exp,asarray,dot,diag,where
        from scipy.sparse import lil_matrix
        from scipy import complex128
        if num_modes == 'analytical' :
            return self.analytic_sigma(ind_contact,E)
        elif num_modes == 'all':
            num_modes =self.maxmode
        Ndim = self.canvas[0]*self.canvas[1]
        #num_modes = where(self.v < E)[0].max()
        print 'Number of Modes used: ', num_modes
        #dd = diag(-self.t0*exp(1j*sqrt((E-self.band_bottom-self.v[:num_modes])/self.t0)))
        dd = diag(-self.t0*exp(1j*sqrt((E-self.v[:num_modes])/self.t0)))
        print 'Energy in Sigma used: ',E-self.v[:num_modes]
        SigmaRet  = asarray(dot(dot(self.d[:,:num_modes],dd),self.d[:,:num_modes].T))
        sigma = lil_matrix((self.multi*Ndim,self.multi*Ndim), dtype=complex128)
        if ind_contact[0].min() == 0:
            self.SigmaRet1 = SigmaRet
            sigma[0:self.wafer.shape[1], 0:self.wafer.shape[1]] = SigmaRet
            if self.multi==2:
                sigma[self.wafer.shape[1]:self.wafer.shape[1]*self.multi, self.wafer.shape[1]:self.wafer.shape[1]*self.multi] = SigmaRet
        elif ind_contact[0].max() == self.wafer.shape[0]-1:
            #import pudb; pudb.set_trace()
            self.SigmaRet2 = SigmaRet
            sigma[-self.block_sizes[-1]:, -self.block_sizes[-1]:] = SigmaRet
            if self.multi == 2:
                sigma[-self.block_sizes[-1]:-self.block_sizes[-1]/2, -self.block_sizes[-1]:-self.block_sizes[-1]/2] = SigmaRet
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
        number_of_nodes = len(self.active_coords)
        #global sigma_l,sigma_r,sigma_in_l,sigma_in_r
        #E_tot=self.Efermi+E_rel
        E_tot=E_rel
        if (not ('graph' in dir(self))):
            print "Using eigenvalue decomp of transfermatrix of blocks to build selfenergy"
            if (not ('lastenergy' in dir(self)) or self.lastenergy != E_rel):
                print 'Generating new selfenergy matrices'
                sigma_l = self.transfersigma(self.contacts[0],E_tot - self.potential_drop[0])
                self.sigma_l = sigma_l
                #sigma_l = self.transfersigma(self.contacts[0], E_tot)
                sigma_r =self.transfersigma(self.contacts[1], E_tot - self.potential_drop[1])
                self.sigma_r = sigma_r
                #sigma_r =self.transfersigma(self.contacts[1], E_tot)
                self.gamma_l = self.gamma(sigma_l)
                self.gamma_r = self.gamma(sigma_r)
                sigma_in_l = -2* sigma_l.imag[0:self.block_sizes[1], 0:self.block_sizes[1]] * self.fermifunction(E_tot, mu=self.mu_l)
                self.sigma_in_l = sigma_in_l
                sigma_in_r = -2* sigma_r.imag[-self.block_sizes[-1]:,-self.block_sizes[-1]:] * self.fermifunction(E_tot, mu=self.mu_r)
                self.sigma_in_r = sigma_in_r
            else:
                sigma_l = self.sigma_l
                sigma_r = self.sigma_r
                sigma_in_l = self.sigma_in_l
                sigma_in_r = self.sigma_in_r
        else:
            print "Using schur decomp of lead_graphs to build selfenergy"
            if (not ('lastenergy' in dir(self)) or self.lastenergy != E_rel):
                print 'Generating new selfenergy matrices'
                sigma_l = self.sigma_from_lead_graph(self.lead_graphs[0],E_tot)
                self.sigma_l = sigma_l
                sigma_r =self.sigma_from_lead_graph(self.lead_graphs[1], E_tot)
                self.sigma_r = sigma_r
                self.gamma_l = self.gamma(sigma_l)
                self.gamma_r = self.gamma(sigma_r)
                sigma_in_l = -2* sigma_l.imag[0:self.block_sizes[1], 0:self.block_sizes[1]] * self.fermifunction(E_tot, mu=self.mu_l)
                self.sigma_in_l = sigma_in_l
                sigma_in_r = -2* sigma_r.imag[-self.block_sizes[-1]:,-self.block_sizes[-1]:] * self.fermifunction(E_tot, mu=self.mu_r)
                self.sigma_in_r = sigma_in_r
            else:
                sigma_l = self.sigma_l
                sigma_r = self.sigma_r
                sigma_in_l = self.sigma_in_l
                sigma_in_r = self.sigma_in_r

        I =eye(number_of_nodes*self.multi,number_of_nodes*self.multi,dtype=complex128, format='lil')
        self.lastenergy = E_rel
        print 'sigmainl',sigma_in_l.shape
        print 'sigmainr',sigma_in_r.shape
        print 'I',I.shape
        print 'H', self.H.shape
        A = (E_tot+self.zplus)*I - self.H - sigma_l  - sigma_r
        return A, sigma_in_l, sigma_in_r

    def adaptiveenergy(self):
        pass

    def dorrgm(self,E_rel):
        from aux import SparseBlocks
        from greensolver import rrgm
        #from io_spinr import writeVTK
        #energy = self.Efermi
        A, sigma_in_l, sigma_in_r = self.spinA(E_rel)
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

    def dolrgm(self,energy):
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
        factor = (self.q**2)/(4* pi* self.eps0 * self.epsr)
        hartree = factor * fftconvolve(target, self.kernel, mode='valid')
        return hartree

    def expand_to_spinspace(self):
        from aux import spingraph_from_graph
        self.expand_contacts_to_spin_space()
        self.graph = spingraph_from_graph(self.graph,self.tso)
        for lead_graph in self.lead_graphs:
            lead_graph[0] = spingraph_from_graph(lead_graph[0], self.tso)

    def setmode(self,mode='normal'):
        if mode == 'normal':
            if ('graph' in dir(self)): del self.graph
            self.multi = 1
            self.order = 'even'
            self.block_sizes=[self.wafer.shape[1]*self.multi]*self.wafer.shape[0]
            self.simpleH()
        elif mode == 'spin':
            if ('graph' in dir(self)): del self.graph
            self.multi = 2
            self.order = 'even'
            self.block_sizes=[self.wafer.shape[1]*self.multi]*self.wafer.shape[0]
            self.spinH()
        elif mode == 'graph':
            self.multi = 1
            self.order = 'even'
            self.generate_graphs()
            self.generate_balanced_levelstructure()
            self.hamiltonian_from_graph()
            self.update_hamil_diag()
        elif mode == 'spin_graph':
            self.multi = 2
            self.order = 'odd'
            self.generate_graphs()
            self.expand_to_spinspace()
            self.generate_balanced_levelstructure()
            self.hamiltonian_from_graph()
            self.update_hamil_diag()
