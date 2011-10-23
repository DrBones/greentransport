class Contact(object):
    import aux

    def __init__(self, contact_temp_tuple_coords, interface_temp_tuple_coords, SO=None,):
        from parameters import p
        self.p = p
        # self.interface_raw_coordinates = None
        # self.contact_raw_coordinates = sc.asarray(contact_temp_coords)
        # the star unpacks the array into lists first and zip zips them to (x,y) tuples
        self.interface_tuple_coordinates = interface_temp_tuple_coords
        self.contact_tuple_coordinates = contact_temp_tuple_coords
        self.length = len(self.contact_tuple_coordinates)
        self.current = 'unpolarized'

    def recreate_graph(self):
        self.graph = Contact.aux.digraph_from_tuple_coords(self.interface_tuple_coordinates)
        print 'Adding Node Names to Interface: ', self.index
        self.add_nodenames()
        print 'Sorting Nodes of Interface: ', self.index
        self.transverse_sorted_nodelist_from_graph()

    def add_nodenames(self):
        self.internal_node_names = set()
        self.node_names = set()
        for initial_tuple in self.contact_tuple_coordinates:
            initial_name = self.p.tuple_canvas_coordinates.index(initial_tuple)
            self.node_names.add(initial_name)
            self.internal_node_names.add(self.interface_tuple_coordinates.index(initial_tuple))
            #for partner_tuple in contact_tuple:
            #    partner_name = self.tuple_of_coords.index(partner_tuple)
            #    self.graph.add_edge(initial_name,partner_name)

    def transverse_sorted_nodelist_from_graph(self):
        self.levelstructure = list(Contact.aux.BreadthFirstLevels(self.graph, self.internal_node_names))
        self.nodelist = [item for level in self.levelstructure for item in sorted(level)]

    def expand_to_spin_space(self):
        list_of_nodenames = list(self.node_names)
        self.node_names=set()
        for contact_node in list_of_nodenames:
            self.node_names.add(contact_node*2)
            self.node_names.add(2*contact_node+1)

        list_of_internal_nodenames = list(self.internal_node_names)
        self.internal_node_names=set()
        for contact_node in list_of_internal_nodenames:
            self.internal_node_names.add(contact_node*2)
            self.internal_node_names.add(2*contact_node+1)

        self.transverse_sorted_nodelist_from_graph()

    def generate_sigma(self):
        import networkx as nx
        from scipy.linalg import inv,schur
        from scipy.sparse import lil_matrix
        from numpy import asarray_chkfinite,isfinite,inf
        from aux import SparseBlocks, eigenvector_from_eigenvalue, all_elements_are_unique
        from scipy import argsort,dot,eye,hstack,vstack,zeros,complex128,split,asarray,diag,array
        # should be able to turn zplus off, but then i need better eigenvalue comparison
        E= self.p.E+self.p.zplus
        Ndim = len(self.p.tuple_canvas_coordinates)
        block=self.graph.order()/2
        I=eye(block)
        Zeros = zeros((block,block))
        Hlead = nx.to_numpy_matrix(self.graph,nodelist=self.nodelist,dtype=complex128)
        self.Hlead = Hlead
        #import pudb; pudb.set_trace()
        # add 4*self.t0 because the matrix from graph lacks diagonal (no self-loops)
        try:
            H00 = asarray_chkfinite(Hlead[:block,:block])+4*self.p.t0 * I
            self.H00 = H00
        except ValueError:
            print 'H00 contains infs or NaNs'
            import pudb; pudb.set_trace()
        try:
            H01 = asarray_chkfinite(Hlead[:block,block:])
            self.H01 = H01
        except ValueError:
            print 'H01 contains infs or NaNs'
            import pudb; pudb.set_trace()
        inv_H01 = inv(H01)
        while not isfinite(inv_H01).all():
            print 'inv_H01 contains infs or NaNs, repeating'
            inv_H01 = inv(H01)
            #import pudb; pudb.set_trace()
        self.inv_H01 = inv_H01
        CompanionMatrix_array =vstack((hstack((Zeros,I)),hstack((dot(-inv_H01,H01.conj().T),dot(inv_H01,E*I-H00)))))
        if not isfinite(CompanionMatrix_array).all():
            print 'CompanionMatrix contains infs or NaNs'
            import pudb; pudb.set_trace()
        #CompanionMatrix_array =vstack((hstack((dot(inv_H01,E*I-H00),dot(-inv_H01,H01.conj().T))),hstack((I,Zeros))))
        # the following 'complex' might be superfluous and only for real input matrices.
        T,Z,number_sorted= schur(CompanionMatrix_array,sort='iuc')
        eigenvalues = diag(T)
        # propagating_eigenvalues = []
        # propagating_eigenvectors = []
        # for eigenvalue in eigenvalues:
        #     if abs(abs(eigenvalue)-1) < 0.01:
        #         propagating_eigenvalues.append(eigenvalue)
        #         eigenvector = eigenvector_from_eigenvalue(CompanionMatrix_array, eigenvalue)
        #         propagating_eigenvectors.append(eigenvector)
        # prop_eig_array = array(propagating_eigenvectors).T
        if not all_elements_are_unique(eigenvalues):
            print "--------------WARNING!!!!!---------------"
            print "One or more eigenvalues are identical, please rotate eigenvectors, I don't know how to do that"
        # sort eigenvalues and Z according to acending abs(eigenvalue), TODO: do better sorting, now it
        # depends on luck. sort using the calulated eigenvectors above
        # sorting_indices = abs(eigenvalues).argsort()
        #T = T[:,sorting_indices][sorting_indices,:]
        #Z = Z[:,sorting_indices][sorting_indices,:]
        Zleft,Zright = split(Z,2,axis=1)
        #S4,S3 = split(Sright,2,axis=0)
        Z11,Z21 = split(Zleft,2,axis=0)
        SigmaRet = dot(H01,dot(Z21,inv(Z11)))
        self.SigmaRet = SigmaRet
        print '- SimgaRet (',self.index,') shape: ',SigmaRet.shape
        sigma = lil_matrix((self.p.multi*Ndim,self.p.multi*Ndim), dtype=complex128)
        print '- sigma shape: ',sigma.shape
#implement polarization like (for spin up) reshape(-1,2) --> [:,1] = 0 --> reshape(shape(SigmaRet))
        if self.index == 0:
            if 'up' in self.current:
                SigmaRet.reshape(-1,2)[:,1] = 0
            elif 'down' in self.current:
                SigmaRet.reshape(-1,2)[:,0] = 0
            sigma[0:SigmaRet.shape[0], 0:SigmaRet.shape[1]] = SigmaRet
        elif self.index == 1:
            sigma[-SigmaRet.shape[0]:, -SigmaRet.shape[1]:] = SigmaRet
        self.sigma=sigma
        return sigma

    def __contact_greensfunction(self, ind_contact, node_i, node_j, E):
        """calculates the value of the transverse mode's square at the point
        of the contact_node, essentially sin(pi)sin(pi) multiplied with
        appropriate amplitude and constats """
        from scipy import arccos, exp, sin, pi, sqrt
        length = (ind_contact.shape[1]-1.0)
        amplitude = 1/sqrt(length)
        mode_energy = self.p.hbar**2 * pi**2/(2*self.mass * (length*self.a)**2)
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
            print "- ValueError probably no modes will fit at that energy"
        if v.max() > self.Efermi-self.band_bottom:
            print 'Some mode energies larger than fermi energy, only up to mode {0} will fit'.format(self.maxmode)
            print 'Argument num_modes="all" takes only modes low enough'
            print ''

    def gamma(self, sigma):
        gamma = 1j*(sigma-sigma.conj())
        return gamma

    def transfersigma(self,ind_contact,E):
        from scipy.linalg import eig,inv
        from scipy.sparse import lil_matrix
        from aux import SparseBlocks
        from scipy import argsort,dot,eye,hstack,vstack,zeros,complex128,asarray,split
        E= E+self.zplus
        Ndim = len(self.p.tuple_canvas_coordinates)
        block=ind_contact.length
        Hblock = SparseBlocks(self.H,self.block_sizes)
        self.Hblock = Hblock
        I=eye(block*self.p.multi)
        Zeros = zeros((block*self.p.multi,block*self.p.multi))
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
        #S2 =S[:block*self.p.multi,:block*self.p.multi]
        #S1= S[block*self.p.multi:,:block*self.p.multi]
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
        sigma = lil_matrix((self.p.multi*Ndim,self.p.multi*Ndim), dtype=complex128)
        print 'sigma shape: ',sigma.shape
        if ind_contact.index == 0:
            sigma[0:SigmaRet.shape[0], 0:SigmaRet.shape[1]] = SigmaRet
            self.sigma1=sigma
        elif ind_contact.index == 1:
            sigma[-SigmaRet.shape[0]:, -SigmaRet.shape[1]:] = SigmaRet
            self.sigma2=sigma
        import pudb; pudb.set_trace()
        return sigma

    def sigma_deprec(self, ind_contact, E, num_modes='all'):
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
        sigma = lil_matrix((self.p.multi*Ndim,self.p.multi*Ndim), dtype=complex128)
        if ind_contact[0].min() == 0:
            self.SigmaRet1 = SigmaRet
            sigma[0:self.wafer.shape[1], 0:self.wafer.shape[1]] = SigmaRet
            if self.p.multi==2:
                sigma[self.wafer.shape[1]:self.wafer.shape[1]*self.p.multi, self.wafer.shape[1]:self.wafer.shape[1]*self.p.multi] = SigmaRet
        elif ind_contact[0].max() == self.wafer.shape[0]-1:
            #import pudb; pudb.set_trace()
            self.SigmaRet2 = SigmaRet
            sigma[-self.block_sizes[-1]:, -self.block_sizes[-1]:] = SigmaRet
            if self.p.multi == 2:
                sigma[-self.block_sizes[-1]:-self.block_sizes[-1]/2, -self.block_sizes[-1]:-self.block_sizes[-1]/2] = SigmaRet
        return sigma

    def analytic_sigma(self, ind_contact, E):
        from scipy.sparse import lil_matrix
        from scipy import complex128
        ystride = self.wafer.shape[1]
        Ndim = self.canvas[0]*self.canvas[1]
        sigma = lil_matrix((self.p.multi*Ndim,self.p.multi*Ndim), dtype=complex128)
        contact_node = 0
        for xypair in ind_contact.T:
            index = xypair[1] + ystride*xypair[0]
            sigma[index, index] = - self.t0 * self.__contact_greensfunction(ind_contact, contact_node, contact_node, E)
            if self.p.multi==1:contact_node+=1; continue
            sigma[index+ystride, index+ystride] = - self.t0 * self.__contact_greensfunction(ind_contact, contact_node, contact_node, E)
            contact_node +=1
        return sigma
