def spindens(self,lrgm_out):
    from scipy import split,pi
    if self.order == 'even':
        Gup, Gdown = split(lrgm_out.reshape(self.wafer.shape[0],self.wafer.shape[1]*2),2,axis=1)
    if self.order == 'odd':
        Gup, Gdown =  split(lrgm_out.reshape(-1,2),2,axis=1)
        Gup, Gdown =  Gup.reshape(self.wafer.shape), Gdown.reshape(self.wafer.shape)
    else:
        print "Please specify order of Nodes, i.e 'even' for allspinup-allspindown per sclice or odd for spinup-spindown-spinup-..."

    Sz = self.hbar/(4*pi*1j*self.a**2)*(Gup-Gdown)
    print 'max Spin Split: ', Sz.imag.max()-Sz.imag.min()
#Realteil scheint wahrscheinlicher, imag oszilliert wie bloed
    return Sz.real

def edens(self,lrgm_out):
    from scipy import split,pi,asarray,sum
# TODO implement sparse and noncoliear structure mapping
    number_of_lattice_points = self.wafer.shape[0]*self.wafer.shape[1]
    number_of_nodes = len(self.active_coords)
    if number_of_nodes  == number_of_lattice_points:
        if self.multi ==1:
                edensity = lrgm_out.reshape(asarray(self.wafer.shape)+[0,0])*2/(2*pi*self.a**2) #times 2 for spin
        if self.multi ==2:
            if self.order == 'even':
                Gup, Gdown = split(lrgm_out.reshape(self.wafer.shape[0],self.wafer.shape[1]*2),2,axis=1)
                edensity = 1/(2*pi*self.a**2)*(Gup+Gdown)
            if self.order == 'odd':
                edensity = 1/(2*pi*self.a**2)*sum(lrgm_out.reshape(-1,2), axis=1).reshape(self.wafer.shape)
            else:
                print "Please specify order of Nodes, i.e 'even' for allspinup-allspindown per sclice or odd for spinup-spindown-spinup-..."
    elif number_of_nodes < number_of_lattice_points:
        from scipy import zeros,complex128
        edensity = zeros(self.wafer.shape,dtype=complex128)
        if self.multi ==1:
            for node_name in self.nodelist:
                edensity[self.tuple_of_coords[node_name]]= lrgm_out[node_name]*2/(2*pi*self.a**2)
        if self.multi ==2:
            print 'multi =2'
            lrgm_out = sum(lrgm_out.reshape(-1,2), axis=1)
            for node_name in self.nodelist:
                edensity[self.tuple_of_coords[node_name]]= lrgm_out[node_name]*1/(2*pi*self.a**2)
    else:
        print 'Number of nodes larger than canvas, something is wrong!'
    print 'max Electron Density: ', edensity.max()
    return edensity.real

def transmission(self,rrgm_out):
     from scipy import matrix,trace
     last_element_index = len(self.block_sizes)-1
     G_pq = matrix(rrgm_out[last_element_index,0])
     T = trace(matrix(-2*self.SigmaRet2.imag)*G_pq*matrix(-2*self.SigmaRet1.imag)*G_pq.getH())
     return T
