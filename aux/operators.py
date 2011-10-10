def spindens(self,lrgm_out):
    from scipy import split,pi,complex128,zeros,repeat
    number_of_lattice_points = self.canvas.shape[0]*self.canvas.shape[1]
    number_of_nodes = len(self.tuple_canvas_coordinates)
    if number_of_nodes  == number_of_lattice_points:
        if self.order == 'even':
            Gup, Gdown = split(lrgm_out.reshape(self.canvas.shape[0],self.canvas.shape[1]*2),2,axis=1)
        if self.order == 'odd':
            Gup, Gdown =  split(lrgm_out.reshape(-1,2),2,axis=1)
            Gup, Gdown =  Gup.reshape(self.canvas.shape), Gdown.reshape(self.canvas.shape)
        else:
            print "Please specify order of Nodes, i.e 'even' for allspinup-allspindown per sclice or odd for spinup-spindown-spinup-..."
        Sz = self.p.hbar/(4*pi*1j*self.p.a**2)*(Gup-Gdown)
    elif number_of_nodes < number_of_lattice_points:
            Sz= zeros(self.canvas.shape,dtype=complex128)
            print 'calculating spin density for sparse structure'
            lrgm_out = self.p.hbar/(4*pi*1j*self.p.a**2)*lrgm_out
            expanded_array_of_coords = repeat(self.tuple_canvas_coordinates,2,axis=0)
            for index,node_name in enumerate(self.nodelist):
                if node_name % 2 == 0:
                    sign = 1
                else:
                    sign = -1
                Sz[tuple(expanded_array_of_coords[node_name])] += sign * lrgm_out[index]
    else:
        print 'Number of nodes larger than canvas, something is wrong!'
    print 'max Spin Split: ', Sz.real.max()-Sz.real.min()
#Realteil scheint wahrscheinlicher, imag oszilliert wie bloed
    return Sz.real

def edens(self,lrgm_out):
    from scipy import split,pi,asarray,sum
# TODO implement sparse and noncoliear structure mapping
    number_of_lattice_points = self.canvas.shape[0]*self.canvas.shape[1]
    number_of_nodes = len(self.tuple_canvas_coordinates)
    if number_of_nodes  == number_of_lattice_points:
        print 'Using stride bases reshape, Attention !!! Probably not what you want!'
        if self.p.multi ==1:
                edensity = lrgm_out.reshape(asarray(self.canvas.shape)+[0,0])*2/(2*pi*self.p.a**2) #times 2 for spin
        if self.p.multi ==2:
            if self.order == 'even':
                Gup, Gdown = split(lrgm_out.reshape(self.canvas.shape[0],self.canvas.shape[1]*2),2,axis=1)
                edensity = 1/(2*pi*self.p.a**2)*(Gup+Gdown)
            if self.order == 'odd':
                edensity = 1/(2*pi*self.p.a**2)*sum(lrgm_out.reshape(-1,2), axis=1).reshape(self.canvas.shape)
            else:
                print "Please specify order of Nodes, i.e 'even' for allspinup-allspindown per sclice or odd for spinup-spindown-spinup-..."
    elif number_of_nodes < number_of_lattice_points:
        from scipy import zeros,complex128,repeat
        edensity_spin_up = zeros(self.canvas.shape,dtype=complex128)
        edensity_spin_down = zeros(self.canvas.shape,dtype=complex128)
        if self.p.multi ==1:
            print 'calculating electron density without SO'
            lrgm_out = lrgm_out *2/(2*pi*self.p.a**2)
            for index,node_name in enumerate(self.nodelist):
                edensity_spin_up[self.tuple_canvas_coordinates[node_name]]= lrgm_out[index]
        if self.p.multi ==2:
            print 'calculating electron density with SO'
            lrgm_out = lrgm_out*1/(2*pi*self.p.a**2)
            expanded_array_of_coords = repeat(self.tuple_canvas_coordinates,2,axis=0)
            for index,node_name in enumerate(self.nodelist):
                if node_name % 2 ==0: 
                    edensity_spin_up[tuple(expanded_array_of_coords[node_name])] = lrgm_out[index]
                else:
                    edensity_spin_down[tuple(expanded_array_of_coords[node_name])] = lrgm_out[index]
        edensity = edensity_spin_up + edensity_spin_down
        edensity = (edensity.real, edensity_spin_up.real , edensity_spin_down.real)
    else:
        print 'Number of nodes larger than canvas, something is wrong!'
    # print 'max Electron Density: ', edensity.max()
    return edensity

def transmission(self,rrgm_out):
     from scipy import matrix,trace
     last_element_index = len(self.block_sizes)-1
     G_pq = matrix(rrgm_out[last_element_index,0])
     T = trace(matrix(-2*self.SigmaRet2.imag)*G_pq*matrix(-2*self.SigmaRet1.imag)*G_pq.getH())
     return T
