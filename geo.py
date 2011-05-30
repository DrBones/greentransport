# -*- coding: utf-8 -*-
"""
Import of geometry to build Hamiltonian

The X-axis is the row Y-Axis the column, X=Down; Y=Right
"""
import scipy
from PIL import Image
#import matplotlib.pyplot as plt
import time
from itertools import tee,izip_longest
from scipy.sparse import lil_matrix,eye
#from scipy.sparse.linalg import splu
#from scipy.sparse.linalg import spsolve
#from scipy.linalg import solve, norm
import scipy.linalg as sl
from collections import OrderedDict

#from itertools import islice
scipy.set_printoptions(precision=3,suppress=True)

def neighbour_zero(iterable):
    iterator = iter(iterable)
    prev = 0
    item = iterator.next()
    for next in iterator:
        yield (prev,item,next)
        prev = item
        item = next
    #yield (prev,item,0)

def pairwise(seq):
    a,b = tee(seq)
    b.next()
    return izip_longest(a,b)

class Parameters:
    """Contains all needed physical constants and parameters"""

    atlas = 'bar500x1onblack.bmp'
    potential_drop = [-0.05, 0.05] #Potential Drop over legth of device
    q = 1.6e-19
    hbar = 1.0545e-34/q                         #eV.s
    a = 2e-10                                  #mesh distance in meter
    m0 = 9.1e-31;                              #kg
    m0 = 0.510e6/((3e8)**2)                     #restmass of electron in eV
    mass = 0.25*m0                                  #effective mass in eV
    t = (hbar**2)/(2*mass*(a**2))
    mu = 0.25                                    #electro-chemical potential in eV
    kT = 0.025                                  #Temperature * k_boltzmann in eV, 0.0025ev~30K
    lambdaf = 10
    BField = 0
    zplus = 1j*1e-12
    Egrid = scipy.linspace(-0.1,0.4,250)+zplus
    dE = Egrid[1]-Egrid[0]

class Device:
    """Contains all methods relating to the device in particular and
    implicit like the defining geometry and the building of the
    Hamiltonian HD"""

    def __init__(self):
        """ Initilizes the device object and executes read_ and compose
        geometry"""
        self.read_geometry()
        self.compose_nodes()
        self.build_HD()

    def read_geometry(self, atlas=Parameters.atlas):
        img = Image.open(atlas)
        arr = scipy.asarray(img)
        contact = []
        contact_shades = [149, 179, 209, 239]
        for shade in contact_shades:
            a = scipy.array(scipy.where(arr == shade))
            if a.shape[1] == 0: continue
            contact.append(a)
        conductor = scipy.transpose(scipy.where(arr > 0))
        header = []
        header.append("# vtk DataFile Version 2.0\nVTK Data of Device\nASCII\nDATASET STRUCTURED_POINTS\n")
        header.append("DIMENSIONS {0} {1} 1".format(arr.shape[0], arr.shape[1]))
        header.append("\nSPACING 1 1 1\nORIGIN 0 0 0\nPOINT_DATA {0}".format(arr.shape[0] * arr.shape[1]))
        header.append("\nSCALARS EDensity double 1\nLOOKUP_TABLE default\n")
        with open('geo.vtk', 'w') as file:
            for line in header:
                file.write(line)
        potential = []
        for i in range(arr.shape[1]):
            potential.append(scipy.vstack(scipy.r_[Parameters.potential_drop[0]:Parameters.potential_drop[1]:arr.shape[0]*1j]))
        potential = scipy.hstack(tuple(potential))

        Parameters.shape = arr.shape
        self.potential_grid = potential
        self.contact = contact
        self.conductor = conductor

    def compose_nodes(self):
        nodes = OrderedDict()
        count = 0
        for item,next_item in pairwise(self.conductor):
            try:
                if item[1] == next_item[1]-1:
                    nodes[tuple(item)] = [count,count+1,None]
                else:
                    nodes[tuple(item)] = [count,None,None]
            except TypeError:
                nodes[tuple(item)] = [count,None,None]
            if item[0]>1:
                try:
                    nodes[tuple(item - [1,0])][2] = count
                except KeyError:
                    pass
            count +=1
        potential_serialized = scipy.array([])
        for key, value in nodes.items():
           potential_serialized = scipy.concatenate((potential_serialized, [self.potential_grid[key[0],key[1]]]))
        self.potential_serialized = potential_serialized
        self.nodes = nodes

    def build_HD(self, t=Parameters.t,BField=Parameters.BField):
        HD = lil_matrix((len(self.nodes),len(self.nodes)),dtype=scipy.complex128)
        for item in self.nodes:
            if self.nodes[item][1] != None:
                HD[self.nodes[item][0],self.nodes[item][1]] = -scipy.exp(1j*BField*self.nodes[item][0])*t
            else: pass
            if self.nodes[item][2] == None: continue
            HD[self.nodes[item][0],self.nodes[item][2]] = -t
        #HD = HD + HD.conjugate().T
        HD = (HD.tocsr() + HD.tocsr().conjugate().T).tolil()# might be faster
        HD.setdiag(self.potential_serialized+2*t+Parameters.zplus)
        self.HD = HD

class Contact:
    """Defines all the methods and attributes needed to evaluate the
    effects of Contacts via the self-energy"""

    def greensfunction_contact(self, ind_contact, contact_node,E):
        """calculates the value of the transverse mode's square at the point
        of the contact_node, essentially sin(pi)sin(pi) multiplied with
        appropriate amplitude and constats """
        #Length = (ind_contact.shape[1]-1.0)
        #Amplitude = scipy.sqrt(1/Length)
        ka = scipy.arccos(1-E/(2*Parameters.t))
        Phase = scipy.exp(1j*ka)
        #xi = Amplitude * scipy.sin(scipy.pi * contact_node/Length)
        xi = 1
        greensfunction = (xi**2) * Phase
        return greensfunction

    def build_SIGMA(self, nodes, contact, E=Parameters.Egrid[0]):
        """Build the self-energy matrix SIGMA by determining the nodes
        adjacent to a Contact and inserting the greensfunction times t**2"""
        contact_index = 0
        SIGMA = []
        for ind_contact in contact:
            SIGMA.append(lil_matrix((len(nodes), len(nodes)), dtype=scipy.complex128))
            for contact_node in range(ind_contact.shape[1]):
                index = nodes[tuple(ind_contact.T[contact_node])][0]
                SIGMA[contact_index][index, index] = - Parameters.t * self.greensfunction_contact(ind_contact, contact_node, E)
            contact_index +=1
        return SIGMA

def calc_Ef(lambdaf=Parameters.lambdaf , t=Parameters.t):
    Ef = 2*t*(scipy.cos(scipy.pi/lambdaf))
    return Ef

def fermifunction(E, kT=Parameters.kT, mu=calc_Ef()):
    """ Simple Fermifunction """
    fermifnc = 1/(scipy.exp((E-mu)/kT)+1)
    return fermifnc

def blocks(cond, matrix):
    block_size = 1
    block_sizes = []
    for i in range(1,len(cond)):
        if cond[i][0] == cond[i-1][0]:
            block_size +=1
        else:
            block_sizes.append(block_size)
            block_size = 1
    block_sizes.append(block_size)
    matrix_block = {}
    current_start = 0
    row_index = 0
    column_index = 0
    matrix_block[row_index, column_index] = matrix[0:block_sizes[0],0:block_sizes[0]]
    for i,current_block, next_block in neighbour_zero(block_sizes):
        next_start = current_start + current_block
        next_end = current_start + current_block +next_block
        matrix_block[row_index+1, column_index+1] = matrix[next_start:next_end,next_start:next_end]
        matrix_block[row_index, column_index+1] = matrix[current_start:next_start, next_start:next_end]
        matrix_block[row_index+1, column_index] = matrix[next_start:next_end, current_start:next_start]
        current_start +=current_block
        row_index +=1
        column_index +=1
    return block_sizes, matrix_block

def RRGM(blocks, Ablock):
    """ Performs recursive algorithm (Svizhenko et. al) to calculate the retarded green's
    function, uses views on A, i.e. the block matrices of A, by Ablock """
    grl = [Ablock[0,0].todense().I]
    prev_greensfnc = grl[0]
    for i in range(1,len(blocks)):
        prev_greensfnc = (Ablock[i,i]-Ablock[i, i-1] * prev_greensfnc * Ablock[i-1,i]).I
        grl.append(prev_greensfnc)
    Gr = {}
    Gr[len(blocks)-1,len(blocks)-1] = grl[-1]
    rev_iterator = reversed(range(1,len(blocks)))
    for i in rev_iterator:
        Gr[i, i-1] = -Gr[i,i] * Ablock[i,i-1] * grl[i-1]
        Gr[i-1, i] = -grl[i-1] * Ablock[i-1,i] * Gr[i,i]
        Gr[i-1, i-1] = grl[i-1]-grl[i-1] * Ablock[i-1,i] * Gr[i, i-1]
    return grl, Gr

def lesserfy(selfenergy_matrix, E):
    """ Calculates the SIGMA-LESSER for both contacts as defined in
    Datta FATT p.227. (caution there SIGMA-IN) """
    lesser_matrix = -(selfenergy_matrix - selfenergy_matrix.conj()) * fermifunction(E)
    return lesser_matrix

def LRGM(blocks, Ablock, grl,sigma_block, E=Parameters.Egrid[0] ):
    """ Performs recursive algorithm (Svizhenko et.al) to calculate
    lesser green's function by the use of the diagonal and offdiagonal
    matrix-elements of the retarded green's function """
    gll = [grl[0] * lesserfy(sigma_block[0,0],E) * grl[0].conj()]
    for i in range(1,len(blocks)-1):
        prev_lesser_greenfnc = grl[i] * (Ablock[i,i-1] * gll[i-1] * Ablock[i-1,i].conj()) * grl[i].conj()
        gll.append(prev_lesser_greenfnc)
    gll.append(grl[i+1] * (lesserfy(sigma_block[i+1,i+1], E) + Ablock[i+1,i] * gll[i] * Ablock[i,i+1].conj()) * grl[i+1].conj())
    Gl = {}
    Gl[len(blocks)-1,len(blocks)-1] = gll[-1]
    rev_iterator = reversed(range(1,len(blocks)))
    for i in rev_iterator:
        Gl[i,i-1] = Gr[i,i] * Ablock[i,i-1] * gll[i-1] + Gl[i,i] * Ablock[i,i-1].conj() * grl[i-1].conj()
        Gl[i-1,i-1] = (gll[i-1] + grl[i-1] * (Ablock[i-1,i] * Gl[i,i] * Ablock[i,i-1].conj()) * grl[i-1].conj() +
                                 (gll[i-1] * Ablock[i-1,i].conj() * Gr[i, i-1].conj() + Gr[i-1,i] * Ablock[i,i-1] * gll[i-1]))
    return Gl


timeittook=time.time()
d = Device()
c = Contact()
#Sb, sigma_block = blocks(cond, sigma[0]+sigma[1])
dos = scipy.array([0]*500)
tic = time.time()
for energy in Parameters.Egrid:
    A = energy*eye(len(d.nodes),len(d.nodes),dtype=scipy.complex128, format='lil') - d.HD - c.build_SIGMA(d.nodes, d.contact,energy - Parameters.potential_drop[0])[0] - c.build_SIGMA(d.nodes, d.contact, energy - Parameters.potential_drop[1])[1]
    block, Ablock = blocks(d.conductor, A)
    grl, Gr = RRGM(block, Ablock)
    green_diag = scipy.array([])
    for i in range(len(block)):
        diag = scipy.array(Gr[i,i].diagonal()).reshape(-1)
        green_diag = scipy.hstack((green_diag, diag))
    dos=dos+(fermifunction(energy)*Parameters.dE*green_diag.imag/(-2*scipy.pi*Parameters.a))
print time.time() -tic

mdos = scipy.array([0]*500)
tic = time.time()
for energy in Parameters.Egrid:
    A = energy*eye(len(d.nodes),len(d.nodes),dtype=scipy.complex128, format='lil') - d.HD - c.build_SIGMA(d.nodes, d.contact,energy - Parameters.potential_drop[0])[0] - c.build_SIGMA(d.nodes, d.contact, energy - Parameters.potential_drop[1])[1]
    green_diag = A.todense().I
    mdos=mdos+(fermifunction(energy)*Parameters.dE*scipy.diagonal(green_diag.imag)/(-2*scipy.pi*Parameters.a))
print time.time() -tic
#lg = LRGM(b, Ablock, grl, sigma_block)
#
#for i in range(len(b)):

#    print i
#    diag = scipy.array(lg[i,i].diagonal()).reshape(-1)
#    lesser_green_diag = scipy.hstack((green_diag, diag))

#lgmatrix = sl.block_diag(lg[0,0], lg[1,1], lg[2,2], lg[3,3], lg[4,4], lg[5,5], lg[6,6], lg[7,7])
#gmatrix = sl.block_diag(Gr[0,0], Gr[1,1], Gr[2,2], Gr[3,3], Gr[4,4], Gr[5,5], Gr[6,6], Gr[7,7])
#plt.imshow(scipy.absolute(scipy.asarray(d.HD().todense())))
#plt.imshow(scipy.asarray(HD.todense()).real)
#imshow(absolute(lgmatrix))
#imshow(-gmatrix.imag/2*scipy.pi)
#plt.show()

#with open('geo.vtk', 'a') as file:
#    for x_pixel in range(Parameters.shape[0]):
#        for y_pixel in range(Parameters.shape[1]):
#            try:
#                file.write(str(gmatrix.imag[nodes[x_pixel,y_pixel][0], nodes[x_pixel,y_pixel][0]]/(2*scipy.pi)) + "\n")
#            except KeyError:
#                file.write('0\n')
print "The Process took", time.time()-timeittook, "seconds"
#Nodes2 = compose_geo2(Cond)
#def main():
#    Cont,Cond = read_geometry(atlas)            
#    Nodes = compose_geo(Cond)  
#
#if __name__ == '__main__':
#    main()

#    
