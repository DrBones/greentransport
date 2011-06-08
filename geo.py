# -*- coding: utf-8 -*-
"""
Import of geometry to build Hamiltonian

The X-axis is the row Y-Axis the column, X=Down; Y=Right
"""
import scipy
from PIL import Image
#import matplotlib.pyplot as plt
import time
from scipy.sparse import lil_matrix,eye
#from scipy.sparse.linalg import splu
#from scipy.sparse.linalg import spsolve
#from scipy.linalg import solve, norm
import scipy.linalg as sl
from collections import OrderedDict
#from itertools import islice
scipy.set_printoptions(precision=3,suppress=True)
import customiterators
from parameters import Parameters
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
        header.append("DIMENSIONS {0} {1} 1".format(arr.shape[1], arr.shape[0]))
        header.append("\nSPACING 1 1 1\nORIGIN 0 0 0\nPOINT_DATA {0}".format(arr.shape[0] * arr.shape[1]))
        header.append("\nSCALARS EDensity double 1\nLOOKUP_TABLE default\n")
        with open(atlas + '.vtk', 'w') as file:
            for line in header:
                file.write(line)
        with open(atlas + 'ldos.vtk', 'w') as fileldos:
            for line in header:
                fileldos.write(line)
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
        for item,next_item in customiterators.pairwise(self.conductor):
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
        HD.setdiag(self.potential_serialized+4*t+Parameters.zplus)
        self.HD = HD

class Contact:
    """Defines all the methods and attributes needed to evaluate the
    effects of Contacts via the self-energy"""

    def greensfunction_contact(self, ind_contact, contact_node,E):
        """calculates the value of the transverse mode's square at the point
        of the contact_node, essentially sin(pi)sin(pi) multiplied with
        appropriate amplitude and constats """
        Length = (ind_contact.shape[1]-1.0)
        #Amplitude = scipy.sqrt(1/Length)
        Amplitude = 1
        ka = scipy.arccos(1-E.real/(2*Parameters.t))
        Phase = scipy.exp(1j*ka)
        xi = Amplitude * scipy.sin(scipy.pi * contact_node/Length)
        #xi = 1
        greensfunction = (xi**2) * Phase
        return greensfunction

    def build_sigma(self, nodes, ind_contact, E=Parameters.Egrid[0]):
        """Build the self-energy matrix SIGMA by determining the nodes
        adjacent to a Contact and inserting the greensfunction times t**2"""
        sigma = lil_matrix((len(nodes), len(nodes)), dtype=scipy.complex128)
        for contact_node in range(ind_contact.shape[1]):
            index = nodes[tuple(ind_contact.T[contact_node])][0]
            sigma[index, index] = - Parameters.t * self.greensfunction_contact(ind_contact, contact_node, E)
        return sigma

def fermifunction(E, kT=Parameters.kT, mu=Parameters.Efermi):
    """ Simple Fermifunction """
    fermifnc = 1/(scipy.exp((E.real-mu)/kT)+1)
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
    for i,current_block, next_block in customiterators.neighbour_zero(block_sizes):
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

def lesserfy(selfenergy_matrix, E, mu):
    """ Calculates the SIGMA-LESSER for both contacts as defined in
    Datta FATT p.227. (caution there SIGMA-IN) """
    gamma = 1j*(selfenergy_matrix - selfenergy_matrix.getH())
    lesser_matrix = 1j*gamma *fermifunction(E,mu)
    return lesser_matrix

def LRGM(block, Ablock, grl,Gr,sigma_in_l, sigma_in_r, E=Parameters.Egrid[0] ):
    """ Performs recursive algorithm (Svizhenko et.al) to calculate
    lesser green's function by the use of the diagonal and offdiagonal
    matrix-elements of the retarded green's function """
    gll = [grl[0] * sigma_in_l * grl[0].getH()]
    for i in range(1,len(block)-1): #len(block)-1 equals N-2 because of pythons way of building the range exclusive the endpoint
        prev_lesser_greenfnc = grl[i] * (Ablock[i,i-1] * gll[i-1] * Ablock[i,i-1].getH()) * grl[i].getH()
        gll.append(prev_lesser_greenfnc)
    i +=1
    gll.append(grl[i] * (sigma_in_r + Ablock[i,i-1] * gll[i-1] * Ablock[i-1,i].conj()) * grl[i].getH()) #N-1 step extra, belongs to for loop
    Gl = {}
    Gl[len(block)-1,len(block)-1] = gll[-1]
    for i in reversed(range(1,len(block))):
        Gl[i,i-1] = -Gr[i,i] * Ablock[i,i-1] * gll[i-1] - Gl[i,i] * Ablock[i-1,i].getH() * grl[i-1].getH()
        Gl[i-1,i-1] = (gll[i-1] + grl[i-1] * (Ablock[i-1,i] * Gl[i,i] * Ablock[i-1,i].getH()) * grl[i-1].getH() - (gll[i-1] * Ablock[i,i-1].getH() * Gr[i-1,i].getH() + Gr[i-1,i] * Ablock[i,i-1] * gll[i-1]))
    return Gl


timeittook=time.time()
d = Device()
c = Contact()
#Sb, sigma_block = blocks(cond, sigma[0]+sigma[1])
def runrrgm():
    ldos = scipy.array([0]*len(d.nodes))
    dos = scipy.array([0]*len(d.nodes))
    for energy in Parameters.Egrid:
        print energy
        tic = time.time()
        sigma_l = c.build_sigma(d.nodes, d.contact[0],energy - Parameters.potential_drop[0])
        sigma_r = c.build_sigma(d.nodes, d.contact[1], energy - Parameters.potential_drop[1])
        A = energy*scipy.sparse.eye(len(d.nodes),len(d.nodes),dtype=scipy.complex128, format='lil') - d.HD - sigma_l  - sigma_r
        block, Ablock = blocks(d.conductor, A)
        grl, Gr = RRGM(block, Ablock)
        print time.time() -tic
        tic = time.time()
        sigma_in_l = -2* sigma_l.imag[0:block[0],0:block[0]] * fermifunction(energy, mu=Parameters.mu_l)
        sigma_in_r = -2* sigma_r.imag[-block[-1]:,-block[0]:] * fermifunction(energy, mu=Parameters.mu_r)
        Gl = LRGM(block, Ablock, grl,Gr, sigma_in_l, sigma_in_r, energy)
        print time.time() -tic
        tic = time.time()
        green_diag = scipy.array([])
        for i in range(len(block)):
            diag = scipy.array(Gr[i,i].diagonal()).reshape(-1)
            green_diag = scipy.hstack((green_diag, diag))
        dos = dos + fermifunction(energy)*green_diag.imag*Parameters.dE/(-2*scipy.pi*Parameters.a)
        less_green_diag = scipy.array([])
        for i in range(len(block)):
            less_diag = scipy.array(Gl[i,i].diagonal()).reshape(-1)
            less_green_diag = scipy.hstack((less_green_diag, less_diag))
        ldos = ldos + less_green_diag.real*Parameters.dE/(scipy.pi*Parameters.a)
        #dos=dos+(fermifunction(energy)*Parameters.dE.real*green_diag.imag/(-2*scipy.pi*parameters.a))
        print time.time() -tic
    return dos,ldos, Gr, Gl

dos,ldos, Gr, Gl = runrrgm()


#lgmatrix = sl.block_diag(lg[0,0], lg[1,1], lg[2,2], lg[3,3], lg[4,4], lg[5,5], lg[6,6], lg[7,7])
#gmatrix = sl.block_diag(Gr[0,0], Gr[1,1], Gr[2,2], Gr[3,3], Gr[4,4], Gr[5,5], Gr[6,6], Gr[7,7])
#plt.imshow(scipy.absolute(scipy.asarray(d.HD().todense())))
#plt.imshow(scipy.asarray(HD.todense()).real)
#imshow(absolute(lgmatrix))
#imshow(-gmatrix.imag/2*scipy.pi)
#plt.show()

def writetofile(dos):
    with open(Parameters.atlas + '.vtk', 'a') as file:
        for x_pixel in range(Parameters.shape[0]):
            for y_pixel in range(Parameters.shape[1]):
                try:
                    file.write(str(dos[d.nodes[x_pixel,y_pixel][0]]) + "\n")
                except KeyError:
                    file.write('0\n')


def lwritetofile(dos):
    with open(Parameters.atlas + 'ldos.vtk', 'a') as file:
        for x_pixel in range(Parameters.shape[0]):
            for y_pixel in range(Parameters.shape[1]):
                try:
                    file.write(str(dos[d.nodes[x_pixel,y_pixel][0]]) + "\n")
                except KeyError:
                    file.write('0\n')
#writetofile(runrrgm())
print "The Process took", time.time()-timeittook, "seconds"
#Nodes2 = compose_geo2(Cond)
#def main():
#    Cont,Cond = read_geometry(atlas)            
#    Nodes = compose_geo(Cond)  
#
#if __name__ == '__main__':
#    main()

#    
