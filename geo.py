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
    """Contains all nedded physical constants and parameters"""

    atlas = 'dot10x10onblack.bmp'
    Confinement = 0
    q = 1.6e-19
    hbar = 1.0545e-34/q                         #eV.s
    a = 2e-10                                  #mesh distance in meter
    m0 = 9.1e-31;                              #kg
    m0 = 0.510e6/((3e8)**2)                     #restmass of electron in eV
    mass = m0                                  #effective mass in eV
    t = (hbar**2)/(2*mass*(a**2))
    mu = 0.25                                    #electro-chemical potential in eV
    kT = 0.0025                                  #Temperature * k_boltzmann in eV, 0.0025ev~30K
    lambdaf = 10
    BField = 0


class Device:
    """Contains all methods relating to the device in particular and
    implicit like the defining geometry and the building of the
    Hamiltonian HD"""

    def read_geometry(self, atlas=Parameters.atlas):
        Img = Image.open(atlas)
        Arr = scipy.asarray(Img)
        Contact = []
        Contact.append(scipy.array(scipy.where(Arr == 149)))
        Contact.append(scipy.array(scipy.where(Arr == 179)))
        Contact.append(scipy.array(scipy.where(Arr == 209)))
        Contact.append(scipy.array(scipy.where(Arr == 239)))
        Conductor = scipy.transpose(scipy.where(Arr > 0))
        return Contact,Conductor



    def compose_geo(self, Cond=read_geometry(Parameters.atlas)[1]):
        Nodes = {}
        Count = 0
        for item,next_item in pairwise(Cond):
            try:
                if item[1] == next_item[1]-1:
                    Nodes[tuple(item)] = [Count,Count+1,None]
                else:
                    Nodes[tuple(item)] = [Count,None,None]
            except TypeError:
                Nodes[tuple(item)] = [Count,None,None]
            if item[0]>1:
                try:
                    Nodes[tuple(item - [1,0])][2] = Count
                except KeyError:
                    pass
            Count +=1
        return Nodes

    def build_HD(self, Nodes, Potential=Parameters.Confinement,t=Parameters.t,BField=Parameters.BField):
        HD = lil_matrix((len(Nodes),len(Nodes)),dtype=scipy.complex128)
        HD.setdiag([Potential+4*t]*len(Nodes))
        for item in Nodes:
            if Nodes[item][1] != None:
                HD[Nodes[item][0],Nodes[item][1]] = -scipy.exp(1j*BField*Nodes[item][0])*t
            else: pass
            if Nodes[item][2] == None: continue
            HD[Nodes[item][0],Nodes[item][2]] = -t
        #HD = HD + HD.conjugate().T
        HD = (HD.tocsr() + HD.tocsr().conjugate().T).tolil()# might be faster
        return HD

    def HD(self):
        """docstring for HD"""
        a, b = self.read_geometry()
        HD = self.build_HD(self.compose_geo(b))
        return HD

def build_EF(Nodes, lambdaf=Parameters.lambdaf , t=Parameters.t):
    EF = 2*t*(scipy.cos(scipy.pi/lambdaf))*eye(len(Nodes),len(Nodes),dtype=scipy.complex128, format='lil')
    return EF

class Contact:
    """Defines all the methods and attributes needed to evaluate the
    effects of Contacts via the self-energy"""

    def build_SIGMA(self, Nodes, contact):
        """Build the self-energy matrix SIGMA by determining the nodes
        adjacent to a Contact and inserting the greensfunction times t**2"""
        SIGMA = lil_matrix((len(Nodes), len(Nodes)), dtype=scipy.complex128)
        for ind_contact in contact:
            for contact_node in range(ind_contact.shape[1]):
                index = Nodes[tuple(ind_contact.T[contact_node])][0]
                SIGMA[index, index] = Parameters.t * self.greensfunction_contact(ind_contact, contact_node)
        return SIGMA

    def greensfunction_contact(self, ind_contact, contact_node):
        """calculates the value of the transverse mode's square at the point
        of the contact_node, essentially sin(pi)sin(pi) multiplied with
        appropriate amplitude and constats """
        Length = ind_contact.shape[1]*Parameters.a
        Amplitude = scipy.sqrt(1/Length)
        Phase = scipy.exp(1j*scipy.pi/Length*Parameters.a)
        xi = Amplitude * scipy.sin(scipy.pi * contact_node/Length)
        greensfunction = (xi**2) * Phase
        return greensfunction

def blocks(cond, A):
    block_size = 1
    block_sizes = []
    for i in range(len(cond)):
        if cond[i][0] == cond[i-1][0]:
            block_size +=1
        else:
            block_sizes.append(block_size)
            block_size = 1
    block_sizes.append(block_size)
    del block_sizes[0]
    Aview = {}
    current_start = 0
    row_index = 0
    column_index = 0
    Aview[row_index, column_index] = A[0:block_sizes[0],0:block_sizes[0]]
    for i,current_block, next_block in neighbour_zero(block_sizes):
        next_start = current_start + current_block
        next_end = current_start + current_block +next_block
        Aview[row_index+1, column_index+1] = A[next_start:next_end,next_start:next_end]
        Aview[row_index, column_index+1] = A[current_start:next_start, next_start:next_end]
        Aview[row_index+1, column_index] = A[next_start:next_end, current_start:next_start]
        current_start +=current_block
        row_index +=1
        column_index +=1
    return block_sizes, Aview

def RRGM(blocks, Aview):
    gl = [Aview[0,0].todense().I]
    iterator = range(1,len(blocks)).__iter__()
    prev_greensfunction = gl[0]
    for i in iterator:
        prev_greensfunction = (Aview[i,i]-Aview[i, i-1] * prev_greensfunction * Aview[i-1,i]).I
        gl.append(prev_greensfunction)
    GR = {}
    GR[len(blocks)-1,len(blocks)-1] = gl[-1]
    rev_iterator = reversed(range(1,len(blocks)))
    for i in rev_iterator:
        print i
        GR[i, i-1] = -GR[i,i] * Aview[i,i-1] * gl[i-1]
        GR[i-1, i] = -gl[i-1] * Aview[i-1,i] * GR[i,i]
        GR[i-1, i-1] = gl[i-1]-gl[i-1] * Aview[i-1,i] * GR[i, i-1]
    return gl, GR

def LRGM(blocks, Aview):
   pass
timeittook=time.time()
d = Device()
c = Contact()

cont, cond = d.read_geometry()
nodes = d.compose_geo()
s = c.build_SIGMA(nodes, cont)
A = build_EF(nodes) - d.HD() - s
b, Aview = blocks(cond, A)
#plt.imshow(scipy.asarray(HD.todense()).real)
#plt.show()

print "The Process took", time.time()-timeittook, "seconds"
#Nodes2 = compose_geo2(Cond)
#def main():
#    Cont,Cond = read_geometry(atlas)            
#    Nodes = compose_geo(Cond)  
#
#if __name__ == '__main__':
#    main()
#    
